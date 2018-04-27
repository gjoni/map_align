/*
 * mp_super.cpp
 *
 *  Created on: Apr 4, 2018
 *      Author: aivan
 */

#include <unistd.h>
#include <ctime>
#include <cassert>
#include <cmath>
#include <cstring>

#include <string>
#include <vector>

#include <omp.h>

#include "Chain.h"
#include "TMalign.h"
#include "kdtree.h"

#define VERSION "V20180405"
#define DMAX 8.0

struct OPTS {
	std::string list; /* list of PDB files */
	std::string cont; /* contacts file */
	std::string bfac; /* B-factors file */
	int dim; /* number of residues */
	unsigned N;
	int nthreads; /* number of threads to use */
	double dmax; /* distance cut-off - for trimming*/
};

bool GetOpts(int argc, char *argv[], OPTS &opts);
void PrintOpts(const OPTS &opts);
void PrintCap(const OPTS &opts);

/* TODO:
 * 		make TMscore(..) function return ALISTAT object
 */
struct ALISTAT {
	std::vector<bool> pair;
	std::vector<double> dist;
	double t[3];
	double u[3][3];
	double tm;
	int i, j;
};

ALISTAT TMscore(const Chain&, const Chain&, const int, bool*);

int main(int argc, char *argv[]) {

	/*
	 * (0) process input parameters
	 */
	OPTS opts = { "", "", "", -1, 0, 1, 0.0 };
	if (!GetOpts(argc, argv, opts)) {
		PrintOpts(opts);
		return 1;
	}

#if defined(_OPENMP)
	omp_set_num_threads(opts.nthreads);
#endif

	/*
	 * (1) read the list
	 */
	std::vector<Chain> chains;
	std::vector<std::string> in_ids;
	std::vector<std::string> out_ids;
	{
		FILE *F = fopen(opts.list.c_str(), "r");
		if (F == NULL) {
			printf("Error: cannot open list file '%s' (-l)\n",
					opts.list.c_str());
			exit(1);
		}

		char buf[1024], in[1024], out[1024];
		while (fgets(buf, 1024, F) != NULL) {
			if (sscanf(buf, "%s %s\n", in, out) == 2) {
				in_ids.push_back(in);
				out_ids.push_back(out);
				chains.push_back(Chain(in));
			} else if (sscanf(buf, "%s\n", in) == 1) {
				in_ids.push_back(in);
				out_ids.push_back("");
				chains.push_back(Chain(in));
			}
		}
		fclose(F);

		opts.N = chains.size();
	}

	/*
	 * (2) guesstimate protein length (if not set)
	 */
	if (opts.dim < 1) {
		for (auto &C : chains) {
			int dim = C.residue[C.nRes - 1].seqNum;
			opts.dim = dim > opts.dim ? dim : opts.dim;
		}
	} else {
		for (auto &C : chains) {
			int dim = C.residue[C.nRes - 1].seqNum;
			assert(dim <= opts.dim); /* user-defined dimension too small */
		}
	}

	PrintCap(opts);

	/* print out processed chains */
	for (unsigned i = 0; i < chains.size(); i++) {
		printf("# %u %s - %d aa", i + 1, in_ids[i].c_str(), chains[i].nRes);
		if (out_ids[i] != "") {
			printf(" - %s", out_ids[i].c_str());
		}
		printf("\n");
	}

	/*
	 * (3) compare structures
	 */

	/* pairwise alignments */
	size_t n = chains.size();
	size_t nij = n * (n + 1) / 2;
	std::vector<ALISTAT> alistat(nij);
	double tm_avg = 0.0;

#if defined(_OPENMP)
#pragma omp parallel for reduction(+:tm_avg) schedule(dynamic)
#endif
	for (size_t ij = 0; ij < nij; ij++) {

		size_t i, j;
		{
			size_t ii = n * (n + 1) / 2 - 1 - ij;
			size_t K = floor((sqrt(8 * ii + 1) - 1) / 2);
			i = n - 1 - K;
			j = ij - n * i + i * (i + 1) / 2;
		}

		/* exclude self-matches */
		if (i == j) {
			alistat[ij].i = alistat[ij].j = 0;
			alistat[ij].tm = -1.0;
			double *t = alistat[ij].t;
			double (*u)[3] = alistat[ij].u;
			memset(t, 0, 3 * sizeof(double));
			memset(u, 0, 9 * sizeof(double));
			u[0][0] = u[1][1] = u[2][2] = 1.0;
			continue;
		}

		/* !!! chain[j] is aligned onto chain[i] !!! */
		alistat[ij] = TMscore(chains[j], chains[i], opts.dim, NULL);
		alistat[ij].i = j;
		alistat[ij].j = i;
		tm_avg += alistat[ij].tm;

	}

	tm_avg /= (nij - n);
	printf("# ----------\n"
			"# average pairwise TM-score: %.6f\n# ----------\n", tm_avg);

	/*
	 * (4) analyze alignments
	 */

	/* matrix for storing combined results*/
	double **mtx = (double**) malloc(opts.dim * sizeof(double*));
	for (int i = 0; i < opts.dim; i++) {
		mtx[i] = (double*) calloc(3, sizeof(double));
	}

	/* residue presence frequency */
	for (auto &C : chains) {
		for (int i = 0; i < C.nRes; i++) {
			mtx[C.residue[i].seqNum - 1][0] += 1.0;
		}
	}

	/* pair statistics */
	for (auto &res : alistat) {
		if (res.i != res.j) {
			for (int i = 0; i < opts.dim; i++) {
				mtx[i][1] += res.pair[i];
				mtx[i][2] += res.dist[i];
			}
		}
	}

	/* normalize */
	for (int i = 0; i < opts.dim; i++) {
		mtx[i][2] = mtx[i][1] > 0.5 ? mtx[i][2] / mtx[i][1] : 0.0;
		if (mtx[i][2] > 20.0) {
			mtx[i][2] = 20.0;
		}
		mtx[i][0] /= n;
		mtx[i][1] /= (nij - n);
	}

	/* trim & recalculate average TM-score */
	if (opts.dmax > 0.0) {

		/* mask residues with deviation > opts.dmax */
		bool *fl = (bool*) malloc(opts.dim * sizeof(bool));
		double cov = 0.0;
		for (int i = 0; i < opts.dim; i++) {
			fl[i] = (mtx[i][2] < opts.dmax);
			cov += fl[i];
		}

		double tm_avg2 = 0.0;

		/* realign */
#if defined(_OPENMP)
#pragma omp parallel for reduction(+:tm_avg2) schedule(dynamic)
#endif
		for (size_t ij = 0; ij < nij; ij++) {

			size_t i, j;
			{
				size_t ii = n * (n + 1) / 2 - 1 - ij;
				size_t K = floor((sqrt(8 * ii + 1) - 1) / 2);
				i = n - 1 - K;
				j = ij - n * i + i * (i + 1) / 2;
			}

			/* exclude self-matches */
			if (i == j) {
				continue;
			}

			ALISTAT res = TMscore(chains[j], chains[i], opts.dim, fl);
			tm_avg2 += res.tm;

		}

		tm_avg2 /= (nij - n);
		printf("# cov= %.3f   tm_new= %.3f   tm_old= %.3f\n", cov / opts.dim,
				tm_avg2, tm_avg);

		free(fl);

	}

//	for (int i = 0; i < opts.dim; i++) {
//		printf("B %4d %7.3f %7.3f %7.3f\n", i + 1, mtx[i][0], mtx[i][1],
//				mtx[i][2]);
//	}

	/*
	 * (5) output coordinates
	 */
	for (size_t i = 0; i < n; i++) {

		ALISTAT &res = alistat[i];

		if (out_ids[res.i] != "") {

			Chain &C = chains[res.i];
			for (int i = 0; i < C.nAtoms; i++) {
				Residue *R = C.atom[i]->residue;
				double temp = 0.0;
				if (mtx[R->seqNum - 1][1] > 0) {
					temp = mtx[R->seqNum - 1][2];
				}
				C.atom[i]->temp = temp;
			}
			printf("--> %s\n", out_ids[res.i].c_str());
			C.Transform(res.t, res.u);
			C.Save(out_ids[res.i].c_str());

		}

	}

	/*
	 * (6) output contacts
	 */

	if (opts.cont != "") {

		/* create file for output */
		FILE *F = fopen(opts.cont.c_str(), "w");
		if (F == NULL) {
			printf("Error: cannot open '%s' file for writing (-c)\n",
					opts.cont.c_str());
			exit(1);
		}

		/* resize mtx[][] */
		for (int i = 0; i < opts.dim; i++) {
			free(mtx[i]);
			mtx[i] = (double*) calloc(opts.dim, sizeof(double));
		}

		/* loop over every chain and save contacts */
		for (auto &C : chains) {

			bool **fl = (bool**) malloc(opts.dim * sizeof(bool*));
			for (int i = 0; i < opts.dim; i++) {
				fl[i] = (bool*) calloc(opts.dim, sizeof(bool));
			}

			/* find contacting residues */
			kdres *res;
			double pos[3];
			for (int i = 0; i < C.nAtoms; i++) {

				Atom *A = C.atom[i];
				int a = A->residue->seqNum - 1;

				if (A->type == 'H') { /* exclude hydrogens */
					continue;
				}

				res = kd_nearest_range3f(C.kd, A->x, A->y, A->z, DMAX);
				while (!kd_res_end(res)) {

					Atom *B = *((Atom**) kd_res_item(res, pos));

					if (B->type == 'H') { /* exclude hydrogens */
						kd_res_next(res);
						continue;
					}

					int b = B->residue->seqNum - 1;
					if (a > b) {
						fl[a][b] += 1.0;
						fl[b][a] += 1.0;
					}
					kd_res_next(res);
				}
				kd_res_free(res);

			}

			/* save residue contacts to mtx[][] */
			for (int i = 0; i < opts.dim; i++) {
				for (int j = 0; j < opts.dim; j++) {
					mtx[i][j] += fl[i][j];
				}
			}

			/* free */
			for (int i = 0; i < opts.dim; i++) {
				free(fl[i]);
			}
			free(fl);

		}

		/* save observed contacts to file */
		for (int i = 0; i < opts.dim; i++) {
			for (int j = i + 1; j < opts.dim; j++) {
				if (mtx[i][j] > 0.5) {
					fprintf(F, "%d %d 0 %.2f %.3f\n", i + 1, j + 1, DMAX,
							mtx[i][j] / chains.size());
				}
			}
		}

		fclose(F);

	}

	/*
	 * free
	 */
	for (int i = 0; i < opts.dim; i++) {
		free(mtx[i]);
	}
	free(mtx);

	return 0;

}

void PrintOpts(const OPTS &opts) {

	printf("\nUsage:   ./mp_super [-option] [argument]\n\n");
	printf("Options:  -l list.txt                    - input, required\n");
	printf("          -n number of residues          - input, optional\n");
	printf("          -c contacts.txt                - output, optional\n");
	printf("          -b bfactors.txt                - output, optional\n");
	printf("          -d distance cut-off              %.1f\n", opts.dmax);
	printf("          -t number of threads             %d\n\n", opts.nthreads);

}

void PrintCap(const OPTS &opts) {

	time_t timer;
	time(&timer);
	struct tm* tm_info = localtime(&timer);
	char buf[100];
	strftime(buf, 26, "%Y:%m:%d / %H:%M:%S", tm_info);

	printf("# %s\n", std::string(70, '-').c_str());
	printf("# mp_super - a program to superimpose protein structures %15s\n",
	VERSION);
	printf("# %s\n", std::string(70, '-').c_str());

	printf("# %20s : %s\n", "start date/time", buf);
	printf("# %20s : %s\n", "list file", opts.list.c_str());
	printf("# %20s : %u\n", "number of structures", opts.N);
	printf("# %20s : %d\n", "number of residues", opts.dim);
	printf("# %20s : %s\n", "contacts file", opts.cont.c_str());
	printf("# %20s : %s\n", "B-factors file", opts.bfac.c_str());
	printf("# %20s : %.1f\n", "distance cut-off", opts.dmax);
	printf("# %20s : %d\n", "threads", opts.nthreads);

	printf("# %s\n", std::string(70, '-').c_str());

}

bool GetOpts(int argc, char *argv[], OPTS &opts) {

	char tmp;
	while ((tmp = getopt(argc, argv, "hl:c:b:n:t:d:")) != -1) {
		switch (tmp) {
		case 'h': /* help */
			printf("!!! HELP !!!\n");
			return false;
			break;
		case 'c': /* contacts file */
			opts.cont = std::string(optarg);
			break;
		case 'b': /* B-factors file */
			opts.bfac = std::string(optarg);
			break;
		case 'l': /* list file of PDBs */
			opts.list = std::string(optarg);
			break;
		case 'n': /* number of residues */
			opts.dim = atoi(optarg);
			break;
		case 't': /* number of threads */
			opts.nthreads = atoi(optarg);
			break;
		case 'd': /* distance cut-off */
			opts.dmax = atof(optarg);
			break;
		default:
			return false;
			break;
		}
	}

	if (opts.list == "") {
		printf("Error: list file not specified ('-l')\n");
		return false;
	}

	return true;

}

ALISTAT TMscore(const Chain &A, const Chain &B, const int dim, bool *fl) {

	ALISTAT res = { std::vector<bool>(dim, 0), std::vector<double>(dim, 0.0), {
			0.0, 0.0, 0.0 }, { { 1.0, 0.0, 0.0 }, { 0.0, 1.0, 0.0 }, { 0.0, 0.0,
			1.0 } }, 0.0, 0, 0 };

	std::vector<int> flag(dim, 0);
	std::vector<int> map;

	for (int i = 0; i < A.nRes; i++) {
		flag[A.residue[i].seqNum - 1]++;
	}

	for (int i = 0; i < B.nRes; i++) {
		flag[B.residue[i].seqNum - 1]++;
	}

	int nali = 0;

	if (fl != NULL) {
		for (int i = 0; i < dim; i++) {
			if (flag[i] == 2 && fl[i] == true) {
				flag[i]++;
				nali++;
				map.push_back(i);
				res.pair[i] = true;
			}
		}
	} else {
		for (int i = 0; i < dim; i++) {
			if (flag[i] == 2) {
				nali++;
				map.push_back(i);
				res.pair[i] = true;
			}
		}

	}

	if (nali < 5) {
		return res;
	}

	/* allocate memory */
	double **x = (double**) malloc(nali * sizeof(double*));
	double **y = (double**) malloc(nali * sizeof(double*));
	for (int i = 0; i < nali; i++) {
		x[i] = (double*) malloc(3 * sizeof(double));
		y[i] = (double*) malloc(3 * sizeof(double));
	}

	/* get aligned residues */
	{
		double **xyz_ptr = x;
		for (int i = 0; i < A.nRes; i++) {
			const Residue &R = A.residue[i];
			if ((flag[R.seqNum - 1] == 2 && fl == NULL)
					|| (flag[R.seqNum - 1] == 3 && fl != NULL)) {
				(*xyz_ptr)[0] = R.CA->x;
				(*xyz_ptr)[1] = R.CA->y;
				(*xyz_ptr)[2] = R.CA->z;
				xyz_ptr++;
			}
		}

		xyz_ptr = y;
		for (int i = 0; i < B.nRes; i++) {
			const Residue &R = B.residue[i];
			if ((flag[R.seqNum - 1] == 2 && fl == NULL)
					|| (flag[R.seqNum - 1] == 3 && fl != NULL)) {
				(*xyz_ptr)[0] = R.CA->x;
				(*xyz_ptr)[1] = R.CA->y;
				(*xyz_ptr)[2] = R.CA->z;
				xyz_ptr++;
			}
		}
	}

	/*
	 * align
	 */
	TMalign TM;
	res.tm = TM.GetTMscore(x, y, nali);

	/* transform y[..][3] coordinates */
	TM.GetTU(res.t, res.u);
	for (int i = 0; i < nali; i++) {

		double x_ = x[i][0];
		double y_ = x[i][1];
		double z_ = x[i][2];

		x[i][0] = res.u[0][0] * x_ + res.u[0][1] * y_ + res.u[0][2] * z_
				+ res.t[0];
		x[i][1] = res.u[1][0] * x_ + res.u[1][1] * y_ + res.u[1][2] * z_
				+ res.t[1];
		x[i][2] = res.u[2][0] * x_ + res.u[2][1] * y_ + res.u[2][2] * z_
				+ res.t[2];

		x_ = x[i][0] - y[i][0];
		y_ = x[i][1] - y[i][1];
		z_ = x[i][2] - y[i][2];

		res.dist[map[i]] = sqrt(x_ * x_ + y_ * y_ + z_ * z_);

	}

	/*
	 * free
	 */
	for (int i = 0; i < nali; i++) {
		free(x[i]);
		free(y[i]);
	}
	free(x);
	free(y);

	return res;

}
