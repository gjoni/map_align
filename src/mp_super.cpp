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

#include <string>
#include <vector>

#include "Chain.h"
#include "TMalign.h"

#define VERSION "V20180404"

struct OPTS {
	std::string list; /* list of PDB files */
	std::string cont; /* contacts file */
	std::string bfac; /* B-factors file */
	int dim; /* number of residues */
	unsigned N;
};

bool GetOpts(int argc, char *argv[], OPTS &opts);
void PrintOpts(const OPTS &opts);
void PrintCap(const OPTS &opts);

double TMscore(const Chain&, const Chain&, const int, double**);

int main(int argc, char *argv[]) {

	/*
	 * (0) process input parameters
	 */
	OPTS opts = { "", "", "", -1, 0 };
	if (!GetOpts(argc, argv, opts)) {
		PrintOpts(opts);
		return 1;
	}

	/*
	 * (1) read the list
	 */
	std::vector<Chain> chains;
	{
		FILE *F = fopen(opts.list.c_str(), "r");
		if (F == NULL) {
			printf("Error: cannot open list file '%s' (-l)\n",
					opts.list.c_str());
			exit(1);
		}

		char buf[1024];
		while (fscanf(F, "%s\n", buf) == 1) {
			printf("--> %s\n", buf);
			chains.push_back(Chain(buf));
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

	/*
	 * (3) compare structures
	 */
	double **mtx = (double**) malloc(opts.dim * sizeof(double*));
	for (int i = 0; i < opts.dim; i++) {
		mtx[i] = (double*) calloc(4, sizeof(double));
	}

	/* residue presence frequency */
	for (auto &C : chains) {
		for (int i = 0; i < C.nRes; i++) {
			mtx[C.residue[i].seqNum - 1][0] += 1.0;
		}
	}

	for (unsigned i = 0; i < chains.size(); i++) {
		for (unsigned j = i + 1; j < chains.size(); j++) {
			double tm = TMscore(chains[i], chains[j], opts.dim, mtx);
			printf("# tm(%u,%u)= %f\n", i + 1, j + 1, tm);
		}
	}

	for (int i = 0; i < opts.dim; i++) {
		printf("%d %f %f %f %f\n", i + 1, mtx[i][0], mtx[i][1], mtx[i][2],
				mtx[i][2] > 0 ? mtx[i][3] / mtx[i][2] : 0.0);
	}

	{
		Chain &C = chains[0];
		for (int i = 0; i < C.nAtoms; i++) {
			Residue *R = C.atom[i]->residue;
			double temp = 0.0;
			if (mtx[R->seqNum - 1][2] > 0) {
				temp = mtx[R->seqNum - 1][3] / mtx[R->seqNum - 1][2];
			}
			C.atom[i]->temp = temp;
		}
		C.Save("model.pdb");
	}

//	double tm = TMscore(chains[0], chains[1], opts.dim, mtx);
//	printf("tm= %f\n", tm);

	/*
	 * TODO:
	 * 		1) residue presence frequency
	 * 		2) TM-aligned frequency
	 * 		3) pair presence frequency
	 * 		4) average distance (B-factors)
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
	printf("          -c contacts.txt                - output, optional\n\n");
	printf("          -b bfactors.txt                - output, optional\n\n");

}

void PrintCap(const OPTS &opts) {

	time_t timer;
	time(&timer);
	struct tm* tm_info = localtime(&timer);
	char buf[100];
	strftime(buf, 26, "%Y:%m:%d / %H:%M:%S", tm_info);

	printf("# %s\n", std::string(70, '-').c_str());
	printf("# mp_super - a program to align protein contact maps %19s\n",
	VERSION);
	printf("# %s\n", std::string(70, '-').c_str());

	printf("# %20s : %s\n", "start date/time", buf);
	printf("# %20s : %s\n", "list file", opts.list.c_str());
	printf("# %20s : %u\n", "number of structures", opts.N);
	printf("# %20s : %d\n", "number of residues", opts.dim);
	printf("# %20s : %s\n", "contacts file", opts.cont.c_str());
	printf("# %20s : %s\n", "B-factors file", opts.bfac.c_str());

	printf("# %s\n", std::string(70, '-').c_str());

}

bool GetOpts(int argc, char *argv[], OPTS &opts) {

	char tmp;
	while ((tmp = getopt(argc, argv, "hl:c:b:n:")) != -1) {
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

double TMscore(const Chain &A, const Chain &B, const int dim, double **mtx) {

	double tm = 0.0;

	std::vector<int> flag(dim, 0);
	std::vector<int> map;

	for (int i = 0; i < A.nRes; i++) {
		flag[A.residue[i].seqNum - 1]++;
	}

	for (int i = 0; i < B.nRes; i++) {
		flag[B.residue[i].seqNum - 1]++;
	}

	int nali = 0;
	for (int i = 0; i < dim; i++) {
		if (flag[i] == 2) {
			nali++;
			map.push_back(i);
		}
	}

	if (nali < 5) {
		return tm;
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
			if (flag[R.seqNum - 1] == 2) {
				(*xyz_ptr)[0] = R.CA->x;
				(*xyz_ptr)[1] = R.CA->y;
				(*xyz_ptr)[2] = R.CA->z;
				xyz_ptr++;
			}
		}

		xyz_ptr = y;
		for (int i = 0; i < B.nRes; i++) {
			const Residue &R = B.residue[i];
			if (flag[R.seqNum - 1] == 2) {
				(*xyz_ptr)[0] = R.CA->x;
				(*xyz_ptr)[1] = R.CA->y;
				(*xyz_ptr)[2] = R.CA->z;
				xyz_ptr++;
			}
		}
	}

//	for (int i = 0; i < nali; i++) {
//		printf("%8.3f %8.3f %8.3f    %8.3f %8.3f %8.3f\n", x[i][0], x[i][1],
//				x[i][2], y[i][0], y[i][1], y[i][2]);
//	}

	/*
	 * align
	 */
	TMalign TM;
	tm = TM.GetTMscore(x, y, nali);

	/* get alignment */
	int *ali = (int*) malloc(nali * sizeof(int));
	TM.GetAliX2Y(ali, nali);

	/* pair presence frequency */
	for (int i = 0; i < nali; i++) {
		mtx[map[i]][1] += 1.0;
		mtx[map[i]][2] += (ali[i] >= 0);
	}

	/* transform y[..][3] coordinates */
	{
		double t[3], u[3][3];
		TM.GetTU(t, u);
		for (int i = 0; i < nali; i++) {

			double x_ = x[i][0];
			double y_ = x[i][1];
			double z_ = x[i][2];

			x[i][0] = u[0][0] * x_ + u[0][1] * y_ + u[0][2] * z_ + t[0];
			x[i][1] = u[1][0] * x_ + u[1][1] * y_ + u[1][2] * z_ + t[1];
			x[i][2] = u[2][0] * x_ + u[2][1] * y_ + u[2][2] * z_ + t[2];

			x_ = x[i][0] - y[i][0];
			y_ = x[i][1] - y[i][1];
			z_ = x[i][2] - y[i][2];

			double d = sqrt(x_ * x_ + y_ * y_ + z_ * z_);

			mtx[map[i]][3] += d;

		}
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
	free(ali);

	return tm;

}
