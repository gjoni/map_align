/*
 Original disclaimer:
 ===============================================================================
 Implementation of TM-align in C/C++ (V20120126)

 This program is written by Jianyi Yang at
 Yang Zhang lab
 Center for Computational Medicine and Bioinformatics
 University of Michigan
 100 Washtenaw Avenue, Ann Arbor, MI 48109-2218


 Please report bugs and questions to yangji@umich.edu or zhng@umich.edu
 ===============================================================================

 [1] Y Zhang, J Skolnick, TM-align: A protein structure alignment algorithm
 based on TM-score. Nucleic Acids Res(2005). 33:2302-9

 */

/*
 Modified into a stand-alone C++ class:
 ******************************************************************************
 *
 * 2013-2016, Ivan Anishchenko <anishchenko.ivan@gmail.com>
 * Vakser Lab, Center for Computational biology, The University of Kansas
 *
 * V20160328
 *
 *****************************************************************************/

#include <cmath>
#include <cstring>
#include <cstdio>
#include <cstdlib>

#include "TMalign.h"
#include "Kabsch.h"

#define D_CA  4.25 /* Max. distance between neighboring CA atoms (to detect broken chains) */
#define LN2   1.4426950409 /* One over natural LOG(2) */
#define sqrt3 1.7320508076 /* Square root of 3 */

double TMalign::dist(double x[], double y[]) {

	double dx = x[0] - y[0];
	double dy = x[1] - y[1];
	double dz = x[2] - y[2];

	return (dx * dx + dy * dy + dz * dz);

}

void TMalign::transform(double t[3], double u[3][3], double x[], double x1[]) {

	x1[0] = t[0] + (u[0][0] * x[0] + u[0][1] * x[1] + u[0][2] * x[2]);
	x1[1] = t[1] + (u[1][0] * x[0] + u[1][1] * x[1] + u[1][2] * x[2]);
	x1[2] = t[2] + (u[2][0] * x[0] + u[2][1] * x[1] + u[2][2] * x[2]);

}

void TMalign::do_rotation(double **x, double **x1, int len, double t[3],
		double u[3][3]) {

	for (int i = 0; i < len; i++) {
		transform(t, u, x[i], x1[i]);
	}

}

void TMalign::set4search(int xlen, int ylen) {

	d0_min = 0.5;
	Lnorm = getmin(xlen, ylen);
	if (Lnorm <= 19) {
		d0 = 0.168;
	} else {
		d0 = (1.24 * pow((Lnorm * 1.0 - 15), 1.0 / 3) - 1.8);
	}

	d0 += 0.8; /* min d0 in search = 0.968 */

	d02 = d0 * d0;

	d0_search = d0;
	if (d0_search > 8.0)
		d0_search = 8.0;
	if (d0_search < 4.5)
		d0_search = 4.5;

	d8_search = 1.5 * pow(Lnorm * 1.0, 0.3) + 3.5;
	d8_search_cut = d8_search * d8_search;

	simplify_step = 40;
	score_sum_method = 8;

	ddcc = 0.4;
	if (Lnorm <= 40) {
		ddcc = 0.1;
	}

}

void TMalign::set4final(int len) {

	d0_min = 0.5;
	Lnorm = len;
	if (Lnorm <= 21) {
		d0 = 0.5;
	} else {
		d0 = (1.24 * pow((Lnorm * 1.0 - 15), 1.0 / 3) - 1.8);
	}

	d02 = d0 * d0;

	d0_search = d0;
	if (d0_search > 8.0)
		d0_search = 8.0;
	if (d0_search < 4.5)
		d0_search = 4.5;

	simplify_step = 1;
	score_sum_method = 0;

}

//1->coil, 2->helix, 3->turn, 4->strand
void TMalign::make_sec(double **x, int len, int sec[]) {
	int j1, j2, j3, j4, j5;
	double d13, d14, d15, d24, d25, d35;
	for (int i = 0; i < len; i++) {
		sec[i] = 1;
		j1 = i - 2;
		j2 = i - 1;
		j3 = i;
		j4 = i + 1;
		j5 = i + 2;

		if (j1 >= 0 && j5 < len) {
			d13 = sqrt(dist(x[j1], x[j3]));
			d14 = sqrt(dist(x[j1], x[j4]));
			d15 = sqrt(dist(x[j1], x[j5]));
			d24 = sqrt(dist(x[j2], x[j4]));
			d25 = sqrt(dist(x[j2], x[j5]));
			d35 = sqrt(dist(x[j3], x[j5]));
			sec[i] = sec_str(d13, d14, d15, d24, d25, d35);
		}
	}
	smooth(sec, len);
}

int TMalign::sec_str(double dis13, double dis14, double dis15, double dis24,
		double dis25, double dis35) {
	int s = 1;

	double delta = 2.1;
	if (fabs(dis15 - 6.37) < delta) {
		if (fabs(dis14 - 5.18) < delta) {
			if (fabs(dis25 - 5.18) < delta) {
				if (fabs(dis13 - 5.45) < delta) {
					if (fabs(dis24 - 5.45) < delta) {
						if (fabs(dis35 - 5.45) < delta) {
							s = 2; //helix
							return s;
						}
					}
				}
			}
		}
	}

	delta = 1.42;
	if (fabs(dis15 - 13) < delta) {
		if (fabs(dis14 - 10.4) < delta) {
			if (fabs(dis25 - 10.4) < delta) {
				if (fabs(dis13 - 6.1) < delta) {
					if (fabs(dis24 - 6.1) < delta) {
						if (fabs(dis35 - 6.1) < delta) {
							s = 4; //strand
							return s;
						}
					}
				}
			}
		}
	}

	if (dis15 < 8) {
		s = 3; //turn
	}

	return s;
}

void TMalign::smooth(int sec[], int len) {
	int i, j;
	//smooth single  --x-- => -----
	for (i = 2; i < len - 2; i++) {
		if (sec[i] == 2 || sec[i] == 4) {
			j = sec[i];
			if (sec[i - 2] != j) {
				if (sec[i - 1] != j) {
					if (sec[i + 1] != j) {
						if (sec[i + 2] != j) {
							sec[i] = 1;
						}
					}
				}
			}
		}
	}

	//   smooth double
	//   --xx-- => ------

	for (i = 0; i < len - 5; i++) {
		//helix
		if (sec[i] != 2) {
			if (sec[i + 1] != 2) {
				if (sec[i + 2] == 2) {
					if (sec[i + 3] == 2) {
						if (sec[i + 4] != 2) {
							if (sec[i + 5] != 2) {
								sec[i + 2] = 1;
								sec[i + 3] = 1;
							}
						}
					}
				}
			}
		}

		//beta
		if (sec[i] != 4) {
			if (sec[i + 1] != 4) {
				if (sec[i + 2] == 4) {
					if (sec[i + 3] == 4) {
						if (sec[i + 4] != 4) {
							if (sec[i + 5] != 4) {
								sec[i + 2] = 1;
								sec[i + 3] = 1;
							}
						}
					}
				}
			}
		}
	}

	//smooth connect
	for (i = 0; i < len - 2; i++) {
		if (sec[i] == 2) {
			if (sec[i + 1] != 2) {
				if (sec[i + 2] == 2) {
					sec[i + 1] = 2;
				}
			}
		} else if (sec[i] == 4) {
			if (sec[i + 1] != 4) {
				if (sec[i + 2] == 4) {
					sec[i + 1] = 4;
				}
			}
		}
	}

}

double TMalign::get_score_fast(double **x, double **y, int ylen, int y2x[]) {

	double rms, tmscore, tmscore1, tmscore2;
	int i, j, k = 0;
	for (j = 0; j < ylen; j++) {
		i = y2x[j];
		if (i >= 0) {
			memcpy(r1[k], x[i], 3 * sizeof(double));
			memcpy(r2[k], y[j], 3 * sizeof(double));
			memcpy(xtm[k], x[i], 3 * sizeof(double));
			memcpy(ytm[k], y[j], 3 * sizeof(double));
			k++;
		}
	}
	if (k < 4) {
		return 0;
	}
	Kabsch(r1, r2, k, 1, &rms, t, u);

	//evaluate score
	double di;
	double *dis = (double*) malloc(k * sizeof(double));
	double d02 = d0 * d0;

	int n_ali = k;
	double xrot[3];
	tmscore = 0.0;
	for (k = 0; k < n_ali; k++) {
		transform(t, u, &xtm[k][0], xrot);
		di = dist(xrot, &ytm[k][0]);
		dis[k] = di;
		tmscore += d02 / (d02 + di);
	}

	//second iteration
	double ddMax_2nd = d0_search * d0_search;
	do {
		j = 0;
		for (k = 0; k < n_ali; k++) {
			if (dis[k] <= ddMax_2nd) {
				memcpy(r1[j], xtm[k], 3 * sizeof(double));
				memcpy(r2[j], ytm[k], 3 * sizeof(double));
				j++;
			}
		}

		//there are not enough feasible pairs, relieve the threshold
		ddMax_2nd += 0.5;
	} while (j < 4);

	if (n_ali > j) {
		Kabsch(r1, r2, j, 1, &rms, t, u);
		tmscore1 = 0;
		for (k = 0; k < n_ali; k++) {
			transform(t, u, &xtm[k][0], xrot);
			di = dist(xrot, &ytm[k][0]);
			dis[k] = di;
			tmscore1 += d02 / (d02 + di);
		}

		//third iteration
		double ddMax_3rd = ddMax_2nd + 1.0;
		do {
			j = 0;
			for (k = 0; k < n_ali; k++) {
				if (dis[k] <= ddMax_3rd) {
					memcpy(r1[j], xtm[k], 3 * sizeof(double));
					memcpy(r2[j], ytm[k], 3 * sizeof(double));
					j++;
				}
			}
			//there are not enough feasible pairs, relieve the threshold
			ddMax_3rd += 0.5;
		} while (j < 4);

		//evaluate the score
		Kabsch(r1, r2, j, 1, &rms, t, u);
		tmscore2 = 0.0;
		for (k = 0; k < n_ali; k++) {
			transform(t, u, &xtm[k][0], xrot);
			di = dist(xrot, &ytm[k][0]);
			tmscore2 += d02 / (d02 + di);
		}
	} else {
		free(dis);
		return tmscore;
	}

	if (tmscore1 >= tmscore)
		tmscore = tmscore1;
	if (tmscore2 >= tmscore)
		tmscore = tmscore2;

	free(dis);

	return tmscore; // no need to normalize this score because it will not be used for latter scoring

}

double TMalign::get_initial(double **x, double **y, int xlen, int ylen,
		int y2x[]) {

	int min_len = getmin(xlen, ylen);

	int min_ali = min_len / 2; //minimum size of considered fragment
	if (min_ali < 5)
		min_ali = 5;
	int n1 = -ylen + min_ali;
	int n2 = xlen - min_ali;

	int i, j, k, k_best;
	double tmscore, tmscore_max = -1;

	k_best = n1;
	for (k = n1; k <= n2; k++) {
		//get the map
		for (j = 0; j < ylen; j++) {
			i = j + k;
			if (i >= 0 && i < xlen) {
				y2x[j] = i;
			} else {
				y2x[j] = -1;
			}
		}

		//evaluate the map quickly in three iterations
		tmscore = get_score_fast(x, y, ylen, y2x);
		if (tmscore >= tmscore_max) {
			tmscore_max = tmscore;
			k_best = k;
		}
	}

	//extract the best map
	for (j = 0; j < ylen; j++) {
		i = j + k_best;
		if (i >= 0 && i < xlen) {
			y2x[j] = i;
		} else {
			y2x[j] = -1;
		}
	}
	return tmscore_max;
}

void TMalign::get_initial_ss(double **x, double **y, int xlen, int ylen,
		int y2x[]) {

	//assign secondary structures
	make_sec(x, xlen, secx);
	make_sec(y, ylen, secy);

	double gap_open = -1.0;
	NWDP_TM(secx, secy, xlen, ylen, gap_open, y2x);

}

bool TMalign::get_initial_local(double **x, double **y, int xlen, int ylen,
		int y2x[]) {

	double GL, rmsd;
	double t[3];
	double u[3][3];

	double d01 = d0 + 1.5;
	double d02 = d01 * d01;

	double GLmax = 0;
	int n_frag; //length of fragment for superposition
	int ns = 20; //tail length to discard

	int *invmap = (int*) malloc((ylen + 1) * sizeof(int));

	int aL = getmin(xlen, ylen);
	if (aL > 250) {
		n_frag = 50;
	} else if (aL > 200) {
		n_frag = 40;
	} else if (aL > 150) {
		n_frag = 30;
	} else {
		n_frag = 20;
	}

	int smallest = aL / 3; // I change here from aL/2 to aL/3

	if (n_frag > smallest)
		n_frag = smallest;
	if (ns > smallest)
		ns = smallest;

	if (n_frag == 0) {
		free(invmap);
		return false;
	}

	int m1 = xlen - n_frag - ns;
	int m2 = ylen - n_frag - ns;

	bool flag = false;
	double gap_open = 0.0;

	memset(y2x, -1, ylen * sizeof(int));

	for (int i = ns - 1; i < m1; i = i + n_frag) //index starts from 0, different from FORTRAN
			{
		for (int j = ns - 1; j < m2; j = j + n_frag) {
			for (int k = 0; k < n_frag; k++) //fragment in y
					{
				memcpy(r1[k], x[k + i], 3 * sizeof(double));
				memcpy(r2[k], y[k + j], 3 * sizeof(double));
			}

			Kabsch(r1, r2, n_frag, 1, &rmsd, t, u);

			NWDP_TM(x, y, xlen, ylen, t, u, d02, gap_open, invmap);
			GL = get_score_fast(x, y, ylen, invmap);

			if (GL > GLmax) {
				GLmax = GL;
				memcpy(y2x, invmap, ylen * sizeof(int));
				flag = true;
			}
		}
	}

	free(invmap);

	return flag;

}

void TMalign::score_matrix_rmsd_sec(double **x, double **y, int xlen, int ylen,
		int y2x[]) {

	double t[3], u[3][3];
	double rmsd, dij;
	double d01 = d0 + 1.5;
	double d02 = d01 * d01;

	double xx[3];
	int i, k = 0;
	for (int j = 0; j < ylen; j++) {
		i = y2x[j];
		if (i >= 0) {
			memcpy(r1[k], x[i], 3 * sizeof(double));
			memcpy(r2[k], y[j], 3 * sizeof(double));
			k++;
		}
	}
	Kabsch(r1, r2, k, 1, &rmsd, t, u);

	for (int ii = 0; ii < xlen; ii++) {
		transform(t, u, &x[ii][0], xx);
		for (int jj = 0; jj < ylen; jj++) {
			dij = dist(xx, &y[jj][0]);
			if (secx[ii] == secy[jj]) {
				score[ii + 1][jj + 1] = 1.0 / (1 + dij / d02) + 0.5;
			} else {
				score[ii + 1][jj + 1] = 1.0 / (1 + dij / d02);
			}
		}
	}
}

void TMalign::get_initial_ssplus(double **x, double **y, int xlen, int ylen,
		int y2x0[], int y2x[]) {

	//create score matrix for DP
	score_matrix_rmsd_sec(x, y, xlen, ylen, y2x0);

	double gap_open = -1.0;
	NWDP_TM(xlen, ylen, gap_open, y2x);

}

void TMalign::find_max_frag(double **x, int len, int *start_max, int *end_max) {

	int fra_min = 4; //minimum fragment for search
	double d;
	int Lfr_max = 0;

	int r_min = len / 3; //minimum fragment, in case too small protein
	if (r_min > fra_min)
		r_min = fra_min;

	int inc = 0;
	double dcu_cut = D_CA * D_CA;

	int iPrev;
	while (Lfr_max < r_min) {
		Lfr_max = 0;
		iPrev = 1;
		for (int i = 1; i < len; ++i) {
			d = dist(x[i - 1], x[i]);
			if (d < dcu_cut) {
				iPrev++;
			} else {
				if (iPrev > Lfr_max) {
					Lfr_max = iPrev;
					*end_max = i - 1;
				}
				iPrev = 1;
			}
		}
		if (iPrev > Lfr_max) { // check last residue
			Lfr_max = iPrev;
			*end_max = len - 1;
		}

		if (Lfr_max < r_min) // need to relieve threshold ?
				{
			inc++;
			dcu_cut = D_CA + 0.1 * inc;
			dcu_cut = dcu_cut * dcu_cut;
		}

	}
	*start_max = *end_max - Lfr_max + 1;

}

double TMalign::get_initial_fgt(double **x, double **y, int xlen, int ylen,
		int y2x[]) {

	int fra_min = 4; //minimum fragment for search
	int fra_min1 = fra_min - 1; //cutoff for shift, save time

	int xstart = 0, ystart = 0, xend = 0, yend = 0;

	find_max_frag(x, xlen, &xstart, &xend);
	find_max_frag(y, ylen, &ystart, &yend);

	int Lx = xend - xstart + 1;
	int Ly = yend - ystart + 1;
	int L_fr = getmin(Lx, Ly);
	int *ifr = (int*) malloc(L_fr * sizeof(int));
	int *y2x_ = (int*) malloc((ylen + 1) * sizeof(int));

	//select what piece will be used (this may araise ansysmetry, but
	//only when L1=L2 and Lfr1=Lfr2 and L1 ne Lfr1
	//if L1=Lfr1 and L2=Lfr2 (normal proteins), it will be the same as initial1

	int start = 0;
	bool scenario1 = false;
	if (Lx < Ly || (Lx == Ly && xlen <= ylen)) {
		start = xstart;
		scenario1 = true;
	} else if (Lx > Ly || (Lx == Ly && xlen > ylen)) {
		start = ystart;
	}
	for (int i = 0; i < L_fr; i++) {
		ifr[i] = start + i;
	}

	int L0 = getmin(xlen, ylen); //non-redundant to get_initial1
	if (L_fr == L0) {
		int n1 = (int) (L0 * 0.1); //my index starts from 0
		int n2 = (int) (L0 * 0.89);

		int j = 0;
		for (int i = n1; i <= n2; i++) {
			ifr[j] = ifr[i];
			j++;
		}
		L_fr = j;
	}

	//gapless threading for the extracted fragment
	double tmscore, tmscore_max = -1;

	if (scenario1) {
		int L1 = L_fr;
		int min_len = getmin(L1, ylen);
		int min_ali = (int) (min_len / 2.5); //minimum size of considered fragment
		if (min_ali <= fra_min1)
			min_ali = fra_min1;
		int n1, n2;
		n1 = -ylen + min_ali;
		n2 = L1 - min_ali;

		int i, j, k;
		for (k = n1; k <= n2; k++) {
			//get the map
			for (j = 0; j < ylen; j++) {
				i = j + k;
				if (i >= 0 && i < L1) {
					y2x_[j] = ifr[i];
				} else {
					y2x_[j] = -1;
				}
			}

			//evaluate the map quickly in three iterations
			tmscore = get_score_fast(x, y, ylen, y2x_);

			if (tmscore >= tmscore_max) {
				tmscore_max = tmscore;
				memcpy(y2x, y2x_, ylen * sizeof(int));
			}
		}
	} else {
		int L2 = L_fr;
		int min_len = getmin(xlen, L2);
		int min_ali = (int) (min_len / 2.5); //minimum size of considered fragment
		if (min_ali <= fra_min1)
			min_ali = fra_min1;
		int n1, n2;
		n1 = -L2 + min_ali;
		n2 = xlen - min_ali;

		int i, j, k;

		for (k = n1; k <= n2; k++) {
			//get the map
			memset(y2x_, -1, ylen * sizeof(int));

			for (j = 0; j < L2; j++) {
				i = j + k;
				if (i >= 0 && i < xlen) {
					y2x_[ifr[j]] = i;
				}
			}

			//evaluate the map quickly in three iterations
			tmscore = get_score_fast(x, y, ylen, y2x_);
			if (tmscore >= tmscore_max) {
				tmscore_max = tmscore;
				memcpy(y2x, y2x_, ylen * sizeof(int));
			}
		}
	}

	free(ifr);
	free(y2x_);

	return tmscore_max;
}

int TMalign::score_fun8(double **x, double **y, int Lali, double d, int y2x[],
		double &score, int score_sum_method) {

	double score_sum, di;
	double dd = d * d;

	int n_cut, inc = 0;

	do {
		n_cut = 0;
		score_sum = 0;
		for (int i = 0; i < Lali; ++i) {
			di = dist(x[i], y[i]);
			if (di < dd) {
				y2x[n_cut] = i;
				n_cut++;
			}
			if (score_sum_method == 8) {
				if (di <= d8_search_cut) {
					score_sum += d02 / (d02 + di);
				}
			} else {
				score_sum += d02 / (d02 + di);
			}
		}

		//there are not enough feasible pairs, relieve the threshold
		inc++;
		di = (d + inc * 0.5);
		dd = di * di;
	} while (n_cut < 3 && Lali > 3);

	score = score_sum / Lnorm;

	return n_cut;
}

double TMalign::TMscore8_search(double **xtm, double **ytm, int Lali,
		double t0[3], double u0[3][3], int simplify_step,
		int score_sum_method) {

	int i, m;
	double score_max, score, rmsd;
	int *k_ali = (int*) malloc(Lali * sizeof(int));
	int ka, k;
	double d;

	double t[3];
	double u[3][3];

	//iterative parameters
	int n_it = 20; //maximum number of iterations
	const int n_init_max = 6; //maximum number of different fragment length
	int L_ini[n_init_max]; //fragment lengths, Lali, Lali/2, Lali/4 ... 4
	int n_init;
	if (Lali <= 4) {
		n_init = 1;
		L_ini[0] = Lali;
	} else {
		double intpart;
		modf(log(Lali) * LN2, &intpart);
		n_init = int(intpart);
		n_init = n_init < n_init_max ? n_init : n_init_max;
		int two_to_the_i = 1;
		for (i = 0; i < n_init; ++i) {
			L_ini[i] = Lali / two_to_the_i;
			two_to_the_i *= 2;
		}
		L_ini[n_init - 1] = 4;
	}

	int i_init;

	score_max = -1;
	//find the maximum score starting from local structures superposition
	int *i_ali = (int*) malloc(Lali * sizeof(int));
	int n_cut;
	int L_frag; //fragment length
	int iL_max; //maximum starting position for the fragment
	int nStep;
	int it;
	bool fl;
	for (i_init = 0; i_init < n_init; i_init++) {
		L_frag = L_ini[i_init];
		iL_max = Lali - L_frag;
		nStep = iL_max / simplify_step + (iL_max % simplify_step != 0) + 1;
		i = 0;
		for (int iStep = 0; iStep < nStep; ++iStep) {
			//extract the fragment starting from position i

			/* MODIFIED */
			for (int ii = 0; ii < L_frag; ii++) {
				memcpy(r1[ii], xtm[i + ii], 3 * sizeof(double));
				memcpy(r2[ii], ytm[i + ii], 3 * sizeof(double));
			}

			//extract rotation matrix based on the fragment
			Kabsch(r1, r2, L_frag, 1, &rmsd, t, u);
			do_rotation(xtm, xt, Lali, t, u);

			//get subsegment of this fragment
			d = d0_search - 1.0;
			n_cut = score_fun8(xt, ytm, Lali, d, i_ali, score,
					score_sum_method);
			if (score > score_max) {
				score_max = score;

				//save the rotation matrix
				memcpy(t0, t, 3 * sizeof(double));
				memcpy(u0, u, 9 * sizeof(double));
			}

			//try to extend the alignment iteratively
			d = d0_search + 1.0;
			it = 0;
			fl = false;
			while (it < n_it && !fl) {
				it++;
				ka = 0;
				for (k = 0; k < n_cut; k++) {
					m = i_ali[k];

					memcpy(r1[k], xtm[m], 3 * sizeof(double));
					memcpy(r2[k], ytm[m], 3 * sizeof(double));

					k_ali[ka] = m;
					ka++;
				}
				//extract rotation matrix based on the fragment
				score = 0.0;
				if (n_cut > 3) {
					Kabsch(r1, r2, n_cut, 1, &rmsd, t, u);
					do_rotation(xtm, xt, Lali, t, u);
					n_cut = score_fun8(xt, ytm, Lali, d, i_ali, score,
							score_sum_method);
				}
				if (score > score_max) {
					score_max = score;

					//save the rotation matrix
					memcpy(t0, t, 3 * sizeof(double));
					memcpy(u0, u, 9 * sizeof(double));
				}

				//check if it converges
				if (n_cut == ka) {
					fl = (memcmp(i_ali, k_ali, n_cut * sizeof(int)) == 0);
				}
			}

			i = i + simplify_step; //shift the fragment
			if (i > iL_max)
				i = iL_max; //do this to use the last missed fragment
		}
	}

	free(i_ali);
	free(k_ali);

	return score_max;

}

double TMalign::detailed_search(double **x, double **y, int xlen, int ylen,
		int y2x[], double t[3], double u[3][3], int simplify_step,
		int score_sum_method) {
	int i, j, k = 0;
	for (i = 0; i < ylen; i++) {
		j = y2x[i];
		if (j >= 0) //aligned
				{
			memcpy(xtm[k], x[j], 3 * sizeof(double));
			memcpy(ytm[k], y[i], 3 * sizeof(double));
			k++;
		}
	}

	//detailed search 40-->1
	double tm = TMscore8_search(xtm, ytm, k, t, u, simplify_step,
			score_sum_method);

	return tm;
}

void TMalign::NWDP_Initialize(int &len1, int &len2, int j2i[]) {

	for (int i = 0; i <= len1; ++i) {
		val[i][0] = 0;
		path[i][0] = false;
	}
	for (int j = 0; j <= len2; ++j) {
		val[0][j] = 0;
	}
	memset(path[0], false, (len2 + 1) * sizeof(bool));
	memset(j2i, -1, (len2 + 1) * sizeof(int));

}

void TMalign::NWDP_TraceBack(int &len1, int &len2, double &gap_open,
		int j2i[]) {

	int i = len1;
	int j = len2;
	double h, v;
	while (i > 0 && j > 0) {
		if (path[i][j]) //from diagonal
		{
			j2i[j - 1] = i - 1;
			i--;
			j--;
		} else {
			h = val[i - 1][j];
			if (path[i - 1][j])
				h += gap_open;

			v = val[i][j - 1];
			if (path[i][j - 1])
				v += gap_open;

			if (v >= h)
				j--;
			else
				i--;
		}
	}
}

void TMalign::NWDP_TM(double **x, double **y, int len1, int len2, double t[3],
		double u[3][3], double d02, double gap_open, int j2i[]) {

	int i, j;
	double h, v, d;

	//initialization
	NWDP_Initialize(len1, len2, j2i);

	double xx[3], dij;

	//decide matrix and path
	for (i = 1; i <= len1; i++) {
		transform(t, u, x[i - 1], xx);
		for (j = 1; j <= len2; j++) {
			//d=val[i-1][j-1]+score[i][j]; //diagonal

			dij = dist(xx, y[j - 1]);
			d = val[i - 1][j - 1] + d02 / (d02 + dij);

			//symbol insertion in horizontal (= a gap in vertical)
			h = val[i - 1][j];
			if (path[i - 1][j]) //aligned in last position
				h += gap_open;

			//symbol insertion in vertical
			v = val[i][j - 1];
			if (path[i][j - 1]) //aligned in last position
				v += gap_open;

			if (d >= h && d >= v) {
				path[i][j] = true; //from diagonal
				val[i][j] = d;
			} else {
				path[i][j] = false; //from horizontal
				if (v >= h)
					val[i][j] = v;
				else
					val[i][j] = h;
			}
		} //for i
	} //for j

	//trace back to extract the alignment
	NWDP_TraceBack(len1, len2, gap_open, j2i);

}

void TMalign::NWDP_TM(int secx[], int secy[], int len1, int len2,
		double gap_open, int j2i[]) {

	int i, j;
	double h, v, d;

	//initialization
	NWDP_Initialize(len1, len2, j2i);

	//decide matrix and path
	for (i = 1; i <= len1; i++) {
		for (j = 1; j <= len2; j++) {
			if (secx[i - 1] == secy[j - 1]) {
				d = val[i - 1][j - 1] + 1.0;
			} else {
				d = val[i - 1][j - 1];
			}

			//symbol insertion in horizontal (= a gap in vertical)
			h = val[i - 1][j];
			if (path[i - 1][j]) //aligned in last position
				h += gap_open;

			//symbol insertion in vertical
			v = val[i][j - 1];
			if (path[i][j - 1]) //aligned in last position
				v += gap_open;

			if (d >= h && d >= v) {
				path[i][j] = true; //from diagonal
				val[i][j] = d;
			} else {
				path[i][j] = false; //from horizontal
				if (v >= h)
					val[i][j] = v;
				else
					val[i][j] = h;
			}
		} //for i
	} //for j

	//trace back to extract the alignment
	NWDP_TraceBack(len1, len2, gap_open, j2i);

}

void TMalign::NWDP_TM(int len1, int len2, double gap_open, int j2i[]) {

	int i, j;
	double h, v, d;

	//initialization
	NWDP_Initialize(len1, len2, j2i);

	//decide matrix and path
	for (i = 1; i <= len1; i++) {
		for (j = 1; j <= len2; j++) {
			d = val[i - 1][j - 1] + score[i][j]; //diagonal

			//symbol insertion in horizontal (= a gap in vertical)
			h = val[i - 1][j];
			if (path[i - 1][j]) //aligned in last position
				h += gap_open;

			//symbol insertion in vertical
			v = val[i][j - 1];
			if (path[i][j - 1]) //aligned in last position
				v += gap_open;

			if (d >= h && d >= v) {
				path[i][j] = true; //from diagonal
				val[i][j] = d;
			} else {
				path[i][j] = false; //from horizontal
				if (v >= h)
					val[i][j] = v;
				else
					val[i][j] = h;
			}
		} //for i
	} //for j

	//trace back to extract the alignment
	NWDP_TraceBack(len1, len2, gap_open, j2i);

}

double TMalign::DP_iter(double **x, double **y, int xlen, int ylen, double t[3],
		double u[3][3], int invmap0[], int g1, int g2, int iteration_max) {

	double gap_open[2] = { -0.6, 0 };
	int *invmap = (int*) malloc((ylen + 1) * sizeof(int));

	int iteration, i, j, k;
	double tmscore, tmscore_max, tmscore_old = 0;
	int score_sum_method = 8, simplify_step = 40;
	tmscore_max = -1;

	for (int g = g1; g < g2; g++) {
		for (iteration = 0; iteration < iteration_max; iteration++) {
			NWDP_TM(x, y, xlen, ylen, t, u, d02, gap_open[g], invmap);

			k = 0;
			for (j = 0; j < ylen; j++) {
				i = invmap[j];

				if (i >= 0) //aligned
						{
					memcpy(xtm[k], x[i], 3 * sizeof(double));
					memcpy(ytm[k], y[j], 3 * sizeof(double));
					k++;
				}
			}

			tmscore = TMscore8_search(xtm, ytm, k, t, u, simplify_step,
					score_sum_method); //, &rmsd);

			if (tmscore > tmscore_max) {
				tmscore_max = tmscore;
				memcpy(invmap0, invmap, ylen * sizeof(int));
			}

			if (iteration > 0) {
				if (fabs(tmscore_old - tmscore) < 0.000001) {
					break;
				}
			}
			tmscore_old = tmscore;
		} // for iteration
	} //for gapopen

	free(invmap);

	return tmscore_max;

}

void TMalign::Allocate(int xlen, int ylen) {

	score = (double**) malloc((xlen + 1) * sizeof(double*));
	path = (bool**) malloc((xlen + 1) * sizeof(bool*));
	val = (double**) malloc((xlen + 1) * sizeof(double*));
	for (int i = 0; i <= xlen; ++i) {
		score[i] = (double*) malloc((ylen + 1) * sizeof(double));
		path[i] = (bool*) malloc((ylen + 1) * sizeof(bool));
		val[i] = (double*) malloc((ylen + 1) * sizeof(double));
	}

	int len = xlen > ylen ? xlen : ylen;
	len++;
	xtm = (double**) malloc(len * sizeof(double*));
	ytm = (double**) malloc(len * sizeof(double*));
	xt = (double**) malloc(len * sizeof(double*));
	r1 = (double**) malloc(len * sizeof(double*));
	r2 = (double**) malloc(len * sizeof(double*));
	for (int i = 0; i < len; i++) {
		xtm[i] = (double*) malloc(3 * sizeof(double));
		ytm[i] = (double*) malloc(3 * sizeof(double));
		xt[i] = (double*) malloc(3 * sizeof(double));
		r1[i] = (double*) malloc(3 * sizeof(double));
		r2[i] = (double*) malloc(3 * sizeof(double));
	}

	secx = (int*) malloc(xlen * sizeof(int));
	secy = (int*) malloc(ylen * sizeof(int));
	invmap = (int*) malloc((ylen + 1) * sizeof(int));
	invmap0 = (int*) malloc((ylen + 1) * sizeof(int));

	memset(invmap, -1, ylen * sizeof(int));
	memset(invmap0, -1, ylen * sizeof(int));

}

void TMalign::Deallocate(int xlen, int ylen) {

	for (int i = 0; i <= xlen; ++i) {
		free(score[i]);
		free(path[i]);
		free(val[i]);
	}
	free(score);
	free(path);
	free(val);

	int len = xlen > ylen ? xlen : ylen;
	len++;
	for (int i = 0; i < len; i++) {
		free(xtm[i]);
		free(ytm[i]);
		free(xt[i]);
		free(r1[i]);
		free(r2[i]);
	}
	free(xtm);
	free(ytm);
	free(xt);
	free(r1);
	free(r2);

	free(secx);
	free(secy);
	free(invmap);
	free(invmap0);

}

double TMalign::GetTMscore(double **x, double **y, int xlen, int ylen,
		int mode) {

	Allocate(xlen, ylen);

	/* initialize parameters for search */
	set4search(xlen, ylen);
	//memset(invmap0, -1, (ylen + 1) * sizeof(int));
	double TM, TMmax = -1.0;

	/*
	 * (1) get initial alignment with gapless threading
	 */
	get_initial(x, y, xlen, ylen, invmap0);

	/* find the max TMscore for this initial alignment
	 * with the simplified search_engine */
	TM = detailed_search(x, y, xlen, ylen, invmap0, t, u, simplify_step,
			score_sum_method);
	if (TM > TMmax) {
		TMmax = TM;
	}

	/* run dynamic programming iteratively to find the best alignment */
	TM = DP_iter(x, y, xlen, ylen, t, u, invmap, 0, 2, 30);
	if (TM > TMmax) {
		TMmax = TM;
		memcpy(invmap0, invmap, ylen * sizeof(int));
	}

	/*
	 * (2) get initial alignment based on secondary structure
	 */
	get_initial_ss(x, y, xlen, ylen, invmap);
	TM = detailed_search(x, y, xlen, ylen, invmap, t, u, simplify_step,
			score_sum_method);
	if (TM > TMmax) {
		TMmax = TM;
		memcpy(invmap0, invmap, ylen * sizeof(int));
	}

	if (TM > TMmax * 0.2) {
		TM = DP_iter(x, y, xlen, ylen, t, u, invmap, 0, 2, 30);
		if (TM > TMmax) {
			TMmax = TM;
			memcpy(invmap0, invmap, ylen * sizeof(int));
		}
	}

	/*
	 * (3) get initial alignment based on local superposition
	 */
	if (get_initial_local(x, y, xlen, ylen, invmap)) {
		TM = detailed_search(x, y, xlen, ylen, invmap, t, u, simplify_step,
				score_sum_method);
		if (TM > TMmax) {
			TMmax = TM;
			memcpy(invmap0, invmap, ylen * sizeof(int));
		}
		if (TM > TMmax * ddcc) {
			TM = DP_iter(x, y, xlen, ylen, t, u, invmap, 0, 2, 2);
			if (TM > TMmax) {
				TMmax = TM;
				memcpy(invmap0, invmap, ylen * sizeof(int));
			}
		}
	}

	/*
	 * (4) get initial alignment based on previous alignment+secondary structure
	 */
	get_initial_ssplus(x, y, xlen, ylen, invmap0, invmap);
	TM = detailed_search(x, y, xlen, ylen, invmap, t, u, simplify_step,
			score_sum_method);
	if (TM > TMmax) {
		TMmax = TM;
		memcpy(invmap0, invmap, ylen * sizeof(int));
	}
	if (TM > TMmax * ddcc) {
		TM = DP_iter(x, y, xlen, ylen, t, u, invmap, 0, 2, 30);
		if (TM > TMmax) {
			TMmax = TM;
			memcpy(invmap0, invmap, ylen * sizeof(int));
		}
	}

	/*
	 * (5) get initial alignment based on fragment gapless threading
	 */
	get_initial_fgt(x, y, xlen, ylen, invmap);
	TM = detailed_search(x, y, xlen, ylen, invmap, t, u, simplify_step,
			score_sum_method);
	if (TM > TMmax) {
		TMmax = TM;
		memcpy(invmap0, invmap, ylen * sizeof(int));
	}
	if (TM > TMmax * ddcc) {
		TM = DP_iter(x, y, xlen, ylen, t, u, invmap, 1, 2, 2);
		if (TM > TMmax) {
			TMmax = TM;
			memcpy(invmap0, invmap, ylen * sizeof(int));
		}
	}

	/*
	 * no alignment was found
	 */
	bool flag = false;
	for (int i = 0; i < ylen; i++) {
		if (invmap0[i] >= 0) {
			flag = true;
			break;
		}
	}
	if (!flag) {
		return -1.0;
	}

	/*
	 * Detailed TMscore search engine  --> prepare for final TMscore
	 */
	simplify_step = 1;
	score_sum_method = 8;
	TM = detailed_search(x, y, xlen, ylen, invmap0, t, u, simplify_step,
			score_sum_method);

	/* select pairs with dis<d8 for final TMscore computation
	 * and output alignment */
	int n_ali8, k = 0, l;
	int n_ali = 0;
	double d;
	do_rotation(x, xt, xlen, t, u);

	for (int j = 0; j < ylen; j++) {
		l = invmap0[j];
		if (l >= 0) //aligned
				{
			n_ali++;
			d = sqrt(dist(&xt[l][0], &y[j][0]));
			if (d <= d8_search) {
				memcpy(xtm[k], x[l], 3 * sizeof(double));
				memcpy(ytm[k], y[j], 3 * sizeof(double));
				k++;
			}
		}
	}
	n_ali8 = k;

	/*
	 * Final TMscore
	 */
	if (mode == 1) {
		set4final(xlen); /* normalized by length of structure X */
	} else if (mode == 2) {
		set4final(ylen); /* normalized by length of structure Y */
	} else if (mode == 3) {
		set4final((xlen + ylen) / 2); /* average of the two */
	} else {
		set4final(xlen);
	}
	TM = TMscore8_search(xtm, ytm, n_ali8, t, u, simplify_step,
			score_sum_method);

	/*
	 * save alignment
	 */
	x2y_map.clear();
	x2y_map.resize(xlen, -1);
	for (int i = 0; i < ylen; ++i) {
		if (invmap0[i] >= 0) {
			x2y_map[invmap0[i]] = i;
		}
	}

	Deallocate(xlen, ylen);

	return TM;

}

double TMalign::GetTMscore(double **x, double **y, int xlen) {

	double ylen = xlen;

	Allocate(xlen, ylen);

	/* initialize parameters for search */
	set4search(xlen, ylen);

	for (int i = 0; i < xlen; ++i) {
		invmap0[i] = i;
	}

	double TM;

	simplify_step = 1;
	score_sum_method = 8;
	TM = detailed_search(x, y, xlen, ylen, invmap0, t, u, simplify_step,
			score_sum_method);

	/* select pairs with dis<d8 for final TMscore computation
	 * and output alignment */
	int n_ali8, k = 0, l;
	int n_ali = 0;
	double d;
	do_rotation(x, xt, xlen, t, u);

	for (int j = 0; j < ylen; j++) {
		l = invmap0[j];
		if (l >= 0) //aligned
				{
			n_ali++;
			d = sqrt(dist(&xt[l][0], &y[j][0]));
			if (d <= d8_search) {
				memcpy(xtm[k], x[l], 3 * sizeof(double));
				memcpy(ytm[k], y[j], 3 * sizeof(double));
				k++;
			}
		}
	}
	n_ali8 = k;

	/*
	 * Final TMscore
	 */
	set4final(xlen);
	TM = TMscore8_search(xtm, ytm, n_ali8, t, u, simplify_step,
			score_sum_method);

	/*
	 * save alignment
	 */
	x2y_map.clear();
	x2y_map.resize(xlen, -1);
	for (int i = 0; i < ylen; ++i) {
		if (invmap0[i] >= 0) {
			x2y_map[invmap0[i]] = i;
		}
	}

	Deallocate(xlen, ylen);

	return TM;

}

void TMalign::GetTU(double T[3], double U[3][3]) {
	memcpy(T, t, 3 * sizeof(double));
	memcpy(U, u, 9 * sizeof(double));
}

void TMalign::GetAliX2Y(int ali[], int xlen) {

	for (int i = 0; i < xlen; i++) {
		ali[i] = x2y_map[i];
	}

}

TMalign::TMalign() :
		d0(0.0), d02(0.0), d0_min(0.0), d0_search(0.0), d8_search(0.0), d8_search_cut(
				0.0), Lnorm(0), ddcc(0.0), simplify_step(0), score_sum_method(
				0), score(
		NULL), path(NULL), val(NULL), xtm(NULL), ytm(NULL), xt(NULL), secx(
		NULL), secy(NULL), r1(NULL), r2(NULL), invmap(NULL), invmap0(NULL) {

	/* empty */

}

TMalign::~TMalign() {

	/* empty */

}
