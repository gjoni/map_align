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

#ifndef TMALIGN_H_
#define TMALIGN_H_

#include <vector>

#define getmax(a,b) a>b?a:b
#define getmin(a,b) a<b?a:b

class TMalign {
private:

	/**************************************************************************
	 *                              VARIABLES                                 *
	 **************************************************************************/
	double d0; /* scale for defining TM-score */
	double d02; /* d02 = d0*d0 */
	double d0_min; /* minimum possible value of d0 */
	double d0_search; /* for quickly calculate TM-score in searching */
	double d8_search; /* remove pairs with dis>d8 during search & final */
	double d8_search_cut; /* d8_search_cut = d8_search*d8_search */
	double Lnorm; /* normalization length */

	double ddcc; /* Don't make refined search if TM-prev > TM-curr*ddcc */

	int simplify_step; /* SS: Check every SSs residue */
	int score_sum_method; /* if = 8 then the sum over pairs with dis<score_d8
	 is calculated */

	double **score; /* Scoring table for dynamic programming (DP) */
	bool **path; /* A table to track path for DP */
	double **val; /* DP matrix */

	/**************************************************************************
	 *   In general, Y is regarded as native structure: superpose X onto Y    *
	 **************************************************************************/
	double **xtm, **ytm; /* for TMscore search engine */
	double **xt; /* superposed version of r1 or xtm */
	int *secx, *secy; /* secondary structure: 1->coil, 2->helix, 3->turn, 4->strand */
	double **r1, **r2; /* for Kabsch rotation */

	double t[3], u[3][3]; /* Kabsch translation vectors and rotation matrixes */

	int *invmap, *invmap0;
	/* Arrays to store the alignments: y2x[j] = i means:
	 the jth element in y is aligned to the ith element in x if i>=0
	 the jth element in y is aligned to a gap in x if i==-1 */

	std::vector<int> x2y_map;

	/**************************************************************************
	 *                        AUXILIARY FUNCTIONS                             *
	 **************************************************************************/

	/* Squared(!!!) euclidean distance between two points */
	double dist(double x[], double y[]);

	/* Translate (t) and rotate (u) a point x[3].
	 * The result is stored into x1[3]. */
	void transform(double t[3], double u[3][3], double x[], double x1[]);

	/* Translate (t) and rotate (u) a set of points x[len][3].
	 * The result is stored into x1[len][3]. */
	void do_rotation(double **x, double **x1, int len, double t[3],
			double u[3][3]);

	/**************************************************************************
	 *                      SETTING INITIAL PARAMETERS                        *
	 **************************************************************************/
	/* Set parameters for searching */
	void set4search(int xlen, int ylen);

	/* Set parameters for final TM-score calculations */
	void set4final(int len);

	/**************************************************************************
	 *                     ASSIGN SECONDARY STRUCTURES                        *
	 **************************************************************************/
	/* Assign secondary structures to a chain of CA atoms x[len][3].
	 * 1->coil, 2->helix, 3->turn, 4->strand */
	void make_sec(double **x, int len, int sec[]);

	/* Define residue SS type based on distances to its neighboring CA atoms */
	int sec_str(double dis13, double dis14, double dis15, double dis24,
			double dis25, double dis35);

	/* Make the SS assignment smooth */
	void smooth(int sec[], int len);

	/**************************************************************************
	 *                        GET INITIAL ALIGNMENTS                          *
	 **************************************************************************/
	/* Perform gapless threading to find the best initial alignment.
	 * Input: x, y, xlen, ylen
	 * Output: y2x stores the best alignment: e.g., y2x[j] = i means:
	 * the jth element in y is aligned to the ith element in x if i>=0
	 * the jth element in y is aligned to a gap in x if i==-1 */
	double get_initial(double **x, double **y, int xlen, int ylen, int y2x[]);

	/* Get initial alignment from secondary structure alignment */
	void get_initial_ss(double **x, double **y, int xlen, int ylen, int y2x[]);

	/* Get initial alignment from local structure superpositions */
	bool get_initial_local(double **x, double **y, int xlen, int ylen,
			int y2x[]);

	/* Fill in the DP score matrix */
	void score_matrix_rmsd_sec(double **x, double **y, int xlen, int ylen,
			int y2x[]);

	/* Get initial alignment from secondary structure and previous alignments */
	void get_initial_ssplus(double **x, double **y, int xlen, int ylen,
			int y2x0[], int y2x[]);

	/* Find the longest continuous (without missing residues) fragment */
	void find_max_frag(double **x, int len, int *start_max, int *end_max);

	/* Perform fragment gapless threading to find the best initial alignment */
	double get_initial_fgt(double **x, double **y, int xlen, int ylen,
			int y2x[]);

	/**************************************************************************
	 *                       TMALIGN CORE FUNCTIONS                           *
	 **************************************************************************/
	/* Calculate TMscore by collecting residues with dis<d */
	int score_fun8(double **x, double **y, int Lali, double d, int y2x[],
			double &score, int score_sum_method);

	/* Compute the score quickly in three iterations */
	double get_score_fast(double **x, double **y, int ylen, int y2x[]);

	/* TMscore search engine
	 * Input:  two aligned vector sets: x[Lali][3], y[Lali][3]
	 *         scale parameter d0
	 *         simplify_step: 1 or 40 or other integers
	 *         score_sum_method: 0 for score over all pairs
	 *                           8 for score over the pairs with dist<score_d8
	 * Output: the best transformation t0, u0 that results in highest TMscore */
	double TMscore8_search(double **xtm, double **ytm, int Lali, double t0[3],
			double u0[3][3], int simplify_step, int score_sum_method);

	/* Comprehensive TMscore search engine
	 * Input:  two vector sets: x[xlen][3], y[ylen][3]
	 *         an alignment y2x[] between x and y
	 *         scale parameter d0
	 *         simplify_step: 1 or 40 or other integers
	 *         score_sum_method: 0 for score over all pairs
	 *                           8 for score over the pairs with dist<score_d8
	 * Output: the best rotation matrix t, u that results in highest TMscore */
	double detailed_search(double **x, double **y, int xlen, int ylen,
			int y2x[], double t[3], double u[3][3], int simplify_step,
			int score_sum_method);

	/**************************************************************************
	 *                         DYNAMIC PROGRAMING                             *
	 **************************************************************************
	 *     Please note the functions below are not a correct implementation   *
	 *     of the N-W dynamic programming because the score tracks back only  *
	 *     one layer of the matrix. This code was exploited in TM-align       *
	 *     because it is about 1.5 times faster than a complete N-W code      *
	 *     and does not influence much the final structure alignment result.  *
	 **************************************************************************/
	/* Initialize DP table (val), scoring matrix (score), and path (path) */
	void NWDP_Initialize(int &len1, int &len2, int j2i[]);

	/* Recreate the path after DP run */
	void NWDP_TraceBack(int &len1, int &len2, double &gap_open, int j2i[]);

	/* DP based on 3D structures */
	void NWDP_TM(double **x, double **y, int len1, int len2, double t[3],
			double u[3][3], double d02, double gap_open, int j2i[]);

	/* DP based on secondary structures */
	void NWDP_TM(int secx[], int secy[], int len1, int len2, double gap_open,
			int j2i[]);

	/* DP based on SS and previous alignments */
	void NWDP_TM(int len1, int len2, double gap_open, int j2i[]);

	/* Heuristic run of DP iteratively to find the best alignment */
	double DP_iter(double **x, double **y, int xlen, int ylen, double t[3],
			double u[3][3], int invmap0[], int g1, int g2, int iteration_max);

	/**************************************************************************
	 *                          MEMORY MANAGEMENT                             *
	 **************************************************************************/
	void Allocate(int xlen, int ylen);
	void Deallocate(int xlen, int ylen);

public:

	/* Structural alignment program for comparing two proteins whose sequences
	 * and size can be different */
	double GetTMscore(double **x, double **y, int xlen, int ylen, int mode = 1);

	/* Structural alignment to compare two models based on their given and
	 * known residue equivalence:
	 * mode = 1 - normalize by xlen (default)
	 *      = 2 - normalize by ylen
	 *      = 3 - normalize by 0.5*(xlen + ylen) */
	double GetTMscore(double **x, double **y, int xlen);

	/* Extract transformation */
	void GetTU(double T[3], double U[3][3]);

	/* Extract alignment */
	void GetAliX2Y(int ali[], int xlen);

	TMalign();
	~TMalign();

};

#endif /* TMALIGN_H_ */
