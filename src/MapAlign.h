/*
 * MapAlign.h
 *
 *  Created on: Jan 26, 2018
 *      Author: ivan
 */

#ifndef MAPALIGN_H_
#define MAPALIGN_H_

#include <vector>

#include "CMap.h"
#include "Chain.h"

using namespace std;

struct MP_RESULT {
	std::string label;
	std::vector<double> sco;
	std::vector<int> len;
	std::vector<int> a2b;
};

class MapAlign {
private:

	MapAlign();
	~MapAlign();

	struct SWDATA {
		const CMap &A;
		const CMap &B;
		const Chain &PA;
		const Chain &PB;
		double **mtx; /* contacts scoring matrix (M x N) */
		double **sco; /* DP scoring matrix (M + 1) x (N + 1) */
		char **label; /* path in the DP matrix (M + 1) x (N + 1) */
		std::vector<double> gap_a; /* gap opening penalties for A */
		std::vector<double> gap_b; /* gap opening penalties for B */
		std::vector<int> a2b;
		std::vector<int> b2a;
		double tot_scoA;
		double tot_scoB;
	};

	static void Alloc(SWDATA*);
	static void Free(SWDATA*);

	static double gaussian(double mean, double std, double x);
	static double sepw(double sep);

	/* step1 Smith-Waterman alignment to calculate score only */
	static double SW1(const NListT&, const NListT&, double, double);

	/* step2 Smith-Waterman alignment */
	static double SW2(SWDATA&, double gap_e);

	/* TODO: current implementation is not very efficient */
	static double Intersect(const NListT&, const NListT&,
			const std::vector<int>&, const std::vector<int>&);

	/* a function to assess current alignment
	 * based on contact/gap scores - returned as a vector */
	static MP_RESULT Assess(const SWDATA&, double);

	static void InitMTX(SWDATA&, double sep_x, double sep_y);

	static void UpdateMTX(SWDATA&, double, int iter);

	/* TMscore of the map_aligned region */
	static double TMscore(const SWDATA&);

	/* MPscore of the tmaligned region */
	static double MPscore(SWDATA&);
	static double TMscore(const Chain&, const Chain&, const std::vector<int>&);

public:

	struct PARAMS {
		double gap_open;
		double gap_ext;
		int kmin;
		int iter;
	};

	static double MaxScore(const CMap&);

	static MP_RESULT Align(const CMap&, const CMap&, const Chain&, const Chain&,
			const PARAMS&);

};

#endif /* MAPALIGN_H_ */
