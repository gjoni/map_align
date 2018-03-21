/*****************************************************************************
 *
 * 2016, Ivan Anishchenko <anishchenko.ivan@gmail.com>
 *
 * Vakser Lab,
 * Center for Computational biology,
 * University of Kansas
 *
 * V20170626 - remove member functions for energy calculations
 * V20170601 - support for copying h,J to external arrays
 * V20160527
 *
 *****************************************************************************/

#ifndef RRCE_H_
#define RRCE_H_

#include <map>
#include <string>
#include <vector>

class RRCE {
private:

	double *h;
	double **J;

	double dmax;

	std::string name;

	int ReadCouplings(const char *name); /* sets h,J */

	void Allocate(int n);
	void Free(int n);

public:

	enum RRCE_TYPE {
		RRCE20RC = 1, /* residue centroids */
		RRCE20SCC = 2, /* side-chain centroids */
		RRCE20CB = 3 /* C-beta atoms */
	};

	RRCE(RRCE_TYPE type, double d, int k);
	RRCE(const RRCE &source);
	RRCE(const char *fname); /* TODO: test this constructor */
	RRCE();

	~RRCE();

	RRCE& operator=(const RRCE &source);

	std::string GetName();

	void GetCouplings(double **J);
	void GetFields(double *h);

	double GetJij(size_t i, size_t j) const {
		return J[i][j];
	}
	;

	double GetHi(size_t i) const {
		return h[i];
	}
	;

};

#endif /* RRCE_H_ */
