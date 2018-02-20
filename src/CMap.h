/*
 * CMap.h
 *
 *  Created on: Jan 26, 2018
 *      Author: ivan
 */

#ifndef CMAP_H_
#define CMAP_H_

#include <vector>
#include <utility>
#include <tuple>
#include <string>
#include <map>

/* neighbors list type: element --> (neighbor, contact_score, separation) */
typedef std::vector<std::tuple<unsigned, double, unsigned> > NListT;

/* adjacency list type */
typedef std::vector<NListT> AListT;

/* edge list type: (i,j) --> (contact_score, separation) */
typedef std::map<std::pair<unsigned, unsigned>, std::pair<double, unsigned> > EListT;

/*
 * contact map class
 */
class CMap {

	friend class MapAlign;

private:

	/* sequence */
	std::string seq;

	/* dimension (aka sequence length) */
	unsigned size;

	/* neighbors to the left and to the right
	 * of the diagonal */
	AListT left, right;

	/* mappings to positions with nonempty
	 * left/right neighbor lists */
	std::vector<unsigned> mleft, mright;

	/* edge list */
	EListT edges;

public:

	/* from file & sequence */
	CMap(const std::string&, const std::string&);

	/* from adjacency list & sequence */
	CMap(const AListT&, const std::string&);

	CMap(const CMap &source);

	CMap();

	~CMap();

	CMap & operator=(const CMap &source);

	unsigned Size() const;

	void Print() const;

	const NListT& GetLeftList(unsigned i) const;
	const NListT& GetRightList(unsigned i) const;

	const std::vector<unsigned>& GetLeftMap() const;
	const std::vector<unsigned>& GetRightMap() const;

};

#endif /* CMAP_H_ */
