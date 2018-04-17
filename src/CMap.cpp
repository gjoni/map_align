/*
 * CMap.cpp
 *
 *  Created on: Jan 26, 2018
 *      Author: ivan
 */

#include <cassert>

#include "CMap.h"

CMap::CMap() :
		seq(""), size(0), size_ali(0), has_cont(0), left(0), right(0), mleft(0), mright(
				0) {

	/* */

}

CMap::~CMap() {
	// TODO Auto-generated destructor stub
}

CMap & CMap::operator=(const CMap &source) {

	seq = source.seq;
	size = source.size;
	size_ali = source.size_ali;
	has_cont = source.has_cont;
	left = source.left;
	right = source.right;
	mleft = source.mleft;
	mright = source.mright;
	edges = source.edges;

	return *this;

}

CMap::CMap(const CMap &source) :
		seq(source.seq), size(source.size), size_ali(source.size_ali), has_cont(
				source.has_cont), left(source.left), right(source.right), mleft(
				source.mleft), mright(source.mright), edges(source.edges) {

	/* */

}

CMap::CMap(const std::string& name, const std::string& sequence) :
		seq(sequence), size(seq.length()), size_ali(0), has_cont(size, 0), left(
				size), right(size) {

	FILE *F = fopen(name.c_str(), "r");
	if (F == NULL) {
		printf("Error: cannot open contacts file '%s'\n", name.c_str());
		exit(1);
	}

	/* temp. adjacency matrix */
	double **mtx = (double**) malloc(size * sizeof(double*));
	for (unsigned i = 0; i < size; i++) {
		mtx[i] = (double*) calloc(size, sizeof(double));
	}

	/* read the file (assuming it is in CASP format)
	 * i  j  d1  d2  p */
	const int SIZE = 1024;
	char buf[SIZE];
	int a, b;
	double p;
	while (fgets(buf, SIZE, F)) {
		if (buf[0] == '#') {
			continue;
		}

		if (sscanf(buf, "%d %d %*s %*s %lf\n", &a, &b, &p) == 3) {

			assert(a > 0 && a <= (int )size); /* index out of range */
			assert(b > 0 && b <= (int )size); /* index out of range */
			assert(p >= 0.0 && p <= 1.0); /* probability out of range */

			/* assume numbering in the input file
			 * starts with 1 -- shift to zero */
			a--;
			b--;

			mtx[a][b] = p;
			mtx[b][a] = p;

			/* mark residues with contacts */
			has_cont[a] = true;
			has_cont[b] = true;

		}
	}

	/* fill in adjacency lists & edge list */
	for (unsigned i = 0; i < size; i++) {

		/* left */
		for (unsigned j = 0; j < i; j++) {
			if (mtx[i][j] > 1.0e-6) {
				left[i].push_back(std::make_tuple(j, mtx[i][j], i - j));
			}
		}

		/* right */
		for (unsigned j = i + 1; j < size; j++) {
			if (mtx[i][j] > 1.0e-6) {
				right[i].push_back(std::make_tuple(j, mtx[i][j], j - i));
				edges[ {i,j}] = {mtx[i][j], j - i};
			}
		}

		/* fill in the maps */
		if (left[i].size()) {
			mleft.push_back(i);
		}

		if (right[i].size()) {
			mright.push_back(i);
		}

		/* count residues with contacts */
		size_ali += has_cont[i];

	}

	/* free */
	for (unsigned i = 0; i < size; i++) {
		free(mtx[i]);
	}
	free(mtx);

}

CMap::CMap(const AListT& adj, const std::string& sequence) :
		seq(sequence), size(sequence.length()), size_ali(0), has_cont(size, 0), left(
				size), right(size) {

	for (unsigned i = 0; i < size; i++) {

		for (auto &n : adj[i]) {
			unsigned j = std::get<0>(n);
			if (j < i) {
				left[i].push_back(n);
			} else if (j > i) {
				right[i].push_back(n);
				edges[ {i, j}] = {std::get<1>(n), std::get<2>(n)};
			}
			has_cont[i] = true;
			has_cont[j] = true;
		}

		if (left[i].size()) {
			mleft.push_back(i);
		}

		if (right[i].size()) {
			mright.push_back(i);
		}

	}

	for (unsigned i = 0; i < size; i++) {
		size_ali += has_cont[i];
	}

}

unsigned CMap::Size() const {

	return size;

}

const std::vector<unsigned>& CMap::GetLeftMap() const {

	return mleft;

}

const std::vector<unsigned>& CMap::GetRightMap() const {

	return mright;

}

const NListT& CMap::GetLeftList(unsigned i) const {

	return left[i];

}

const NListT& CMap::GetRightList(unsigned i) const {

	return right[i];

}

void CMap::Print() const {

	printf("LEN %lu\n", seq.length());
	for (auto &e : edges) {
		printf("CON %u %u %f\n", e.first.first, e.first.second, e.second.first);
	}

	/*
	 printf("%s\n", seq.c_str());

	 for (unsigned i = 0; i < size; i++) {
	 printf("%d:", i);
	 for (auto &c : left[i]) {
	 printf(" L(%u,%.3f,%u)", std::get<0>(c), std::get<1>(c),
	 std::get<2>(c));
	 }
	 for (auto &c : right[i]) {
	 printf(" R(%u,%.3f,%u)", std::get<0>(c), std::get<1>(c),
	 std::get<2>(c));
	 }
	 printf("\n");
	 }
	 */

}

const std::vector<bool>& CMap::GetContFl() const {

	return has_cont;

}
