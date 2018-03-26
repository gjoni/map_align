/*
 * MSAclass.cpp
 *
 *  Created on: Jul 8, 2016
 *      Author: ivan
 */

#include "MSAclass.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cassert>
#include <cmath>
#include <cfloat>

#include <algorithm>

#include <omp.h>

/* A B C D E F G H I J K L M N O P Q R S T U V W X Y Z */
const unsigned char MSAclass::AMINO_INDICES[26] = { 0, 20, 4, 3, 6, 13, 7, 8, 9,
		20, 11, 10, 12, 2, 20, 14, 5, 1, 15, 16, 20, 19, 17, 20, 18, 20 };

/* 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 *
 * A R N D C Q E G H I  L  K  M  F  P  S  T  W  Y  V  - */
const unsigned char MSAclass::CHAR_INDICES[NAA] = { 65, 82, 78, 68, 67, 81, 69,
		71, 72, 73, 76, 75, 77, 70, 80, 83, 84, 87, 89, 86, 45 };

unsigned char MSAclass::aatoi(unsigned char aa) {

	if (!isalpha(aa)) {
		return 20;
	}

	aa = toupper(aa);
	if (aa < 65 || aa > 90) {
		return 20;
	}

	return MSAclass::AMINO_INDICES[aa - 65];

}

unsigned char MSAclass::itoaa(unsigned char i) {

	if (i < 0 || i > 20) {
		return '-';
	} else {
		return MSAclass::CHAR_INDICES[i];
	}

}

void MSAclass::aatoi(unsigned char *str, size_t len) {

	for (size_t i = 0; i < len; i++) {
		*str = aatoi(*str);
		str++;
	}

}

MSAclass::MSAclass() :
		len_ref(0), nrow(0), ncol(0) {

	/* nothing to be done */

}

MSAclass::MSAclass(const char *name) :
		len_ref(0), nrow(0), ncol(0) {

	FILE *F = fopen(name, "r");
	if (F == NULL) {
		printf("Error: cannot open %s file\n", name);
		exit(1);
	}

	const size_t STRLEN = 65536;

	char buf[STRLEN];

	fgets(buf, STRLEN, F);
	if (buf[0] != '>') {
		printf("Error: MSA '%s' is not in A3M/FASTA format\n", name);
		exit(1);
	}

	TrimRight(buf);
	a3m.push_back(std::make_pair(std::string(buf), std::string("")));

	while (fgets(buf, STRLEN, F)) {

		TrimRight(buf);

		if (buf[0] == '>') /* create a new entry */{
			a3m.push_back(std::make_pair(std::string(buf), std::string("")));
		} else /* update existing */{
			a3m.back().second += std::string(buf);
		}
	}
	fclose(F);

	/* clean lowercase */
	len_ref = a3m[0].second.size(); /* reference sequence length */
	for (size_t i = 1; i < a3m.size(); i++) {
		size_t len = CleanLowercase(a3m[i].second);
		if (len != len_ref) {
			printf("Error: sequence length mismatch at a3m[%ld]\n", i + 1);
			exit(1);
		}
	}

	/* create row & col maps (no gaps will be cleaned here) */
	nrow = a3m.size();
	row_map.resize(nrow);
	std::generate(row_map.begin(), row_map.end(),
			[] {static int i {0}; return i++;});

	ncol = len_ref;
	col_map.resize(ncol);
	std::generate(col_map.begin(), col_map.end(),
			[] {static int i {0}; return i++;});

	a3m_to_msa = col_map;

	/* populate msa */
	Allocate();

	/* set default weights:
	 * all sequences get a weight of 1 */
	weight.assign(nrow, 1.0);

	/* set frequencies (with pseudocounts) */
	SetFreq();

}

MSAclass::MSAclass(const std::string &name) :
		MSAclass(name.c_str()) {

}

MSAclass::MSAclass(const MSAclass &source) :
		a3m(source.a3m), len_ref(source.len_ref), row_map(source.row_map), col_map(
				source.col_map), a3m_to_msa(source.a3m_to_msa), nrow(
				source.nrow), ncol(source.ncol), weight(source.weight), fi(
				source.fi) {

	Allocate();

}

MSAclass & MSAclass::operator =(const MSAclass & source) {

	assert(this != &source); /* an attempt to assign Residue to itself */

	Free();

	a3m = source.a3m;
	len_ref = source.len_ref;

	row_map = source.row_map;
	col_map = source.col_map;

	a3m_to_msa = source.a3m_to_msa;

	nrow = source.nrow;
	ncol = source.ncol;

	weight = source.weight;

	fi = source.fi;

	Allocate();

	return *this;

}

void MSAclass::Allocate() {

	/* */

}

void MSAclass::Free() {

	/* */

}

MSAclass::~MSAclass() {

	Free();

}

void MSAclass::TrimRight(char *str) {

	/*
	 * code from CCMpred
	 */

// find first non-whitespace character from right
	char *end = str + strlen(str) - 1;
	while (end > str && isspace(*end))
		end--;

// add new null terminator
	*(end + 1) = 0;

}

size_t MSAclass::CleanLowercase(std::string &str) {

	size_t pos = 0;

	for (size_t i = 0; i < str.size(); i++) {
		if (!islower(str[i])) {
			str[pos] = str[i];
			pos++;
		}
	}
	str.resize(pos);

	return pos;

}

char MSAclass::GetA3Mres(size_t i, size_t j) const {

	return a3m[i].second[j];

}

size_t MSAclass::GetNrow() const {

	return row_map.size();

}

size_t MSAclass::GetNcol() const {

	return col_map.size();

}

size_t MSAclass::GetMsaIdx(size_t idx) const {

	if (idx < 0 || idx >= len_ref) {
		printf("Error: A3M index (%ld) out of range\n", idx);
		exit(1);
	}

	return a3m_to_msa[idx];

}

size_t MSAclass::GetA3MIdx(size_t idx) const {

	if (idx < 0 || idx >= ncol) {
		printf("Error: MSA index (%ld) out of range\n", idx);
		exit(1);
	}

	return col_map[idx];

}

unsigned char * MSAclass::GetMsa() const {

	unsigned char *msa_ = (unsigned char*) malloc(
			nrow * ncol * sizeof(unsigned char));

	size_t idx = 0;

	for (size_t i = 0; i < row_map.size(); i++) {

		size_t irow = row_map[i];
		const std::string &seq = a3m[irow].second;

		for (size_t j = 0; j < col_map.size(); j++) {
			size_t icol = col_map[j];
			msa_[idx] = seq[icol];
			idx++;
		}
	}

	return msa_;

}

size_t MSAclass::CleanRows(double gaps_frac) {

	row_map.clear();

	for (size_t i = 0; i < a3m.size(); i++) {

		std::string &seq = a3m[i].second;
		size_t ngaps = std::count(seq.begin(), seq.end(), '-');
		if ((double) ngaps / seq.size() < gaps_frac) {
			row_map.push_back(i);
		}

	}

	nrow = row_map.size();

	weight.clear();
	weight.assign(nrow, 1.0);

	SetFreq();

	return nrow;

}

size_t MSAclass::CleanCols(double cols_frac) {

	col_map.clear();
	a3m_to_msa.resize(len_ref);

	/*
	 * calculate numbers of gaps for each column
	 * (only for sequences in row_map[])
	 */
	size_t *counts = (size_t*) malloc(len_ref * sizeof(size_t));
	memset(counts, 0, len_ref * sizeof(size_t));

	for (size_t i = 0; i < row_map.size(); i++) {

		std::string &seq = a3m[row_map[i]].second;
		for (size_t j = 0; j < len_ref; j++) {
			if (seq[j] == '-') {
				counts[j]++;
			}
		}

	}

	/*
	 * populate col_map[] and a3m_to_msa[]
	 */
	size_t count = 0;
	for (size_t i = 0; i < len_ref; i++) {

		if ((double) counts[i] / nrow < cols_frac) {
			col_map.push_back(i);
			a3m_to_msa[i] = count;
			count++;
		} else {
			a3m_to_msa[i] = SIZE_MAX;
		}

	}

	free(counts);

	ncol = col_map.size();

	return ncol;

}

void MSAclass::CleanMsa(double rgaps, double cgaps) {

	CleanRows(rgaps);
	CleanCols(cgaps);

	printf("# Initial MSA: %ld x %ld\n", len_ref, a3m.size());
	printf("# Cleaned MSA: %ld x %ld (%.1f%% x %.1f%% gaps)\n", ncol, nrow,
			rgaps * 100, cgaps * 100);

	const char *seq = a3m.front().second.c_str();
	printf("# %s\n# ", seq);
	for (size_t i = 0; i < len_ref; i++) {
		if (a3m_to_msa[i] < SIZE_MAX) {
			printf("%c", seq[i]);
		} else {
			printf("-");
		}
	}
	printf("\n");

}

void MSAclass::SaveMSA(const std::string &name) const {

	FILE *F = fopen(name.c_str(), "w");

	if (F == NULL) {
		printf("Error: cannot save MSA '%s'\n", name.c_str());
		exit(1);
	}

	for (size_t i = 0; i < row_map.size(); i++) {

		size_t irow = row_map[i];
		const std::string &seq = a3m[irow].second;

		for (size_t j = 0; j < col_map.size(); j++) {
			size_t icol = col_map[j];
			fprintf(F, "%c", seq[icol]);
		}

		fprintf(F, "\n");

	}

	fclose(F);

}

void MSAclass::PrintMSA() const {

	for (size_t i = 0; i < row_map.size(); i++) {

		size_t irow = row_map[i];
		const std::string &seq = a3m[irow].second;

		for (size_t j = 0; j < col_map.size(); j++) {
			size_t icol = col_map[j];
			printf("%c", seq[icol]);
		}

		printf("\n");

	}

}

std::vector<std::pair<size_t, size_t> > MSAclass::CastToMsa(
		const std::vector<std::pair<size_t, size_t> > &contacts) const {

	std::vector<std::pair<size_t, size_t> > contacts_new;

	for (const auto& c : contacts) {
		size_t i = GetMsaIdx(c.first);
		size_t j = GetMsaIdx(c.second);
		if (i < SIZE_MAX && j < SIZE_MAX) {
			contacts_new.push_back(std::make_pair(i, j));
		}
	}

	return contacts_new;

}

size_t MSAclass::GetLen() const {

	return len_ref;

}

size_t MSAclass::GetLen(size_t frag) const {

	size_t len = 0;

	for (size_t i = 0; i < frag; i++) {
		if (a3m_to_msa[i] < SIZE_MAX) {
			len++;
		}
	}

	return len;

}

void MSAclass::Hx(double *hx) const {

	memset(hx, 0, ncol * sizeof(double));

	double **px = (double**) malloc(ncol * sizeof(double**));
	for (size_t i = 0; i < ncol; i++) {
		px[i] = (double*) calloc(NAA, sizeof(double));
	}

	/*
	 * (1) collect aa counts for every sequence position
	 */
	double Beff = 0.0;
	for (size_t i = 0; i < nrow; i++) {
		size_t row = row_map[i];
		double w = weight[i];
		for (size_t j = 0; j < ncol; j++) {
			size_t col = col_map[j];
			unsigned char c = aatoi(a3m[row].second[col]);
			if (c >= 0 && c < NAA) {
				px[j][c] += w;
			}
		}
		Beff += w;
	}

	/*
	 * (2) convert counts into entropies
	 */
	for (size_t i = 0; i < ncol; i++) {
		for (size_t c = 0; c < NAA; c++) {
			if (px[i][c] > 0.0) {
				px[i][c] /= Beff;
				hx[i] -= px[i][c] * log(px[i][c]);
			}
		}
	}

	/*
	 * free
	 */
	for (size_t i = 0; i < ncol; i++) {
		free(px[i]);
	}
	free(px);

}

void MSAclass::HxHy(double **hxhy) const {

	for (size_t i = 0; i < ncol; i++) {
		memset(hxhy[i], 0, ncol * sizeof(double));
	}

	double *hx = (double*) malloc(ncol * sizeof(double));

	Hx(hx);

	for (size_t i = 0; i < ncol; i++) {
		for (size_t j = i + 1; j < ncol; j++) {
			hxhy[i][j] = hx[i] + hx[j];
			hxhy[j][i] = hxhy[i][j];
		}
	}

	free(hx);

}

void MSAclass::GxGy(double **gxgy) {

	for (size_t i = 0; i < ncol; i++) {
		memset(gxgy[i], 0, ncol * sizeof(double));
	}

	double *gx = (double*) calloc(ncol, sizeof(double));

	/*
	 * (1) collect gap counts for every sequence position
	 */
	double Beff = 0.0;
	for (size_t i = 0; i < nrow; i++) {
		size_t row = row_map[i];
		double w = weight[i];
		for (size_t j = 0; j < ncol; j++) {
			size_t col = col_map[j];
			char c = a3m[row].second[col];
			if (c == '-') {
				gx[j] += w;
			}
		}
		Beff += w;
	}

	/*
	 * (2) convert counts into frequencies
	 */
	for (size_t j = 0; j < ncol; j++) {
		gx[j] /= Beff;
	}

	/*
	 * (3) 
	 */
	for (size_t i = 0; i < ncol; i++) {
		double ei = gx[i] > 0 ? log(gx[i]) : 0.0;
		for (size_t j = i + 1; j < ncol; j++) {
			double ej = gx[j] > 0 ? log(gx[j]) : 0.0;
			gxgy[i][j] = ei + ej;
			gxgy[j][i] = gxgy[i][j];
		}
	}

}

void MSAclass::Gxy(double **gxy) {

	for (size_t i = 0; i < ncol; i++) {
		memset(gxy[i], 0, ncol * sizeof(double));
	}

	double Beff = 0.0;
	for (size_t i = 0; i < nrow; i++) {

		size_t row = row_map[i];
		double w = weight[i];

		for (size_t p = 0; p < ncol; p++) {

			char cp = a3m[row].second[col_map[p]];
			for (size_t q = p + 1; q < ncol; q++) {
				char cq = a3m[row].second[col_map[q]];
				if (cp == '-' && cq == '-') {
					gxy[p][q] += w;
				}
			}
		}
		Beff += w;
	}

	for (size_t p = 0; p < ncol; p++) {
		for (size_t q = p + 1; q < ncol; q++) {
			gxy[p][q] = gxy[p][q] > 0.0 ? log(gxy[p][q] / Beff) : 0.0;
			gxy[q][p] = gxy[p][q];
		}
	}

}

void MSAclass::Hxy(double **hxy) const {

	for (size_t i = 0; i < ncol; i++) {
		memset(hxy[i], 0, ncol * sizeof(double));
	}

	/* loop over all i<j pairs */
	for (size_t i = 0; i < ncol; i++) {

		for (size_t j = i + 1; j < ncol; j++) {

			/* store aa pair counts here */
			double *pxy = (double*) calloc(NAA * NAA, sizeof(double));

			/* loop over all sequences */
			double Beff = 0.0;
			for (size_t k = 0; k < nrow; k++) {
				double w = weight[k];
				const std::string &seq = a3m[row_map[k]].second;
				unsigned char a = aatoi(seq[col_map[i]]);
				unsigned char b = aatoi(seq[col_map[j]]);
				if (a >= 0 && a < NAA && b >= 0 && b < NAA) {
					pxy[a * NAA + b] += w;
				}
				Beff += w;
			}

			/* convert counts into entropies */
			for (size_t ab = 0; ab < NAA * NAA; ab++) {
				if (pxy[ab] > 0.0) {
					pxy[ab] /= Beff;
					hxy[i][j] -= pxy[ab] * log(pxy[ab]);
				}
			}

			/* make the hxy[][] matrix symmetric */
			hxy[j][i] = hxy[i][j];

			free(pxy);

		}
	}

}

void MSAclass::MI(double **mi) const {

	for (size_t i = 0; i < ncol; i++) {
		memset(mi[i], 0, ncol * sizeof(double));
	}

	double *hx = (double*) malloc(ncol * sizeof(double));
	double **hxy = (double**) malloc(ncol * sizeof(double*));
	for (size_t i = 0; i < ncol; i++) {
		hxy[i] = (double*) malloc(ncol * sizeof(double));
	}

	Hx(hx);
	Hxy(hxy);

	for (size_t i = 0; i < ncol; i++) {
		for (size_t j = 0; j < ncol; j++) {
			mi[i][j] = hx[i] + hx[j] - hxy[i][j];
//				printf("%ld %ld %.4f %.4f %.4f %.4f\n", i + 1, j + 1, hx[i], hx[j],
//						hxy[i][j], mi[i][j]);
		}
	}

	/*
	 * free
	 */
	for (size_t i = 0; i < ncol; i++) {
		free(hxy[i]);
	}
	free(hxy);
	free(hx);

}

double MSAclass::Reweight(double t) {

	assert(t > 0.0 && t < 1.0); /* out of range */

	size_t idthres = (size_t) ceil(t * ncol);

	weight.assign(nrow, 0.0);

	size_t nij = nrow * (nrow + 1) / 2;

	unsigned char *msa = GetMsa();

#if defined(_OPENMP)
#pragma omp parallel for
#endif
	for (size_t ij = 0; ij < nij; ij++) {

		// compute i and j from ij
		// http://stackoverflow.com/a/244550/1181102
		size_t i, j;
		{
			size_t ii = nrow * (nrow + 1) / 2 - 1 - ij;
			size_t K = floor((sqrt(8 * ii + 1) - 1) / 2);
			i = nrow - 1 - K;
			j = ij - nrow * i + i * (i + 1) / 2;
		}

		size_t ids = 0;

		unsigned char *seqi = msa + i * ncol;
		unsigned char *seqj = msa + j * ncol;

		for (size_t k = 0; k < ncol; k++) {
			ids += (*seqi++ == *seqj++);
		}

		if (ids > idthres) {
			weight[i]++;
			weight[j]++;
		}
	}

	free(msa);

	double wsum = 0, wmin = DBL_MAX, wmax = DBL_MIN;
	for (size_t i = 0; i < nrow; i++) {
		weight[i] = 1. / (weight[i] - 1);
		wsum += weight[i];
		wmin = weight[i] < wmin ? weight[i] : wmin;
		wmax = weight[i] > wmax ? weight[i] : wmax;
	}

	printf("# threshold= %.1f Beff= %g mean= %g min= %g max= %g\n", t, wsum,
			wsum / nrow, wmin, wmax);

	SetFreq();

	return wsum;

}

double MSAclass::GetWeight(size_t i) const {

	assert(i >= 0 && i < nrow); /* out of range */
	return weight[i];

}

double MSAclass::GetNeff() const {

	double wsum = 0;
	for (size_t i = 0; i < nrow; i++) {
		wsum += weight[i];
	}

	return wsum;

}

const std::string& MSAclass::GetSequence(size_t i) const {

	assert(i >= 0 && i < a3m.size()); /* out of range */

	return a3m[i].second;

}

void MSAclass::SetFreq() {

	fi.clear();

	double Beff = GetNeff();

	/* pseudocount */
	double lpseudo = 1.0 * Beff;

	/* resize fi[ncol][NAA] and
	 * initialize with pseudocounts */
	fi.resize(ncol, std::vector<double>(NAA, lpseudo / NAA));

	/* fill in fi[ncol][NAA] woth counts */
	for (size_t i = 0; i < nrow; i++) {
		const char *seqi = a3m[row_map[i]].second.c_str();
		double w = weight[i];
		for (size_t j = 0; j < ncol; j++) {
			unsigned char a = aatoi(seqi[col_map[j]]);
			fi[j][a] += w;
		}
	}

	/* convert counts into frequencies */
	for (size_t j = 0; j < ncol; j++) {
		for (size_t a = 0; a < NAA; a++) {
			fi[j][a] /= (Beff + lpseudo);
		}
	}

}

double MSAclass::GetFi(size_t i, unsigned char a) const {

	assert(i >= 0 && i < ncol); /* residue index out of range */
	assert(a >= 0 && a < NAA); /* aa identity out of range */

	return fi[i][a];

}
