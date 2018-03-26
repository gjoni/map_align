/*
 * MSAclass.h
 *
 *  Created on: Jul 8, 2016
 *      Author: ivan
 */

#ifndef MSACLASS_H_
#define MSACLASS_H_

#include <string>
#include <vector>

class MSAclass {

	friend class ProblemBase;
	friend class ProblemFull;
	friend class ProblemFullOMP;
	friend class ProblemFullAsym;
	friend class ProblemPNAS;
	friend class MRFclass;

private:

	/* raw MSA in A3M/FASTA */
	std::vector<std::pair<std::string, std::string> > a3m;

	/* reference sequence length (1st sequence in the MSA) */
	size_t len_ref;

	size_t CleanLowercase(std::string &str); /* in-place removal of lowercase letters */

	/* cleaned rows and columns */
	std::vector<size_t> row_map; /* cleaned to reference */
	std::vector<size_t> col_map; /* cleaned to reference */
	std::vector<size_t> a3m_to_msa;

	size_t CleanRows(double gaps_frac); /* sets row_map */
	size_t CleanCols(double gaps_frac); /* sets col_map, ref_to_cleaned */

	/* size of the CLEANED alignment msa[] */
	size_t nrow; // sequence length - 65536 MAX
	size_t ncol; // number of sequences in MSA - 65536 MAX

	/* sequence weights (CLEANED alignment) */
	std::vector<double> weight;

	/* 1-site frequencies (CLEANED alignment) */
	std::vector<std::vector<double> > fi; /* ncol * NAA */
	void SetFreq();

	/* alloc/free msa*/
	void Allocate();
	void Free();

	void TrimRight(char *str);

public:

	/*
	 * conversions between AA letters and indices
	 */
	static const unsigned char NAA = 21;
	static const unsigned char AMINO_INDICES[26];
	static const unsigned char CHAR_INDICES[NAA];
	static unsigned char aatoi(unsigned char aa);
	static unsigned char itoaa(unsigned char aa);
	static void aatoi(unsigned char *str, size_t len);

	MSAclass(const char *name);
	MSAclass(const std::string&);
	MSAclass(const MSAclass &source);
	MSAclass();

	~MSAclass();

	MSAclass& operator=(const MSAclass &source);

	/* cleaned MSA dimensions */
	size_t GetNrow() const;
	size_t GetNcol() const;

	/* cleaned length of frag */
	size_t GetLen(size_t frag) const;

	/* reference length */
	size_t GetLen() const;

	/* sequence reweighting */
	double Reweight(double t = 0.8);
	double GetWeight(size_t i) const;
	double GetNeff() const; /* sum over weights */

	/* A3M->MSA mapping*/
	size_t GetMsaIdx(size_t idx) const;

	/* MSA->A3M mapping*/
	size_t GetA3MIdx(size_t idx) const;

	/* change residue indices in the array of contacts
	 * to match the cleaned msa[]
	 * (contacts with residues missing in the msa[] are omitted) */
	std::vector<std::pair<size_t, size_t> > CastToMsa(
			const std::vector<std::pair<size_t, size_t> > &contacts) const;

	/* for a complex */
//	std::vector<std::pair<size_t, size_t> > CastToMsa(
//			const std::vector<std::pair<int, int> > &contacts,
//			size_t len_rec) const;
	char GetA3Mres(size_t i, size_t j) const;
//	char GetMSAres(size_t i, size_t j) const;

	/* rows are cleaned first, then columns */
	void CleanMsa(double rgaps, double cgaps);

	/* 1d array [nrow x ncol] with the cleaned alignment
	 * !!! allocates memory with malloc !!! */
	unsigned char * GetMsa() const;

	const std::string& GetSequence(size_t i) const;

	void SaveMSA(const std::string &name) const;
	void PrintMSA() const;

	void Hx(double *hx) const;
	void Hxy(double **hxy) const;
	void MI(double **mi) const;
	void HxHy(double **hxhy) const; // H(x) + H(y)

	void GxGy(double **gxgy); // Gap(x) + Gap(y)
	void Gxy(double **gxy); // Gap(x && y);

	double GetFi(size_t i, unsigned char a) const;

};

#endif /* MSACLASS_H_ */
