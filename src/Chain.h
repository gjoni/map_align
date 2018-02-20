/*
 * Chain.h
 *
 *  Created on: Jun 3, 2014
 *      Author: ivan
 */

#ifndef CHAIN_H_
#define CHAIN_H_

#include <string>
#include <map>
#include <utility>

#include "AtomRecord.h"
#include "Residue.h"
#include "kdtree.h"

class Chain {

private:

	/* create AtomRecord from C-style string */
	AtomRecord ReadAtomRecord(char *str); /* 'str' will be modified !!! */
	AtomRecord ReadAtomRecordPQR(char *str); /* PQR support for APBS */

	/* initialize Chain from an array of AtomRecords */
	void InitChain(std::vector<AtomRecord> &atomRecordVector);

	/* initialize Atom **atoms */
	void SetAtomPointers();

	/* A function to detect N- and C- termini based on
	 * CA-CA distance between 2 sequential residues */
	void SetTermini();

	/* a map to quickly access residues by their sequence number and
	 * insertion code in PDB */
	std::map<std::pair<int /* number */, char /* insCode */>, Residue*> resmap;
	void SetResMap();

public:

	/* MEMBER VARIABLES */
	Atom **atom; /* array of pointers to all atoms in the chain */
	int nAtoms;
	Residue *residue; /* array of residues */
	int nRes;

	kdtree *kd; /* kd-tree for all atoms */
	kdtree *kdCA; /* kd-tree for C-alphas */
	kdtree *kdCent; /* kd-tree for residues' centroids */
	void SetKD(); /* set all 3 kd-trees */

	double **ca_trace; /* array of CAs' Cartesian coordinates (for TMalign) */

	/* CONSTRUCTORS */
	Chain(); /* default constructor */
	Chain(const char *name); /* constructs Chain by reading a PDB-file */
	Chain(const std::string &name); /* constructs Chain by reading a PDB-file */
	Chain(const Chain &source); /* copy constructor */
	Chain(const std::vector<std::string> &atoms_str);

	/* DESTRUCTOR */
	~Chain();

	/* METHODS */
	Chain& operator=(const Chain &source); /* assignment operator */
	void Transform(double t[3], double u[3][3]); /* r' = ur + t */
	void Load(const char *name);
	void Save(const char *name);

	void Renumber(const std::vector<int> &num);

	/* TODO: get residue by number in PDB */
	Residue* GetResidue(int n, char ins = ' ');

	std::string GetSequence();

	int IfProtein();
	/* TODO:
	 * 1) void SetChainID();
	 * 2) Chain(const std::vector<std::string> &atoms);
	 *      or
	 *    Chain(std::vector<AtomRecord>::iterator begin_it,
	 *          std::vector<AtomRecord>::iterator end_it);
	 *  3) int IfProtein();  number of standard AA */

};

#endif /* CHAIN_H_ */
