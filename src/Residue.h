/*****************************************************************************
 *
 * 2015-2016, Ivan Anishchenko <anishchenko.ivan@gmail.com>
 * Vakser Lab, Center for Computational biology, The University of Kansas
 *
 * V20160810
 *
 *****************************************************************************/

#ifndef RESIDUE_H_
#define RESIDUE_H_

#include "Atom.h"
#include "AtomRecord.h"

#include <vector>

class Residue {
private:

	int SetBBAtoms();

	char SetType(const char *name); /* sets 'type' variable */

public:

	/* MEMBER VARIABLES */
	char type; /* residue type: 0, 1, 2, ..., 24; 23 = UNKNOWN (-: unexpectedly :-) */
	char name[4];
	int seqNum;
	char insCode;
	Atom *atom; /* array of atoms */
	int nAtoms;

	/* TODO: test */
	char chainId; /* chain ID from PDB */

	double centroid[3]; /* residue's center of mass (Ca for GLY) */
	void SetCentroid();

	Atom dummy[3]; /* dummy atoms for correspondence with TopResidue */

	bool ntFlag; /* N-terminal */
	bool ctFlag; /* C-terminal */

	Atom *N, *CA, *C, *O; /* pointers to backbone N, CA, C, O */

	/* CONSTRUCTORS */
	Residue(); /* default constructor */
	Residue(const Residue &source); /* copy constructor */

	/* TODO: test these two constructors */
	Residue(const std::vector<AtomRecord> &atoms);
	Residue(std::vector<AtomRecord>::iterator begin_it,
			std::vector<AtomRecord>::iterator end_it);

	/* DESTRUCTOR */
	~Residue();

	/* METHODS */
	Residue& operator=(const Residue &source); /* assignment operator */

	/* TODO: implement saving in PDB */
	void Save(const char *name);

	/* minimal distance between a pair of residues:
	 * mode = 0  -  SC + BB
	 *      = 1  -  SC
	 *      = 2  -  BB */
	static double MinDist(const Residue &R1, const Residue &R2, int mode = 0);

};

#endif /* RESIDUE_H_ */
