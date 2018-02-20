/*****************************************************************************
 *
 * 2013-2016, Ivan Anishchenko <anishchenko.ivan@gmail.com>
 * Vakser Lab, Center for Computational biology, The University of Kansas
 *
 * V20160216
 *
 *****************************************************************************/

#ifndef ATOM_H_
#define ATOM_H_

#include "AtomRecord.h"

class Residue;

class Atom {
private:

	char GetType(char name[]);

public:

	/* MEMBER VARIABLES */
	int atomNum; /* Atom serial number (PDB) */
	char name[5]; /* Atom name (PDB) */
	char altLoc; /* Alternate location indicator (PDB) */

	double x, y, z; /* Orthogonal coordinates for X, Y, Z (PDB) */

	double occup; /* Occupancy (PDB) */
	double temp; /* Temperature factor (PDB) */
	char element[3]; /* Element symbol, right-justified (PDB) */
	char charge[3]; /* Charge on the atom (PDB) */

	Residue *residue; /* pointer to a Residue this atom belongs to */

	double q; /* atomic charge */
	double r; /* atomic radius */
	/* TODO: double p; - some parameter
	 * TODO: or remove r,q,p ??? */

	char type; /* C, N, O, H, S, ... */

	/* CONSTRUCTORS */
	Atom(); /* default constructor */
	Atom(const Atom &source); /* copy constructor */
	Atom(const AtomRecord &source); /* constructs Atom from AtomRecord */

	/* DESTRUCTOR */
	~Atom();

	/* METHODS */
	Atom& operator=(const Atom &source); /* assignment operator */

	static double Dist(const Atom &A, const Atom &B);

};

#endif /* ATOM_H_ */
