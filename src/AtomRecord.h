
#ifndef ATOM_RECORD_H_
#define ATOM_RECORD_H_

/*
 * A structure to fully describe the ATOM record in PDB
 */

/* TODO: add a hetero flag */

struct AtomRecord
{

	int atomNum; /* Atom serial number */
	char atomName[5]; /* Atom name */
	char altLoc; /* Alternate location indicator */
	char resName[4]; /* Residue name */
	char chainId; /* Chain identifier */
	int resNum; /* Residue sequence number */
	char insCode; /* Code for insertion of residues */
	double x, y, z; /* Orthogonal coordinates for X, Y, Z */
	double occup; /* Occupancy */
	double temp; /* Temperature factor */
	char segId[5]; /* Segment identifier, left-justified */
	char element[3]; /* Element symbol, right-justified */
	char charge[3]; /* Charge on the atom */

	/* for handling PQR files */
	double q; /* atomic charge */
	double r; /* atomic radius */

};

#endif /* ATOM_RECORD_H_ */
