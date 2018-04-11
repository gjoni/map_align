/*
 * Residue.cpp
 *
 *  Created on: Jun 3, 2014
 *      Author: ivan
 */

#include "Residue.h"
#include "Info.h"

#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <new>
#include <map>

Residue::Residue() :
		type(23), seqNum(-1), insCode(' '), atom(NULL), nAtoms(0), chainId(' '), ntFlag(
				false), ctFlag(false), N(NULL), CA(NULL), C(NULL), O(NULL), CB(
		NULL) {

	strcpy(name, "UNK");

	centroid[0] = centroid[1] = centroid[2] = 0;

}

Residue::Residue(const Residue &source) :
		type(source.type), seqNum(source.seqNum), insCode(source.insCode), atom(
		NULL), nAtoms(source.nAtoms), chainId(source.chainId), ntFlag(
				source.ntFlag), ctFlag(source.ctFlag), N(NULL), CA(NULL), O(
		NULL), CB(NULL) {

	assert(nAtoms > 0); /* pedantic check */

	strcpy(name, source.name);

	atom = new Atom[nAtoms];
	for (int i = 0; i < nAtoms; i++) {
		atom[i] = source.atom[i];
		atom[i].residue = this;
	}

	SetBBAtoms(); /* N, CA, C, O, CB are pointers - they need  to be reinitialized */

	memcpy(centroid, source.centroid, 3 * sizeof(double));

}

Residue::Residue(const std::vector<AtomRecord>& source) :
		seqNum(source[0].resNum), insCode(source[0].insCode), atom(NULL), nAtoms(
				source.size()), chainId(source[0].chainId), ntFlag(false), ctFlag(
				false), N(NULL), CA(NULL), C(
		NULL), O(NULL), CB(NULL) {

	assert(nAtoms > 0); /* pedantic check */

	strcpy(name, source[0].resName); /* first atom is used to name the residue*/

	atom = new Atom[nAtoms];
	int idx = 0;
	for (int i = 0; i < nAtoms; i++) {
		atom[i] = Atom(source[i]);
		atom[i].residue = this;
		idx++;
	}

	SetBBAtoms(); /* N, CA, C, O, CB are pointers - they need to be reinitialized */

	type = SetType(name);

	SetCentroid();

}

Residue::Residue(std::vector<AtomRecord>::iterator begin_it,
		std::vector<AtomRecord>::iterator end_it) :
		seqNum(begin_it->resNum), insCode(begin_it->insCode), atom(NULL), nAtoms(
				end_it - begin_it + 1), chainId(begin_it->chainId), ntFlag(
				false), ctFlag(false), N(NULL), CA(
		NULL), C(NULL), O(NULL), CB(NULL) {

	assert(nAtoms > 0); /* pedantic check */

	strcpy(name, begin_it->resName); /* residue name is read from the 1st atom */

	/* check for alternative conformations */
	std::vector<AtomRecord>::iterator it;
	bool fl = false;
	for (it = begin_it; it <= end_it; it++) {
		if (it->altLoc != ' ') {
			fl = true;
			break;
		}
	}

	if (!fl) { /* no alternative conformations */

		/* use all atoms */
		atom = new Atom[nAtoms];
		int idx = 0;
		for (it = begin_it; it <= end_it; it++) {
			atom[idx] = Atom(*it);
			atom[idx].residue = this;
			idx++;
		}

	} else { /* there are alternative conformations */

		/* save top-populated conformation for each atom */
		std::map<std::string, std::vector<AtomRecord>::iterator> atoms_map;
		std::map<std::string, std::vector<AtomRecord>::iterator>::iterator itmap;

		atoms_map[begin_it->atomName] = begin_it;
		for (it = begin_it + 1; it <= end_it; it++) {

			/* save only atoms with the same residue name as the first atom */
			if (strcmp(begin_it->resName, it->resName) == 0) {
				itmap = atoms_map.find(it->atomName);
				if (itmap == atoms_map.end()) {
					atoms_map[it->atomName] = it;
				} else {
					if (itmap->second->occup < it->occup) {
						itmap->second = it;
					}
				}
			}
		}

		/* save marked atoms */
		atom = new Atom[atoms_map.size()];
		nAtoms = 0;
		for (itmap = atoms_map.begin(); itmap != atoms_map.end(); itmap++) {
			atom[nAtoms] = Atom(*(itmap->second));
			atom[nAtoms].residue = this;
			nAtoms++;
		}

	}

	SetBBAtoms(); /* N, CA, C, O, CB are pointers - they need to be reinitialized */

	type = SetType(name);

	SetCentroid();

}

Residue::~Residue() {

	delete[] atom;
	atom = NULL;

}

Residue & Residue::operator =(const Residue & source) {

	assert(this != &source); /* an attempt to assign Residue to itself */

	type = source.type;
	strcpy(name, source.name);
	seqNum = source.seqNum;
	insCode = source.insCode;
	nAtoms = source.nAtoms;
	chainId = source.chainId;

	delete[] atom;

	atom = new Atom[nAtoms];
	for (int i = 0; i < nAtoms; i++) {
		atom[i] = source.atom[i];
		atom[i].residue = this;
	}

	SetBBAtoms();

	memcpy(centroid, source.centroid, 3 * sizeof(double));

	return *this;

}

int Residue::SetBBAtoms() {

	bool flN = false, flCA = false, flC = false, flO = false, flCB = false;

	for (int i = 0; i < nAtoms; i++) {
		if (strcmp(atom[i].name, "N") == 0) {
			N = &(atom[i]);
			flN = true;
		} else if (strcmp(atom[i].name, "CA") == 0) {
			CA = &(atom[i]);
			flCA = true;
		} else if (strcmp(atom[i].name, "CB") == 0) {
			CB = &(atom[i]);
			flCB = true;
		} else if (strcmp(atom[i].name, "C") == 0) {
			C = &(atom[i]);
			flC = true;
		} else if (strcmp(atom[i].name, "O") == 0
				|| strcmp(atom[i].name, "OXT") == 0) {
			O = &(atom[i]);
			flO = true;
		} else if (strcmp(atom[i].name, "H1") == 0) {
			ntFlag = true;
		} else if (strcmp(atom[i].name, "H2") == 0) {
			ntFlag = true;
		} else if (strcmp(atom[i].name, "H3") == 0) {
			ntFlag = true;
		} else if (strcmp(atom[i].name, "OXT") == 0) {
			ctFlag = true;
		}
	}

	int n = flN + flCA + flC + flO + flCB;

	return n;

}

char Residue::SetType(const char* name) {

	for (int i = 0; i < 24; i++) {
		if (strcmp(name, AAA3[i]) == 0) {
			return i;
		}
	}
	return 23;
}

void Residue::SetCentroid() {

	centroid[0] = centroid[1] = centroid[2] = 0.0;

	if (type == 7 /* GLY */) {
		centroid[0] = CA->x;
		centroid[1] = CA->y;
		centroid[2] = CA->z;
	} else /* other residues */{
		for (int j = 0; j < nAtoms; j++) {
			centroid[0] += atom[j].x;
			centroid[1] += atom[j].y;
			centroid[2] += atom[j].z;
		}
		centroid[0] /= nAtoms;
		centroid[1] /= nAtoms;
		centroid[2] /= nAtoms;
	}

}

double Residue::MinDist(const Residue &R1, const Residue &R2, int mode) {

	double dmin = 9999.99;

	switch (mode) {

	case 0: /* SC + BB */
		for (int i = 0; i < R1.nAtoms; i++) {
			if (R1.atom[i].type == 'H') {
				continue;
			}
			for (int j = 0; j < R2.nAtoms; j++) {
				if (R2.atom[j].type == 'H') {
					continue;
				}
				double d = Atom::Dist(R1.atom[i], R2.atom[j]);
				dmin = d < dmin ? d : dmin;
			}
		}
		break;

	case 1: /* SC */
		for (int i = 0; i < R1.nAtoms; i++) {
			Atom *A = &(R1.atom[i]);
			if (A == R1.C || A == R1.CA || A == R1.N || A == R1.O) {
				continue;
			}
			if (R1.atom[i].type == 'H') {
				continue;
			}
			for (int j = 0; j < R2.nAtoms; j++) {
				Atom *B = &(R2.atom[j]);
				if (B == R2.C || B == R2.CA || B == R2.N || B == R2.O) {
					continue;
				}
				if (R2.atom[j].type == 'H') {
					continue;
				}
				double d = Atom::Dist(R1.atom[i], R2.atom[j]);
				dmin = d < dmin ? d : dmin;
			}
		}
		break;

	case 2: /* BB */
		for (int i = 0; i < R1.nAtoms; i++) {
			Atom *A = &(R1.atom[i]);
			if (A != R1.C && A != R1.CA && A != R1.N && A != R1.O) {
				continue;
			}
			if (R1.atom[i].type == 'H') {
				continue;
			}
			for (int j = 0; j < R2.nAtoms; j++) {
				Atom *B = &(R2.atom[j]);
				if (B != R2.C && B != R2.CA && B != R2.N && B != R2.O) {
					continue;
				}
				if (R2.atom[j].type == 'H') {
					continue;
				}
				double d = Atom::Dist(R1.atom[i], R2.atom[j]);
				dmin = d < dmin ? d : dmin;
			}
		}
		break;

	default:
		return dmin;
		break;

	}

	return dmin;

}
