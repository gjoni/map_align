/*
 * Chain.cpp
 *
 *  Created on: Jun 3, 2014
 *      Author: ivan
 */

#include "Chain.h"
#include "Info.h"

#include <cmath>
#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <cassert>
#include <vector>

Chain::Chain() :
		atom(NULL), nAtoms(0), residue(NULL), nRes(0), kd(NULL), kdCA(
		NULL), kdCent(NULL), ca_trace(NULL) {

}

Chain::Chain(const std::string &name) :
		Chain(name.c_str()) {

	/* */

}

Chain::Chain(const char *name) :
		atom(NULL), nAtoms(0), residue(NULL), nRes(0), kd(NULL), kdCA(
		NULL), kdCent(NULL), ca_trace(NULL) {

	/* open file */
	FILE *F;
	F = fopen(name, "r");
	if (F == NULL) {
		fprintf(stderr, "ERROR: Can't open PDB file '%s'\n", name);
		exit(1);
	}

	/* read file */
	char buf[81];
	std::vector<AtomRecord> atomRecordVector;

	/* read according to type */
	if (strncmp(name + strlen(name) - 4, ".pqr", 4) == 0) {

		/* read as PQR */
		while (fgets(buf, 81, F) != NULL) {
			if ((buf[0] == 'A') && (buf[1] == 'T')) {
				atomRecordVector.push_back(ReadAtomRecordPQR(buf));
			}
		}

	} else {

		/* read as PDB (in all other cases) */
		while (fgets(buf, 81, F) != NULL) {
			if ((buf[0] == 'A') && (buf[1] == 'T')) {
				atomRecordVector.push_back(ReadAtomRecord(buf));
			}
		}

	}
	fclose(F);

	if (atomRecordVector.size() == 0) {
		fprintf(stderr, "ERROR: No atoms were read from %s\n", name);
		exit(1);
	}

	InitChain(atomRecordVector);

	SetResMap();

}

Chain::Chain(const std::vector<std::string> &atoms_str) :
		atom(NULL), nAtoms(0), residue(NULL), nRes(0), kd(NULL), kdCA(
		NULL), kdCent(NULL), ca_trace(NULL) {

	std::vector<AtomRecord> atomRecordVector(atoms_str.size());
	char buf[81];
	for (unsigned i = 0; i < atoms_str.size(); i++) {
		strcpy(buf, atoms_str[i].c_str());
		atomRecordVector[i] = ReadAtomRecord(buf);
	}

	InitChain(atomRecordVector);

	SetResMap();

}

void Chain::InitChain(std::vector<AtomRecord> &atomRecordVector) {

	/* define boundaries for all residues : insCode + resNum */
	std::vector<int> first, last; /* residues boundaries */
	first.push_back(0); /* first atom of the first residue */
	nAtoms = atomRecordVector.size();
	for (int i = 1; i < nAtoms; i++) {
		if (atomRecordVector[i].insCode != atomRecordVector[i - 1].insCode
				|| atomRecordVector[i].resNum
						!= atomRecordVector[i - 1].resNum) {
			first.push_back(i);
			last.push_back(i - 1);
		}
	}
	last.push_back(nAtoms - 1); /* last atom of the last residue */

	assert(first.size() == last.size()); /* unexpectedly :-) */

	/*
	 * TODO: temp. solution for eliminating incomplete residues
	 */
	bool *fl = (bool*) calloc(first.size(), sizeof(bool));
	nRes = 0;
	for (unsigned i = 0; i < first.size(); i++) {
		for (int j = first[i]; j <= last[i]; j++) {
			if (strcmp(atomRecordVector[j].atomName, "CA") == 0) {
				fl[i] = true;
			}
		}
		nRes += fl[i];
	}

	/* initialize residues */
	residue = new Residue[nRes];
	int idx = 0;
	nAtoms = 0;
	for (unsigned i = 0; i < first.size(); i++) {
		if (fl[i]) {
			residue[idx] = Residue(atomRecordVector.begin() + first[i],
					atomRecordVector.begin() + last[i]);
			//nAtoms += last[i] - first[i] + 1;
			nAtoms += residue[idx].nAtoms;
			idx++;
		}
	}
	free(fl);

	/* initialize pointers to atoms */
	atom = new Atom*[nAtoms]; /* array of pointers to atoms */
	SetAtomPointers();

	SetKD();

	/* initialize CA trace */
	ca_trace = (double**) malloc(nRes * sizeof(double*));
	for (int i = 0; i < nRes; i++) {
		ca_trace[i] = (double*) malloc(3 * sizeof(double));
		ca_trace[i][0] = residue[i].CA->x;
		ca_trace[i][1] = residue[i].CA->y;
		ca_trace[i][2] = residue[i].CA->z;
	}

}

Chain::Chain(const Chain & source) :
		resmap(source.resmap), atom(NULL), nAtoms(source.nAtoms), residue(NULL), nRes(
				source.nRes), kd(
		NULL), kdCA(NULL), kdCent(NULL), ca_trace(
		NULL) {

	/* copy residues */
	residue = new Residue[nRes];
	for (int i = 0; i < nRes; i++) {
		residue[i] = source.residue[i];
	}

	/* set pointers to atoms */
	atom = new Atom*[nAtoms];

	SetAtomPointers();

	SetKD();

	/* initialize CA trace */
	ca_trace = (double**) malloc(nRes * sizeof(double*));
	for (int i = 0; i < nRes; i++) {
		ca_trace[i] = (double*) malloc(3 * sizeof(double));
		memcpy(ca_trace[i], source.ca_trace[i], 3 * sizeof(double));
	}

	SetResMap();

}

Chain::~Chain() {

	delete[] atom;
	atom = NULL;

	delete[] residue;
	residue = NULL;

	if (kd != NULL) {
		kd_free(kd);
	}

	if (kdCA != NULL) {
		kd_free(kdCA);
	}

	if (kdCent != NULL) {
		kd_free(kdCent);
	}

	if (ca_trace != NULL) {
		for (int i = 0; i < nRes; i++) {
			free(ca_trace[i]);
		}
		free(ca_trace);
	}

}

Chain & Chain::operator =(const Chain & source) {

	delete[] atom;
	delete[] residue;
	for (int i = 0; i < nRes; i++) {
		free(ca_trace[i]);
	}
	free(ca_trace);

	/* copy residues */
	nRes = source.nRes;
	residue = new Residue[nRes];
	for (int i = 0; i < nRes; i++) {
		residue[i] = source.residue[i];
	}

	/* set pointers to atoms */
	nAtoms = source.nAtoms;
	atom = new Atom*[nAtoms];

	SetAtomPointers();

	SetKD();

	//resmap = source.resmap;
	SetResMap();

	/* initialize CA trace */
	ca_trace = (double**) malloc(nRes * sizeof(double*));
	for (int i = 0; i < nRes; i++) {
		ca_trace[i] = (double*) malloc(3 * sizeof(double));
		ca_trace[i][0] = residue[i].CA->x;
		ca_trace[i][1] = residue[i].CA->y;
		ca_trace[i][2] = residue[i].CA->z;
	}

	return *this;

}

AtomRecord Chain::ReadAtomRecord(char* str) {

	AtomRecord A;

	/* complete string to 80 characters */
	char *pch = strpbrk(str, "\n");
	if (pch != NULL && pch - str < 80) {
		while (pch != str + 80) {
			*pch++ = ' ';
		}
		pch = '\0';
	}

	/*
	 * reading PDB ATOM string in the reverse direction
	 */

	/* Charge on the atom */
	str[80] = '\0';
	strncpy(A.charge, str + 78, 3);

	/* Element symbol, right-justified */
	str[78] = '\0';
	strncpy(A.element, str + 76, 3);

	/* Segment identifier, left-justified */
	str[76] = '\0';
	strncpy(A.segId, str + 72, 5);

	/* Temperature factor */
	str[67] = '\0';
	A.temp = atof(str + 60);

	/* Occupancy */
	str[60] = '\0';
	A.occup = atof(str + 54);

	/* Z coordinate */
	str[54] = '\0';
	A.z = atof(str + 46);

	/* Y coordinate */
	str[46] = '\0';
	A.y = atof(str + 38);

	/* X coordinate */
	str[38] = '\0';
	A.x = atof(str + 30);

	/* Code for insertion of residues */
	A.insCode = str[26];

	/* Residue sequence number */
	str[26] = '\0';
	A.resNum = atoi(str + 22);

	/* Chain identifier */
	A.chainId = str[21];

	/* Residue name */
	str[20] = '\0';
	//strncpy(A.resName, str + 17, 4);
	sscanf(str + 17, "%s", A.resName); /* no spaces in residue name */

	/* Alternate location indicator */
	A.altLoc = str[16];

	/* Atom name */
	str[16] = '\0';
	//strncpy(A.atomName, str + 12, 5);
	sscanf(str + 12, "%s", A.atomName);

	/* Atom serial number */
	str[11] = '\0';
	A.atomNum = atoi(str + 6);

	return A;

}

AtomRecord Chain::ReadAtomRecordPQR(char* str) {

	AtomRecord A;

	/*
	 * reading PDB ATOM string in the reverse direction
	 */

	/* fields missing from PQR files */
	A.segId[0] = '\0';
	A.element[0] = '\0';
	A.charge[0] = '\0';

	/* save atom charge and radius in Temperature factor
	 * and Occupancy respectively */
	str[80] = '\0';
	sscanf(str + 54, "%lf %lf", &(A.q), &(A.r));
	A.temp = A.q;
	A.occup = A.r;

	/* Z coordinate */
	str[54] = '\0';
	A.z = atof(str + 46);

	/* Y coordinate */
	str[46] = '\0';
	A.y = atof(str + 38);

	/* X coordinate */
	str[38] = '\0';
	A.x = atof(str + 30);

	/* Code for insertion of residues */
	A.insCode = str[26];

	/* Residue sequence number */
	str[26] = '\0';
	A.resNum = atoi(str + 22);

	/* Chain identifier */
	A.chainId = str[21];

	/* Residue name */
	str[20] = '\0';
	//strncpy(A.resName, str + 17, 4);
	sscanf(str + 17, "%s", A.resName); /* no spaces in residue name */

	/* Alternate location indicator */
	A.altLoc = str[16];

	/* Atom name */
	str[16] = '\0';
	//strncpy(A.atomName, str + 12, 5);
	sscanf(str + 12, "%s", A.atomName);

	/* Atom serial number */
	str[11] = '\0';
	A.atomNum = atoi(str + 6);

	A.charge[0] = A.element[0] = '\0';

	return A;

}

void Chain::SetAtomPointers() {

	int idx = 0;
	for (int i = 0; i < nRes; i++) {
		for (int j = 0; j < residue[i].nAtoms; j++) {
			atom[idx] = &(residue[i].atom[j]);
			idx++;
		}
	}

}

void Chain::SetTermini() {

	/* first residue is N-terminal */
	residue[0].ntFlag = true;

	/* last residue is C-terminal */
	residue[nRes - 1].ctFlag = true;

	/* check residues in between */
	double x, y, z, d;
	for (int i = 1; i < nRes; i++) {
		if (residue[i - 1].CA != NULL && residue[i].CA != NULL) {
			x = residue[i - 1].CA->x - residue[i].CA->x;
			y = residue[i - 1].CA->y - residue[i].CA->y;
			z = residue[i - 1].CA->z - residue[i].CA->z;
			d = sqrt(x * x + y * y + z * z);
			if (d > 4.25) { /* as in TM-align */
				residue[i - 1].ctFlag = true;
				residue[i].ntFlag = true;
			}
		}
	}

}

void Chain::SetResMap() {

	resmap.clear();
	for (int i = 0; i < nRes; i++) {
		Residue *R = residue + i;
		resmap[ { R->seqNum, R->insCode }] = R;
	}

}

void Chain::SetKD() {

	if (kd != NULL) {
		kd_free(kd);
	}

	if (kdCA != NULL) {
		kd_free(kdCA);
	}

	if (kdCent != NULL) {
		kd_free(kdCent);
	}

	kd = kd_create(3);
	kdCA = kd_create(3);
	kdCent = kd_create(3);

	Atom *A;
	int idx = 0;
	for (int i = 0; i < nRes; i++) {

		/* kd */
		for (int j = 0; j < residue[i].nAtoms; j++) {
			A = &(residue[i].atom[j]);

			/* exclude hydrogens from the tree */
			if (A->type == 'H') {
				continue;
			}

			kd_insert3(kd, A->x, A->y, A->z, atom + idx);
			idx++;
		}

		A = residue[i].CA;
		if (A != NULL) {

			/* kdCA */
			kd_insert3(kdCA, A->x, A->y, A->z, A);

			/* kdCent - GLY */
			if (residue[i].type == 7 /* GLY */) {
				kd_insert3(kdCent, A->x, A->y, A->z, residue[i].atom);
			}

		} else {
			//printf("WARNING: Missing CA atom in %s-%d\n", residues[i].name,
			//		residues[i].seqNum);
		}

		/* kdCent - all other residues */
		if (residue[i].type != 7) {
			kd_insert(kdCent, residue[i].centroid, residue[i].atom);
		}

	}

}

void Chain::Transform(double t[3], double u[3][3]) {

	for (int i = 0; i < nRes; i++) {

		for (int j = 0; j < residue[i].nAtoms; j++) {

			Atom *A = &(residue[i].atom[j]);

			double x = A->x;
			double y = A->y;
			double z = A->z;

			A->x = u[0][0] * x + u[0][1] * y + u[0][2] * z + t[0];
			A->y = u[1][0] * x + u[1][1] * y + u[1][2] * z + t[1];
			A->z = u[2][0] * x + u[2][1] * y + u[2][2] * z + t[2];

		}

		residue[i].SetCentroid();

	}

	SetKD();

	for (int i = 0; i < nRes; i++) {
		ca_trace[i][0] = residue[i].CA->x;
		ca_trace[i][1] = residue[i].CA->y;
		ca_trace[i][2] = residue[i].CA->z;
	}

}

void Chain::Load(const char *name) {

	Chain tmp(name);
	*this = tmp;

}

void Chain::Save(const char *name) {

	FILE *F = fopen(name, "w");
	if (F == NULL) {
		fprintf(stderr, "ERROR: Cannot write to file '%s'\n", name);
		return;
	}
	Atom *A;
	for (int i = 0; i < nRes; ++i) {
		for (int j = 0; j < residue[i].nAtoms; j++) {
			A = &(residue[i].atom[j]);
			if (strlen(A->name) == 4 || A->name[0] == '1' || A->name[0] == '2'
					|| A->name[0] == '3' || A->name[0] == '4') {
				fprintf(F,
						"ATOM  %5d %-4s%c%3s %c%4d%c   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2s\n",
						A->atomNum, A->name, A->altLoc, residue[i].name,
						residue[i].chainId, residue[i].seqNum,
						residue[i].insCode, A->x, A->y, A->z, A->occup, A->temp,
						A->element, A->charge);
			} else {
				fprintf(F,
						"ATOM  %5d  %-3s%c%3s %c%4d%c   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2s\n",
						A->atomNum, A->name, A->altLoc, residue[i].name,
						residue[i].chainId, residue[i].seqNum,
						residue[i].insCode, A->x, A->y, A->z, A->occup, A->temp,
						A->element, A->charge);
			}
		}
	}
	fprintf(F, "TER\n");
	fclose(F);

}

int Chain::IfProtein() {

	int n = 0;
	for (int i = 0; i < nRes; i++) {
		n += (residue[i].type >= 0 && residue[i].type < 20);
	}
	return n;

}

std::string Chain::GetSequence() {

	std::string seq;
	for (int i = 0; i < nRes; i++) {
		seq += AAA1[(int) residue[i].type];
	}

	return seq;

}

void Chain::Renumber(const std::vector<int> &num) {

	assert((int )num.size() == nRes);

	for (int i = 0; i < nRes; i++) {
		residue[i].seqNum = num[i];
		residue[i].insCode = ' ';
	}

	SetResMap();

}

Residue* Chain::GetResidue(int n, char ins) {

	std::map<std::pair<int, char>, Residue*>::iterator it;
	it = resmap.find( { n, ins });
	if (it == resmap.end()) {
		return NULL;
	} else {
		return it->second;
	}

}

