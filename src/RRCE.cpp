/*
 * RRCE.cpp
 *
 *  Created on: May 27, 2016
 *      Author: ivan
 */

#include <cstdlib>
#include <cstring>
#include <cassert>
#include <cstdio>

#include "RRCE.h"

RRCE::RRCE(RRCE_TYPE type, double d, int k) :
		h(NULL), J(NULL), dmax(d) {

	char *DATADIR = getenv("TMDOCKDAT");
	if (DATADIR == NULL) {
		printf("Error: environment variable 'TMDOCKDAT' not set\n");
		exit(1);
	}

	const int STRMAX = 1000;
	char tname[STRMAX], rrce_name[STRMAX];

	switch (type) {

	case RRCE20RC:
		sprintf(tname, "%s/RRCE20RC/table.%.1fA_k%d", DATADIR, dmax, k);
		sprintf(rrce_name, "RRCE20RC_%.1fA_k%d", dmax, k);
		break;

	case RRCE20SCC:
		sprintf(tname, "%s/RRCE20SCC/table.%.1fA_k%d", DATADIR, dmax, k);
		sprintf(rrce_name, "RRCE20SCC_%.1fA_k%d", dmax, k);
		break;

	case RRCE20CB:
		sprintf(tname, "%s/RRCE20CB/table.%.1fA_k%d", DATADIR, dmax, k);
		sprintf(rrce_name, "RRCE20CB_%.1fA_k%d", dmax, k);
		break;

	default:
		printf("Error: wrong RRCE type (%d)\n", type);
		exit(1);
		break;

	}

	Allocate(20);

	ReadCouplings(tname);

	name = rrce_name;

}

RRCE::RRCE(const RRCE &source) :
		h(NULL), J(NULL), dmax(source.dmax), name(source.name) {

	const int N = 20;

	Allocate(N);

	memcpy(h, source.h, N * sizeof(double));
	for (int i = 0; i < N; i++) {
		memcpy(J[i], source.J[i], N * sizeof(double));
	}

}

RRCE::RRCE(const char *fname) {

	Allocate(20);
	ReadCouplings(fname);
	name = fname;

}

RRCE::RRCE() :
		h(NULL), J(NULL), dmax(0.0), name("") {

}

RRCE::~RRCE() {

	Free(20);

}

RRCE & RRCE::operator =(const RRCE & source) {

	assert(this != &source);

	const int N = 20;

	Free(N);

	dmax = source.dmax;
	name = source.name;

	Allocate(N);

	memcpy(h, source.h, N * sizeof(double));
	for (int i = 0; i < N; i++) {
		memcpy(J[i], source.J[i], N * sizeof(double));
	}

	return *this;

}

void RRCE::Allocate(int N) {

	h = (double*) malloc(N * sizeof(double));
	J = (double**) malloc(N * sizeof(double*));
	for (int i = 0; i < N; i++) {
		J[i] = (double*) malloc(N * sizeof(double));
	}

}

void RRCE::Free(int N) {

	if (J != NULL) {
		for (int i = 0; i < N; i++) {
			free(J[i]);
		}
		free(J);
	}
	free(h);

}

int RRCE::ReadCouplings(const char *name) {

	FILE *F = fopen(name, "r");
	if (F == NULL) {
		printf("Error: cannot open %s file\n", name);
		exit(1);
	}

	const int STRMAX = 10000;
	char buf[STRMAX];

	const int N = 20;

	for (int i = 0; i < N; i++) {
		fgets(buf, STRMAX, F);
		char * pch;
		pch = strtok(buf, " "); /* label */
		pch = strtok(NULL, " "); /* local field */
		h[i] = atof(pch);
		for (int j = 0; j < N; j++) {
			pch = strtok(NULL, " ");
			J[i][j] = atof(pch);
		}
	}
	fclose(F);

	return 0;

}

std::string RRCE::GetName() {

	return name;

}

void RRCE::GetCouplings(double **J_) {

	for (int i = 0; i < 20; i++) {
		for (int j = 0; j < 20; j++) {
			J_[i][j] = J[i][j];
		}
	}

}

void RRCE::GetFields(double *h_) {

	for (int i = 0; i < 20; i++) {
		h_[i] = h[i];
	}

}
