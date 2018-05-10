/*
 * mp_trim.cpp
 *
 *  Created on: May 10, 2018
 *      Author: aivan
 */

#include <unistd.h>
#include <cassert>
#include <cstring>

#include "Chain.h"
#include "Info.h"

#define VERSION "V20180510"

struct OPTS {
	std::string in; /* input PDB to trim */
	std::string out; /* output (trimmed) PDB */
	std::string tmplt; /* template PDB with B-factors */
	double dmax; /* distance cut-off - for trimming*/
};

bool GetOpts(int argc, char *argv[], OPTS &opts);
void PrintOpts(const OPTS &opts);
void PrintCap(const OPTS &opts);

int main(int argc, char *argv[]) {

	/*
	 * (0) process input parameters
	 */
	OPTS opts = { "", "", "", 999.9 };
	if (!GetOpts(argc, argv, opts)) {
		PrintOpts(opts);
		return 1;
	}

	/*
	 * (1) read PDBs
	 */
	Chain IN(opts.in.c_str());
	Chain REF(IN);
	if (opts.tmplt != "") {
		REF = Chain(opts.tmplt.c_str());
	}

	assert(IN.nRes == REF.nRes);

	/*
	 * (2) process
	 */
	FILE *F = fopen(opts.out.c_str(), "w");
	if (F == NULL) {
		printf("Error: cannot open %s for saving (-o)\n", opts.out.c_str());
		exit(1);
	}

	char *seq_in = (char*) malloc((REF.nRes + 1) * sizeof(char));
	char *seq_out = (char*) malloc((REF.nRes + 1) * sizeof(char));
	seq_in[REF.nRes] = seq_out[REF.nRes] = '\0';

	int remained = 0;

	for (int i = 0; i < REF.nRes; i++) {

		Residue *R = &(IN.residue[i]);
		seq_in[i] = AAA1[(int) R->type];

		/* skip residues with large B-factors */
		if (REF.residue[i].CA->temp > opts.dmax) {
			seq_out[i] = '-';
			continue;
		} else {
			seq_out[i] = seq_in[i];
			remained++;
		}

		/* save otherwise */
		for (int j = 0; j < R->nAtoms; j++) {
			Atom *A = &(R->atom[j]);
			if (strlen(A->name) == 4 || A->name[0] == '1' || A->name[0] == '2'
					|| A->name[0] == '3' || A->name[0] == '4') {
				fprintf(F,
						"ATOM  %5d %-4s%c%3s %c%4d%c   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2s\n",
						A->atomNum, A->name, A->altLoc, R->name, R->chainId,
						R->seqNum, R->insCode, A->x, A->y, A->z, A->occup,
						A->temp, A->element, A->charge);
			} else {
				fprintf(F,
						"ATOM  %5d  %-3s%c%3s %c%4d%c   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2s\n",
						A->atomNum, A->name, A->altLoc, R->name, R->chainId,
						R->seqNum, R->insCode, A->x, A->y, A->z, A->occup,
						A->temp, A->element, A->charge);
			}
		}

	}
	fprintf(F, "TER\n");
	fclose(F);

	printf("Coverage= %.3f\n", 1.0 * remained / REF.nRes);
	printf("%s\n", seq_in);
	printf("%s\n", seq_out);

	free(seq_in);
	free(seq_out);

	return 0;

}

void PrintOpts(const OPTS &opts) {

	printf("\nUsage:   ./mp_trim [-option] [argument]\n\n");
	printf("Options:  -i input.pdb                   - input, required\n");
	printf("          -o output.pdb                  - output, required\n");
	printf("          -t template.pdb                - input, optional\n");
	printf("          -d distance cut-off              %.1f\n", opts.dmax);

}

void PrintCap(const OPTS &opts) {

	time_t timer;
	time(&timer);
	struct tm* tm_info = localtime(&timer);
	char buf[100];
	strftime(buf, 26, "%Y:%m:%d / %H:%M:%S", tm_info);

	printf("# %s\n", std::string(70, '-').c_str());
	printf("# mp_trim - a program to trim PDBs based on B-factors  %15s\n",
	VERSION);
	printf("# %s\n", std::string(70, '-').c_str());

	printf("# %20s : %s\n", "input PDB", opts.in.c_str());
	printf("# %20s : %s\n", "output PDB", opts.out.c_str());
	printf("# %20s : %s\n", "template PDB", opts.tmplt.c_str());
	printf("# %20s : %.1f\n", "distance cut-off", opts.dmax);

	printf("# %s\n", std::string(70, '-').c_str());

}

bool GetOpts(int argc, char *argv[], OPTS &opts) {

	char tmp;
	while ((tmp = getopt(argc, argv, "hl:i:o:t:d:")) != -1) {
		switch (tmp) {
		case 'h': /* help */
			printf("!!! HELP !!!\n");
			return false;
			break;
		case 'i': /* input PDB */
			opts.in = std::string(optarg);
			break;
		case 'o': /* output PDB */
			opts.out = std::string(optarg);
			break;
		case 't': /* template PDB */
			opts.tmplt = std::string(optarg);
			break;
		case 'd': /* distance cut-off */
			opts.dmax = atof(optarg);
			break;
		default:
			return false;
			break;
		}
	}

	if (opts.in == "") {
		printf("Error: input PDB not specified ('-i')\n");
		return false;
	}

	if (opts.out == "") {
		printf("Error: output PDB not specified ('-i')\n");
		return false;
	}

	return true;

}
