/*
 * Atom.cpp
 *
 *  Created on: Jun 3, 2014
 *      Author: ivan
 */

#include "Atom.h"

#include <cstring>
#include <cstdlib>
#include <cmath>

Atom::Atom() :
		atomNum(1), altLoc(' '), x(0.0), y(0.0), z(0.0), occup(0.0), temp(0.0), residue(
		NULL), q(0.0), r(0.0), type(' ') {

	name[0] = '\n';
	x = y = z = 0.0;
	element[0] = '\n';
	charge[0] = '\n';

}

Atom::Atom(const Atom &source) :
		atomNum(source.atomNum), altLoc(source.altLoc), x(source.x), y(
				source.y), z(source.z), occup(source.occup), temp(source.temp), residue(
				source.residue), q(source.q), r(source.r), type(source.type) {

	strcpy(name, source.name);
	strcpy(element, source.element);
	strcpy(charge, source.charge);

}

Atom::Atom(const AtomRecord &source) :
		atomNum(source.atomNum), altLoc(source.altLoc), x(source.x), y(
				source.y), z(source.z), occup(source.occup), temp(source.temp), residue(
		NULL), q(source.q), r(source.r), type(' ') {

	strcpy(name, source.atomName);
	strcpy(element, source.element);
	strcpy(charge, source.charge);
	type = GetType(name);

}

Atom::~Atom() {

	/* nothing to do */

}

Atom& Atom::operator =(const Atom& source) {

	atomNum = source.atomNum;
	altLoc = source.altLoc;
	strcpy(name, source.name);
	x = source.x;
	y = source.y;
	z = source.z;
	occup = source.occup;
	temp = source.temp;
	strcpy(element, source.element);
	strcpy(charge, source.charge);
	residue = source.residue;
	q = source.q;
	r = source.r;
	type = source.type;

	return *this;

}

char Atom::GetType(char name[]) {

	if (name[0] == '0' || name[0] == '1' || name[0] == '2' || name[0] == '3'
			|| name[0] == '4' || name[0] == '5' || name[0] == '6'
			|| name[0] == '7' || name[0] == '8' || name[0] == '9') {

		/* first letter is a number */
		return name[1];

	} else {

		/* all other cases */
		return name[0];

	}

}

double Atom::Dist(const Atom &A, const Atom &B) {

	double x = A.x - B.x;
	double y = A.y - B.y;
	double z = A.z - B.z;

	return sqrt(x * x + y * y + z * z);

}
