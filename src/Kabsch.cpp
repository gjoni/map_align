/*
 * Kabsch.cpp
 *
 *  Created on: Jun 2, 2015
 *      Author: ivan
 */

#include <cmath>
#include <cstring>
#include <cstdlib>
#include <cstdio>

#define SWAP_KABSCH(a,b) do { double tmp = b ; b = a ; a = tmp ; } while(0)

double Normalize(double *a) {

	double norm = sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);

	a[0] /= norm;
	a[1] /= norm;
	a[2] /= norm;

	return norm;

}

void Cross(double *a, double *b, double *c) {

	c[0] = a[1] * b[2] - b[1] * a[2];
	c[1] = a[2] * b[0] - a[0] * b[2];
	c[2] = a[0] * b[1] - a[1] * b[0];

}

double Dot(double *a, double *b) {

	return (a[0] * b[0] + a[1] * b[1] + a[2] * b[2]);

}

bool Kabsch(double **x, double **y, int n, int mode, double *rms, double t[3],
		double u[3][3]) {

	/* centers of mass for x(i) and y(i) */
	double comx[3] = { 0.0, 0.0, 0.0 };
	double comy[3] = { 0.0, 0.0, 0.0 };
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < 3; j++) {
			comx[j] += x[i][j];
			comy[j] += y[i][j];
		}
	}
	for (int i = 0; i < 3; i++) {
		comx[i] /= n;
		comy[i] /= n;
	}

	/* calculate correlation matrix R and determine E0 */
	double E0 = 0.0;
	double R[3][3] = { { 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0 } };

	for (int i = 0; i < n; i++) {

		/* shift x and y to their centers of mass */
		double x_shift[3] = { x[i][0] - comx[0], x[i][1] - comx[1], x[i][2]
				- comx[2] };

		double y_shift[3] = { y[i][0] - comy[0], y[i][1] - comy[1], y[i][2]
				- comy[2] };

		for (int j = 0; j < 3; j++) {

			/*
			 * E0 = 1/2 * sum(over j): y(j)*y(j) + x(j)*x(j)
			 */
			E0 += x_shift[j] * x_shift[j] + y_shift[j] * y_shift[j];

			for (int k = 0; k < 3; k++) {
				/*
				 * correlation matrix R:
				 *   R[j,k] = sum(over i): y(i,j) * x(i,k)
				 */
				R[j][k] += y_shift[j] * x_shift[k];

			}
		}

	}
	E0 *= 0.5;

	/*
	 * form symmetric matrix M = trans(R) x R (only upper triangle is saved):
	 *
	 *  M = d0 e0 f0
	 *         d1 e1
	 *            d2
	 */
	double d0, d1, d2, e0, e1, f0;
	{
		d0 = R[0][0] * R[0][0] + R[1][0] * R[1][0] + R[2][0] * R[2][0];

		/* divide matrix elements by d0, so that the cubic root algorithm
		 * can be applied */
		double _d0 = 1.0 / d0;

		d1 = (R[0][1] * R[0][1] + R[1][1] * R[1][1] + R[2][1] * R[2][1]) * _d0;
		d2 = (R[0][2] * R[0][2] + R[1][2] * R[1][2] + R[2][2] * R[2][2]) * _d0;

		e0 = (R[0][0] * R[0][1] + R[1][0] * R[1][1] + R[2][0] * R[2][1]) * _d0;
		e1 = (R[0][1] * R[0][2] + R[1][1] * R[1][2] + R[2][1] * R[2][2]) * _d0;

		f0 = (R[0][0] * R[0][2] + R[1][0] * R[1][2] + R[2][0] * R[2][2]) * _d0;

	}

	/* find cubic roots */
	double e[3] = { 0.0, 0.0, 0.0 }; /* eigenvalues */
	{
		double B, C, D, q, r;

		/*
		 * solving for eigenvalues of det(M - I * lambda) = 0  (1)
		 */

		/* the above eigenvalue problem (1) is represented in the form of
		 * a cubic equation: x^3 + B x^2 + C x + D = 0  (2) */
		B = -1.0 - d1 - d2;
		C = d1 + d2 + d1 * d2 - e0 * e0 - f0 * f0 - e1 * e1;
		D = e0 * e0 * d2 + e1 * e1 + f0 * f0 * d1 - d1 * d2 - 2 * e0 * f0 * e1;

		/*
		 **************************************************
		 * the code below is adapted from the GSL library *
		 * (poly/solve_cubic.c)                           *
		 **************************************************
		 */

		/* reduction of (2) to a depressed cubic:
		 * x^3 + p x + q = 0  (3) */
		q = B * B - 3.0 * C;
		r = 2.0 * B * B * B - 9.0 * B * C + 27.0 * D;

		double Q = q / 9;
		double R = r / 54;

		double Q3 = Q * Q * Q;
		double R2 = R * R;

		double CR2 = 729 * r * r;
		double CQ3 = 2916 * q * q * q;

		if (R == 0 && Q == 0) {

			e[0] = -B / 3;
			e[1] = -B / 3;
			e[2] = -B / 3;

		} else if (CR2 == CQ3) {

			/* this test is actually R2 == Q3, written in a form suitable
			 for exact computation with integers */

			/* Due to finite precision some double roots may be missed, and
			 considered to be a pair of complex roots z = x +/- epsilon i
			 close to the real axis. */

			double sqrtQ = sqrt(Q);

			if (R > 0) {
				e[0] = sqrtQ - B / 3;
				e[1] = sqrtQ - B / 3;
				e[2] = -2 * sqrtQ - B / 3;
			} else {
				e[0] = 2 * sqrtQ - B / 3;
				e[1] = -sqrtQ - B / 3;
				e[2] = -sqrtQ - B / 3;
			}

		} else if (R2 < Q3) {

			double sgnR = (R >= 0 ? 1 : -1);
			double ratio = sgnR * sqrt(R2 / Q3);
			double theta = acos(ratio);
			double norm = -2 * sqrt(Q);
			e[0] = norm * cos(theta / 3) - B / 3;
			e[1] = norm * cos((theta + 2.0 * M_PI) / 3) - B / 3;
			e[2] = norm * cos((theta - 2.0 * M_PI) / 3) - B / 3;

			/* Sort e[0], e[1], e[2] into DECREASING order */
			if (e[2] > e[1])
				SWAP_KABSCH(e[2], e[1]);

			if (e[1] > e[0]) {
				SWAP_KABSCH(e[1], e[0]);

				if (e[2] > e[1])
					SWAP_KABSCH(e[2], e[1]);

			}

		} else {

			/* one root */
			double sgnR = (R >= 0 ? 1 : -1);
			double C1 = -sgnR * pow(fabs(R) + sqrt(R2 - Q3), 1.0 / 3.0);
			double C2 = Q / C1;
			e[0] = C1 + C2 - B / 3;

			//return false;

		}

		/*
		 * end GSL
		 */

	}

	/* due to finite precision some roots may appear to be negative -
	 * make them equal to zero */
	int nzero = 0; /* number of zero roots */
	for (int i = 0; i < 3; i++) {
		if (e[i] < 1.0e-10) {
			e[i] = 0.0;
			nzero++;
		}
	}

	/* undo the d0 norm to get eigenvalues */
	e[0] *= d0;
	e[1] *= d0;
	e[2] *= d0;

	d1 *= d0;
	d2 *= d0;
	e0 *= d0;
	e1 *= d0;
	f0 *= d0;

	/*
	 if (mode == 3) {

	 double omega;
	 double det = d0 * d1 * d2 + 2.0 * e0 * e1 * f0 - d1 * f0 * f0
	 - e1 * e1 * d0 - e0 * e0 * d2;
	 if (det > 0.0) {
	 omega = 1.0;
	 } else {
	 omega = -1.0;
	 }

	 E0 = E0 - sqrt(e[0]) - sqrt(e[1]) - omega * sqrt(e[2]);
	 *rms = sqrt(2.0 * E0 / n);

	 return true;

	 }
	 */

	/*
	 * compute eigenvectors of trans(R) x R
	 * and store them in a(i)
	 */
	double a[3][3];

	int imax = 2;
	if (nzero == 2) {

		imax = 1;

		/* arbitrary set a(2) = [ 1.0, 0.0, 0.0 ] */
		a[1][0] = 1.0;
		a[1][1] = a[1][2] = 0.0;
	}

	/* find first and second eigenvectors a(1) and a(2) */
	for (int i = 0; i < imax; i++) {

		/*
		 * solve linear system Cx = 0 by Gaussian elimination
		 */

		/* matrix C of the linear system */
		double c[3][3] = { { d0 - e[i], e0, f0 }, { e0, d1 - e[i], e1 }, { f0,
				e1, d2 - e[i] } };

		int idx[3];

		/* find the largest coefficient (by the absolute value)
		 * in the first column of C and renumber rows such that
		 * it appears in the first row (new numbering is stored in idx[]) */
		if (fabs(c[0][0]) >= fabs(c[1][0]) && fabs(c[0][0]) >= fabs(c[2][0])) {
			idx[0] = 0;
			idx[1] = 1;
			idx[2] = 2;
		} else if (fabs(c[1][0]) >= fabs(c[0][0])
				&& fabs(c[1][0]) >= fabs(c[2][0])) {
			idx[0] = 1;
			idx[1] = 0;
			idx[2] = 2;
		} else {
			idx[0] = 2;
			idx[1] = 0;
			idx[2] = 1;

		}

		c[idx[0]][1] /= c[idx[0]][0];
		c[idx[0]][2] /= c[idx[0]][0];
		c[idx[1]][1] -= c[idx[1]][0] * c[idx[0]][1];
		c[idx[1]][2] -= c[idx[1]][0] * c[idx[0]][2];
		c[idx[2]][1] -= c[idx[2]][0] * c[idx[0]][1];
		c[idx[2]][2] -= c[idx[2]][0] * c[idx[0]][2];

		if (fabs(c[idx[1]][1]) < fabs(c[idx[2]][1])) {
			idx[1] = idx[2];
		}

		a[i][2] = -c[idx[1]][1];
		a[i][1] = c[idx[1]][2];

		a[i][0] = -a[i][1] * c[idx[0]][1] - a[i][2] * c[idx[0]][2];

		Normalize(a[i]);

	}

	/* set a(3) = a(1) x a(2)
	 * a(1) and a(2) are already normalized, so there is no need to
	 * normalize a(3) */
	Cross(a[0], a[1], a[2]);

	/* b(i) = Ra(i) */
	double b[3][3];

	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			b[i][j] = R[j][0] * a[i][0] + R[j][1] * a[i][1] + R[j][2] * a[i][2];
		}
	}

	double b3[3];
	Normalize(b[0]);
	Normalize(b[1]);
	Cross(b[0], b[1], b3);

	int omega = -1;
	if (Dot(b3, b[2]) > 0.0) {
		omega = 1;
	}

	E0 = E0 - sqrt(e[0]) - sqrt(e[1]) - omega * sqrt(e[2]);
	if (E0 < 0.0) {
		E0 = 0.0;
	}
	*rms = sqrt(2.0 * E0 / n);

	memcpy(b[2], b3, 3 * sizeof(double));

	memset(u, 0, 9 * sizeof(double));

	memcpy(t, comy, 3 * sizeof(double));

	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {

			/*
			 * u[i,j] = sum(over k): b(k,i) * a(k,j)
			 */
			for (int k = 0; k < 3; k++) {
				u[i][j] += b[k][i] * a[k][j];
			}

			/*
			 * t[i] = comy[i] - sum(over j): u[i][j] * comx[j]
			 */
			t[i] -= u[i][j] * comx[j];
		}
	}

	return true;

}
