/*****************************************************************************
 *
 * 2015-2016, Ivan Anishchenko <anishchenko.ivan@gmail.com>
 * Vakser Lab, Center for Computational biology, The University of Kansas
 *
 * V20160215
 *
 *****************************************************************************
 *
 * Implementation of the Kabsch algorithm [1,2] to find RMSD, translation
 * vector and the least-squares rotation matrix to superpose two vectors of
 * same dimension.
 *
 * [1] W Kabsch. A solution for the best rotation to relate two sets of
 *     vectors. Acta Cryst (1976). A32:922-3
 * [2] W Kabsch. A discussion of the solution for the best rotation to relate
 *     two sets of vectors. Acta Cryst (1978). A34:827-8
 *
 * This implementation was inspired by the
 * 'rmsd.c' (c) 2005 Bosco K Ho routine.
 *
 *****************************************************************************/

#ifndef KABSCH_H_
#define KABSCH_H_

/*
 * x    - x(m,i) are coordinates of atom m in set x           (input)
 * y    - y(m,i) are coordinates of atom m in set y           (input)
 * n    - n is number of atom pairs                           (input)
 * mode - 0:calculate rms only                                (input)
 *        1:calculate rms,u,t                                 (takes longer)
 * rms  - sum of (ux+t-y)**2 over all atom pairs              (output)
 * u    - u(i,j) is  rotation  matrix for best superposition  (output)
 * t    - t(i)  is translation vector for best superposition  (output)
 */

bool Kabsch(double **x, double **y, int n, int mode, double *rms, double t[3],
		double u[3][3]);

#endif /* KABSCH_H_ */
