/* src/ctl.h.  Generated from ctl.h.in by configure.  */
/* libctl: flexible Guile-based control files for scientific software 
 * Copyright (C) 1998-2014 Massachusetts Institute of Technology and Steven G. Johnson
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA  02111-1307, USA.
 *
 * Steven G. Johnson can be contacted at stevenj@alum.mit.edu.
 */

/* 
 * file ctl.h from libctl/src directory filtered to retain
 *            only what is necessary for libctl-geom
 */

#ifndef CTL_NOSCHEME_H
#define CTL_NOSCHEME_H

#ifdef __cplusplus
extern "C" {
#endif                          /* __cplusplus */

  /* Basic types: */

typedef int integer;
typedef double number;
typedef struct { number re, im; } cnumber; /* complex numbers! */
typedef short boolean;
typedef char *string;

  /* define vector3 as a structure, not an array, so that it can
     be a function return value and so that simple assignment works. */
typedef struct { number x,y,z; } vector3;

  /* similarly for matrix3x3 */
typedef struct { vector3 c0, c1, c2; /* the columns */ } matrix3x3;

/* define complex equivalents: */
typedef struct { cnumber x,y,z; } cvector3;
typedef struct { cvector3 c0, c1, c2; /* the columns */ } cmatrix3x3;

/**************************************************************************/

  /* vector3 and matrix3x3 utilities: */

extern number vector3_dot(vector3 v1,vector3 v2);
extern number vector3_norm(vector3 v);
extern vector3 vector3_scale(number s, vector3 v);
extern vector3 unit_vector3(vector3 v);
extern vector3 vector3_cross(vector3 v1,vector3 v2);
extern vector3 vector3_plus(vector3 v1,vector3 v2);
extern vector3 vector3_minus(vector3 v1,vector3 v2);
extern int vector3_equal(vector3 v1, vector3 v2);

extern vector3 matrix3x3_vector3_mult(matrix3x3 m, vector3 v);
extern vector3 matrix3x3_transpose_vector3_mult(matrix3x3 m, vector3 v);
extern matrix3x3 matrix3x3_mult(matrix3x3 m1, matrix3x3 m2);
extern matrix3x3 matrix3x3_transpose(matrix3x3 m);
extern number matrix3x3_determinant(matrix3x3 m);
extern matrix3x3 matrix3x3_inverse(matrix3x3 m);
extern int matrix3x3_equal(matrix3x3 m1, matrix3x3 m2);

extern vector3 matrix3x3_row1(matrix3x3 m);
extern vector3 matrix3x3_row2(matrix3x3 m);
extern vector3 matrix3x3_row3(matrix3x3 m);

/**************************************************************************/

  /* complex number utilities */


extern cnumber make_cnumber(number r, number i);
extern cnumber cnumber_conj(cnumber c);
extern int cnumber_equal(cnumber c1, cnumber c2);
#define cnumber_re(c) ((c).re)
#define cnumber_im(c) ((c).im)

extern vector3 cvector3_re(cvector3 cv);
extern vector3 cvector3_im(cvector3 cv);
extern cvector3 make_cvector3(vector3 vr, vector3 vi);
extern int cvector3_equal(cvector3 v1, cvector3 v2);

extern matrix3x3 cmatrix3x3_re(cmatrix3x3 cm);
extern matrix3x3 cmatrix3x3_im(cmatrix3x3 cm);
extern cmatrix3x3 make_cmatrix3x3(matrix3x3 mr, matrix3x3 mi);
cmatrix3x3 make_hermitian_cmatrix3x3(number m00, number m11, number m22,
                                     cnumber m01, cnumber m02, cnumber m12);
extern int cmatrix3x3_equal(cmatrix3x3 m1, cmatrix3x3 m2);

/**************************************************************************/

typedef double (*integrand) (unsigned ndim, const double *x, void *);
double adaptive_integration(integrand f, double *xmin, double *xmax,
			    int n, void *fdata,
			    double abstol, double reltol, int maxnfe,
			    double *esterr, int *errflag);

/**************************************************************************/

#ifdef __cplusplus
	   }                               /* extern "C" */
#endif                          /* __cplusplus */

#endif /* CTL_NOSCHEME_H */
