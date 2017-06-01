
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
 * file ctl.c from libctl/src directory filtered to retain
 *            only what is necessary for libctl-geom
 */


#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "ctl-noscheme.h"

/**************************************************************************/

/* vector3 and matrix3x3 utilities: */

number vector3_dot(vector3 v1,vector3 v2)
{
  return (v1.x * v2.x + v1.y * v2.y + v1.z * v2.z);
}

number vector3_norm(vector3 v)
{
  return (sqrt(vector3_dot(v,v)));
}

vector3 vector3_scale(number s, vector3 v)
{
  vector3 vnew;

  vnew.x = s * v.x;
  vnew.y = s * v.y;
  vnew.z = s * v.z;
  return vnew;
}

vector3 unit_vector3(vector3 v)
{
  number norm = vector3_norm(v);
  if (norm == 0.0)
    return v;
  else
    return vector3_scale(1.0/norm, v);
}

vector3 vector3_plus(vector3 v1,vector3 v2)
{
  vector3 vnew;

  vnew.x = v1.x + v2.x;
  vnew.y = v1.y + v2.y;
  vnew.z = v1.z + v2.z;
  return vnew;
}

vector3 vector3_minus(vector3 v1,vector3 v2)
{
  vector3 vnew;

  vnew.x = v1.x - v2.x;
  vnew.y = v1.y - v2.y;
  vnew.z = v1.z - v2.z;
  return vnew;
}

vector3 vector3_cross(vector3 v1,vector3 v2)
{
  vector3 vnew;

  vnew.x = v1.y * v2.z - v2.y * v1.z;
  vnew.y = v1.z * v2.x - v2.z * v1.x;
  vnew.z = v1.x * v2.y - v2.x * v1.y;
  return vnew;
}

int vector3_equal(vector3 v1, vector3 v2)
{
     return (v1.x == v2.x && v1.y == v2.y && v1.z == v2.z);
}

vector3 matrix3x3_vector3_mult(matrix3x3 m, vector3 v)
{
  vector3 vnew;

  vnew.x = m.c0.x * v.x + m.c1.x * v.y + m.c2.x * v.z;
  vnew.y = m.c0.y * v.x + m.c1.y * v.y + m.c2.y * v.z;
  vnew.z = m.c0.z * v.x + m.c1.z * v.y + m.c2.z * v.z;
  return vnew;
}

vector3 matrix3x3_transpose_vector3_mult(matrix3x3 m, vector3 v)
{
  vector3 vnew;

  vnew.x = m.c0.x * v.x + m.c0.y * v.y + m.c0.z * v.z;
  vnew.y = m.c1.x * v.x + m.c1.y * v.y + m.c1.z * v.z;
  vnew.z = m.c2.x * v.x + m.c2.y * v.y + m.c2.z * v.z;
  return vnew;
}

matrix3x3 matrix3x3_mult(matrix3x3 m1, matrix3x3 m2)
{
  matrix3x3 m;

  m.c0.x = m1.c0.x * m2.c0.x + m1.c1.x * m2.c0.y + m1.c2.x * m2.c0.z;
  m.c0.y = m1.c0.y * m2.c0.x + m1.c1.y * m2.c0.y + m1.c2.y * m2.c0.z;
  m.c0.z = m1.c0.z * m2.c0.x + m1.c1.z * m2.c0.y + m1.c2.z * m2.c0.z;

  m.c1.x = m1.c0.x * m2.c1.x + m1.c1.x * m2.c1.y + m1.c2.x * m2.c1.z;
  m.c1.y = m1.c0.y * m2.c1.x + m1.c1.y * m2.c1.y + m1.c2.y * m2.c1.z;
  m.c1.z = m1.c0.z * m2.c1.x + m1.c1.z * m2.c1.y + m1.c2.z * m2.c1.z;

  m.c2.x = m1.c0.x * m2.c2.x + m1.c1.x * m2.c2.y + m1.c2.x * m2.c2.z;
  m.c2.y = m1.c0.y * m2.c2.x + m1.c1.y * m2.c2.y + m1.c2.y * m2.c2.z;
  m.c2.z = m1.c0.z * m2.c2.x + m1.c1.z * m2.c2.y + m1.c2.z * m2.c2.z;

  return m;
}

matrix3x3 matrix3x3_transpose(matrix3x3 m)
{
     matrix3x3 mt;
    
     mt.c0.x = m.c0.x;
     mt.c1.x = m.c0.y;
     mt.c2.x = m.c0.z;
     mt.c0.y = m.c1.x;
     mt.c1.y = m.c1.y;
     mt.c2.y = m.c1.z;
     mt.c0.z = m.c2.x;
     mt.c1.z = m.c2.y;
     mt.c2.z = m.c2.z;
     return mt;
}

number matrix3x3_determinant(matrix3x3 m)
{
     return(m.c0.x*m.c1.y*m.c2.z - m.c2.x*m.c1.y*m.c0.z +
	    m.c1.x*m.c2.y*m.c0.z + m.c0.y*m.c1.z*m.c2.x -
	    m.c1.x*m.c0.y*m.c2.z - m.c2.y*m.c1.z*m.c0.x);
}

matrix3x3 matrix3x3_inverse(matrix3x3 m)
{
     matrix3x3 minv;
     number detinv = matrix3x3_determinant(m);

     if (detinv == 0.0) {
	  fprintf(stderr, "error: singular matrix in matrix3x3_inverse!\n");
	  exit(EXIT_FAILURE);
     }
     detinv = 1.0/detinv;

     minv.c0.x = detinv * (m.c1.y * m.c2.z - m.c2.y * m.c1.z);
     minv.c1.y = detinv * (m.c0.x * m.c2.z - m.c2.x * m.c0.z);
     minv.c2.z = detinv * (m.c1.y * m.c0.x - m.c0.y * m.c1.x);
     
     minv.c0.z = detinv * (m.c0.y * m.c1.z - m.c1.y * m.c0.z);
     minv.c0.y = -detinv * (m.c0.y * m.c2.z - m.c2.y * m.c0.z);
     minv.c1.z = -detinv * (m.c0.x * m.c1.z - m.c1.x * m.c0.z);
     
     minv.c2.x = detinv * (m.c1.x * m.c2.y - m.c1.y * m.c2.x);
     minv.c1.x = -detinv * (m.c1.x * m.c2.z - m.c1.z * m.c2.x);
     minv.c2.y = -detinv * (m.c0.x * m.c2.y - m.c0.y * m.c2.x);

     return minv;
}

int matrix3x3_equal(matrix3x3 m1, matrix3x3 m2)
{
     return (vector3_equal(m1.c0, m2.c0)
	     && vector3_equal(m1.c1, m2.c1)
	     && vector3_equal(m1.c2, m2.c2));
}

vector3 matrix3x3_row1(matrix3x3 m)
{
     vector3 v;
     v.x = m.c0.x;
     v.y = m.c1.x;
     v.z = m.c2.x;
     return v;
}

vector3 matrix3x3_row2(matrix3x3 m)
{
     vector3 v;
     v.x = m.c0.y;
     v.y = m.c1.y;
     v.z = m.c2.y;
     return v;
}

vector3 matrix3x3_row3(matrix3x3 m)
{
     vector3 v;
     v.x = m.c0.z;
     v.y = m.c1.z;
     v.z = m.c2.z;
     return v;
}

/**************************************************************************/

/* complex number utilities */

cnumber make_cnumber(number r, number i)
{
     cnumber c;
     c.re = r; c.im = i;
     return c;
}

cnumber cnumber_conj(cnumber c)
{
     return make_cnumber(c.re, -c.im);
}

int cnumber_equal(cnumber c1, cnumber c2)
{
     return (c1.re == c2.re && c1.im == c2.im);
}

vector3 cvector3_re(cvector3 cv)
{
     vector3 v;
     v.x = cv.x.re; v.y = cv.y.re; v.z = cv.z.re;
     return v;
}

vector3 cvector3_im(cvector3 cv)
{
     vector3 v;
     v.x = cv.x.im; v.y = cv.y.im; v.z = cv.z.im;
     return v;
}

cvector3 make_cvector3(vector3 vr, vector3 vi)
{
     cvector3 cv;
     cv.x = make_cnumber(vr.x, vi.x);
     cv.y = make_cnumber(vr.y, vi.y);
     cv.z = make_cnumber(vr.z, vi.z);
     return cv;
}

int cvector3_equal(cvector3 v1, cvector3 v2)
{
     return (vector3_equal(cvector3_re(v1), cvector3_re(v2)) &&
	     vector3_equal(cvector3_im(v1), cvector3_im(v2)));
}

matrix3x3 cmatrix3x3_re(cmatrix3x3 cm)
{
     matrix3x3 m;
     m.c0 = cvector3_re(cm.c0);
     m.c1 = cvector3_re(cm.c1);
     m.c2 = cvector3_re(cm.c2);
     return m;
}

matrix3x3 cmatrix3x3_im(cmatrix3x3 cm)
{
     matrix3x3 m;
     m.c0 = cvector3_im(cm.c0);
     m.c1 = cvector3_im(cm.c1);
     m.c2 = cvector3_im(cm.c2);
     return m;
}

cmatrix3x3 make_cmatrix3x3(matrix3x3 mr, matrix3x3 mi)
{
     cmatrix3x3 cm;
     cm.c0 = make_cvector3(mr.c0, mi.c0);
     cm.c1 = make_cvector3(mr.c1, mi.c1);
     cm.c2 = make_cvector3(mr.c2, mi.c2);
     return cm;
}

cmatrix3x3 make_hermitian_cmatrix3x3(number m00, number m11, number m22,
				     cnumber m01, cnumber m02, cnumber m12)
{
     cmatrix3x3 cm;
     cm.c0.x = make_cnumber(m00, 0);
     cm.c1.y = make_cnumber(m11, 0);
     cm.c2.z = make_cnumber(m22, 0);
     cm.c1.x = m01; cm.c0.y = cnumber_conj(m01);
     cm.c2.x = m02; cm.c0.z = cnumber_conj(m02);
     cm.c2.y = m12; cm.c1.z = cnumber_conj(m12);
     return cm;
}

int cmatrix3x3_equal(cmatrix3x3 m1, cmatrix3x3 m2)
{
     return (matrix3x3_equal(cmatrix3x3_re(m1), cmatrix3x3_re(m2)) &&
	     matrix3x3_equal(cmatrix3x3_im(m1), cmatrix3x3_im(m2)));
}
