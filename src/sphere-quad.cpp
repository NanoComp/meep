/* Copyright (C) 2003 Massachusetts Institute of Technology.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define SHIFT3(x,y,z) {double SHIFT3_dummy = z; z = y; y = x; x = SHIFT3_dummy;}

#define CHECK(condition, message) do { \
     if (!(condition))  { \
          fprintf(stderr, "CHECK failure on line %d of " __FILE__ ": " \
                  message "\n", __LINE__);  exit(EXIT_FAILURE); \
     } \
} while (0)

/* Compute quadrature points and weights for integrating on the
   unit sphere.  x, y, z, and weight should be arrays of num_sq_pts
   values to hold the coordinates and weights of the quadrature points
   on output.   Currently, num_sq_pts = 12, 50, and 72 are supported. */
void spherical_quadrature_points(double *x, double *y, double *z,
				 double *weight, int num_sq_pts)
{
     int i,j,k,l, n = 0;
     double x0, y0, z0, w;

     if (num_sq_pts == 50) {
	  /* Computes quadrature points and weights for 50-point 11th degree
	     integration formula on a unit sphere.  This particular quadrature
	     formula has the advantage, for our purposes, of preserving the
	     symmetry group of an octahedraon (i.e. simple cubic symmetry, with
	     respect to the Cartesian xyz axes).  References:
	
	     A. D. McLaren, "Optimal Numerical Integration on a Sphere,"
	     Math. Comp. 17, pp. 361-383 (1963).
	     
	     Also in: Arthur H. Stroud, "Approximate Calculation of Multiple
	     Integrals" (Prentice Hall, 1971) (formula number U3:11-1). 

	     This code was written with the help of example code by
	     John Burkardt: 
	           http://www.psc.edu/~burkardt/src/stroud/stroud.html */
     
	  x0 = 1; y0 = z0 = 0;
	  w = 9216 / 725760.0;
	  for (i = 0; i < 2; ++i) {
	       x0 = -x0;
	       for (j = 0; j < 3; ++j) {
		    SHIFT3(x0,y0,z0);
		    x[n] = x0; y[n] = y0; z[n] = z0; weight[n++] = w;
	       }
	  }
	  
	  x0 = y0 = sqrt(0.5); z0 = 0;
	  w = 16384 / 725760.0;
	  for (i = 0; i < 2; ++i) {
	       x0 = -x0;
	       for (j = 0; j < 2; ++j) {
		    y0 = -y0;
		    for (k = 0; k < 3; ++k) {
			 SHIFT3(x0,y0,z0);
			 x[n] = x0; y[n] = y0; z[n] = z0; weight[n++] = w;
		    }
	       }
	  }
	  
	  x0 = y0 = z0 = sqrt(1.0 / 3.0);
	  w = 15309 / 725760.0;
	  for (i = 0; i < 2; ++i) {
	       x0 = -x0;
	       for (j = 0; j < 2; ++j) {
		    y0 = -y0;
		    for (k = 0; k < 2; ++k) {
			 z0 = -z0;
			 x[n] = x0; y[n] = y0; z[n] = z0; weight[n++] = w;
		    }
	       }
	  }
	  
	  x0 = y0 = sqrt(1.0 / 11.0); z0 = 3 * x0;
	  w = 14641 / 725760.0;
	  for (i = 0; i < 2; ++i) {
	       x0 = -x0;
	       for (j = 0; j < 2; ++j) {
		    y0 = -y0;
		    for (k = 0; k < 2; ++k) {
			 z0 = -z0;
			 for (l = 0; l < 3; ++l) {
			      SHIFT3(x0,y0,z0);
			      x[n] = x0; y[n] = y0; z[n] = z0; weight[n++] = w;
			 }
		    }
	       }
	  }
     }  
     else if (num_sq_pts == 72 || num_sq_pts == 12) {
	  /* As above (same references), but with a 72-point 14th degree
	     formula, this time with the symmetry group of an icosohedron.
	     (Stroud formula number U3:14-1.)  Alternatively, just use
	     the 12-point 5th degree formula consisting of the vertices
	     of a regular icosohedron. */
	  
	  /* first, the vertices of an icosohedron: */
	  x0 = sqrt(0.5 - sqrt(0.05));
	  y0 = sqrt(0.5 + sqrt(0.05));
	  z0 = 0;
	  if (num_sq_pts == 72)
	       w = 125 / 10080.0;
	  else
	       w = 1 / 12.0;
	  for (i = 0; i < 2; ++i) {
	       x0 = -x0;
	       for (j = 0; j < 2; ++j) {
		    y0 = -y0;
		    for (k = 0; k < 3; ++k) {
			 SHIFT3(x0,y0,z0);
			 x[n] = x0; y[n] = y0; z[n] = z0; weight[n++] = w;
		    }
	       }
	  }
	  
	  if (num_sq_pts == 72) {
	       /* it would be nice, for completeness, to have more
                  digits here: */
	       double coords[3][5] = {
		    { -0.151108275, 0.315838353, 0.346307112, -0.101808787, -0.409228403 },
		    { 0.155240600, 0.257049387, 0.666277790,  0.817386065, 0.501547712 },
		    { 0.976251323, 0.913330032, 0.660412970,  0.567022920, 0.762221757 }
	       };
	       
	       w = 143 / 10080.0;
	       for (l = 0; l < 5; ++l) {
		    x0 = coords[0][l]; y0 = coords[1][l]; z0 = coords[2][l];
		    for (i = 0; i < 3; ++i) {
			 double dummy = x0;
			 x0 = z0;
			 z0 = -y0;
			 y0 = -dummy;
			 for (j = 0; j < 3; ++j) {
			      SHIFT3(x0,y0,z0);
			      x[n] = x0; y[n] = y0; z[n] = z0; weight[n++] = w;
			 }
			 y0 = -y0;
			 z0 = -z0;
			 x[n] = x0; y[n] = y0; z[n] = z0; weight[n++] = w;
		    }
	       }
	  }
     }
     else
	  CHECK(0, "spherical_quadrature_points: passed unknown # points!");

     CHECK(n == num_sq_pts,
	   "bug in spherical_quadrature_points: wrong number of points!");
}

#define NQUAD3 50 /* use 50-point quadrature formula by default */

/**********************************************************************/
#define K_PI 3.141592653589793238462643383279502884197
#define NQUAD2 12
/**********************************************************************/

int main(void)
{
     int i;
     double x[NQUAD3], y[NQUAD3], z[NQUAD3], weight[NQUAD3];

     printf(
"/* For 1d, 2d, and 3d, quadrature points and weights on a unit sphere.\n"
"   There are num_sphere_quad[dim-1] points i, with the i-th point at\n"
"   (x,y,z) = (sphere_quad[dim-1][i][ 0, 1, 2 ]), and with a quadrature\n"
"   weight sphere_quad[dim-1][i][3]. */\n\n");


     printf("static const int num_sphere_quad[3] = { %d, %d, %d };\n\n",
	    2, NQUAD2, NQUAD3);

     printf("static const double sphere_quad[3][%d][4] = {\n", NQUAD3);

     printf("    { {0,0,1,0.5}, {0,0,-1,0.5} },\n");

     printf("    {\n");
     for (i = 0; i < NQUAD2; ++i) {
	  printf("        { %0.20g, %0.20g, %0.20g, %0.20g },\n",
		 cos(2*i * K_PI / NQUAD2),
		 sin(2*i * K_PI / NQUAD2),
		 0.0,
		 1.0 / NQUAD2);
     }
     printf("    },\n");

     printf("    {\n");
     spherical_quadrature_points(x,y,z, weight, NQUAD3);
     for (i = 0; i < NQUAD3; ++i) {
	  printf("        { %0.20g, %0.20g, %0.20g, %0.20g },\n",
		 x[i], y[i], z[i], weight[i]);
     }
     printf("    }\n");

     printf("};\n");
}
