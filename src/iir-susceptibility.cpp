/* Copyright (C) 2005-2019 Massachusetts Institute of Technology.
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

typedef struct {
  size_t sz_data;
  size_t ntot;
  realnum *P;
  realnum *P_prev;
  realnum *W_prev;
  realnum data[1];
} iir_data;

// constructor must take the s domain coefficients and transform
// them to the z domain and store them in the class.
iir_susceptibility::iir_susceptibility(double *num, int num_N, double *den, int den_N, double dt){
    numz_N = num_N;
    denz_N = den_N;
    tustins_method(num,num_N,den,den_N, dt);
    numz = num;
    denz = den;
}