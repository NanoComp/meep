/* Copyright (C) 2005-2015 Massachusetts Institute of Technology
%
%  This program is free software; you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation; either version 2, or (at your option)
%  any later version.
%
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%
%  You should have received a copy of the GNU General Public License
%  along with this program; if not, write to the Free Software Foundation,
%  Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
*/

/* given an arbitrary field configuration F, compute the coefficients
   in an expansion of F as a linear combination of normal modes
   of the geometry as computed by mpb.
*/


#include <stdio.h>
#include <string.h>
#include <math.h>

#include "meep_internals.hpp"

using namespace std;

typedef complex<double> cdouble;

namespace meep {

// prototype for optional user-supplied function to provide an
// initial estimate of the wavevector of band #band at frequency freq
typedef vec (*kpoint_func)(void *user_data, double freq, int band);

typedef struct {

} mode_projection_data;

#define UNUSED(x) (void) x // silence compiler warnings

/***************************************************************/
/* callback function passed to loop_in_chunks to evaluate the  */
/* projection of the user-specified fields onto one or more    */
/* normal modes as computed by mpb                             */
/***************************************************************/
static void mode_projection_chunkloop(fields_chunk *fc, 
                                      int ichnk, component cgrid,
				      ivec is, ivec ie,
				      vec s0, vec s1, vec e0, vec e1,
				      double dV0, double dV1,
				      ivec shift, complex<double> shift_phase,
				      const symmetry &S, int sn,
				      void *data_)
{
  UNUSED(ichnk);UNUSED(cgrid);UNUSED(s0);UNUSED(s1);UNUSED(e0);UNUSED(e1);
  UNUSED(dV0);UNUSED(dV1);UNUSED(shift_phase); UNUSED(fc);
  mode_projection_data *data = (mode_projection_data *) data_;

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
std::vec<cdouble> fields::get_mode_coefficients(dft_flux flux,
                                                direction d,
                                                const volume where,
                                                std::vec<int> bands,
                                                kpoint_func user_func,
                                                void *user_data)
{
  am_now_working_on(ModeExpansion);

  // create output array
  std::vec<cdouble> coefficients(bands.size());

  // some default inputs for add_eigenmode_source
  component DefaultComponent = Dielectric;
  vec kpoint(0,0,0);
  int parity = 0; /* NO_PARITY */
  bool match_frequency = true
  double eig_resolution = a;
  double eigensolver_tol = 1.0e-7;
  cdouble amp=1.0;

  // loop over all frequencies in the dft_flux
  for(int nfreq=0; nfreq<flux.NFreq; nfreq++)
   { 
     double freq = flux.freq_min + ((double)nfreq)*flux.dfreq;
     continuous_src_time src(freq);

     // loop over caller's list of requested bands
     for(int nband=0; nband<bands.size(); nband++)
      { 
        int band = bands[nband];

        // query user's function (if present) for initial k-point guess
        if (user_func)
         kpoint = user_func(user_data, freq, nband);
       
        // add source
        add_eigenmode_source(DefaultComponent, src, d, where, where,
                             band, kpoint, match_frequency, parity,
                             eig_resolution, eigensolver_tol, amp);
       
        // call loop_in_chunks to evaluate projection and normalization 
        // integrals
      };
   };

}

} // namespace meep
