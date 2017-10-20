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

/***************************************************************/
/* precisely one of fun, rfun should be non-NULL               */
/***************************************************************/
void *fields::do_get_array_slice(const volume &where,
                                 std::vector<component> components,
                                 field_function fun,
                                 field_rfunction rfun,
                                 void *fun_data,
                                 void *vslice) {

  am_now_working_on(FieldOutput);

  /***************************************************************/
  /* call get_array_slice_dimensions to get slice dimensions and */
  /* partially initialze an array_slice_data struct              */
  /***************************************************************/
  int dims[3];
  array_slice_data data;
  int rank=get_array_slice_dimensions(where, dims, &data);
  int slice_size=data.slice_size;
  if (rank==0 || slice_size==0) return 0; // no data to write

  bool complex_data = (rfun==0);
  cdouble *zslice;
  double *slice;
  if (vslice==0)
   { if (complex_data)
      { zslice = new cdouble[slice_size];
        vslice = (void *)zslice;
      } 
     else
      { slice  = new double[slice_size];
        vslice = (void *)slice;
      };
   };
   
  data.vslice     = vslice;
  data.fun        = fun;
  data.rfun       = rfun;
  data.fun_data   = fun_data;
  data.components = components;

  int num_components = components.size();

  data.cS      = new component[num_components];
  data.ph      = new cdouble[num_components];
  data.fields  = new cdouble[num_components];
  
  data.offsets = new int[2 * num_components];
  for (int i = 0; i < 2 * num_components; ++i)
    data.offsets[i] = 0;

  /* compute inverse-epsilon directions for computing Dielectric fields */
  data.ninveps = 0;
  bool needs_dielectric = false;
  for (int i = 0; i < num_components; ++i)
    if (components[i] == Dielectric) { needs_dielectric = true; break; }
  if (needs_dielectric) 
    FOR_ELECTRIC_COMPONENTS(c) if (gv.has_field(c)) {
      if (data.ninveps == 3) abort("more than 3 field components??");
      data.inveps_cs[data.ninveps] = c;
      data.inveps_ds[data.ninveps] = component_direction(c);
      ++data.ninveps;
    }
  
  /* compute inverse-mu directions for computing Permeability fields */
  data.ninvmu = 0;
  bool needs_permeability = false;
  for (int i = 0; i < num_components; ++i)
    if (components[i] == Permeability) { needs_permeability = true; break; }
  if (needs_permeability) 
    FOR_MAGNETIC_COMPONENTS(c) if (gv.has_field(c)) {
      if (data.ninvmu == 3) abort("more than 3 field components??");
      data.invmu_cs[data.ninvmu] = c;
      data.invmu_ds[data.ninvmu] = component_direction(c);
      ++data.ninvmu;
    }
 
  loop_in_chunks(get_array_slice_chunkloop, (void *) &data,
		 where, Centered, true, true);

  /***************************************************************/
  /* repeatedly call sum_to_all to consolidate full array slice  */
  /* on all cores                                                */
  /***************************************************************/
#define BUFSIZE 1<<16 // use 64k buffer
  if (complex_data)
   { cdouble *buffer = new cdouble[BUFSIZE];
     cdouble *slice = (cdouble *)vslice;
     int offset=0, remaining=slice_size;
     while(remaining!=0)
      { 
        int size = (remaining > BUFSIZE ? BUFSIZE : remaining);
        sum_to_all(slice + offset, buffer, size);
        memcpy(slice+offset, buffer, size*sizeof(cdouble));
        remaining-=size;
        offset+=size;
      };
     delete[] buffer;
   }
  else
   { double *buffer = new double[BUFSIZE];
     double *slice = (double *)vslice;
     int offset=0, remaining=slice_size;
     while(remaining!=0)
      { int size = (remaining > BUFSIZE ? BUFSIZE : remaining);
        sum_to_all(slice + offset, buffer, size);
        memcpy(slice+offset, buffer, size*sizeof(double));
        remaining-=size;
        offset+=size;
      };
     delete[] buffer;
   };

  delete[] data.offsets;
  delete[] data.fields;
  delete[] data.ph;
  delete[] data.cS;
  finished_working();

  return vslice;
}

/***************************************************************/
/* entry points to get_array_slice                             */
/***************************************************************/
double *fields::get_array_slice(const volume &where,
                                std::vector<component> components,
                                field_rfunction rfun, void *fun_data,
                                double *slice)
{
  return (double *)do_get_array_slice(where, components,
                                      0, rfun, fun_data,
                                      (void *)slice);
} 

cdouble *fields::get_complex_array_slice(const volume &where,
                                         std::vector<component> components,
                                         field_function fun, void *fun_data,
                                         cdouble *slice)
{
  return (cdouble *)do_get_array_slice(where, components,
                                       fun, 0, fun_data,
                                       (void *)slice);
}

double *fields::get_array_slice(const volume &where, component c,
                                double *slice, int slice_length)
{
  (void) slice_length;
  std::vector<component> components(1);
  components[0]=c;
  return (double *)do_get_array_slice(where, components,
                                      0, default_field_rfunc, 0,
                                      (void *)slice);
}

double *fields::get_array_slice(const volume &where,
                                derived_component c,
                                double *slice, int slice_length)
{
  (void) slice_length;
  int nfields;
  component carray[12];
  field_rfunction rfun = derived_component_func(c, gv, nfields, carray);
  std::vector<component> cs(carray, carray+nfields);
  return (double *)do_get_array_slice(where, cs,
                                      0, rfun, &nfields,
                                      (void *)slice);
}

cdouble *fields::get_complex_array_slice(const volume &where, component c,
                                         cdouble *slice, int slice_length)
{
  (void) slice_length;
  std::vector<component> components(1);
  components[0]=c;
  return (cdouble *)do_get_array_slice(where, components,
                                       default_field_func, 0, 0,
                                       (void *)slice);
}

} // namespace meep
