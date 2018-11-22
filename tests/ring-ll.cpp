/***************************************************************/
/* C++ port of meep/examples/ring.ctl, using the               */
/* "low-level" meep C++ interface stack, which consists of     */
/* libmeep_geom + libctlgeom + libmeep                         */
/***************************************************************/
/*
; Calculating 2d ring-resonator modes, from the Meep tutorial.
*/

#include <stdio.h>
#include <stdlib.h>
#include <complex>
#include <vector>

#include "meep.hpp"

#include "ctl-math.h"
#include "ctlgeom.h"

#include "meepgeom.hpp"

using namespace meep;

typedef std::complex<double> cdouble;

vector3 v3(double x, double y=0.0, double z=0.0)
{
  vector3 v;
  v.x=x; v.y=y; v.z=z;
  return v;
}

/***************************************************************/
/* dummy material function needed to pass to structure( )      */
/* constructor as a placeholder before we can call             */
/* set_materials_from_geometry                                 */
/***************************************************************/
double dummy_eps(const vec &) { return 1.0; }

/***************************************************************/
/***************************************************************/
/***************************************************************/
int main(int argc, char *argv[])
{
  initialize mpi(argc, argv);

  double n=3.4;     // index of waveguide
  double w=1.0;     // width of waveguide
  double r=1.0;     // inner radius of ring

  double pad=4;     // padding between waveguide and edge of PML
  double dpml=2;    // thickness of PML

  double sxy        = 2.0*(r+w+pad+dpml);  // cell size
  double resolution = 10.0;

  // (set-param! resolution 10)
  // (set! geometry-lattice (make lattice (size sxy sxy no-size)))
  geometry_lattice.size.x=sxy;
  geometry_lattice.size.y=sxy;
  geometry_lattice.size.z=0.0;
  grid_volume gv = voltwo(sxy, sxy, resolution);
  gv.center_origin();

  // (set! symmetries (list (make mirror-sym (direction Y))))
  //symmetry sym=mirror(Y, gv);
  symmetry sym=identity();

  // (set! pml-layers (list (make pml (thickness dpml))))
  // ; exploit the mirror symmetry in structure+source:
  structure the_structure(gv, dummy_eps, pml(dpml), sym);

  // ; Create a ring waveguide by two overlapping cylinders - later objects
  // ; take precedence over earlier objects, so we put the outer cylinder first.
  // ; and the inner (air) cylinder second.
  // (set! geometry (list
  //	(make cylinder (center 0 0) (height infinity)
  //		(radius (+ r w)) (material (make dielectric (index n))))
  // 	(make cylinder (center 0 0) (height infinity)
  // 		(radius r) (material air))))
  meep_geom::material_type dielectric = meep_geom::make_dielectric(n*n);
  geometric_object objects[2];
  vector3 v3zero = {0.0, 0.0, 0.0};
  vector3 zaxis  = {0.0, 0.0, 1.0};
  objects[0] = make_cylinder(dielectric, v3zero, r+w, ENORMOUS, zaxis);
  objects[1] = make_cylinder(meep_geom::vacuum,  v3zero, r, ENORMOUS, zaxis);
  geometric_object_list g={ 2, objects };
  meep_geom::set_materials_from_geometry(&the_structure, g);
  fields f(&the_structure);

  // ; If we don't want to excite a specific mode symmetry, we can just
  // ; put a single point source at some arbitrary place, pointing in some
  // ; arbitrary direction.  We will only look for TM modes (E out of the plane).
  // (set! sources (list
  //              (make source
  //                (src (make gaussian-src (frequency fcen) (fwidth df)))
  //                (component Ez) (center (+ r 0.1) 0))))
  double fcen = 0.15;  // ; pulse center frequency
  double df   = 0.1;   // ; df
  gaussian_src_time src(fcen, df);
  volume v( vec(r+0.1,0.0), vec(0.0, 0.0) );
  f.add_volume_source(Ez, src, v);

  // (run-sources+ 300
  // 	(at-beginning output-epsilon)
  // 	(after-sources (harminv Ez (vector3 (+ r 0.1)) fcen df)))
  while( f.round_time() < f.last_source_time() )
   f.step();

  double T = 300.0;
  double stop_time = f.round_time() + T;
  std::vector<cdouble> fieldData;
  vec eval_pt(r+0.1,0.0);
  while( f.round_time() < stop_time )
   { f.step();
     fieldData.push_back( f.get_field(Ez, eval_pt) );
   };

  #define MAXBANDS 100
  cdouble amp[MAXBANDS];
  double freq_re[MAXBANDS];
  double freq_im[MAXBANDS];
  double err[MAXBANDS];
  master_printf("starting do_harminv...\n");
  int bands=do_harminv( &fieldData[0], fieldData.size(), f.dt,
                        fcen-0.5*df, fcen+0.5*df, MAXBANDS,
                        amp, freq_re, freq_im, err);
  master_printf("harminv0: | real(freq) | imag(freq)  |     Q      |  abs(amp)  |          amp          | err\n");
  master_printf("--------------------------------------------------------------------------------------------\n");
  for(int nb=0; nb<bands; nb++)
   master_printf("harminv0: | %.4e | %+.4e | %.4e | %.4e | {%+.2e,%+.2e} | %.1e \n",
                  freq_re[nb], freq_im[nb], 0.5*freq_re[nb]/freq_im[nb],
                  abs(amp[nb]), real(amp[nb]), imag(amp[nb]), err[nb]);

  // test comparison with expected values
  int ref_bands=3;
  double ref_freq_re[3] = {  1.1807e-01, 1.4716e-01,  1.7525e-01  };
  double ref_freq_im[3] = { -7.6133e-04, -2.1156e-04, -5.2215e-05 };
  cdouble ref_amp[3]    = { cdouble(-8.28e-04,-1.34e-03),
    cdouble(1.23e-03,-1.25e-02), cdouble(2.83e-03,-6.52e-04) };
  if (bands!=3)
   abort("harminv found only %i/%i bands\n",bands,ref_bands);
  for(int nb=0; nb<bands; nb++)
   if (    fabs(freq_re[nb]-ref_freq_re[nb]) > 1.0e-2*fabs(ref_freq_re[nb])
        || fabs(freq_im[nb]-ref_freq_im[nb]) > 1.0e-2*fabs(ref_freq_im[nb])
        ||  abs(    amp[nb]-    ref_amp[nb]) > 1.0e-2* abs(    ref_amp[nb])

      )
    abort("harminv band %i disagrees with ref: {re f, im f, re A, im A}={%e,%e,%e,%e}!= {%e,%e,%e,%e}\n",
           nb, freq_re[nb],     freq_im[nb],     real(amp[nb]),     imag(amp[nb]),
               ref_freq_re[nb], ref_freq_im[nb], real(ref_amp[nb]), imag(ref_amp[nb]));

  master_printf("all harminv results match reference values\n");

  // ; Output fields for one period at the end.  (If we output
  // ; at a single time, we might accidentally catch the Ez field
  // ; when it is almost zero and get a distorted view.)
  // (run-until (/ 1 fcen) (at-every (/ 1 fcen 20) output-efield-z))
  double DeltaT = 1.0/(20*fcen);
  double NextOutputTime = f.round_time() + DeltaT;
  while( f.round_time() < 1.0/fcen )
   { f.step();
     if ( f.round_time() >= NextOutputTime  )
      { f.output_hdf5(Ez, f.total_volume());
        NextOutputTime += DeltaT;
      };
   };

  // this seems to be necessary to prevent failures
  all_wait();

  // success if we made it here
  return 0;

}
