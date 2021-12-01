/***************************************************************/
/* C++ port of meep/examples/pw-source.ctl, using the          */
/* "low-level" meep C++ interface stack, which consists of     */
/* libmeep_geom + libctlgeom + libmeep                         */
/***************************************************************/
/*
; This example creates an approximate TM planewave in vacuum
; propagating at a 45-degree angle, by using a couple of current sources
; with amplitude exp(ikx) corresponding to the desired planewave.
*/

#include <stdio.h>
#include <stdlib.h>
#include <complex>

#include <meep.hpp>

#include "ctl-math.h"
#include "ctlgeom.h"

#include "meepgeom.hpp"

using namespace meep;

/***************************************************************/
/*
; pw-amp is a function that returns the amplitude exp(ik(x+x0)) at a
; given point x.  (We need the x0 because current amplitude functions
; in Meep are defined relative to the center of the current source,
; whereas we want a fixed origin.)  Actually, it is a function of k
; and x0 that returns a function of x ...
(define ((pw-amp k x0) x)
  (exp (* 0+1i (vector3-dot k (vector3+ x x0)))))
*/
/***************************************************************/
typedef struct pw_amp_data {
  vec k;
  vec x0;
} pw_amp_data;

std::complex<double> pw_amp(vec x, void *UserData) {
  pw_amp_data *data = (pw_amp_data *)UserData;
  vec k = data->k;
  vec x0 = data->x0;

  const std::complex<double> II(0.0, 1.0);
  return exp(II * (k & (x + x0)));
}

/***************************************************************/
/* note: meep::fields::add_volume_source needs a more flexible */
/* interface so that the amplitude function can accept user    */
/* data! without this we must have multiple hard-coded         */
/* amplitude functions and global variables.                   */
/***************************************************************/
pw_amp_data pw_amp_data_left;
std::complex<double> pw_amp_left(const vec &x) { return pw_amp(x, (void *)&pw_amp_data_left); }

pw_amp_data pw_amp_data_bottom;
std::complex<double> pw_amp_bottom(const vec &x) { return pw_amp(x, (void *)&pw_amp_data_bottom); }

/***************************************************************/
/* dummy material function needed to pass to structure( )      */
/* constructor as a placeholder before we can call             */
/* set_materials_from_geometry                                 */
/***************************************************************/
double dummy_eps(const vec &) { return 1.0; }

/***************************************************************/
/***************************************************************/
/***************************************************************/
int main(int argc, char **argv) {
  meep::initialize mpi(argc, argv);

  int s = 11;             // size of computational cell, excluding PML
  int dpml = 1;           // thickness of PML layers
  int sxy = s + 2 * dpml; // cell size, including PML
  int resolution = 10;    // pixel spacing

  // (set! geometry-lattice (make lattice (size sxy sxy no-size)))
  // (set! pml-layers (list (make pml (thickness dpml))))
  geometry_lattice.size.x = sxy;
  geometry_lattice.size.y = sxy;
  geometry_lattice.size.z = 0.0;
  grid_volume gv = voltwo(sxy, sxy, resolution);
  gv.center_origin();
  structure the_structure(gv, dummy_eps, pml(dpml));

  geometric_object_list g = {0, 0};
  meep_geom::set_materials_from_geometry(&the_structure, g);

  fields f(&the_structure);

  // (set! sources (list
  //   ; left
  //   (make source
  //     (src (make continuous-src (frequency fcen) (fwidth df)))
  //      (component Ez) (center (* -0.5 s) 0) (size 0 s)
  //      (amp-func (pw-amp k (vector3 (* -0.5 s) 0))))
  //   ; bottom
  //   (make source
  //     (src (make continuous-src (frequency fcen) (fwidth df)))
  //      (component Ez) (center 0 (* -0.5 s)) (size s 0)
  //      (amp-func (pw-amp k (vector3 0 (* -0.5 s)))))
  double fcen = 0.8; // pulse center frequency
  double df = 0.02;  // turn-on bandwidth
  continuous_src_time src(fcen, df);

  vec kdir(1.0, 1.0); // k direction (length is irrelevant)
  vec k = kdir * 2.0 * pi * fcen / abs(kdir);

  vec x0_left(-0.5 * s, 0.0);
  vec size_left(0.0, s);

  vec x0_bottom(0.0, -0.5 * s);
  vec size_bottom(s, 0.0);

  pw_amp_data_left.k = k;
  pw_amp_data_left.x0 = x0_left;
  meep::volume vleft(x0_left, size_left);
  f.add_volume_source(Ez, src, vleft, pw_amp_left);

  pw_amp_data_bottom.k = k;
  pw_amp_data_bottom.x0 = x0_bottom;
  meep::volume vbottom(x0_bottom, size_bottom);
  f.add_volume_source(Ez, src, vbottom, pw_amp_bottom);

  // (run-until T (at-end output-efield-z))
  double T = 400.0;
  while (f.time() < T)
    f.step();

  f.output_hdf5(Ez, f.total_volume());

  // success if we made it here
  return 0;
}
