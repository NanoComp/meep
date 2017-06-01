/***************************************************************/
/* C++ port of meep/examples/pw-source.ctl                     */
/***************************************************************/
/*
; This example creates an approximate TM planewave in vacuum
; propagating at a 45-degree angle, by using a couple of current sources
; with amplitude exp(ikx) corresponding to the desired planewave.
*/

#include <stdio.h>
#include <stdlib.h>
#include <complex>
#include "mpSession.hpp"

using namespace meep;
using namespace meepSession;

typedef std::complex<double> cdouble;

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
typedef struct pw_amp_data
 { vector3 k;
   vector3 x0;
 } pw_amp_data;

cdouble pw_amp(vector3 x, void *UserData)
{
  pw_amp_data *data = (pw_amp_data *)UserData;
  vector3 k  = data->k;
  vector3 x0 = data->x0;

  const cdouble II(0.0,1.0);
  return exp( II * vector3_dot(k, vector3_plus(x, x0 ) ) );
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
int main(int argv, char *argc[])
{
  mpSession S;

  int s=11;           // size of computational cell, excluding PML
  int dpml=1;         // thickness of PML layers
  int sxy = s+2*dpml; // cell size, including PML

  //(set! geometry-lattice (make lattice (size sxy sxy no-size)))
  S.set_geometry_lattice(sxy,sxy);

  // (set-param! resolution 10)
  S.set_resolution(10);

  // (set! pml-layers (list (make pml (thickness dpml))))
  S.add_pml_layer(dpml);
  
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
  double fcen  = 0.8;  // pulse center frequency
  double df    = 0.02; // turn-on bandwidth
  vector3 kdir =       // k direction (length is irrelevant)
   {1.0, 1.0, 0.0};
  vector3 k =          // k with correct length
   vector3_scale(2.0*pi*fcen,unit_vector3(kdir));

  vector3 size = vector3_scale(s, v3_yaxis);

  vector3 x0Left = vector3_scale(-0.5*s, v3_xaxis);
  pw_amp_data leftData = { .k=k, .x0 = x0Left };
  S.add_continuous_src(fcen, df, component::Ez, x0Left, size,
                       pw-amp, (void *)&leftData);

  //vector3 x0Bottom = vector3_scale(-0.5*s, v3_yaxis);
 // pw_amp_data bottomData = { .k=k, .x0 = x0Bottom };
 // S.add_continuous_src(fcen, df, component::Ez, x0Bottom, size,
 //                      pw-amp, (void *)&bottomData);

  // (run-until T (at-end output-efield-z))
  double T=400.0;
 // S.add_output(component::Ez);
 // S.run_until(T);
}
