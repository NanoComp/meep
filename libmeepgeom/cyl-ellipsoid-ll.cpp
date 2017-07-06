/***************************************************************/
/* simple test for libmeepgeom, modeled after meep_test.ctl    */
/***************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <complex>

#include "meep.hpp"

#include "ctl-math.h"
#include "ctlgeom.h"

#include "meepgeom.hpp"

using namespace meep;

typedef std::complex<double> cdouble;

/***************************************************************/
/* dummy material function needed to pass to structure( )      */
/* constructor as a placeholder before we can call             */
/* set_materials_from_geometry                                 */
/***************************************************************/
double dummy_eps(const vec &) { return 1.0; }

/***************************************************************/
/* usage: cyl-ellipsoid [ --polarization xx ]                  */
/*  where xx = S for TE polarization (default)                 */
/*             P for TM polarization                           */
/***************************************************************/
int main(int argc, char *argv[])
{
  initialize mpi(argc, argv);

  // simple argument parsing 
  meep::component src_cmpt=Ez;
  for(int narg=1; narg<argc; narg++)
   if ( argv[narg] && !strcmp(argv[narg],"--polarization") )
    { if (narg+1 == argc) 
       meep::abort("no option specified for --polarization");
      else if (!strcasecmp(argv[narg+1],"S"))
       printf("Using S-polarization\n");
      else if (!strcasecmp(argv[narg+1],"P"))
       { src_cmpt=Hz;
         printf("Using P-polarization\n");
       }
      else 
       meep::abort("invalid --polarization %s",argv[narg+1]);
    };

  //(set-param! resolution 100)
  double resolution = 100.0; 

  // (set! geometry-lattice (make lattice (size 10 10 no-size)))
  // (set! pml-layers (list (make pml (thickness 1))))
  // (if (= src-cmpt Ez)
  //  (set! symmetries (list (make mirror-sym (direction X))
  //                         (make mirror-sym (direction Y)))))
  // (if (= src-cmpt Hz) 
  //  (set! symmetries (list (make mirror-sym (direction X) (phase -1))
  //  (set! symmetries (list (make mirror-sym (direction Y) (phase -1))
  grid_volume gv = voltwo(10.0, 10.0, resolution);
  gv.center_origin();
  symmetry sym = (src_cmpt==Ez) ?  mirror(X,gv) + mirror(Y,gv)
                                : -mirror(X,gv) - mirror(Y,gv);
  structure the_structure(gv, dummy_eps, pml(1.0), sym);

  // (set! geometry (list 
  //    (make cylinder (center 0 0 0) (radius 3) (height infinity)
  //                   (material (make medium (index 3.5))))
  //    (make ellipsoid (center 0 0 0) (size 1 2 infinity)
  //                   (material air))))
  double n=3.5; // index of refraction
  material_type dielectric = meep_geom::make_dielectric(n*n);
  geometric_object objects[2];
  vector3 center = {0.0, 0.0, 0.0};
  double radius  = 3.0;
  double height  = HUGE_VAL;
  vector3 xhat   = {1.0, 0.0, 0.0};
  vector3 yhat   = {0.0, 1.0, 0.0};
  vector3 zhat   = {0.0, 0.0, 1.0};
  //vector3 size   = {1.0, 2.0, HUGE_VAL};
  vector3 size   = {1.0, 2.0, 2.0/1.0e20};
  objects[0] = make_cylinder(dielectric, center, radius, height, zhat);
  //objects[1] = make_ellipsoid(meep_geom::vacuum, center, xhat, yhat, zhat, size);
  //geometric_object_list g={ 2, objects };
geometric_object_list g={ 1, objects };
  meep_geom::set_materials_from_geometry(&the_structure, g);

  // (set! sources (list (make source (src (make gaussian-src (frequency 1) (fwidth 0.1)))
                     //  (center 0 0 0) (component src-cmpt))))
  fields f(&the_structure);
  double fcen = 1.0;
  double df   = 0.1;
  gaussian_src_time src(fcen, df);
  vec src_point = vec(0.0, 0.0);
 // vec src_size  = vec(10.0, 10.0);
  vec src_size  = vec(2.0/1e20, 2.0/1e20);
  f.add_volume_source(src_cmpt, src, volume(src_point, src_size));

  //(define print-stuff (lambda () (print "field:, " (get-field-point src-cmpt (vector3 4.13 3.75 0)) "\n")))
  //
  //(run-until 23 (at-beginning output-epsilon)
  //	      (at-end print-stuff))
  f.output_hdf5(Dielectric, f.total_volume());
  double stop_time=23.0;
  while( f.round_time() < stop_time)
   f.step();

  meep::vec eval_pt=vec(4.13, 3.75);
  std::complex<double> out_field=f.get_field(src_cmpt, eval_pt);
  printf("field: %e + i%e = %e\n",real(out_field), imag(out_field), abs(out_field));

}
