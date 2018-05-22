/***************************************************************/
/* C++ port of meep/examples/bend-flux.ctl, using the          */
/* "low-level" meep C++ interface stack, which consists of     */
/* libmeep_geom + libctlgeom + libmeep                         */
/***************************************************************/

/*
; From the Meep tutorial: transmission around a 90-degree waveguide
; bend in 2d.
*/

#include <stdio.h>
#include <stdlib.h>
#include <complex>

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
void bend_flux(bool no_bend, bool use_prisms)
{
  double sx=16.0;       // size of cell in X direction
  double sy=32.0;       // size of cell in Y direction
  double pad=4.0;       // padding distance between waveguide and cell edge
  double w=1.0;         // width of waveguide
  double resolution=10; // (set-param! resolution 10)

  // (set! geometry-lattice (make lattice (size sx sy no-size)))
  // (set! pml-layers (list (make pml (thickness 1.0))))
  geometry_lattice.size.x=16.0;
  geometry_lattice.size.y=32.0;
  geometry_lattice.size.z=0.0;
  grid_volume gv = voltwo(sx, sy, resolution);
  gv.center_origin();
  structure the_structure(gv, dummy_eps, pml(1.0));

  double wvg_ycen = -0.5*(sy - w - 2.0*pad); //y center of horiz. wvg
  double wvg_xcen =  0.5*(sx - w - 2.0*pad); //x center of vert. wvg

   //  (set! geometry
   //      (if no-bend?
   //	  (list
   //	  (list
   //	   (make block
   //	     (center wvg-xcen (* 0.5 pad))
   //	     (size w (- sy pad) infinity)
   //	     (material (make dielectric (epsilon 12)))))))
   //	   (make block
   //	     (center 0 wvg-ycen)
   //	     (size infinity w infinity)
   //	     (material (make dielectric (epsilon 12)))))
   //	   (make block
   //	     (center (* -0.5 pad) wvg-ycen)
   //	     (size (- sx pad) w infinity)
   //	     (material (make dielectric (epsilon 12))))
   vector3 e1=v3(1.0, 0.0, 0.0);
   vector3 e2=v3(0.0, 1.0, 0.0);
   vector3 e3=v3(0.0, 0.0, 1.0);

   meep_geom::material_type dielectric = meep_geom::make_dielectric(12.0);
   if (no_bend)
    {
      if (use_prisms)
       { vector3 vertices[4];
         vertices[0]=v3(-0.5*sx, wvg_ycen - 0.5*w, 0.0);
         vertices[1]=v3(+0.5*sx, wvg_ycen - 0.5*w, 0.0);
         vertices[2]=v3(+0.5*sx, wvg_ycen + 0.5*w, 0.0);
         vertices[3]=v3(-0.5*sx, wvg_ycen + 0.5*w, 0.0);
         double height=0.0;
         vector3 axis=v3(0.0,0.0,1.0);
         geometric_object objects[1];
         objects[0] = make_prism( dielectric, vertices, 4, height, axis);
         geometric_object_list g={ 1, objects };
         meep_geom::set_materials_from_geometry(&the_structure, g);
       }
      else 
       {
         vector3 center = v3(0.0, wvg_ycen, 0.0);
         vector3 size   = v3(sx,  w,        0.0);
         geometric_object objects[1];
         objects[0] = make_block( dielectric, center, e1, e2, e3, size );
         geometric_object_list g={ 1, objects };
         meep_geom::set_materials_from_geometry(&the_structure, g);
       }
    }
   else
    {
      if (use_prisms)
       { vector3 vertices[6];
         vertices[0]=v3(        -0.5*sx, wvg_ycen -0.5*w);
         vertices[1]=v3(wvg_xcen+0.5*w,  wvg_ycen -0.5*w);
         vertices[2]=v3(wvg_xcen+0.5*w,           +0.5*sy);
         vertices[3]=v3(wvg_xcen-0.5*w,           +0.5*sy);
         vertices[4]=v3(wvg_xcen-0.5*w,  wvg_ycen +0.5*w);
         vertices[5]=v3(        -0.5*sx, wvg_ycen +0.5*w);
         double height=0.0;
         vector3 axis=v3(0.0,0.0,1.0);
         geometric_object objects[1];
         objects[0] = make_prism( dielectric, vertices, 6, height, axis);
         geometric_object_list g={ 1, objects };
         meep_geom::set_materials_from_geometry(&the_structure, g);
       }
      else 
       { 
         geometric_object objects[2];
         vector3 hcenter = v3(-0.5*pad, wvg_ycen), hsize = v3(sx-pad, w);
         vector3 vcenter = v3(wvg_xcen, 0.5*pad),  vsize = v3(w, sy-pad);
         objects[0] = make_block(dielectric, hcenter, e1, e2, e3, hsize);
         objects[1] = make_block(dielectric, vcenter, e1, e2, e3, vsize);
         geometric_object_list g={ 2, objects };
         meep_geom::set_materials_from_geometry(&the_structure, g);
       }
    }

  fields f(&the_structure);

  //  (set! sources (list
  //	       (make source
  //		 (src (make gaussian-src (frequency fcen) (fwidth df)))
  //		 (component Ez)
  //		 (center (+ 1 (* -0.5 sx)) wvg-ycen)
  //		 (size 0 w))))
  //
  double fcen = 0.15;  // ; pulse center frequency
  double df   = 0.1;   // ; df
  gaussian_src_time src(fcen, df);
  volume v( vec(1.0-0.5*sx, wvg_ycen), vec(0.0,w) );
  f.add_volume_source(Ez, src, v);

  //(define trans ; transmitted flux
  //      (add-flux fcen df nfreq
  //		(if no-bend?
  //		    (make flux-region
  //		      (center (- (/ sx 2) 1.5) wvg-ycen) (size 0 (* w 2)))
  //		    (make flux-region
  //		      (center wvg-xcen (- (/ sy 2) 1.5)) (size (* w 2) 0)))))
  double f_start = fcen-0.5*df;
  double f_end   = fcen+0.5*df;
  int nfreq      = 100; //  number of frequencies at which to compute flux

  volume *trans_volume=
   no_bend ? new volume(vec(0.5*sx-1.5, wvg_ycen), vec(0.0, 2.0*w))
           : new volume(vec(wvg_xcen, 0.5*sy-1.5), vec(2.0*w, 0.0));
  volume_list trans_vl = volume_list(*trans_volume, Sz);
  dft_flux trans=f.add_dft_flux(&trans_vl, f_start, f_end, nfreq);

  //(define refl ; reflected flux
  //      (add-flux fcen df nfreq
  //		(make flux-region
  //		  (center (+ (* -0.5 sx) 1.5) wvg-ycen) (size 0 (* w 2)))))
  //
  volume refl_volume( vec(-0.5*sx+1.5, wvg_ycen), vec(0.0,2.0*w));
  volume_list refl_vl= volume_list(refl_volume, Sz);
  dft_flux refl=f.add_dft_flux(&refl_vl, f_start, f_end, nfreq);

  //; for normal run, load negated fields to subtract incident from refl. fields
  //(if (not no-bend?) (load-minus-flux "refl-flux" refl))
  const char *dataname="refl-flux";
  if (!no_bend)
   { refl.load_hdf5(f, dataname);
     refl.scale_dfts(-1.0);
   };

  //(run-sources+
  // (stop-when-fields-decayed 50 Ez
  //			   (if no-bend?
  //			       (vector3 (- (/ sx 2) 1.5) wvg-ycen)
  //			       (vector3 wvg-xcen (- (/ sy 2) 1.5)))
  //			   1e-3)
  // (at-beginning output-epsilon))
  char *prefix = const_cast<char *>(no_bend ? "straight" : "bent");
  f.output_hdf5(Dielectric, f.total_volume(), 0, false, true, prefix);
  vec eval_point = no_bend ? vec(0.5*sx-1.5, wvg_ycen)
                           : vec(wvg_xcen, 0.5*sy - 1.5);
  double DeltaT=50.0, NextCheckTime = f.round_time() + DeltaT;
  double Tol=1.0e-3;
  double max_abs=0.0, cur_max=0.0;
  bool Done=false;
  do
   {
     f.step();

     // manually check fields-decayed condition
     double absEz = abs(f.get_field(Ez, eval_point));
     cur_max = fmax(cur_max, absEz);
     if ( f.round_time() >= NextCheckTime )
      { NextCheckTime += DeltaT;
        max_abs = fmax(max_abs, cur_max);
        if ( (max_abs>0.0) && cur_max < Tol*max_abs)
         Done=true;
        cur_max=0.0;
      };
     //printf("%.2e %.2e %.2e %.2e\n",f.round_time(),absEz,max_abs,cur_max);
   } while(!Done);

  //; for normalization run, save flux fields for refl. plane
  //(if no-bend? (save-flux "refl-flux" refl))
  if (no_bend)
   refl.save_hdf5(f, dataname);

  //(display-fluxes trans refl)
  printf("%11s | %12s | %12s\n", "   Time    ", " trans flux", "  refl flux");
  double f0=fcen-0.5*df, fstep=df/(nfreq-1);
  double *trans_flux=trans.flux();
  double *refl_flux=refl.flux();
  for(int nf=0; nf<nfreq; nf++)
   printf("%.4e | %+.4e | %+.4e\n",f0+nf*fstep,trans_flux[nf],refl_flux[nf]);

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
int main(int argc, char *argv[])
{
  initialize mpi(argc, argv);

  bool use_prisms=false;
  for(int narg=1; narg<argc; narg++)
   if (!strcasecmp(argv[narg],"--use-prisms"))
    use_prisms=true;

  bend_flux(false, use_prisms);
  bend_flux(true, use_prisms);

  // success if we made it here
  return 0;
}
