/***************************************************************/
/***************************************************************/
/***************************************************************/

/*
; Meep Tutorial: TE transmission and reflection through a cavity
; formed by a periodic sequence of holes in a dielectric waveguide,
; with a defect formed by a larger spacing between one pair of holes.

; This structure is based on one analyzed in:
;    S. Fan, J. N. Winn, A. Devenyi, J. C. Chen, R. D. Meade, and
;    J. D. Joannopoulos, "Guided and defect modes in periodic dielectric
;    waveguides," J. Opt. Soc. Am. B, 12 (7), 1267-1272 (1995).
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
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
void usage(char *progname)
{ master_printf("usage: %s [options]\n",progname);
  master_printf("options: \n");
  master_printf(" --use-symmetry         use geometric symmetries\n");
  abort();
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
int main(int argc, char *argv[])
{
  initialize mpi(argc, argv);

  bool use_symmetry=false;
  for(int narg=1; narg<argc; narg++)
   { if ( argv[narg]==0 )
      continue;
     if (!strcasecmp(argv[narg],"--use-symmetry") )
      { use_symmetry=true;
        master_printf("Using symmetry.\n");
      }
     else
      { master_printf("unknown command-line option %s (aborting)",argv[narg]);
        usage(argv[0]);
      }; 
   };

  // Some parameters to describe the geometry:
  double eps=13.0;   // dielectric constant of waveguide
  double w=1.2;      // width of waveguide
  double r=0.36;     // radius of holes
  double d=1.4;      // defect spacing (ordinary spacing = 1)
  int N=3;           // number of holes on either side of defect

  // The cell dimensions
  double sy=6.0;     // size of cell in y direction (perpendicular to wvg.)
  double pad=2.0;    // padding between last hole and PML edge
  double dpml=1.0;   // PML thickness

  // (define sx (+ (* 2 (+ pad dpml N)) d -1))
  // (set! geometry-lattice (make lattice (size sx sy no-size)))
  // (set! pml-layers (list (make pml (thickness dpml))))
  // (set-param! resolution 20)
  double sx = 2.0*(pad + dpml + N) + d -1.0; // size of cell in x dir
  double resolution=20.0;
  geometry_lattice.size.x=sx;
  geometry_lattice.size.y=sy;
  geometry_lattice.size.z=0.0;
  grid_volume gv = voltwo(sx, sy, resolution);
  gv.center_origin();

  symmetry sym = use_symmetry ? -mirror(Y,gv) : identity();
  structure the_structure(gv, dummy_eps, pml(dpml), sym);

  // (set! geometry
  //    (append ; combine lists of objects:
  //     (list (make block (center 0 0) (size infinity w infinity)
  // (material (make dielectric (epsilon eps)))))
  //      (geometric-object-duplicates (vector3 1 0) 0 (- N 1)
  // (make cylinder (center (/ d 2) 0) (radius r) (height infinity)
  // (material air)))
  // (geometric-object-duplicates (vector3 -1 0) 0 (- N 1)
  // (make cylinder (center (/ d -2) 0) (radius r) (height infinity)
  // (material air)))))
  material_type vacuum     = meep_geom::vacuum;
  material_type dielectric = meep_geom::make_dielectric(eps);
  geometric_object objects[7];
  vector3 origin = v3(0.0,   0.0,  0.0);
  vector3 xhat   = v3(1.0,   0.0,  0.0);
  vector3 yhat   = v3(0.0,   1.0,  0.0);
  vector3 zhat   = v3(0.0,   0.0,  1.0);
  vector3 size   = v3(ENORMOUS, w, ENORMOUS);
  double x0      = 0.5*d;
  double deltax  = 1.0;
  double height  = ENORMOUS;

  objects[0] = make_block(dielectric, origin, xhat, yhat, zhat, size);
  int no=1;
  for(int n=0; n<N; n++)
   { vector3 center=v3(x0 + n*deltax,  0.0, 0.0);
     objects[no++] = make_cylinder(vacuum, center, r, height, zhat);
   };
  for(int n=0; n<N; n++)
   { vector3 center=v3(-x0 - n*deltax,  0.0, 0.0);
     objects[no++] = make_cylinder(vacuum, center, r, height, zhat);
   };
  geometric_object_list g={ no, objects };
  meep_geom::set_materials_from_geometry(&the_structure, g);
  fields f(&the_structure);

  //
  double fcen = 0.25;  // pulse center frequency
  double df   = 0.2;   // pulse width (in frequency)
  //int nfreq   = 500;   // number of frequencies at which to compute flux
  
  //  (set! sources (list
  //  	(make source
  //		(src (make gaussian-src (frequency fcen) (fwidth df)))
  // 		(component Hz) (center 0 0))))
  gaussian_src_time src(fcen, df);
  component src_cmpt = Hz;
  f.add_point_source(src_cmpt, src, vec(0.0, 0.0));

  // (run-sources+ 400
  //    (at-beginning output-epsilon)
  //    (after-sources (harminv Hz (vector3 0) fcen df)))
  f.output_hdf5(Dielectric, f.total_volume());
  std::vector<cdouble> field_data;
  vec eval_pt(0.0,0.0);
  while( f.round_time() < f.last_source_time())
   f.step();

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  int dims[2], rank;
    
  volume v1d( vec(-0.25*sx, 0.0), vec(+0.25*sx, 0.0) );
  rank=f.get_array_slice_dimensions(v1d, dims);
  printf("1d slice: rank=%i, dim=%i\n",rank,dims[0]);
  double *slice1d=f.get_array_slice(v1d, Hz);
  FILE *f1=fopen("slice1d.out","w");
  for(int nx=0; nx<dims[0]; nx++)
   fprintf(f1,"%e %e\n",-0.5*sx + nx/resolution, slice1d[nx]);
  fclose(f1);

  volume v2d( vec(-0.25*sx, -0.15*sy), vec(+0.25*sx, +0.15*sy) );
  rank=f.get_array_slice_dimensions(v2d, dims);
  printf("2d slice: rank=%i, dims={%i,%i}\n",rank,dims[0],dims[1]);
  double *slice2d=f.get_array_slice(v2d, Hz);
  FILE *f2=fopen("slice2d.out","w");
  for(int nx=0; nx<dims[0]; nx++)
   for(int ny=0; ny<dims[1]; ny++)
   fprintf(f2,"%e %e %e\n",-0.5*sx + nx/resolution, 
                          -0.5*sy + ny/resolution, 
                          slice2d[nx*dims[1] + ny]);
  fclose(f2);
  master_printf("Thank you for your support.\n");

}
