/***************************************************************/
/* test of user-defined materials in libmeepgeom.              */
/* this code creates a user-defined material (defined by the   */
/* function my_material_func below) that reproduces the        */
/* material geometry of the cyl_ellipsoid example.             */
/***************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <complex>

#include "meep.hpp"

#include "ctl-math.h"
#include "ctlgeom.h"

#include "meepgeom.hpp"

#ifndef DATADIR
  #define DATADIR "./"
#endif

using namespace meep;

typedef std::complex<double> cdouble;

/***************************************************************/
/* dummy material function needed to pass to structure( )      */
/* constructor as a placeholder before we can call             */
/* set_materials_from_geometry                                 */
/***************************************************************/
double dummy_eps(const vec &) { return 1.0; }

/***************************************************************/
/* material function that recreates the ellipsoid-in-cylinder  */
/* configuration of the cyl-ellipsoid sample code              */
/***************************************************************/
void my_material_func(vector3 p, void *user_data, meep_geom::medium_struct *m)
{ 
  (void) user_data;
#define R1X 0.5
#define R1Y 1.0
#define R2  3.0

  double x=p.x, y=p.y;

  // test for point inside inner ellipsoid
  double nn;
  if ( (x*x/(R1X*R1X) + y*y/(R1Y*R1Y)) < 1.0 )
   nn=1.0;
  else if ( ( x*x/(R2*R2) + y*y/(R2*R2) ) < 1.0 )
   nn=3.5*3.5;
  else 
   nn=1.0;

  m->epsilon_diag.x = m->epsilon_diag.y = m->epsilon_diag.z = nn*nn;

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
int main(int argc, char *argv[])
{
  initialize mpi(argc, argv);

  // simple argument parsing
  meep::component src_cmpt=Ez;
  std::string eps_ref_file = "cyl-ellipsoid-eps-ref.h5";
  for(int narg=1; narg<argc; narg++)
   {
     if (argv[narg] && !strcmp(argv[narg],"--eps_ref_file"))
      { if (narg+1 == argc)
         meep::abort("no option specified for --eps_ref_file");
        eps_ref_file=argv[++narg];
      }
     else
      meep::abort("unrecognized command-line option %s",argv[narg]);
   };
   std::string eps_ref_path = DATADIR + eps_ref_file;

  double resolution = 100.0;

  geometry_lattice.size.x=10.0;
  geometry_lattice.size.y=10.0;
  geometry_lattice.size.z=0.0;
  grid_volume gv = voltwo(10.0, 10.0, resolution);
  gv.center_origin();
  symmetry sym = (src_cmpt==Ez) ?  mirror(X,gv) + mirror(Y,gv)
                                : -mirror(X,gv) - mirror(Y,gv);
  structure the_structure(gv, dummy_eps, pml(1.0), sym);

  meep_geom::material_type my_material
   = meep_geom::make_user_material(my_material_func, 0);

  geometric_object_list g={0, 0};
  bool use_anisotropic_averaging=true;
  bool ensure_periodicity=true;
  bool verbose=false;
  meep_geom::set_materials_from_geometry(&the_structure, g,
                                         use_anisotropic_averaging,
                                         DEFAULT_SUBPIXEL_TOL,
                                         DEFAULT_SUBPIXEL_MAXEVAL,
                                         ensure_periodicity,
                                         verbose,
                                         my_material);

  fields f(&the_structure);

  // extract dielectric constant on a line of points
  // lying on the x-axis from the origin to the right cell wall
  volume v1d( vec(0.0, 0.0), vec(10.0, 0.0) );
  int dims1D[1];
  int rank=f.get_array_slice_dimensions(v1d, dims1D);
  double *slice1d=f.get_array_slice(v1d,Dielectric);

  FILE *ff=fopen("/tmp/doom.out","w");
  for(int n=0; n<dims1D[0]; n++)
   fprintf(ff,"%i %e\n",n,slice1d[n]);
  fclose(ff);

  return 0;

}
