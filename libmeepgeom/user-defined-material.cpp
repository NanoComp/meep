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
/* return true if the datasets match, false if not             */
/***************************************************************/
bool compare_hdf5_datasets(const char *file1, const char *name1,
                           const char *file2, const char *name2,
                           int expected_rank=2,
                           double rel_tol=1.0e-2)
{
  h5file f1(file1, h5file::READONLY, false);
  int rank1;
  int *dims1=new int[expected_rank];
  double *data1 = f1.read(name1, &rank1, dims1, expected_rank);
  if (!data1) return false;

  h5file f2(file2, h5file::READONLY, false);
  int rank2;
  int *dims2=new int[expected_rank];
  double *data2 = f2.read(name2, &rank2, dims2, expected_rank);
  if (!data2) return false;

  if ( rank1!=expected_rank || rank2!=expected_rank ) return false;

  size_t size = 1;
  for(int r=0; r<expected_rank; r++)
   { if (dims1[r]!=dims2[r])
      return false;
     size *= dims1[r];
   };

  double norm1=0.0, norm2=0.0, normDelta=0.0;
  for(size_t n=0; n<size; n++)
   { double d1=data1[n], d2=data2[n], delta=d1-d2; 
     norm1     += d1*d1;
     norm2     += d2*d2;
     normDelta += delta*delta;
   };
  norm1     = sqrt(norm1)/ size;
  norm2     = sqrt(norm2) / size;
  normDelta = sqrt(normDelta) / size;

  double norm12 = fmax( sqrt(norm1), sqrt(norm2) );
printf("Norm{1,2,d} = {%e,%e,%e} (%e)\n",norm1,norm2,normDelta,normDelta/norm12);
  if (normDelta > rel_tol*norm12)
   return false;
  return true;
}

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
   nn=3.5;
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
         abort("no option specified for --eps_ref_file");
        eps_ref_file=argv[++narg];
      }
     else
      abort("unrecognized command-line option %s",argv[narg]);
   };
   std::string eps_ref_path = DATADIR + eps_ref_file;

  double resolution = 100.0;

  vector3 lattice_size = {10.0, 10.0, 0.0};
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
  meep_geom::set_materials_from_geometry(lattice_size,
                                         &the_structure, g,
                                         use_anisotropic_averaging,
                                         DEFAULT_SUBPIXEL_TOL,
                                         DEFAULT_SUBPIXEL_MAXEVAL,
                                         ensure_periodicity,
                                         verbose,
                                         my_material);

  fields f(&the_structure);

  f.output_hdf5(Dielectric, f.total_volume());
  bool status=compare_hdf5_datasets("eps-000000000.h5", "eps",
                                     eps_ref_path.c_str(), "eps");
  if (status)
   master_printf("user-defined-material test successful.\n");
  else
   abort("user-defined-material test failed.\n");
  return 0;

}
