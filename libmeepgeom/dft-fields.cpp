/***************************************************************/
/* unit test for dft_fields and output_dft functionality       */
/***************************************************************/
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
void Run(bool Pulse, double resolution)
{
  /***************************************************************/
  /* initialize geometry                                         */
  /***************************************************************/
  double n=3.4;     // index of waveguide
  double w=1.0;     // width of waveguide
  double r=1.0;     // inner radius of ring
  double pad=4;     // padding between waveguide and edge of PML
  double dpml=2;    // thickness of PML

  double sxy        = 2.0*(r+w+pad+dpml);  // cell size
  geometry_lattice.size.x=sxy;
  geometry_lattice.size.y=sxy;
  geometry_lattice.size.z=0.0;
  grid_volume gv = voltwo(sxy, sxy, resolution);
  gv.center_origin();
  symmetry sym=identity(); // mirror(Y, gv);
  structure the_structure(gv, dummy_eps, pml(dpml), sym);

  /***************************************************************/
  /* add objects                                                 */
  /***************************************************************/
  meep_geom::material_type dielectric = meep_geom::make_dielectric(n*n);
  geometric_object objects[2];
  vector3 v3zero = {0.0, 0.0, 0.0};
  vector3 zaxis  = {0.0, 0.0, 1.0};
  objects[0] = make_cylinder(dielectric, v3zero, r+w, ENORMOUS, zaxis);
  objects[1] = make_cylinder(meep_geom::vacuum,  v3zero, r, ENORMOUS, zaxis);
  geometric_object_list g={ 2, objects };
  meep_geom::set_materials_from_geometry(&the_structure, g);
  fields f(&the_structure);
  f.step(); // single timestep to trigger internal initialization
  f.output_hdf5(Dielectric, f.v);

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  double fcen = 0.118; // ; pulse center frequency
  double df   = 0.1;   // ; df
  vec x0(r+0.1,0.0);   // ; source location
  if(Pulse)
   { 
     f.add_point_source(Ez, gaussian_src_time(fcen,df), x0);

     component components[6] = {Ex, Ey, Ez, Hx, Hy, Hz};
     dft_fields dftFields    = f.add_dft_fields(components, 6, f.v, fcen, fcen, 1);
     dft_flux dftFlux        = f.add_dft_flux(X, f.v, fcen, fcen, 1);

     while( f.round_time() < f.last_source_time() + 100.0 )
      f.step();

     f.output_dft(dftFlux,   "dft-flux");
     f.output_dft(dftFields, "dft-fields");
   }
  else
   { 
     f.add_point_source(Ez, continuous_src_time(fcen,df), x0);
     f.solve_cw(1e-8, 10000, 10);
     f.output_hdf5(Ez, f.v);
     f.output_hdf5(Hx, f.v);
     f.output_hdf5(Hy, f.v);
   }

}

/***************************************************************/
/* compute the L2 norm of two complex-valued HDF5 datasets,    */
/* after (a) normalizing each data set by its maximum amplitude*/
/* and (b) compensating for a constant overall phase factor    */
/* between the datasets.                                       */
/***************************************************************/
double compare_complex_hdf5_datasets(const char *file1, const char *name1,
                                     const char *file2, const char *name2,
                                     int expected_rank=2)
{
  char dataname[100];

  // read dataset 1
  h5file f1(file1, h5file::READONLY, false);
  int rank1, *dims1=new int[expected_rank];
  snprintf(dataname,100,"%s.r",name1);
  double *rdata1 = f1.read(dataname, &rank1, dims1, expected_rank);
  snprintf(dataname,100,"%s.i",name1);
  double *idata1 = f1.read(dataname, &rank1, dims1, expected_rank);
  if (!rdata1 || !idata1) return -1.0;

  // read dataset 2
  h5file f2(file2, h5file::READONLY, false);
  int rank2, *dims2=new int[expected_rank];
  snprintf(dataname,100,"%s.r",name2);
  double *rdata2 = f2.read(dataname, &rank2, dims2, expected_rank);
  snprintf(dataname,100,"%s.i",name2);
  double *idata2 = f2.read(dataname, &rank2, dims2, expected_rank);
  if (!rdata2 || !idata2) return -1.0;

  // check same size
  bool same_size = (rank1==rank2);
  for(int d=0; same_size && d<rank1; d++)
   if(dims1[d]!=dims2[d])
    same_size = false;
  if (!same_size)
   meep::abort("data size mismatch\n");

  // first pass to normalize each dataset to its maximum absolute magnitude;
  // we also note the phase difference between the datasets at their points
  // of maximum magnitude so we can compensate for this in the comparison below.
  int length = dims1[0];
  for(int d=1; d<rank1; d++)
   length*=dims1[d];

  double max_abs1=0.0, max_abs2=0.0;
  double max_arg1=0.0, max_arg2=0.0;
  for(int n=0; n<length; n++)
   { cdouble z1 = cdouble(rdata1[n], idata1[n]);
     if (abs(z1) > max_abs1)
      { max_abs1 = abs(z1);
        max_arg1 = arg(z1);
      }
     cdouble z2 = cdouble(rdata2[n], idata2[n]);
     if (abs(z2) > max_abs2)
      { max_abs2 = abs(z2);
        max_arg2 = arg(z2);
      }
   }

  // second pass to get L2 norm of difference between normalized data sets
  double norm1=0.0, norm2=0.0, normdiff=0.0;
  cdouble phase_fac = exp( cdouble(0,1)*(max_arg1 - max_arg2) );
  for(int n=0; n<length; n++)
   { cdouble z1 = cdouble(rdata1[n], idata1[n]) / max_abs1;
     cdouble z2 = phase_fac*cdouble(rdata2[n], idata2[n]) / max_abs2;
     norm1    += norm(z1);
     norm2    += norm(z2);
     normdiff += norm(z1-z2);
   }
  norm1    = sqrt(norm1)    / ((double)length);
  norm2    = sqrt(norm2)    / ((double)length);
  normdiff = sqrt(normdiff) / ((double)length);

  if (am_master())
   { FILE *ff=fopen("/tmp/dft-fields.log","a");
     fprintf(ff,"length %i:\n",length);
     fprintf(ff,"norm1    : %e\n",norm1);
     fprintf(ff,"norm2    : %e\n",norm2);
     fprintf(ff,"normdiff : %e\n",normdiff);
     fprintf(ff,"\n");
     fclose(ff);
   }

  return normdiff / fmax(norm1,norm2);
}


/***************************************************************/
/***************************************************************/
/***************************************************************/
int main(int argc, char *argv[])
{
  initialize mpi(argc, argv);

#define NUM_RESOLUTIONS 3    
  double resolutions[NUM_RESOLUTIONS]={10.0, 20.0, 40.0};
  double L2Errors[NUM_RESOLUTIONS];
  for(int nr=0; nr<NUM_RESOLUTIONS; nr++)
   {
     double res = resolutions[nr];
     Run(true,  res);
     Run(false, res);

     L2Errors[nr] 
      = compare_complex_hdf5_datasets("dft-fields.h5","ez_0",
                                      "ez-000000.05.h5","ez",2);
     if (L2Errors[nr]==-1.0) // files couldn't be read or datasets had different sizes
      return 1;
   }

  // condition for test to pass is simply that L2 norm decreases
  // with increasing resolution; note that this is also implicitly
  // testing the success of the dft_fields and output_dft 
  // implementations
  if ( (L2Errors[0] > L2Errors[1]) && (L2Errors[1] > L2Errors[2]) )
   return 0;

  return 1;
}
