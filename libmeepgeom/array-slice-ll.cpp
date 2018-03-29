/***************************************************************/
/***************************************************************/
/***************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
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

vector3 v3(double x, double y=0.0, double z=0.0)
{
  vector3 v;
  v.x=x; v.y=y; v.z=z;
  return v;
}

// passthrough field function
cdouble default_field_function(const cdouble *fields,
			const vec &loc, void *data_)
{
  (void) loc; // unused
  (void) data_; // unused
  return fields[0];
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
#define RELTOL 1.0e-6
double Compare(double *d1, double *d2, int N, const char *Name)
{
  double Norm1=0.0, Norm2=0.0, NormDelta=0.0;
  for(int n=0; n<N; n++)
   { Norm1     += d1[n]*d1[n];
     Norm2     += d2[n]*d2[n];
     NormDelta += (d1[n]-d2[n])*(d1[n]-d2[n]);
   };
  Norm1=sqrt(Norm1);
  Norm2=sqrt(Norm2);
  NormDelta=sqrt(NormDelta);
  double RelErr = NormDelta / (0.5*(Norm1+Norm2));
  if (RelErr > RELTOL)
   abort("fail: rel error in %s data = %e\n",Name,RelErr);
  return RelErr;
}

double Compare(cdouble *d1, cdouble *d2, int N, const char *Name)
{
  double Norm1=0.0, Norm2=0.0, NormDelta=0.0;
  for(int n=0; n<N; n++)
   { Norm1     += norm(d1[n]);
     Norm2     += norm(d2[n]);
     NormDelta += norm(d1[n]-d2[n]);
   };
  Norm1=sqrt(Norm1);
  Norm2=sqrt(Norm2);
  NormDelta=sqrt(NormDelta);
  double RelErr = NormDelta / (0.5*(Norm1+Norm2));
  if (RelErr > RELTOL)
   abort("fail: rel error in %s data = %e\n",Name,RelErr);
  return RelErr;
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
  master_printf(" --use-symmetry    use geometric symmetries\n");
  master_printf(" --write-files     write reference data files\n");
  abort();
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
int main(int argc, char *argv[])
{
  initialize mpi(argc, argv);

  /***************************************************************/
  /* parse command-line options **********************************/
  /***************************************************************/
  bool use_symmetry=false;
  bool write_files=false;
  for(int narg=1; narg<argc; narg++)
   { if ( argv[narg]==0 )
      continue;
     if (!strcasecmp(argv[narg],"--use-symmetry") )
      { use_symmetry=true;
        master_printf("Using symmetry.\n");
      }
     else if (!strcasecmp(argv[narg],"--write-files") )
      { write_files=true;
        master_printf("writing HDF5 data files");
      }
     else
      { master_printf("unknown command-line option %s (aborting)",argv[narg]);
        usage(argv[0]);
      };
   };

  /***************************************************************/
  /* initialize geometry, similar to holey_wvg_cavity   **********/
  /***************************************************************/
  double eps=13.0;   // dielectric constant of waveguide
  double w=1.2;      // width of waveguide
  double r=0.36;     // radius of holes
  double d=1.4;      // defect spacing (ordinary spacing = 1)
  int N=3;           // number of holes on either side of defect
  double sy=6.0;     // size of cell in y direction (perpendicular to wvg.)
  double pad=2.0;    // padding between last hole and PML edge
  double dpml=1.0;   // PML thickness
  double sx = 2.0*(pad + dpml + N) + d -1.0; // size of cell in x dir
  double resolution=20.0;
  geometry_lattice.size.x=sx;
  geometry_lattice.size.y=sy;
  geometry_lattice.size.z=0.0;
  grid_volume gv = voltwo(sx, sy, resolution);
  gv.center_origin();
  symmetry sym = use_symmetry ? -mirror(Y,gv) : identity();
  structure the_structure(gv, dummy_eps, pml(dpml), sym);
  meep_geom::material_type vacuum     = meep_geom::vacuum;
  meep_geom::material_type dielectric = meep_geom::make_dielectric(eps);
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

  /***************************************************************/
  /* add source and timestep until source has finished (no later)*/
  /***************************************************************/
  double fcen = 0.25;  // pulse center frequency
  double df   = 0.2;   // pulse width (in frequency)
  gaussian_src_time src(fcen, df);
  component src_cmpt = Hz;
  f.add_point_source(src_cmpt, src, vec(0.0, 0.0));
  while( f.round_time() < f.last_source_time())
   f.step();

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  double xMin = -0.25*sx, xMax=+0.25*sx;
  double yMin = -0.15*sy, yMax=+0.15*sy;
  volume v1d( vec(xMin, 0.0),  vec(xMax, 0.0) );
  volume v2d( vec(xMin, yMin), vec(xMax, yMax) );

  int rank;
  size_t dims1D[1], dims2D[2];
  cdouble *file_slice1d=0;
  double *file_slice2d=0;

#define H5FILENAME DATADIR"array-slice-ll"
#define NX 126
#define NY 38
  if (write_files)
   {
     h5file *file = f.open_h5file(H5FILENAME);
     f.output_hdf5(Hz, v1d, file);
     f.output_hdf5(Sy, v2d, file);
     master_printf("Wrote binary data to file %s.h5\n",H5FILENAME);
     delete file;
     exit(0);
   }
  else
   {
     //
     // read 1D and 2D array-slice data from HDF5 file
     //
     h5file *file = f.open_h5file(H5FILENAME, h5file::READONLY);
     double *rdata = file->read("hz.r", &rank, dims1D, 1);
     if (rank!=1 || dims1D[0]!=NX)
      abort("failed to read 1D data(hz.r) from file %s.h5",H5FILENAME);
     double *idata = file->read("hz.i", &rank, dims1D, 1);
     if (rank!=1 || dims1D[0]!=NX)
      abort("failed to read 1D data(hz.i) from file %s.h5",H5FILENAME);
     file_slice1d = new cdouble[dims1D[0]];
     for(size_t n=0; n<dims1D[0]; n++)
      file_slice1d[n] = cdouble(rdata[n], idata[n]);
     delete[] rdata;
     delete[] idata;

     file_slice2d = file->read("sy", &rank, dims2D, 2);
     if (rank!=2 || dims2D[0]!=NX || dims2D[1]!=NY)
      abort("failed to read 2D reference data from file %s.h5",H5FILENAME);
     delete file;

     //
     // generate 1D and 2D array slices and compare to
     // data read from file
     //
     rank=f.get_array_slice_dimensions(v1d, dims1D);
     if (rank!=1 || dims1D[0]!=NX)
      abort("incorrect dimensions for 1D slice");
     cdouble *slice1d=f.get_complex_array_slice(v1d, Hz);
     double RelErr1D=Compare(slice1d, file_slice1d, NX, "Hz_1d");
     master_printf("1D: rel error %e\n",RelErr1D);

     rank=f.get_array_slice_dimensions(v2d, dims2D);
     if (rank!=2 || dims2D[0]!=NX || dims2D[1]!=NY)
      abort("incorrect dimensions for 2D slice");
     double *slice2d=f.get_array_slice(v2d, Sy);
     double RelErr2D=Compare(slice2d, file_slice2d, NX*NY, "Sy_2d");
     master_printf("2D: rel error %e\n",RelErr2D);

   }; // if (write_files) ... else ...

  return 0;
}
