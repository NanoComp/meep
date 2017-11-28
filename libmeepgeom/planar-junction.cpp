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

/***************************************************************/
/* dummy material function needed to pass to structure( )      */
/* constructor as a placeholder before we can call             */
/* set_materials_from_geometry                                 */
/***************************************************************/
double dummy_eps(const vec &) { return 1.0; }

/***************************************************************/
/***************************************************************/
/***************************************************************/
void usage(char *progname, const char *msg=0)
{ 
  if (msg) master_printf("%s\n",msg);
  master_printf("usage: %s [options]\n",progname);
  master_printf("options: \n");
  master_printf(" --use-symmetry    use geometric symmetries\n");
  master_printf(" --write-files     write reference data files\n");
  abort();
}

cdouble DefaultAmplitude(const vec &v) 
{ (void) v; // unused
  return 1.0;
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
  bool use_symmetry = false;
  bool point_source = false;
  int  band_num     = 1;
  for(int narg=1; narg<argc; narg++)
   { if ( argv[narg]==0 )
      continue;
     if (!strcasecmp(argv[narg],"--use-symmetry") )
      { use_symmetry=true;
        master_printf("Using symmetry.\n");
      }
     else if (!strcasecmp(argv[narg],"--point-source") )
      { point_source=true;
        master_printf("using point-source excitation\n");
      }
     else if (!strcasecmp(argv[narg],"--band-num"))
      { if ( ++narg >= argc )
         usage(argv[0], "error: no argument given for --band-num");
        sscanf(argv[narg], "%i", &band_num);
        master_printf("setting band-num=%i\n",band_num);
      }
     else if (!strcasecmp(argv[narg],"--ratio"))
      { if ( ++narg >= argc )
         usage(argv[0], "error: no argument given for --ratio");
        sscanf(argv[narg], "%le", &ratio);
        master_printf("setting ratio=%e\n",ratio);
      }
     else
      { master_printf("unknown command-line option %s (aborting)",argv[narg]);
        usage(argv[0]);
      }; 
   };

  /***************************************************************/
  /* initialize geometry, similar to holey_wvg_cavity   **********/
  /****************** *********************************************/
  double ncore=2.5;   // index of waveguide core
  double H=4.0;       // half-height of computational cell (Y direction)
  double L=8.0;       // half-length of computational cell (Z direction)
  double h1=1.0;      // half-height of narrow waveguide
  double h2=h1*ratio; // half-height of wider waveguide
  double a=10.0;      // resolution
  double dpml=1.0;    // PML thickness

  geometry_lattice.size.x=0.0;
  geometry_lattice.size.y=2.0*H;
  geometry_lattice.size.z=2.0*L;
  grid_volume gv=vol3d(0, 2.0*H, 2.0*L, a);
  gv.center_origin();
  symmetry sym = use_symmetry ? -mirror(Y,gv) : identity();
  structure the_structure(gv, dummy_eps); //, pml(dpml)); //, pml(dpml), sym);

  material_type dielectric = meep_geom::make_dielectric(ncore*ncore);
  vector3 xhat    = v3(1.0,   0.0,  0.0);
  vector3 yhat    = v3(0.0,   1.0,  0.0);
  vector3 zhat    = v3(0.0,   0.0,  1.0);
  geometric_object objects[2];
  objects[0] = make_block( dielectric,
                           v3(0.0, 0.0, -0.5*L),     // center
                           xhat, yhat, zhat,
                           v3(ENORMOUS, 2.0*h1, L)   // size;
                         );
  objects[1] = make_block( dielectric,
                           v3(0.0, 0.0, +0.5*L),     // center
                           xhat, yhat, zhat,
                           v3(ENORMOUS, 2.0*h2, L)   // size;
                         );
  geometric_object_list g={ 2, objects };
  meep_geom::set_materials_from_geometry(&the_structure, g);
  fields f(&the_structure);
  f.output_hdf5(Dielectric,f.total_volume());

  /***************************************************************/
  /* add source                                                  */
  /***************************************************************/
  double fcen = 0.15;  // center frequency
  double df   = 0.1;   // bandwidth
  int nfreq   = 10;    // number of frequency points
  gaussian_src_time src(fcen, df);
  component src_cmpt = Ex;

  if (point_source)
   { 
     f.add_point_source(src_cmpt, src, vec(0.0, 0.0, -0.5*L));
   }
  else
   { 
     volume eigenslice( vec(0,-h1,-0.5*L), vec(0,h1,-0.5*L) );
     src_cmpt=Dielectric;
     f.add_eigenmode_source(src_cmpt,   // component c
                            src,        // src_time src
                            Z,          // direction d
                            eigenslice, // const volume &where
                            eigenslice, // const volume &eig_vol
                            band_num,   // int band_num
                            vec(0,0,0.25), // const vec &kpoint,
                            true,       // match_frequency
                            2,          // no-parity,
                            a,          // resolution
                            1.e-7,      // eigensolver_tol 
                            1.0        // amp
                            //DefaultAmplitude   //A
                           );
   };

  /***************************************************************/
  /* add flux planes                                             */
  /***************************************************************/
  volume fv1 = volume( vec(0,-0.5*H,-0.25*L), vec(0,+0.5*H,-0.25*L) );
  volume fv2 = volume( vec(0,-0.5*H, 0.00*L), vec(0,+0.5*H, 0.00*L) );
  volume fv3 = volume( vec(0,-0.5*H,+0.25*L), vec(0,+0.5*H,+0.25*L) );
  volume fv4 = volume( vec(0,-0.5*H,+0.50*L), vec(0,+0.5*H,+0.50*L) );
  dft_flux flux1=f.add_dft_flux_plane(fv1, fcen-0.5*df, fcen+0.5*df, nfreq);
  dft_flux flux2=f.add_dft_flux_plane(fv2, fcen-0.5*df, fcen+0.5*df, nfreq);
  dft_flux flux3=f.add_dft_flux_plane(fv3, fcen-0.5*df, fcen+0.5*df, nfreq);
  dft_flux flux4=f.add_dft_flux_plane(fv4, fcen-0.5*df, fcen+0.5*df, nfreq);

  /***************************************************************/
  /* timestep until fields decayed + 50 time unit ****************/
  /***************************************************************/
  vec eval_point(0.0, 0.0, 0.0);
  double DeltaT=500.0, NextCheckTime = f.round_time() + DeltaT;
  double max_abs=0.0, cur_max=0.0;
  double Tol=1.0e-3;
  bool Done=false;
  do
   {  
     f.step();

     // manually check fields-decayed condition
     double absEx = abs(f.get_field(Ex, eval_point));
     cur_max = fmax(cur_max, absEx);
     if ( f.round_time() >= NextCheckTime )
      { NextCheckTime += DeltaT;
        max_abs = fmax(max_abs, cur_max);
        if ( (max_abs>0.0) && cur_max < Tol*max_abs)
         Done=true;
        cur_max=0.0;
      };
   } while(!Done);
  printf("Done timestepping at t=%e.\n",f.round_time());
 
  flux1.save_hdf5( f, "flux1" );
  flux2.save_hdf5( f, "flux2" );
  flux3.save_hdf5( f, "flux3" );
  flux4.save_hdf5( f, "flux4" );

  f.output_hdf5(Ex,f.total_volume());
  f.output_hdf5(Hy,f.total_volume());

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  volume vSlice( vec(0, -H, -L), vec(0, H, L) );

  int dims[2];
  int rank=f.get_array_slice_dimensions(vSlice, dims);
  master_printf("Rank=%i, dims={",rank);
  for(int r=0; r<rank; r++) master_printf("%i ",dims[r]);
  master_printf("\n");
  if (rank!=2)
   abort("incorrect dimensions for 2D slice");

  double *Eps=f.get_array_slice(vSlice, Dielectric);
  cdouble *Ex_slice=f.get_complex_array_slice(vSlice, Ex);
  cdouble *Hy_slice=f.get_complex_array_slice(vSlice, Hy);
  double *Sz_slice=f.get_array_slice(vSlice, Sz);
  int NY=dims[0], NZ=dims[1];
  FILE *ff=fopen("waveguide-junction.out","w");
  for(int ny=0; ny<NY; ny++)
   for(int nz=0; nz<NZ; nz++)
    fprintf(ff,"%e %e %e %e %e %e %e %e\n",-H + ny/a, -L +nz/a,
               Eps[ny + nz*NY], Sz_slice[ny+nz*NY],
               real(Ex_slice[ny+nz*NY]), imag(Ex_slice[ny + nz*NY]),
               real(Hy_slice[ny+nz*NY]), imag(Hy_slice[ny + nz*NY])
           );
  fclose(ff);

  ff=fopen("waveguide-junction.out2","w");
  for(int ny=0; ny<NY; ny++)
   for(int nz=0; nz<NZ; nz++)
    fprintf(ff,"%e %e %e %e %e %e %e %e\n",-H + ny/a, -L +nz/a,
               Eps[ny*NZ + nz], Sz_slice[ny*NZ+ny],
               real(Ex_slice[ny*NZ+nz]), imag(Ex_slice[ny*NZ+ny]),
               real(Hy_slice[ny*NZ+nz]), imag(Hy_slice[ny*NZ+nz])
           );
  fclose(ff);

  return 0;
}
