/***************************************************************/
/* simple demonstration of mode expansion in a z-invariant     */
/* cylindrical waveguide                                       */
/***************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <complex>
#include <vector>

#include "meep.hpp"

#include "ctl-math.h"
#include "ctlgeom.h"
#include "mpb.h"

#include "meepgeom.hpp"

using namespace meep;

typedef std::complex<double> cdouble;

namespace meep{
double get_group_velocity(void *vedata);
void destroy_eigenmode_data(void *vedata);
}

/***************************************************************/
/* function to provide initial guess for the wavevector of the */
/* eigenmode with frequency freq and band index band           */
/***************************************************************/
bool equal_float(double x, double y) { return (float)x == (float)y; }

vec k_guess(void *user_data, double freq, int band_num)
{
  (void) user_data;
  (void) freq;
  (void) band_num;

// hard-coded dispersion relations
  if ( equal_float(freq,0.25) )
   { if (band_num==1) return vec(0.0, 0.0, 0.736917);
     if (band_num==2) return vec(0.0, 0.0, 0.622793);
     if (band_num==3) return vec(0.0, 0.0, 0.591875);
     if (band_num==4) return vec(0.0, 0.0, 0.521444);
   }
  else if ( equal_float(freq,0.50) )
   { if (band_num==1) return vec(0.0, 0.0, 1.62266);
     if (band_num==2) return vec(0.0, 0.0, 1.58266);
     if (band_num==3) return vec(0.0, 0.0, 1.55956);
     if (band_num==4) return vec(0.0, 0.0, 1.52044);
   }
  else if ( equal_float(freq,1.00) )
   { if (band_num==1) return vec(0.0, 0.0, 3.37072);
     if (band_num==2) return vec(0.0, 0.0, 1.58266);
     if (band_num==3) return vec(0.0, 0.0, 1.55956);
     if (band_num==4) return vec(0.0, 0.0, 3.31898);
   };

  return vec(0.0, 0.0, 0.736917);
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
int main(int argc, char *argv[])
{
  initialize mpi(argc, argv);

  /***************************************************************/
  /* parse command-line options **********************************/
  /***************************************************************/
  bool use_symmetry = false;
  bool point_source = false;
  int  band_num     = 1;
  double ratio      = 1.0;
  bool CW           = false;
  char *filebase    = "dj";
  for(int narg=1; narg<argc; narg++)
   { if ( argv[narg]==0 )
      continue;
     if (!strcasecmp(argv[narg],"--use-symmetry") )
      { use_symmetry=true;
        master_printf("Using symmetry.\n");
      }
     if (!strcasecmp(argv[narg],"--CW") )
      { CW=true;
        master_printf("Using CW source.\n");
      }
     else if (!strcasecmp(argv[narg],"--point-source") )
      { point_source=true;
        master_printf("using point-source excitation\n");
      }
     else if (!strcasecmp(argv[narg],"--band-num"))
      { if ( ++narg >= argc )
         abort("error: no argument given for --band-num");
        sscanf(argv[narg], "%i", &band_num);
        master_printf("setting band-num=%i\n",band_num);
      }
     else if (!strcasecmp(argv[narg],"--ratio"))
      { if ( ++narg >= argc )
         abort("error: no argument given for --ratio");
        sscanf(argv[narg], "%le", &ratio);
        master_printf("setting ratio=%e\n",ratio);
      }
     else if (!strcasecmp(argv[narg],"--filebase"))
      { if ( ++narg >= argc )
         abort("error: no argument given for --filebase");
        filebase=strdup(argv[narg]);
      }
     else
      { master_printf("unknown command-line option %s (aborting)",argv[narg]);
        //usage(argv[0]);
      }; 
   };

  /***************************************************************/
  /* set up geometry: cylinder of radius r in a 2D computational */
  /* cell of size LxL (z-invariant)                              */
  /***************************************************************/
  double Eps=11.7;       // permittivity of waveguide

  double SX=1.0;         // duct size in X direction
  double SYA=2.0;        // smaller duct size in Y direction
  double SYB=ratio*SYA;  // larger duct size in Y direction

  double WXY =  5.0;     // size of computational cell in transverse (XY) directions
  double WZ  = 10.0;     // size of computational cell in propagation direction

  double dpml = 1.0;     // thickness of PML

  double resolution = 5.0;

  geometry_lattice.size.x=WXY;
  geometry_lattice.size.y=WXY;
  geometry_lattice.size.z=WZ;
  grid_volume gv = vol3d(WXY, WXY, WZ, resolution);
  gv.center_origin();

  symmetry sym=identity();

  structure the_structure(gv, dummy_eps, pml(dpml), sym);

  meep_geom::material_type dielectric = meep_geom::make_dielectric(Eps);
  geometric_object objects[2];
  vector3 xhat    = {1.0,0.0,0.0};
  vector3 yhat    = {0.0,1.0,0.0};
  vector3 zhat    = {0.0,0.0,1.0};
  double zcenter  = 0.25*WZ;
  vector3 centerA = {0.0,0.0,-1.0*zcenter};
  vector3 sizeA   = {SX,SYA,+0.5*WZ};
  vector3 centerB = {0.0,0.0,+1.0*zcenter};
  vector3 sizeB   = {SX,SYB,+0.5*WZ};
  objects[0] = make_block(dielectric, centerA, xhat, yhat, zhat, sizeA);
  objects[1] = make_block(dielectric, centerB, xhat, yhat, zhat, sizeB);
  geometric_object_list g={ 2, objects }; 
  meep_geom::set_materials_from_geometry(&the_structure, g);
  fields f(&the_structure);

  f.output_hdf5(Dielectric, f.total_volume(), 0, false, true);

  /***************************************************************/
  /* add sources: point source (if --point_source option present)*/
  /* or eigenmode source of band band_num                        */
  /***************************************************************/
  double fcen = 0.25;  // ; pulse center frequency
  double df   = 0.25;  // ; df
  int nfreq   = 1;
  gaussian_src_time gsrc(fcen, df);
  continuous_src_time csrc(fcen);

  vec minA(-0.5*WXY, -0.5*WXY, -1.0*zcenter);
  vec maxA(+0.5*WXY, +0.5*WXY, -1.0*zcenter);
  vec minB(-0.5*WXY, -0.5*WXY, +1.0*zcenter);
  vec maxB(+0.5*WXY, +0.5*WXY, +1.0*zcenter);
  vec minC(-0.5*WXY,  0.0,     -WZ);
  vec maxC(+0.5*WXY,  0.0,     +WZ);

  volume fvA(minA, maxA);
  volume fvB(minB, maxB);
  volume fvC(minC, maxC);

  if (point_source)
   { 
     f.add_point_source(Ez, gsrc, vec(0.1, 0.2, -1.0*zcenter) ); 
   }
  else
   {
     vec kpoint=k_guess(0, fcen, band_num);
     bool match_frequency = true;
     int parity = 0; // NO_PARITY
     double eigensolver_tol=1.0e-3;
     if (CW)
      f.add_eigenmode_source(Dielectric, csrc,
                             Z, fvA, fvA, band_num,
                             kpoint, match_frequency, parity,
                             resolution, eigensolver_tol, 1.0);
     else                       
      f.add_eigenmode_source(Dielectric, gsrc,
                             Z, fvA, fvA, band_num,
                             kpoint, match_frequency, parity,
                             resolution, eigensolver_tol, 1.0);
   };

  /***************************************************************/
  /* if the --CW option was specified, just export the           */
  /* instantaneous fields on the upper and lower flux planes     */
  /***************************************************************/
  if (CW)
   { while(f.round_time() < 50.0)
      f.step();
     f.output_hdf5(Ex, fvA, 0, false, false, "A");
     f.output_hdf5(Ey, fvA, 0, false, false, "A");
     f.output_hdf5(Hx, fvA, 0, false, false, "A");
     f.output_hdf5(Hy, fvA, 0, false, false, "A");
     f.output_hdf5(Ex, fvB, 0, false, false, "B");
     f.output_hdf5(Ey, fvB, 0, false, false, "B");
     f.output_hdf5(Hx, fvB, 0, false, false, "B");
     f.output_hdf5(Hy, fvB, 0, false, false, "B");
     f.output_hdf5(Sz, fvC, 0, false, false, "Sz");
     return 0;
   };

  /***************************************************************/
  /* add flux plane and timestep to accumulate frequency-domain  */
  /* fields at nfreq frequencies                                 */
  /***************************************************************/
  dft_flux fluxA=f.add_dft_flux_plane(fvA, fcen, fcen, nfreq);
  dft_flux fluxB=f.add_dft_flux_plane(fvB, fcen, fcen, nfreq);
  dft_flux fluxC=f.add_dft_flux_plane(fvC, fcen, fcen, nfreq);

  // (run-sources+ 10
  while( f.round_time() < (f.last_source_time() + 1.0/fcen) )
   f.step();

  /***************************************************************/
  /* write HDF5 files for visualization **************************/
  /***************************************************************/
  char filename[100];
  snprintf(filename,100,"%s_fluxA",filebase);
  f.output_flux_fields(&fluxA, fvA, filename);

  snprintf(filename,100,"%s_fluxB",filebase);
  f.output_flux_fields(&fluxB, fvB, filename);

  snprintf(filename,100,"%s_fluxC",filebase);
  f.output_flux_fields(&fluxC, fvC, filename);

  snprintf(filename,100,"%s.dat",filebase);
  FILE *ff=fopen(filename,"w");

  double eigensolver_tol=1.0e-3;
  void *vedata[4];
  double vgrp[4];
  for(int nb=0; nb<4; nb++)
   { int band_num=nb+1;
     vedata[nb]=f.get_eigenmode(fcen, Z, fvB, fvB,
                                band_num, k_guess(0,fcen,band_num),
                                true, 0, resolution,
                                eigensolver_tol);
     vgrp[nb]=get_group_velocity(vedata[nb]);
     char filename[100];
     snprintf(filename,100,"%s_mode%i",filebase,band_num);
     f.output_mode_fields(vedata[nb], &fluxB, fvB, filename);
   };

  /***************************************************************/
  /* compute mode expansion coefficients *************************/
  /***************************************************************/
  std::vector<int> bands(4);
  bands[0]=1;
  bands[1]=2;
  bands[2]=3;
  bands[3]=4;

  int num_bands = bands.size();
  int num_freqs = fluxB.Nfreq;
  std::vector<cdouble> coeffs =
   f.get_eigenmode_coefficients(&fluxB, Z, fvB, bands, k_guess, 0);

  cdouble cnorm=0.0;
  for(unsigned n=0; n<bands.size(); n++)
   cnorm+=norm(coeffs[n]);
  cnorm=sqrt(cnorm);

  snprintf(filename,100,"%s.coefficients",filebase);
  ff=fopen(filename,"w");
  printf("Mode  | coefficient \n");
  printf("------------------------------------------------\n");
  for(unsigned n=0; n<bands.size(); n++)
   { cdouble coeff = coeffs[n] / cnorm;
     double cmag=abs(coeff); 
     double cang=atan2(imag(coeff),real(coeff)) * (180.0/M_PI);
     fprintf(ff,"%2i   %2e (@ %2e degrees)\n",n,cmag,cang);
     printf("%2i   %2e (@ %2e degrees)\n",n,cmag,cang);
   };

  return 0;
}
