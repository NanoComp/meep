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

#include "meepgeom.hpp"

using namespace meep;

typedef std::complex<double> cdouble;

namespace meep{
void output_hdf5_flux(fields *f, dft_flux *flux, const volume where,
                      const char *HDF5FileName,
                      bool retain_integration_weight=false);
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void nudge_onto_dielectric_grid(grid_volume gv, vec &v)
{ 
  ivec iv = gv.round_vec(v);
  if ( (iv.x() % 2) == 0 )
   iv += ivec(1,0,0);
  if ( (iv.y() % 2) == 0 )
   iv += ivec(0,1,0);
  if ( (iv.z() % 2) == 0 )
   iv += ivec(0,0,1);

  v=gv[iv];
}

/***************************************************************/
/* function to provide initial guess for the wavevector of the */
/* eigenmode with frequency freq and band index band           */
/***************************************************************/
vec k_guess(void *user_data, double freq, int band_num)
{
  (void) user_data;
  (void) freq;
  (void) band_num;

// hard-coded dispersion relations
  if (band_num==0)
   return vec(0.0, 0.0, 3.8978*freq - 0.273203);
  else 
   return vec(0.0, 0.0, 4.23559*freq - 0.44189);
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
  double ratio      = 2.0;
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

  double WXY = 10.0;     // size of computational cell in transverse (XY) directions
  double WZ  = 20.0;     // size of computational cell in propagation direction

  double dpml=0.5;       // thickness of PML

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
  gaussian_src_time src(fcen, df);

  vec minA(-0.5*WXY, -0.5*WXY, -1.0*zcenter);
  vec maxA(+0.5*WXY, +0.5*WXY, -1.0*zcenter);
  vec minB(-0.5*WXY, -0.5*WXY, +1.0*zcenter);
  vec maxB(+0.5*WXY, +0.5*WXY, +1.0*zcenter);

/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
printf("before: min,max A={%+f, %+f, %+f} {%+f, %+f, %+f}\n",
        minA.x(), minA.y(), minA.z(),
        maxA.x(), maxA.y(), maxA.z());
printf("before: min,max B={%+f, %+f, %+f} {%+f, %+f, %+f}\n",
        minB.x(), minB.y(), minB.z(),
        maxB.x(), maxB.y(), maxB.z());
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
  nudge_onto_dielectric_grid(gv, minA);
  nudge_onto_dielectric_grid(gv, maxA);
  nudge_onto_dielectric_grid(gv, minB);
  nudge_onto_dielectric_grid(gv, maxB);
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
printf("after: min,max A={%+f, %+f, %+f} {%+f, %+f, %+f}\n",
        minA.x(), minA.y(), minA.z(),
        maxA.x(), maxA.y(), maxA.z());
printf("after: min,max B={%+f, %+f, %+f} {%+f, %+f, %+f}\n",
        minB.x(), minB.y(), minB.z(),
        maxB.x(), maxB.y(), maxB.z());
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/

  volume fvA(minA, maxA), fvB(minB, maxB);
  if (point_source)
   { 
     f.add_point_source(Ez, src, vec(0.1, 0.2, -1.0*zcenter) ); 
   }
  else
   {
     vec kpoint=k_guess(0, fcen, band_num);
     bool match_frequency = true;
     int parity = 0; // NO_PARITY
     double eigensolver_tol=1.0e-4;
     f.add_eigenmode_source(Dielectric, src, Z, fvA, fvA, band_num,
                            kpoint, match_frequency, parity, 
                            resolution, eigensolver_tol, 1.0);
   };

  /***************************************************************/
  /* add flux plane and timestep to accumulate frequency-domain  */
  /* fields at nfreq frequencies                                 */
  /***************************************************************/
  dft_flux flux=f.add_dft_flux_plane(fvB, fcen, fcen, nfreq);

  // (run-sources+ 100
  while( f.round_time() < (f.last_source_time() + 10.0) )
   f.step();

  f.output_hdf5_flux(&flux, fvB, "flux");
             
  /***************************************************************/
  /* compute mode expansion coefficients *************************/
  /***************************************************************/
  std::vector<int> bands(4);
  bands[0]=1;
  bands[1]=2;
  bands[2]=3;
  bands[3]=4;

  int num_bands = bands.size();
  int num_freqs = flux.Nfreq;

  std::vector<cdouble> coeffs =
   f.get_eigenmode_coefficients(&flux, Z, fvB, bands, k_guess, 0);
   
  if (am_master())
   {
     FILE *ff=fopen("duct-junction.out","a");
     fprintf(ff,"\n\n");
     for(int nf=0; nf<num_freqs; nf++)
      for(int nb=0; nb<num_bands; nb++)
       { 
         cdouble alpha = coeffs[nf*num_bands + nb];
         double alphaMag = abs(alpha);
         double alphaArg = atan2(imag(alpha),real(alpha))*180.0/M_PI;

         printf("(nf,f,nb)=(%2i,%4f,%i): {%.2e,%.2e} (%.4f @ %3.0f)\n",
                 nf,flux.freq_min + nf*flux.dfreq,bands[nb],
                 real(alpha), imag(alpha), alphaMag, alphaArg);

         fprintf(ff,"(nf,f,nb)=(%2i,%4f,%i): {%.2e,%.2e} (%.4f @ %3.0f)\n",
                     nf,flux.freq_min + nf*flux.dfreq,bands[nb],
                     real(alpha), imag(alpha), alphaMag, alphaArg);
       };
     fclose(ff);
   };
  
  return 0;

}
