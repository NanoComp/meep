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

/***************************************************************/
/* function to provide initial guess for the wavevector of the */
/* eigenmode with frequency freq and band index band           */
/***************************************************************/
vec k_guess(void *user_data, double freq, int band_num)
{
  (void) user_data;
  (void) freq;
  (void) band_num;
 
  return vec(0.0, 0.0, 0.303278);
} 
vector3 v3(double x, double y=0.0, double z=0.0)
{
  (void) user_data;
  (void) freq;
  (void) band_num;
 
  return vec(0.0, 0.0, 0.303278);
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
      { //if ( ++narg >= argc )
        // usage(argv[0], "error: no argument given for --band-num");
        sscanf(argv[narg], "%i", &band_num);
        master_printf("setting band-num=%i\n",band_num);
      }
     else if (!strcasecmp(argv[narg],"--ratio"))
      { // if ( ++narg >= argc )
        // usage(argv[0], "error: no argument given for --ratio");
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

  /***************************************************************/
  /* set up geometry: cylinder of radius r in a 2D computational */
  /* cell of size LxL (z-invariant)                              */

  double n=3.0;     // index of waveguide
  double r=1.0;     // cylinder radius
  double L=5.0;     // size of computational cell
  double dpml=1;    // thickness of PML

  double resolution = 10.0; 

  geometry_lattice.size.x=L;
  geometry_lattice.size.y=L;
  geometry_lattice.size.z=0.0;
  grid_volume gv = voltwo(L, L, resolution);
  gv.center_origin();

  symmetry sym=identity();

  structure the_structure(gv, dummy_eps, pml(dpml), sym);

  meep_geom::material_type dielectric = meep_geom::make_dielectric(n*n);
  geometric_object objects[1];
  vector3 v3zero = {0.0,0.0,0.0};
  vector3 zaxis  = {0.0,0.0,1.0};
  objects[0] = make_cylinder(dielectric, v3zero, r, ENORMOUS, zaxis);
  geometric_object_list g={ 1, objects }; 
  meep_geom::set_materials_from_geometry(&the_structure, g);
  fields f(&the_structure);
  f.output_hdf5(Dielectric,f.total_volume());

  /***************************************************************/
  /* add sources: point source (if --point_source option present */
  /* or eigenmode source of band band_num                        */
  /***************************************************************/

  f.output_hdf5(Dielectric,f.total_volume());

  /***************************************************************/
  /* add sources: point source (if --point_source option present */
  /* or eigenmode source of band band_num                        */
  /***************************************************************/
  double fcen = 0.15;  // ; pulse center frequency
  double df   = 0.1;   // ; df
  int nfreq   = 10;
nfreq=1;
  int nfreq   = 1;
  gaussian_src_time src(fcen, df);
  volume fv = volume( vec(-0.5*L, -0.5*L), vec(+0.5*L, +0.5*L));
  if (point_source)
   { 
     f.add_point_source(Ez, src, vec(0.1, 0.2) ); }
     f.add_point_source(Ez, src, vec(0.1, 0.2) );
   }
  else
   {
     vec kpoint(0,0,0.303278);
     bool match_frequency = true;
     int parity = 0; // NO_PARITY
     double eigensolver_tol=1.0e-7;
     f.add_eigenmode_source(Dielectric, src, Z, fv, fv, band_num,
                            kpoint, match_frequency, parity, 
                            resolution, eigensolver_tol, 1.0);
   };

  /***************************************************************/
  /* add flux plane and timestep to accumulate frequency-domain  */
  /* fields at nfreq frequencies                                 */
  /***************************************************************/
  dft_flux flux=f.add_dft_flux(Z, fv, fcen-0.5*df, fcen+0.5*df, nfreq);
<<<<<<< HEAD

  // (run-sources+ 100 
  while( f.round_time() < (f.last_source_time() + 100.0) )
   f.step();
             
  /***************************************************************/
  /* compute mode expansion coefficients *************************/
  /***************************************************************/
  std::vector<int> bands(2);
  bands[0]=1;
  bands[1]=2;

  int num_bands = bands.size();
  int num_freqs = flux.Nfreq;

  std::vector<cdouble> coeffs =
   f.get_eigenmode_coefficients(&flux, Z, fv, bands, k_guess, 0);
   
  if (am_master())
   {
     FILE *ff=fopen("fiber-junction.out","a");
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
       };
  //dft_flux flux=f.add_dft_flux(Z, fv, fcen-0.5*df, fcen+0.5*df, nfreq);
  dft_flux flux=f.add_dft_flux(Z, fv, fcen, fcen, 1);

  // (run-sources+ 100 
  while( f.round_time() < (f.last_source_time() + 100.0) )
   f.step();
             
  /***************************************************************/
  /* compute mode expansion coefficients *************************/
  /***************************************************************/
  std::vector<int> bands(2);
  bands[0]=1;
  bands[1]=2;

  int num_bands = bands.size();
  int num_freqs = flux.Nfreq;

  std::vector<cdouble> coeffs =
   f.get_eigenmode_coefficients(&flux, Z, fv, bands, k_guess, 0);
   
  if (am_master())
   {
     FILE *ff=fopen("fiber-junction.out","a");
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
       };
  //dft_flux flux=f.add_dft_flux(Z, fv, fcen-0.5*df, fcen+0.5*df, nfreq);
  dft_flux flux=f.add_dft_flux(Z, fv, fcen, fcen, 1);

  // (run-sources+ 100 
  while( f.round_time() < (f.last_source_time() + 100.0) )
   f.step();
             
  /***************************************************************/
  /* compute mode expansion coefficients *************************/
  /***************************************************************/
  std::vector<int> bands(2);
  bands[0]=1;
  bands[1]=2;

  int num_bands = bands.size();
  int num_freqs = flux.Nfreq;

  std::vector<cdouble> coeffs =
   f.get_eigenmode_coefficients(&flux, Z, fv, bands, k_guess, 0);
   
  if (am_master())
   {
     FILE *ff=fopen("fiber-junction.out","a");
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
       };
   };
  
  return 0;

}
