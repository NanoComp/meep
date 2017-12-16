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

  volume fvA(minA, maxA), fvB(minB, maxB);
  if (point_source)
   { 
     f.add_point_source(Ez, csrc, vec(0.1, 0.2, -1.0*zcenter) ); 
   }
  else
   {
     vec kpoint=k_guess(0, fcen, band_num);
     bool match_frequency = true;
     int parity = 0; // NO_PARITY
     double eigensolver_tol=1.0e-3;
     f.add_eigenmode_source(Dielectric, csrc,
                            Z, fvA, fvA, band_num,
                            kpoint, match_frequency, parity,
                            resolution, eigensolver_tol, 1.0);
   };

  /***************************************************************/
  /* if the --CW option was specified, just export the           */
  /* instantaneous fields on the upper and lower flux planes     */
  /***************************************************************/
  if (CW)
   { while(f.round_time() < 10.0)
      f.step();
     f.output_hdf5(Ex, fvA, 0, false, false, "A");
     f.output_hdf5(Ey, fvA, 0, false, false, "A");
     f.output_hdf5(Hx, fvA, 0, false, false, "A");
     f.output_hdf5(Hy, fvA, 0, false, false, "A");
     f.output_hdf5(Ex, fvB, 0, false, false, "B");
     f.output_hdf5(Ey, fvB, 0, false, false, "B");
     f.output_hdf5(Hx, fvB, 0, false, false, "B");
     f.output_hdf5(Hy, fvB, 0, false, false, "B");
     return 0;
   };

  /***************************************************************/
  /* add flux plane and timestep to accumulate frequency-domain  */
  /* fields at nfreq frequencies                                 */
  /***************************************************************/
  dft_flux flux=f.add_dft_flux_plane(fvB, fcen, fcen, nfreq);

  // (run-sources+ 10
  while( f.round_time() < (f.last_source_time() + 10.0) )
   f.step();

  f.output_flux_fields(&flux, fvB, "flux");

/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
  FILE *ff=fopen("duct-junction.dat","w");
  cdouble mode_flux[4];
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/

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
     snprintf(filename,100,"mode%i",band_num);
     f.output_mode_fields(vedata[nb], &flux, fvB, filename);
     
     cdouble mm[8], mf[8];

     setenv("MEEP_OVERLAP_ALGORITHM","x",1);
     mm[0]=f.get_mode_mode_overlap(vedata[nb],vedata[nb], &flux, fvB);
     mf[0]=f.get_mode_flux_overlap(vedata[nb],&flux,0,fvB);
     setenv("MEEP_OVERLAP_ALGORITHM","5",1);
     mm[1]=f.get_mode_mode_overlap(vedata[nb],vedata[nb], &flux, fvB);
     mf[1]=f.get_mode_flux_overlap(vedata[nb],&flux,0,fvB);
     setenv("MEEP_OVERLAP_ALGORITHM","4",1);
     mm[2]=f.get_mode_mode_overlap(vedata[nb],vedata[nb], &flux, fvB);
     mf[2]=f.get_mode_flux_overlap(vedata[nb],&flux,0,fvB);
     setenv("MEEP_OVERLAP_ALGORITHM","0",1);
     mm[3]=f.get_mode_mode_overlap(vedata[nb],vedata[nb], &flux, fvB);
     mf[3]=f.get_mode_flux_overlap(vedata[nb],&flux,0,fvB);
     setenv("MEEP_OVERLAP_ALGORITHM","1",1);
     mm[4]=f.get_mode_mode_overlap(vedata[nb],vedata[nb], &flux, fvB);
     mf[4]=f.get_mode_flux_overlap(vedata[nb],&flux,0,fvB);
     setenv("MEEP_OVERLAP_ALGORITHM","2",1);
     mm[5]=f.get_mode_mode_overlap(vedata[nb],vedata[nb], &flux, fvB);
     mf[5]=f.get_mode_flux_overlap(vedata[nb],&flux,0,fvB);
     setenv("MEEP_OVERLAP_ALGORITHM","3",1);
     mm[6]=f.get_mode_mode_overlap(vedata[nb],vedata[nb], &flux, fvB);
     mf[6]=f.get_mode_flux_overlap(vedata[nb],&flux,0,fvB);

#if 0
     char s[3]="0";
     int ns=1;
     for(s[0]='0'; s[0]<='6'; s[0]++)
      { setenv("MEEP_OVERLAP_ALGORITHM",s,1);
        mm[ns]=f.get_mode_mode_overlap(vedata[nb],vedata[nb], &flux, fvB);
        mf[ns]=f.get_mode_flux_overlap(vedata[nb],&flux,0,fvB);
        ns++;
      };
#endif

     setenv("MEEP_OVERLAP_ALGORITHM","x",1);
     fprintf(ff,"mode %i, gv=%e\n",band_num,get_group_velocity(vedata[nb]));
     for(int n=0; n<7; n++)
      fprintf(ff,"%i    mm={%+.4e,%+.4e}  mf={%+.4e,%+.4e}\n",n,real(mm[n]),imag(mm[n]),real(mf[n]),imag(mf[n]));
     fprintf(ff,"\n\n");
     fflush(ff);
   };

  cdouble mode_mode[4][4];
  for(int nb=0; nb<4; nb++)
   for(int nbb=0; nbb<4; nbb++)
    mode_mode[nb][nbb]=f.get_mode_mode_overlap(vedata[nb],vedata[nbb], &flux, fvB);

  fprintf(ff,"\n\n");
  for(int nb=0; nb<4; nb++)
   fprintf(ff,"  %.2f   %.2f   %.2f   %.2f\n",
            abs(mode_mode[nb][0]/vgrp[nb]),
            abs(mode_mode[nb][1]/vgrp[nb]),  
            abs(mode_mode[nb][2]/vgrp[nb]),  
            abs(mode_mode[nb][3]/vgrp[nb]));
  fclose(ff);

  /***************************************************************/
  /* compute mode expansion coefficients *************************/
  /***************************************************************/
  std::vector<int> bands(3);
  bands[0]=1;
  bands[1]=2;
  bands[2]=3;
  //bands[3]=4;

  int num_bands = bands.size();
  int num_freqs = flux.Nfreq;

 // std::vector<cdouble> coeffs =
 //  f.get_eigenmode_coefficients(&flux, Z, fvB, bands, k_guess, 0);

#if 0
  for(int nb=0; nb<num_bands; nb++)
   { 
     int band_num=bands[nb];
     cdouble components1[4], components2[4];
     double vgrp;
     //cdouble flux_mode = f.get_eigenmode_coefficient(&flux, 0, Z, fvB, band_num, k_guess, 0, &vgrp);

     if (am_master())
      { FILE *ff=fopen("duct-junction.log",nb==0 ? "w" : "a");
        fprintf(ff,"mode %i (vgrp=%e): \n",band_num,vgrp);
        fprintf(ff," flux-mode {%+e,%+e}\n", real(global_flux_dot_mode_components[0]), imag(global_flux_dot_mode_components[0]));
        fprintf(ff,"           {%+e,%+e}\n", real(global_flux_dot_mode_components[1]), imag(global_flux_dot_mode_components[1]));
        fprintf(ff,"           {%+e,%+e}\n", real(global_flux_dot_mode_components[2]), imag(global_flux_dot_mode_components[2]));
        fprintf(ff,"           {%+e,%+e}\n", real(global_flux_dot_mode_components[3]), imag(global_flux_dot_mode_components[3]));
        fprintf(ff,"-----------------------------------------\n");
        fprintf(ff,"           {%+e,%+e}\n", real(flux_mode), imag(flux_mode));
        fprintf(ff,"\n");
        fflush(ff);

        for(int nbp=nb; nbp<num_bands; nbp++)
         { int band_nump=bands[nbp]; 
           cdouble mode_mode_components[4];
           printf("Getting mode_mode overlap(%i,%i)\n",band_nump,band_num);
           cdouble mode_mode=get_mode_mode_overlap(&f, &flux, fvB, fcen, band_nump, band_num, k_guess, 0, mode_mode_components);
   
           fprintf(ff,"mm%i-%i:   {%+e,%+e}\n", band_nump, band_num, real(mode_mode_components[0]), imag(mode_mode_components[0]));
           fprintf(ff,"           {%+e,%+e}\n", real(mode_mode_components[1]), imag(mode_mode_components[1]));
           fprintf(ff,"           {%+e,%+e}\n", real(mode_mode_components[2]), imag(mode_mode_components[2]));
           fprintf(ff,"           {%+e,%+e}\n", real(mode_mode_components[3]), imag(mode_mode_components[3]));
           fprintf(ff,"-----------------------------------------\n");
           fprintf(ff,"           {%+e,%+e}\n", real(mode_mode), imag(mode_mode));
           fprintf(ff,"\n");
           fflush(ff);
         };
        fclose(ff);
      };
   };
#endif

#if 0   
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
#endif

  return 0;

}
