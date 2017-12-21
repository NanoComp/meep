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

/***************************************************************/
/* function to provide initial guess for the wavevector of the */
/* eigenmode with frequency freq and band index band           */
/***************************************************************/
bool equal_float(double x, double y) { return (float)x == (float)y; }

vec k_guess(void *user_data, double freq, int band_num)
{
  (void) freq;
  (void) band_num;
 
  double h=*((double *) user_data);

  // hard-coded dispersion relations for h=1.5
  if ( equal_float(h,1.5) && equal_float(freq,0.15) )
   { if (band_num==1)  return vec(0.0, 0.0, 0.501600);
     if (band_num==2)  return vec(0.0, 0.0, 0.492976);
     if (band_num==3)  return vec(0.0, 0.0, 0.442391);
     if (band_num==4)  return vec(0.0, 0.0, 0.404212);
     if (band_num==5)  return vec(0.0, 0.0, 0.330163);
     if (band_num==6)  return vec(0.0, 0.0, 0.213681);
     if (band_num==7)  return vec(0.0, 0.0, 0.139897);
     if (band_num==8)  return vec(0.0, 0.0, 0.139368);
     if (band_num==9)  return vec(0.0, 0.0, 0.139407);
     if (band_num==10) return vec(0.0, 0.0, 0.139322);
     if (band_num==11) return vec(0.0, 0.0, 0.0316779);
     if (band_num>=12) return vec(0.0, 0.0, 0.0297807);
   }
  else if ( equal_float(h,1.0) && equal_float(freq,0.15) )
   { if (band_num==1)  return vec(0.0, 0.0, 0.425568);
     if (band_num==2)  return vec(0.0, 0.0, 0.263221);
     if (band_num==3)  return vec(0.0, 0.0, 0.133216);
     if (band_num==4)  return vec(0.0, 0.0, 0.131931);
     if (band_num>=5)  return vec(0.0, 0.0, 0.100000);
     if (band_num==6)  return vec(0.0, 0.0, 0.213681);
     if (band_num==7)  return vec(0.0, 0.0, 0.139897);
     if (band_num==8)  return vec(0.0, 0.0, 0.139368);
     if (band_num==9)  return vec(0.0, 0.0, 0.139407);
     if (band_num==10) return vec(0.0, 0.0, 0.139322);
     if (band_num==11) return vec(0.0, 0.0, 0.0316779);
     if (band_num>=12) return vec(0.0, 0.0, 0.0297807);
   }
  else // if ( equal_float(h,0.5) && equal_float(freq,0.15) )
   { if (band_num==1)  return vec(0.0, 0.0, 0.425568);
     if (band_num==2)  return vec(0.0, 0.0, 0.263354);
     if (band_num==3)  return vec(0.0, 0.0, 0.143212);
     if (band_num==4)  return vec(0.0, 0.0, 0.142938);
     if (band_num==5)  return vec(0.0, 0.0, 0.142935);
     if (band_num==6)  return vec(0.0, 0.0, 0.142860);
     if (band_num==7)  return vec(0.0, 0.0, 0.0725943);
     if (band_num>=8)  return vec(0.0, 0.0, 0.0710456);
   };

  return vec(0.0, 0.0, 0.425568);
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
  bool use_symmetry     = false;
  bool point_source     = false;
  int  band_num         = 1;
  int  num_bands        = 4;
  double ratio          = 1.0;
  char *filebase        = const_cast<char *>("pj");
  bool plot_modes       = false;
  bool plot_flux        = false;
  double frame_interval = 0.0;
  for(int narg=1; narg<argc; narg++)
   { if ( argv[narg]==0 )
      continue;
     if (!strcasecmp(argv[narg],"--use-symmetry") )
      { use_symmetry=true;
        master_printf("Using symmetry.\n");
      }
     else if (!strcasecmp(argv[narg],"--plot-modes") )
      plot_modes=true;
     else if (!strcasecmp(argv[narg],"--plot-flux") )
      plot_flux=true;
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
     else if (!strcasecmp(argv[narg],"--num-bands"))
      { if ( ++narg >= argc )
         usage(argv[0], "error: no argument given for --num-bands");
        sscanf(argv[narg], "%i", &num_bands);
        master_printf("setting num-bands=%i\n",num_bands);
      }
     else if (!strcasecmp(argv[narg],"--ratio"))
      { if ( ++narg >= argc )
         usage(argv[0], "error: no argument given for --ratio");
        sscanf(argv[narg], "%le", &ratio);
        master_printf("setting ratio=%e\n",ratio);
      }
     else if (!strcasecmp(argv[narg],"--frame-interval"))
      { if ( ++narg >= argc )
         usage(argv[0], "error: no argument given for --frame-interval");
        sscanf(argv[narg], "%le", &frame_interval);
        master_printf("setting frame-interval=%e\n",frame_interval);
      }
     else if (!strcasecmp(argv[narg],"--filebase"))
      { if ( ++narg >= argc )
         usage(argv[0], "error: no argument given for --filebase");
        filebase=strdup(argv[narg]);
      }
     else
      { master_printf("unknown command-line option %s (aborting)",argv[narg]);
        usage(argv[0]);
      }; 
   };

  /***************************************************************/
  /* initialize geometry, similar to holey_wvg_cavity   **********/
  /****************** ********************************************/
  double epsilon=12.0; // permeability of waveguide core
  double H=6.0;        // half-height of computational cell (Y direction)
  double L=8.0;        // half-length of computational cell (Z direction)
  double hA=0.5;       // half-height of narrow waveguide
  double hB=hA*ratio;  // half-height of wider waveguide
  double dpml=1.0;     // PML thickness
  double a= 10.0;      // resolution
  geometry_lattice.size.x=0.0;
  geometry_lattice.size.y=2.0*H;
  geometry_lattice.size.z=2.0*L;
  grid_volume gv=vol3d(0, 2.0*H, 2.0*L, a);
  gv.center_origin();
  symmetry sym = use_symmetry ? -mirror(Y,gv) : identity();
  structure the_structure(gv, dummy_eps, pml(dpml), sym);

  meep_geom::material_type dielectric = meep_geom::make_dielectric(epsilon);
  vector3 xhat  = {1.0, 0.0, 0.0};
  vector3 yhat  = {0.0, 1.0, 0.0};
  vector3 zhat  = {0.0, 0.0, 1.0};
  geometric_object objects[2];
  objects[0] = make_block( dielectric, 
                           v3(0.0, 0.0, -0.5*L),     // center
                           xhat, yhat, zhat,
                           v3(ENORMOUS, 2.0*hA, L)   // size;
                         );
  objects[1] = make_block( dielectric,
                           v3(0.0, 0.0, +0.5*L),     // center
                           xhat, yhat, zhat,
                           v3(ENORMOUS, 2.0*hB, L)   // size;
                         );
  geometric_object_list g={ 2, objects };
  meep_geom::set_materials_from_geometry(&the_structure, g);
  fields f(&the_structure);
  f.output_hdf5(Dielectric,f.total_volume(),0,false,true,filebase);

  /***************************************************************/
  /* add source                                                  */
  /***************************************************************/
  double fcen = 0.15;  // center frequency
  double df   = 0.15;  // bandwidth
  int nfreq   = 1;     // number of frequency points
  gaussian_src_time gsrc(fcen, df);

  volume fvA( vec(0,-H,-0.5*L), vec(0,H,-0.5*L) );
  volume fvB( vec(0,-H,+0.5*L), vec(0,H,+0.5*L) );
  volume fvC( vec(0,-H,-L),     vec(0,H,L) );

  if (point_source)
   { 
      f.add_point_source(Ex, gsrc, vec(0.0, 0.0, -0.5*L));
   }
  else
   { 
     volume fvA( vec(0,-H,-0.5*L), vec(0,H,-0.5*L) );
     volume fvB( vec(0,-H,+0.5*L), vec(0,H,+0.5*L) );
     volume cs(  vec(0,-H,-L),     vec(0,H,+L) );
     bool match_frequency = true;
     int parity = 0; // NO_PARITY
     double resolution=f.a;
     double eigensolver_tol=1.0e-4;
     vec kpoint=k_guess((void *)&hA, fcen, band_num);
     f.add_eigenmode_source(Dielectric, gsrc,
                            Z, fvA, fvA, band_num,
                            kpoint, match_frequency, parity,
                            resolution, eigensolver_tol, 1.0);
   };

  /***************************************************************/
  /* add flux planes                                             */
  /***************************************************************/
  dft_flux fluxA=f.add_dft_flux_plane(fvA, fcen-0.5*df, fcen+0.5*df, nfreq);
  dft_flux fluxB=f.add_dft_flux_plane(fvB, fcen-0.5*df, fcen+0.5*df, nfreq);
  dft_flux fluxC=f.add_dft_flux_plane(fvC, fcen-0.5*df, fcen+0.5*df, nfreq);

  /***************************************************************/
  /* timestep until fields decayed + 50 time unit ****************/
  /***************************************************************/
  double NextFileTime=f.round_time();
  while( f.round_time() < (f.last_source_time() + 5.0/df) )
   { f.step();
     if (frame_interval>0.0 && f.round_time()>NextFileTime)
      { NextFileTime += frame_interval;
        f.output_hdf5(Hz,f.total_volume(),0,false,false,filebase);
      };
   };
 
  /***************************************************************/
  /* write output files ******************************************/
  /***************************************************************/
  if (plot_flux)
   { 
     char filename[100];
     snprintf(filename,100,"%s_fluxA",filebase);
     f.output_flux_fields(&fluxA, fvA, filename);
     snprintf(filename,100,"%s_fluxB",filebase);
     f.output_flux_fields(&fluxB, fvB, filename);
     snprintf(filename,100,"%s_fluxC",filebase);
     f.output_flux_fields(&fluxC, fvC, filename);
   };

  if (plot_modes)
   { 
     void **vedata = (void **)malloc(num_bands * sizeof(void *));
     double *vgrp  = new double[num_bands];
     char filename[100];
     for(int nb=0; nb<num_bands; nb++)
      { int band_num=nb+1;
        vedata[nb]=f.get_eigenmode(fcen, Z, fvB, fvB,
                                   band_num, k_guess((void *)&hB,fcen,band_num),
                                   true, 0, f.a, 1.0e-4);
        vgrp[nb]=get_group_velocity(vedata[nb]);
        snprintf(filename,100,"%s_mode%i",filebase,band_num);
        f.output_mode_fields(vedata[nb], &fluxB, fvB, filename);
      };
   
     FILE *ff1=0, *ff2=0, *ff3=0, *ff4=0, *ff5=0;
     if (am_master()) 
      { snprintf(filename,100,"%s_exh_mag.dat",filebase);
        ff1=fopen(filename,"w");
        snprintf(filename,100,"%s_exh_reim.dat",filebase);
        ff2=fopen(filename,"w");
        snprintf(filename,100,"%s_hxe_mag.dat",filebase);
        ff3=fopen(filename,"w");
        snprintf(filename,100,"%s_hxe_reim.dat",filebase);
        ff4=fopen(filename,"w");
        snprintf(filename,100,"%s.modeData",filebase);
        ff5=fopen(filename,"w");
      };
     double vol=2.0*H; // 1-dimensional "volume" of flux plane
     for(int nb=0; nb<num_bands; nb++)
      for(int nbb=0; nbb<num_bands; nbb++)
       { cdouble mm[2];
         f.get_mode_mode_overlap(vedata[nb],vedata[nbb], &fluxB, fvB, mm);
         double normfac=vol*vgrp[nb];
         mm[0]/=normfac;
         mm[1]/=normfac;
         if (am_master())
          { char c= ( (nbb==num_bands-1) ? '\n' : ' ');
            fprintf(ff1,"%.4f  %c",abs(mm[0]),c);
            fprintf(ff2,"{%+.4f,%+.4f}   %c",real(mm[0]),imag(mm[0]),c);
            fprintf(ff3,"%.4f  %c",abs(mm[1]),c);
            fprintf(ff4,"{%+.4f,%+.4f}   %c",real(mm[1]),imag(mm[1]),c);
            if(nbb==0)
             { double vgrp=get_group_velocity(vedata[nb]);
               vec k=get_k(vedata[nb]);
               fprintf(ff5,"nb=%2i   vgrp=%e   ",nb,vgrp);
               fprintf(ff5,"kpoint={%e,%e,%e}\n",k.x(),k.y(),k.z());
             };
          };
       }; 
     if (am_master()) 
      { fclose(ff1);
        fclose(ff2);
        fclose(ff3);
        fclose(ff4);
        fclose(ff5);
      };
   };

  /***************************************************************/
  /* compute mode-expansion coefficients *************************/
  /***************************************************************/
  std::vector<int> bands(num_bands);
  for(int n=0; n<num_bands; n++)
   bands[n] = n+1;

  int num_freqs = fluxB.Nfreq;
  std::vector<cdouble> coeffs =
   f.get_eigenmode_coefficients(&fluxB, Z, fvB, bands, k_guess, (void *)&hB);

  if (am_master())
   { 
     char filename[100];
     snprintf(filename,100,"%s.coefficients",filebase);
     FILE *ff=fopen(filename,"w");
     printf("freq | band | alpha^+ | alpha^-\n");
     printf("------------------------------------------------\n");
     for(int nf=0; nf<num_freqs; nf++)
      for(unsigned nb=0; nb<bands.size(); nb++)
       { 
         double anorm=0.0;
         for(int nbb=0; nbb<num_bands; nbb++)
          for(int pm=0; pm<2; pm++)
           anorm += norm( coeffs[2*nbb*num_freqs + 2*nf + pm] );
         anorm=sqrt(anorm);
   
         cdouble aP = coeffs[2*nb*num_freqs + 2*nf + 0];
         cdouble aM = coeffs[2*nb*num_freqs + 2*nf + 1];
         printf("%2i  %2i  (+)  %e {%+e,%+e} (%e %%)\n",nf,nb,abs(aP),real(aP),imag(aP),100.0*abs(aP/anorm));
         printf("%2i  %2i  (-)  %e {%+e,%+e} (%e %%)\n",nf,nb,abs(aM),real(aM),imag(aM),100.0*abs(aM/anorm));
         fprintf(ff,"%2i  %2i   %e {%+e,%+e} (%e %%)  %e  {%+e, %+e} (%e %%)\n",nf,nb,
                     abs(aP),real(aP),imag(aP),abs(aP/anorm),
                     abs(aM),real(aM),imag(aM),abs(aM/anorm));
      };
     fclose(ff);
   };
   
  return 0;
}
