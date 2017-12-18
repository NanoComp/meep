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

bool equal_float(double x, double y) { return (float)x == (float)y; }

vec k_guess(void *user_data, double freq, int band_num)
{
  (void) user_data;
  (void) freq;
  (void) band_num;

  // hard-coded dispersion relation
  if ( equal_float(freq,0.25) )
   { if (band_num==1) return vec(0.0, 0.0, 0.593618);
     if (band_num==2) return vec(0.0, 0.0, 0.577254);
     if (band_num==3) return vec(0.0, 0.0, 0.494127);
     if (band_num==4) return vec(0.0, 0.0, 0.416103);
   }
  else if ( equal_float(freq,0.5) )
   { if (band_num==1) return vec(0.0, 0.0, 1.23024);
     if (band_num==2) return vec(0.0, 0.0, 1.22527);
     if (band_num==3) return vec(0.0, 0.0, 1.16948);
     if (band_num==4) return vec(0.0, 0.0, 0.507376);
   }
 

  return vec(0.0, 0.0, 0.736917);
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
  double ratio      = 1.0;
  bool CW           = false;
  char *filebase    = "pj";
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
  /****************** *********************************************/
  double ncore=2.5;   // index of waveguide core
  double H=6.0;       // half-height of computational cell (Y direction)
  double L=10.0;      // half-length of computational cell (Z direction)
  double h1=1.0;      // half-height of narrow waveguide
  double h2=h1*ratio; // half-height of wider waveguide
  double dpml=1.0;    // PML thickness
  double a= 10.0;     // resolution

  geometry_lattice.size.x=0.0;
  geometry_lattice.size.y=2.0*H;
  geometry_lattice.size.z=2.0*L;
  grid_volume gv=vol3d(0, 2.0*H, 2.0*L, a);
  gv.center_origin();
  symmetry sym = use_symmetry ? -mirror(Y,gv) : identity();
  structure the_structure(gv, dummy_eps, pml(dpml), sym);

  meep_geom::material_type dielectric = meep_geom::make_dielectric(ncore*ncore);
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
  double fcen = 0.25;  // center frequency
  double df   = 0.25;  // bandwidth
  int nfreq   = 1;     // number of frequency points
  gaussian_src_time gsrc(fcen, df);
  continuous_src_time csrc(fcen);

  volume fvA( vec(0,-H,-0.5*L), vec(0,H,-0.5*L) );
  volume fvB( vec(0,-H,+0.5*L), vec(0,H,+0.5*L) );

  double vol=2.0*H; // 1-dimensional "volume" of flux plane

  volume fvC( vec(0,-H,-L),     vec(0,H,L) );

  if (point_source)
   { 
     if (CW)
      f.add_point_source(Ex, csrc, vec(0.0, 0.0, -0.5*L));
     else
      f.add_point_source(Ex, gsrc, vec(0.0, 0.0, -0.5*L));
   }
  else
   { 
     bool match_frequency = true;
     int parity = 0; // NO_PARITY
     double resolution=f.a;
     double eigensolver_tol=1.0e-4;
     vec kpoint=k_guess(0, fcen, band_num);
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
  /* add flux planes                                             */
  /***************************************************************/
  dft_flux fluxA=f.add_dft_flux_plane(fvA, fcen-0.5*df, fcen+0.5*df, nfreq);
  dft_flux fluxB=f.add_dft_flux_plane(fvB, fcen-0.5*df, fcen+0.5*df, nfreq);
  dft_flux fluxC=f.add_dft_flux_plane(fvC, fcen-0.5*df, fcen+0.5*df, nfreq);

  /***************************************************************/
  /* timestep until fields decayed + 50 time unit ****************/
  /***************************************************************/
  while( f.round_time() < (f.last_source_time() + 1.0/df) )
   f.step();
 
  /***************************************************************/
  /* write output files ******************************************/
  /***************************************************************/
  char filename[100];
  snprintf(filename,100,"%s_fluxA",filebase);
  f.output_flux_fields(&fluxA, fvA, filename);
  snprintf(filename,100,"%s_fluxB",filebase);
  f.output_flux_fields(&fluxB, fvB, filename);
  snprintf(filename,100,"%s_fluxC",filebase);
  f.output_flux_fields(&fluxC, fvC, filename);

  f.output_hdf5(Ex,f.total_volume());
  f.output_hdf5(Ey,f.total_volume());
  f.output_hdf5(Hx,f.total_volume());
  f.output_hdf5(Hy,f.total_volume());

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  void *vedata[4];
  double vgrp[4];
  snprintf(filename,100,"%s.dat",filebase);
  FILE *ff=fopen(filename,"w");
  for(int nb=0; nb<4; nb++)
   { int band_num=nb+1;
     vedata[nb]=f.get_eigenmode(fcen, Z, fvB, fvB,
                                band_num, k_guess(0,fcen,band_num),
                                true, 0, f.a, 1.0e-4);
     vgrp[nb]=get_group_velocity(vedata[nb]);
     snprintf(filename,100,"%s_mode%i",filebase,band_num);
     f.output_mode_fields(vedata[nb], &fluxB, fvB, filename);
     
     cdouble mm[8], mf[8];

     setenv("MEEP_OVERLAP_ALGORITHM","x",1);
     mm[0]=f.get_mode_mode_overlap(vedata[nb],vedata[nb], &fluxB, fvB);
     mf[0]=f.get_mode_flux_overlap(vedata[nb],&fluxB,0,fvB);

     setenv("MEEP_OVERLAP_ALGORITHM","5",1);
     mm[1]=f.get_mode_mode_overlap(vedata[nb],vedata[nb], &fluxB, fvB);
     mf[1]=f.get_mode_flux_overlap(vedata[nb],&fluxB,0,fvB);

     setenv("MEEP_OVERLAP_ALGORITHM","4",1);
     mm[2]=f.get_mode_mode_overlap(vedata[nb],vedata[nb], &fluxB, fvB);
     mf[2]=f.get_mode_flux_overlap(vedata[nb],&fluxB,0,fvB);

     setenv("MEEP_OVERLAP_ALGORITHM","0",1);
     mm[3]=f.get_mode_mode_overlap(vedata[nb],vedata[nb], &fluxB, fvB);
     mf[3]=f.get_mode_flux_overlap(vedata[nb],&fluxB,0,fvB);

     setenv("MEEP_OVERLAP_ALGORITHM","1",1);
     mm[4]=f.get_mode_mode_overlap(vedata[nb],vedata[nb], &fluxB, fvB);
     mf[4]=f.get_mode_flux_overlap(vedata[nb],&fluxB,0,fvB);

     setenv("MEEP_OVERLAP_ALGORITHM","2",1);
     mm[5]=f.get_mode_mode_overlap(vedata[nb],vedata[nb], &fluxB, fvB);
     mf[5]=f.get_mode_flux_overlap(vedata[nb],&fluxB,0,fvB);

     setenv("MEEP_OVERLAP_ALGORITHM","3",1);
     mm[6]=f.get_mode_mode_overlap(vedata[nb],vedata[nb], &fluxB, fvB);
     mf[6]=f.get_mode_flux_overlap(vedata[nb],&fluxB,0,fvB);

     double normfac=vol*vgrp[nb];
     for(int n=0; n<7; n++)
      { mm[n] /= normfac;
        mf[n] /= normfac;
      };

     setenv("MEEP_OVERLAP_ALGORITHM","x",1);
     fprintf(ff,"mode %i, gv=%e\n",band_num,vgrp[nb]);
     for(int n=0; n<7; n++)
      fprintf(ff,"%i    mm=%.2e {%+.4e,%+.4e}  mf=%.2e {%+.4e,%+.4e}\n",n,
                        abs(mm[n]),real(mm[n]),imag(mm[n]),
                        abs(mf[n]),real(mf[n]),imag(mf[n]));
     fprintf(ff,"\n\n");
     fflush(ff);
   };

  cdouble mode_mode[4][4];
  for(int nb=0; nb<4; nb++)
   for(int nbb=0; nbb<4; nbb++)
    mode_mode[nb][nbb]=f.get_mode_mode_overlap(vedata[nb],vedata[nbb], &fluxB, fvB);

  fprintf(ff,"\n\n");
  for(int nb=0; nb<4; nb++)
   fprintf(ff,"  %.2f   %.2f   %.2f   %.2f\n",
            abs(mode_mode[nb][0]/(vol*vgrp[nb])),
            abs(mode_mode[nb][1]/(vol*vgrp[nb])),
            abs(mode_mode[nb][2]/(vol*vgrp[nb])),
            abs(mode_mode[nb][3]/(vol*vgrp[nb])));
  fclose(ff);

  /***************************************************************/
  /* compute mode-expansion coefficients *************************/
  /***************************************************************/
  std::vector<int> bands(8);
  bands[0]=1;
  bands[1]=2;
  bands[2]=3;
  bands[3]=4;
  bands[4]=5;
  bands[5]=6;
  bands[6]=7;
  bands[7]=8;

  int num_bands = bands.size();
  int num_freqs = fluxB.Nfreq;
  std::vector<cdouble> coeffs =
   f.get_eigenmode_coefficients(&fluxB, Z, fvB, bands, k_guess, 0);

  cdouble cnorm=0.0;
  for(unsigned n=0; n<bands.size(); n++)
   cnorm+=abs(coeffs[n])*abs(coeffs[n]);
  cnorm=sqrt(cnorm);

  snprintf(filename,100,"%s.coefficients",filebase);
  ff=fopen(filename,"w");
  printf("Mode  | coefficient \n");
  printf("------------------------------------------------\n");
  for(unsigned n=0; n<bands.size(); n++)
   { cdouble coeff = coeffs[n] / cnorm;
     double cmag=abs(coeff); 
     double cang=atan2(imag(coeff),real(coeff)) * (180.0/M_PI);
     fprintf(ff,"%2i   %6f (@ %2e degrees)\n",n,cmag,cang);
     printf("%2i   %6f (@ %2e degrees)\n",n,cmag,cang);
   };

}
