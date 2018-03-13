/***************************************************************/
/* demonstration of mode expansion in a strip-waveguide taper  */
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

extern std::vector<double> mode_group_velocities;

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
/* user-supplied routine for estimating dispersion relations   */
/* to accelerate MPB calculations                              */
/***************************************************************/
bool equal_float(double x, double y) { return (float)x == (float)y; }

vec k_guess(void *user_data, double freq, int band_num)
{
  double w     = *((double *)user_data);

  if ( equal_float(freq,0.41) && equal_float(w,0.5))
   { if (band_num==1)  return vec(0.9892331, 0.0, 0.0);
     if (band_num==2)  return vec(0.6175083, 0.0, 0.0);
     if (band_num==3)  return vec(0.5469879, 0.0, 0.0);
     if (band_num==4)  return vec(0.5245156, 0.0, 0.0);
     if (band_num==5)  return vec(0.4267270, 0.0, 0.0);
     if (band_num>=6)  return vec(0.4245740, 0.0, 0.0);
   }
  else if ( equal_float(freq,0.41) && equal_float(w,1.0))
   { if (band_num==1)  return vec(1.073627,  0.0, 0.0);
     if (band_num==2)  return vec(0.9856316, 0.0, 0.0);
     if (band_num==3)  return vec(0.8233921, 0.0, 0.0);
     if (band_num==4)  return vec(0.5844210, 0.0, 0.0);
     if (band_num==5)  return vec(0.5497692, 0.0, 0.0);
     if (band_num>=6)  return vec(0.5116020, 0.0, 0.0);
   }
  else if ( equal_float(freq,0.15) && equal_float(w,3.0))
   { if (band_num==1)  return vec(0.494476, 0.0, 0.0);
     if (band_num==2)  return vec(0.486399, 0.0, 0.0);
     if (band_num==3)  return vec(0.435861, 0.0, 0.0);
     if (band_num==4)  return vec(0.397068, 0.0, 0.0);
     if (band_num==5)  return vec(0.322812, 0.0, 0.0);
     if (band_num>=6)  return vec(0.211186, 0.0, 0.0);
   }
  else if ( equal_float(freq,0.15) && equal_float(w,2.0))
   { if (band_num==1)  return vec(0.426302,  0.0, 0.0);
     if (band_num==2)  return vec(0.4534589, 0.0, 0.0);
     if (band_num==3)  return vec(0.3446421, 0.0, 0.0);
     if (band_num==4)  return vec(0.1860106, 0.0, 0.0);
     if (band_num>=5)  return vec(0.1475703, 0.0, 0.0);
   }
  else // if ( equal_float(freq,0.15) && equal_float(w,1.0))
   { if (band_num==1)  return vec(0.419984, 0.0, 0.0);
     if (band_num>=2)  return vec(0.426302, 0.0, 0.0);
   };

  return vec(0.0, 0.0, 0.0);
}

/***************************************************************/
/* user-defined material function for 2D/3D waveguide geometry */
/***************************************************************/
typedef struct wvg_data
 { double wA, wB;        // width of left, right waveguides
   double taper_length;  // taper length (=0 for no taper)
   int taper_order;      // index of first discontinuous derivative - 1 (0,1,2)
   double eps_wvg;       // permittivity of waveguide material
   double eps_ambient;   // permittivity of surrounding medium
 } wvg_data;

void wvg_material(vector3 loc, void *user_data, meep_geom::medium_struct *m)
{
  wvg_data *wdata=(wvg_data *)user_data;

  double eps = wdata->eps_ambient; // assume we are in ambient medium
  double x0 = loc.x / wdata->taper_length;
  double wA  = wdata->wA, wB = wdata->wB, w;
  if (x0 <= -0.5 )
   w = wA;
  else if ( x0 >= 0.5)
   w = wB;
  else if (wdata->taper_order==2)
   w = 0.5*(wA+wB) + (wB-wA)*x0*(15.0 + x0*x0*(-40.0 + x0*x0*48.0))/8.0;
  else if (wdata->taper_order==1)
   w = 0.5*(wA+wB) + (wB-wA)*x0*(1.5 - 2.0*x0*x0);
  else // p=0, i.e. linear taper
   w = 0.5*(wA+wB) + (wB-wA)*x0;

  eps = ( fabs(loc.y)<=0.5*w ) ? wdata->eps_wvg : wdata->eps_ambient;

  m->epsilon_diag.x=m->epsilon_diag.y=m->epsilon_diag.z=eps;

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
  double wA             = 3.0;
  double wB             = 1.0;
  double taper_length   = 0.0;
  int taper_order       = 0;
  double wvg_length     = 5.0;
  int band_num          = 1;
  int  num_bands        = 2;
  double freq           = 0.25;
  bool use_symmetry     = false;
  bool plot_modes       = false;
  bool plot_flux        = false;
  bool plot_structure   = false;
  double frame_interval = 0.0;
  double res            = 25.0; // resolution
  char *filebase        = const_cast<char *>("wt");
  for(int narg=1; narg<argc; narg++)
   { if ( argv[narg]==0 )
      continue;
     if (!strcasecmp(argv[narg],"--wA"))
      { if ( ++narg >= argc )
         usage(argv[0], "error: no argument given for --wA");
        sscanf(argv[narg], "%le", &wA);
        master_printf("setting wA=%e\n",wA);
      }
     else if (!strcasecmp(argv[narg],"--wB"))
      { if ( ++narg >= argc )
         usage(argv[0], "error: no argument given for --wB");
        sscanf(argv[narg], "%le", &wB);
        master_printf("setting wB=%e\n",wB);
      }
     else if (!strcasecmp(argv[narg],"--taper-length"))
      { if ( ++narg >= argc )
         usage(argv[0], "error: no argument given for --taper-length");
        sscanf(argv[narg], "%le", &taper_length);
        master_printf("setting taper_length=%e\n",taper_length);
      }
     else if (!strcasecmp(argv[narg],"--taper-order"))
      { if ( ++narg >= argc )
         usage(argv[0], "error: no argument given for --taper-order");
        sscanf(argv[narg], "%i", &taper_order);
        master_printf("setting taper_order=%i\n",taper_order);
      }
     else if (!strcasecmp(argv[narg],"--wvg-length"))
      { if ( ++narg >= argc )
         usage(argv[0], "error: no argument given for --wvg-length");
        sscanf(argv[narg], "%le", &wvg_length);
        master_printf("setting wvg_length=%e\n",wvg_length);
      }
     else if (!strcasecmp(argv[narg],"--plot-modes") )
      plot_modes=true;
     else if (!strcasecmp(argv[narg],"--plot-flux") )
      plot_flux=true;
     else if (!strcasecmp(argv[narg],"--plot-structure") )
      plot_structure=true;
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
     else if (!strcasecmp(argv[narg],"--frame-interval"))
      { if ( ++narg >= argc )
         usage(argv[0], "error: no argument given for --frame-interval");
        sscanf(argv[narg], "%le", &frame_interval);
        master_printf("setting frame-interval=%e\n",frame_interval);
      }
     else if (!strcasecmp(argv[narg],"--freq"))
      { if ( ++narg >= argc )
         usage(argv[0], "error: no argument given for --freq");
        sscanf(argv[narg], "%le", &freq);
        master_printf("setting freq=%e\n",freq);
      }
     else if (!strcasecmp(argv[narg],"--resolution"))
      { if ( ++narg >= argc )
         usage(argv[0], "error: no argument given for --resolution");
        sscanf(argv[narg], "%le", &res);
        master_printf("setting resolution=%e\n",res);
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

  if (taper_length==0.0)
   taper_length=0.1/res;

  /***************************************************************/
  /* initialize computational cell                               */
  /****************** ********************************************/
  double dpml       = 0.5;
  double wmax       = fmax(wA, wB);
  double air_gap    = 0.25*wmax;
  double LX         = dpml + wvg_length + taper_length + wvg_length + dpml;
  double LY         = dpml + air_gap + wmax + air_gap + dpml;
  geometry_lattice.size.x = LX;
  geometry_lattice.size.y = LY;
  geometry_lattice.size.z = 0.0;
  grid_volume gv = voltwo(LX, LY, res);
  gv.center_origin();
  symmetry sym = use_symmetry ? -mirror(Y,gv) : identity();
  structure the_structure(gv, dummy_eps, pml(dpml), sym);

  /***************************************************************/
  /* specify user-defined material function                      */
  /***************************************************************/
  wvg_data wdata;
  wdata.wA            = wA;
  wdata.wB            = wB;
  wdata.taper_length  = taper_length;
  wdata.taper_order   = taper_order;
  wdata.eps_wvg       = 11.7;
  wdata.eps_ambient   = 1.0;  // ambient medium is vacuum
  meep_geom::material_type my_material
   = meep_geom::make_user_material(wvg_material, (void *)&wdata);
  bool use_anisotropic_averaging = true;
  double sbtol                   = DEFAULT_SUBPIXEL_TOL;
  int maxeval                    = DEFAULT_SUBPIXEL_MAXEVAL;
  bool ensure_periodicity        = false;
  bool verbose                   = false;
  geometric_object_list g={0,0};
  meep_geom::set_materials_from_geometry(&the_structure, g,
                                         use_anisotropic_averaging,
                                         sbtol, maxeval, ensure_periodicity,
                                         verbose, my_material);
  fields f(&the_structure);

  /***************************************************************/
  /* plot structure if requested *********************************/
  /***************************************************************/
  char filename[100];
  if (plot_structure)
   { snprintf(filename,100,"%s_L%g_p%i",filebase,taper_length,taper_order);
     h5file *eps_file=f.open_h5file("eps", h5file::WRITE, filename, false);
     f.output_hdf5(Dielectric,f.total_volume(),eps_file,false,false,0);
     delete eps_file;
   }

  /***************************************************************/
  /* add eigenmode source at left end of left waveguide          */
  /***************************************************************/
  double fcen = freq;
  double df   = 0.5*fcen;  // bandwidth
  int nfreq   = 1;         // number of frequency points
  gaussian_src_time gsrc(fcen, df);

  double x0  = -0.5*LX + dpml + 0.25*wvg_length; // inlet of left waveguide
  double xA  = -0.5*LX + dpml + 0.75*wvg_length; // midpoint of left waveguide
  double xB  = +0.5*LX - dpml - 0.5*wvg_length; // midpoint of right waveguide
  volume fv0 ( vec(x0,  -0.5*LY), vec(x0,  +0.5*LY) );
  volume fvA ( vec(xA,  -0.5*LY), vec(xA,  +0.5*LY) );
  volume fvB ( vec(xB,  -0.5*LY), vec(xB,  +0.5*LY) );
  double vol = LY;  // volume of flux planes
  direction d = X; // f.normal_direction(*fvA);

  bool match_frequency = true;
  int parity = 0; // NO_PARITY
  double tol=1.0e-4;
  vec kpoint=k_guess((void *)&wA, fcen, band_num);
  f.add_eigenmode_source(Dielectric, gsrc, d, fv0, fv0, band_num,
                         kpoint, match_frequency, parity, res, tol, 1.0);

  /***************************************************************/
  /* add flux planes                                             */
  /***************************************************************/
  dft_flux fluxA=f.add_dft_flux(d, fvA, fcen-0.5*df, fcen+0.5*df, nfreq);
  dft_flux fluxB=f.add_dft_flux(d, fvB, fcen-0.5*df, fcen+0.5*df, nfreq);

  /***************************************************************/
  /* timestep until all conditions for stopping are met.         */
  /*  conditions for stopping:                                   */
  /*   (1) sources have expired                                  */
  /*   (2) poynting flux through destination flux plane has      */
  /*       decayed below 1% of its maximum value                 */
  /***************************************************************/
  double PVCheckInterval=1.0, PVTol=1.0e-7;
  double NextPVCheckTime=f.round_time() + PVCheckInterval;
  double NextFileTime=f.round_time();
  double MaxPV=0.0;
  bool Stop=false;
  while(!Stop)
   {
     f.step();

     // do poynting-flux magnitude check at regular time intervals
     bool FieldsDecayed=false;
     if ( f.round_time() > NextPVCheckTime )
      { NextPVCheckTime += PVCheckInterval;
        double ThisPV = f.flux_in_box(X,fvA);
        if ( fabs(ThisPV) > MaxPV)
         MaxPV = fabs(ThisPV);
        else if ( fabs(ThisPV) < PVTol*MaxPV )
         FieldsDecayed=true;
      }

     // output HDF5 data at regular intervals if user requested that
     if (frame_interval>0.0 && f.round_time()>NextFileTime)
      { NextFileTime += frame_interval;
        f.output_hdf5(Sx,f.total_volume(),0,false,false,filebase);
      }

     bool SourcesFinished = ( f.round_time() > f.last_source_time() );

     Stop = SourcesFinished && FieldsDecayed;
   }

  /***************************************************************/
  /* plot eigenmode field patterns if requested ******************/
  /***************************************************************/
  void *mode_data_A, **mode_data_B = new void*[num_bands];
  if (plot_modes)
   for(int nb=-1; nb<num_bands; nb++)
    {
      int band_num   = (nb==-1) ? 1       : nb+1;
      volume *fv     = (nb==-1) ? &fvA    : &fvB;
      double ww      = (nb==-1) ? wA      : wB;
      char AB        = (nb==-1) ? 'A'     : 'B';
      dft_flux *flux = (nb==-1) ? &fluxA  : &fluxB;

      void *mode_data=f.get_eigenmode(fcen, d, *fv, *fv, band_num,
                                      k_guess((void *)&ww,fcen,band_num));
      if (nb==-1)
       mode_data_A=mode_data;
      else
       mode_data_B[nb]=mode_data;

      double vgrp=get_group_velocity(mode_data);

      snprintf(filename,100,"%s_mode%c%i",filebase,AB,band_num);
      f.output_mode_fields(mode_data, *flux, filename);
      cdouble mfOverlap[2], mmOverlap[2];
      f.get_mode_flux_overlap(mode_data, *flux, 0, d, mfOverlap);
      f.get_mode_mode_overlap(mode_data, mode_data, *flux, d, mmOverlap);
      master_printf("...is {%e,%e} {%e,%e}, should be %e\n",
                    real(mmOverlap[0]),imag(mmOverlap[0]),
                    real(mmOverlap[1]),imag(mmOverlap[1]),
                    vgrp*vol);

      master_printf("Computing <mode|flux> %c%i...\n",AB,band_num);
      if (am_master())
       {
         double vgrp=get_group_velocity(mode_data);
         vec k=get_k(mode_data);
         snprintf(filename,100,"%s.modeData",filebase);
         FILE *ff = fopen(filename,"a");
         fprintf(ff,"nb=%c%i: \n",AB,band_num);
         fprintf(ff,"  vgrp=%e, k=%e\n",vgrp,k.x());
         fprintf(ff,"  vgrp*area=%e\n",vgrp*vol);
         fprintf(ff,"  mf0={%e,%e}\n",real(mfOverlap[0]),imag(mfOverlap[0]));
         fprintf(ff,"  mf1={%e,%e}\n",real(mfOverlap[1]),imag(mfOverlap[1]));

         cdouble alpha[2];
         alpha[0]=0.5*(mfOverlap[0]+mfOverlap[1])/(vgrp*vol);
         alpha[1]=0.5*(mfOverlap[0]-mfOverlap[1])/(vgrp*vol);

         fprintf(ff,"   aP={%e,%e}\n",real(alpha[0]),imag(alpha[0]));
         fprintf(ff,"   aM={%e,%e}\n",real(alpha[1]),imag(alpha[1]));

         fclose(ff);
       };

    }; // if (plot_modes...) ... for(int nb=...)

  double *Aflux=fluxA.flux();
  double *Bflux=fluxB.flux();
  double Aflux0=0.0;

  /***************************************************************/
  /* for wA==wB, write flux to reference file and exit           */
  /***************************************************************/
  snprintf(filename,100,"%s_w%g_freq%g_res%g",filebase,wA,freq,res);
  if (wA==wB)
   { 
     if (am_master())
      { FILE *f=fopen(filename,"w");
        fprintf(f,"%e\n",Aflux[0]);
        fclose(f);
      }
   }
  else
   { FILE *f=fopen(filename,"r");
     if (!f)
      fprintf(stderr,"warning: could not open reference file %s\n",filename);
     else
      { fscanf(f,"%le",&Aflux0);
        fclose(f);
      };
   };

  /***************************************************************/
  /* write output files ******************************************/
  /***************************************************************/
  if (plot_flux)
   { char filename[100];
     snprintf(filename,100,"%s_fluxA",filebase);
     f.output_dft(fluxA, filename);
     snprintf(filename,100,"%s_fluxB",filebase);
     f.output_dft(fluxB, filename);
     snprintf(filename,100,"%s.fluxData",filebase);
     if (am_master())
      { FILE *ff=fopen(filename,"a");
        fprintf(ff,"flux A = %e\n",fluxA.flux()[0]);
        fprintf(ff,"flux B  = %e\n",fluxB.flux()[0]);
        fclose(ff);
      };
   };

  /***************************************************************/
  /* compute mode-expansion coefficients *************************/
  /***************************************************************/
  int *bands = new int[num_bands];
  for(int nb=0; nb<num_bands; nb++)
   bands[nb]=nb+1;
  int num_freqs = fluxB.Nfreq;
  cdouble *alphaA = new cdouble[2*num_freqs*num_bands];
  double *vgrpA   = new double[num_freqs*num_bands];
  cdouble *alphaB = new cdouble[2*num_freqs*num_bands];
  double *vgrpB   = new double[num_freqs*num_bands];
  f.get_eigenmode_coefficients(fluxA, d, fvA, bands, num_bands, alphaA, vgrpA);
  f.get_eigenmode_coefficients(fluxB, d, fvB, bands, num_bands, alphaB, vgrpB);

  if (am_master())
   {
     char filename[100];
     snprintf(filename,100,"%s.out",filebase);
     FILE *ff=fopen(filename,"a");
     for(int nf=0; nf<num_freqs; nf++)
      for(int nb=0; nb<num_bands; nb++)
       {
         cdouble aP = alphaA[2*nb*num_freqs + 2*nf + 0];
         cdouble aM = alphaA[2*nb*num_freqs + 2*nf + 1];
         double vgA = vgrpA[nb*num_freqs + nf];

         cdouble bP = alphaB[2*nb*num_freqs + 2*nf + 0];
         cdouble bM = alphaB[2*nb*num_freqs + 2*nf + 1];
         double vgB = vgrpB[nb*num_freqs + nf];

         fprintf(ff,"%e %e %e %i %e %i ",wA,wB,taper_length,taper_order,res,nb);
         fprintf(ff,"%e ",Aflux0);
         fprintf(ff,"%e %e %e %e ",Aflux[0], vgA, abs(aP), abs(aM));
         fprintf(ff,"%e %e %e %e ",Bflux[0], vgB, abs(bP), abs(bM));
         fprintf(ff,"\n");
      };
     fclose(ff);
   };

  return 0;
}
