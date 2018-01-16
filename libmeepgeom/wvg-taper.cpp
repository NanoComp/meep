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

  return vec(0.419984, 0.0, 0.0);
} 

/***************************************************************/
/* user-defined material function for 2D/3D waveguide geometry */
/***************************************************************/
typedef struct wvg_data
 { double wA;            // width of smaller waveguide
   double wB;            // width of larger waveguide
   double taper_length;  // taper length (=0 for no taper)
   int taper_order;      // index of first discontinuous derivative (1,2)
   double eps_wvg;       // permittivity of waveguide material
   double eps_ambient;   // permittivity of surrounding medium
   bool three_d;         // true if we are in the 3D case

   // these fields used only for 3D case
   double z_substrate;   // z-coordinate of substrate-oxide interface
   double z_oxide;       // z-coordinate of oxide-wvg interface
   double z_wvg;         // z-coordinate of wvg upper surface
   double eps_substrate; // dielectric constant of substrate
   double eps_oxide;     // dielectric constant of oxide layer
 } wvg_data;

void wvg_material(vector3 loc, void *user_data, meep_geom::medium_struct *m)
{
  wvg_data *wdata=(wvg_data *)user_data;

  double z = loc.z;
  double eps = wdata->eps_ambient; // assume we are in ambient medium
  if (wdata->three_d && z<wdata->z_substrate)  // 3D, in substrate
   { 
     eps = wdata->eps_substrate;
   }
  else if (wdata->three_d && z<wdata->z_oxide) // 3D, in oxide layer 
   {
     eps = wdata->eps_oxide;
   }
  else if (wdata->three_d && z>wdata->z_wvg )  // 3D, above waveguide
   {
     eps = wdata->eps_ambient;
   }
  else // 2D or 3D, inside waveguide
   {
     double x0 = loc.x / wdata->taper_length;
     double wA  = wdata->wA, wB = wdata->wB, w;
     if (x0 <= -0.5 )
      w = wA;
     else if ( x0 >= 0.5)
      w = wB;
     else if (wdata->taper_order==0)   // linear taper
      w = 0.5*(wA+wB) + (wB-wA)*x0;
     else // (wdata->taper_order==1)   // second-order taper
      w = 0.5*(wA+wB) + (wB-wA)*x0*(1.5 - 2.0*x0*x0);

     eps = ( fabs(loc.y)<=0.5*w ) ? wdata->eps_wvg : wdata->eps_ambient;
   };

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
  bool three_d          = false;
  double taper_length   = 0.0;
  double ratio          = 3.0;
  int taper_order       = 0;
  int  band_num         = 1;
  int  num_bands        = 6;
  bool use_symmetry     = false;
  bool plot_modes       = false;
  bool plot_flux        = false;
  bool plot_structure   = false;
  double frame_interval = 0.0;
  char *filebase        = const_cast<char *>("wt");
  double LY             = 3.0;  // half-width of cell in transverse direction double LZ=1.5;
  double LZ             = 2.0;  // half-width of cell in transverse direction double LZ=1.5;
  double dpml           = 0.50; // PML thickness
  double res            = 25.0; // resolution
  for(int narg=1; narg<argc; narg++)
   { if ( argv[narg]==0 )
      continue;
     if (!strcasecmp(argv[narg],"--three-d") )
      { three_d=true;
        master_printf("Using three-dimensional taper.\n");
      }
     else if (!strcasecmp(argv[narg],"--ratio"))
      { if ( ++narg >= argc )
         usage(argv[0], "error: no argument given for --ratio");
        sscanf(argv[narg], "%le", &ratio);
        master_printf("setting ratio=%e\n",ratio);
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
     else if (!strcasecmp(argv[narg],"--use-symmetry") )
      { use_symmetry=true;
        master_printf("Using symmetry.\n");
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
     else if (!strcasecmp(argv[narg],"--resolution"))
      { if ( ++narg >= argc )
         usage(argv[0], "error: no argument given for --resolution");
        sscanf(argv[narg], "%le", &res);
        master_printf("setting resolution=%e\n",res);
      }
     else if (!strcasecmp(argv[narg],"--LY"))
      { if ( ++narg >= argc )
         usage(argv[0], "error: no argument given for --LY");
        sscanf(argv[narg], "%le", &LY);
        master_printf("setting LY=%e\n",LY);
      }
     else if (!strcasecmp(argv[narg],"--LZ"))
      { if ( ++narg >= argc )
         usage(argv[0], "error: no argument given for --LZ");
        sscanf(argv[narg], "%le", &LZ);
        master_printf("setting LZ=%e\n",LZ);
      }
     else if (!strcasecmp(argv[narg],"--dpml"))
      { if ( ++narg >= argc )
         usage(argv[0], "error: no argument given for --dpml");
        sscanf(argv[narg], "%le", &dpml);
        master_printf("setting dpml=%e\n",dpml);
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
  double LX = dpml + 1.0 + 0.5*taper_length;
  geometry_lattice.size.x = 2*LX;
  geometry_lattice.size.y = 2*LY;
  geometry_lattice.size.z = ( (three_d) ? 2*LZ : 0.0 );
  grid_volume gv = three_d ? vol3d(2*LX, 2*LY, 2*LZ, res) : voltwo(2*LX, 2*LY, res);
  gv.center_origin();
  symmetry sym = use_symmetry ? -mirror(Y,gv) : identity();
  structure the_structure(gv, dummy_eps, pml(dpml), sym);

  /***************************************************************/
  /* specify user-defined material function                      */
  /***************************************************************/
  double wA            = 1.0;
  double wB            = wA*ratio;
  double h_substrate   = 0.0;    // no substrate by default
  double h_oxide       = 1.5;    // oxide layer thickness (3D case)
  double h_wvg         = 0.22;   // waveguide thickness (3D case)
  double eps_Si        = 11.7;   // dielectric constant of silicon
  double eps_SiO2      = 2.1;    // dielectric constant of SiO2

  wvg_data wdata;
  wdata.wA            = wA;
  wdata.wB            = wB;
  wdata.taper_length  = taper_length;
  wdata.taper_order   = taper_order;
  wdata.eps_wvg       = eps_Si;
  wdata.eps_ambient   = 1.0;  // ambient medium is vacuum
  wdata.three_d       = three_d;
  if (three_d)
   { wdata.z_substrate   = -LZ + h_substrate;
     wdata.z_oxide       = -LZ + h_substrate + h_oxide;
     wdata.z_wvg         = -LZ + h_substrate + h_oxide + h_wvg;
     wdata.eps_substrate = eps_Si;
     wdata.eps_oxide     = eps_SiO2;
   };
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
  if (plot_structure)
   { char filename[100];
     snprintf(filename,100,"%s_L%g_p%i",filebase,taper_length,taper_order);
     h5file *eps_file=f.open_h5file("eps", h5file::WRITE, filename, false);
     f.output_hdf5(Dielectric,f.total_volume(),eps_file,false,false,0);
     delete eps_file;
   }

  /***************************************************************/
  /* add source                                                  */
  /***************************************************************/
  double fcen = three_d ? 0.41 : 0.15;  // center frequency
  double df   = 0.5*fcen;  // bandwidth
  int nfreq   = 1;         // number of frequency points
  gaussian_src_time gsrc(fcen, df);

  double LXP = LX-dpml;
  double LYP = LY-dpml;
  double LZP = three_d ? LZ-dpml : 0.0;
  volume *fvA, *fvB, *fvC;
  if (three_d)
   { fvA = new volume( vec(-0.5*LX, -LYP, -LZP), vec(-0.5*LX, +LYP, +LZP) );
     fvB = new volume( vec(+0.5*LX, -LYP, -LZP), vec(+0.5*LX, +LYP, +LZP) );
     fvC = new volume( vec(   -LXP,    0, -LZP), vec(    LXP,    0,  LZP) );
   }
  else
   { fvA = new volume( vec(-0.5*LX, -LYP), vec(-0.5*LX, +LYP) );
     fvB = new volume( vec(+0.5*LX, -LYP), vec(+0.5*LX, +LYP) );
     fvC = new volume( vec(   -LXP,    0), vec(    LXP,    0) );
   };
  direction dA = f.normal_direction(*fvA);
  direction dB = f.normal_direction(*fvB);

  bool match_frequency = true;
  int parity = 0; // NO_PARITY
  double tol=1.0e-4;
  vec kpoint=k_guess((void *)&wA, fcen, band_num);
  f.add_eigenmode_source(Dielectric, gsrc, dA, *fvA, *fvA, band_num,
                         kpoint, match_frequency, parity, res, tol, 1.0);

  /***************************************************************/
  /* add flux planes                                             */
  /***************************************************************/
  dft_flux fluxA=f.add_dft_flux_plane(*fvA, fcen-0.5*df, fcen+0.5*df, nfreq);
  dft_flux fluxB=f.add_dft_flux_plane(*fvB, fcen-0.5*df, fcen+0.5*df, nfreq);
  dft_flux fluxC=f.add_dft_flux_plane(*fvC, fcen-0.5*df, fcen+0.5*df, nfreq);

  /***************************************************************/
  /* plot eigenmode field patterns if requested ******************/
  /***************************************************************/
  if (plot_modes)
   { for(int nb=0; nb<num_bands; nb++)
      { int band_num=nb+1;
        void *vedata=f.get_eigenmode(fcen, dB, *fvB, *fvB,
                                     band_num, k_guess((void *)&wB,
                                     fcen,band_num),
                                     true, 0, f.a, 1.0e-4);
        char filename[100];
        snprintf(filename,100,"%s_mode%i",filebase,band_num);
        f.output_mode_fields(vedata, fluxB, *fvB, filename);
        double vgrp=get_group_velocity(vedata);
        vec k=get_k(vedata);
        if (am_master()) 
         { 
           snprintf(filename,100,"%s.modeData",filebase);
           static char mode[2]="w";
           FILE *ff = fopen(filename,mode);
           mode[0]='a';
           fprintf(ff,"nb=%2i vgrp=%e k=%e \n",nb,vgrp,k.x());
           fclose(ff);
         };
        destroy_eigenmode_data(vedata);
      };
   };

  /***************************************************************/
  /* timestep until all conditions for stopping are met.         */
  /*  conditions for stopping:                                   */
  /*   (1) sources have expired                                  */
  /*   (2) poynting flux through destination flux plane has      */
  /*       decayed below 1% of its maximum value                 */
  /***************************************************************/
  double PVCheckInterval=1.0, PVTol=0.01;
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
        double ThisPV = f.flux_in_box(X,*fvB);
        if (ThisPV > MaxPV)
         MaxPV = ThisPV;
        else if ( ThisPV < PVTol*MaxPV )
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
  /* write output files ******************************************/
  /***************************************************************/
  if (plot_flux)
   { char filename[100];
     snprintf(filename,100,"%s_fluxA",filebase);
     f.output_flux_fields(fluxA, *fvA, filename);
     snprintf(filename,100,"%s_fluxB",filebase);
     f.output_flux_fields(fluxB, *fvB, filename);
     snprintf(filename,100,"%s_fluxC",filebase);
     f.output_flux_fields(fluxC, *fvC, filename);
   };

  /***************************************************************/
  /* compute mode-expansion coefficients *************************/
  /***************************************************************/
  std::vector<int> bands(num_bands);
  for(int n=0; n<num_bands; n++)
   bands[n] = n+1;

  int num_freqs = fluxB.Nfreq;
  std::vector<cdouble> coeffs =
   f.get_eigenmode_coefficients(fluxB, dB, *fvB, bands, k_guess, (void *)&wB);

  double *bflux=fluxB.flux();
  double bvol1=fvB->computational_volume();
  double bvol2=fvB->integral_volume();
  double bvol3=fvB->full_volume();

  if (am_master())
   {

     char filename[100];
     snprintf(filename,100,"%s.coefficients",filebase);
     FILE *ff=fopen(filename,"a");
     printf("freq | band | alpha^+ | alpha^-\n");
     printf("------------------------------------------------\n");
     for(int nf=0; nf<num_freqs; nf++)
      for(unsigned nb=0; nb<bands.size(); nb++)
       { 
         double atot=0.0;
         for(int nbb=0; nbb<num_bands; nbb++)
          for(int pm=0; pm<2; pm++)
           atot += norm( coeffs[2*nbb*num_freqs + 2*nf + pm] );
   
         cdouble aP = coeffs[2*nb*num_freqs + 2*nf + 0];
         cdouble aM = coeffs[2*nb*num_freqs + 2*nf + 1];
         printf("%2i  %2i  (+)  %e {%+e,%+e} (%e %%)\n",nf,nb,abs(aP),real(aP),imag(aP),100.0*norm(aP)/atot);
         printf("%2i  %2i  (-)  %e {%+e,%+e} (%e %%)\n",nf,nb,abs(aM),real(aM),imag(aM),100.0*norm(aM)/atot);
         fprintf(ff,"%g %.2f %i %g %2i %2i  %e %e %e  %e %e %e %e %e \n",ratio,taper_length,taper_order,res,nb,nf,norm(aP), arg(aP), norm(aP)/atot, bflux[nf], bvol1, bvol2, bvol3, mode_group_velocities[nb*num_freqs + nf]);
      };
     fclose(ff);
   };
   
  return 0;
}
