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
  (void) freq;
  (void) band_num;
 
  double w=*((double *) user_data);

  // hard-coded dispersion relations for w=0.5, 1
  if ( equal_float(w,0.5) && equal_float(freq,0.41) )
   { if (band_num==1)  return vec(0.9892331, 0.0, 0.0);
     if (band_num==2)  return vec(0.6175083, 0.0, 0.0);
     if (band_num==3)  return vec(0.5469879, 0.0, 0.0);
     if (band_num==4)  return vec(0.5245156, 0.0, 0.0);
     if (band_num==5)  return vec(0.4267270, 0.0, 0.0);
     if (band_num>=6)  return vec(0.4245740, 0.0, 0.0);
   }
  else // if ( equal_float(w,1.0) && equal_float(freq,0.41) )
   { if (band_num==1)  return vec(1.073627,  0.0, 0.0);
     if (band_num==2)  return vec(0.9856316, 0.0, 0.0);
     if (band_num==3)  return vec(0.8233921, 0.0, 0.0);
     if (band_num==4)  return vec(0.5844210, 0.0, 0.0);
     if (band_num==5)  return vec(0.5497692, 0.0, 0.0);
     if (band_num>=6)  return vec(0.5116020, 0.0, 0.0);
   }

  return vec(0.2, 0.0, 0.0);
} 

/***************************************************************/
/***************************************************************/
/***************************************************************/
typedef struct material_func_data
 { double wA;            // width of smaller waveguide
   double wB;            // width of larger waveguide
   double z_substrate;   // z-coordinate of substrate-oxide interface
   double z_oxide;       // z-coordinate of oxide-strip interface
   double z_strip;       // z-coordinate of strip-air interface
   double taper_length;  // taper length (=0 for no taper)
   double taper_power;   // taper exponent (1,2 for linear,quadratic taper)
   double eps_substrate; // dielectric constant of substrate
   double eps_oxide;     // dielectric constant of oxide layer
   double eps_strip;     // dielectric constant of waveguide (strip)
   double eps_ambient;   // dielectric constant of ambient medium
 } material_func_data;

void material_func(vector3 loc, void *user_data, meep_geom::medium_struct *m)
{
  material_func_data *mfdata=(material_func_data *)user_data;

  double z = loc.z;
  double eps = mfdata->eps_ambient; // assume we are in air above waveguide
  if (z<mfdata->z_substrate)
   { 
     eps = mfdata->eps_substrate;
   }
  else if (z<mfdata->z_oxide)
   {
     eps = mfdata->eps_oxide;
   }
  else if (z<mfdata->z_strip) // in waveguide
   {
     double x = loc.x; 
     double y = loc.y;
     double x_min = -0.5*mfdata->taper_length;
     double x_max = +0.5*mfdata->taper_length;
     double wA = mfdata->wA, wB = mfdata->wB, w;
     if (x <= x_min)
      w = wA;
     else if (x>=x_max)
      w = wB;
     else
      w = wA + (wB-wA)*pow( (x-x_min)/(x_max-x_min), mfdata->taper_power);

     eps = ( fabs(y)<=w ) ? mfdata->eps_strip : mfdata->eps_ambient;
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
  bool use_symmetry     = false;
  bool point_source     = false;
  int  band_num         = 1;
  int  num_bands        = 6;
  double ratio          = 1.0;
  char *filebase        = const_cast<char *>("fj");
  bool plot_modes       = false;
  bool plot_flux        = false;
  double frame_interval = 0.0;
  double taper_length   = 0.0;
  double taper_power    = 1.0;
  double resolution     = 10.0; // resolution
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
     else if (!strcasecmp(argv[narg],"--taper-length"))
      { if ( ++narg >= argc )
         usage(argv[0], "error: no argument given for --taper-length");
        sscanf(argv[narg], "%le", &taper_length);
        master_printf("setting taper_length=%e\n",taper_length);
      }
     else if (!strcasecmp(argv[narg],"--taper-power"))
      { if ( ++narg >= argc )
         usage(argv[0], "error: no argument given for --taper-power");
        sscanf(argv[narg], "%le", &taper_power);
        master_printf("setting taper_length=%e\n",taper_power);
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
        sscanf(argv[narg], "%le", &resolution);
        master_printf("setting resolution=%e\n",resolution);
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
  /* initialize computational cell                               */
  /****************** ********************************************/
  double LX=3.0;       // half-length (in X) of computational cell
  double LY=2.0;       // half-width (in Y)
  double LZ=1.5;       // half-width (in Z)
  double dpml=0.50;    // PML thickness
  geometry_lattice.size.x=0;
  geometry_lattice.size.y=2*LY;
  geometry_lattice.size.z=2*LZ;
  grid_volume gv=vol3d(2.0*LX, 2.0*LY, 2.0*LZ, resolution);
  gv.center_origin();
  symmetry sym = use_symmetry ? -mirror(Y,gv) : identity();
  structure the_structure(gv, dummy_eps, pml(dpml), sym);

  /***************************************************************/
  /* specify user-defined material function             **********/
  /***************************************************************/
  double wA            = 0.5;
  double wB            = wA*ratio;
  double h_substrate   = 0.0;    // no substrate by default
  double h_oxide       = 1.5;    // oxide layer thickness
  double h_strip       = 0.22;   // strip thickness
  double eps_Si        = 11.7;   // dielectric constant of silicon
  double eps_SiO2      = 2.1;    // dielectric constant of SiO2

  material_func_data mfdata;
  mfdata.wA            = wA;
  mfdata.wB            = wB;
  mfdata.z_substrate   = -LZ + h_substrate;
  mfdata.z_oxide       = -LZ + h_substrate + h_oxide;
  mfdata.z_strip       = -LZ + h_substrate + h_oxide + h_strip;
  mfdata.taper_length  = taper_length;
  mfdata.taper_power   = taper_power;
  mfdata.eps_substrate = eps_Si;
  mfdata.eps_oxide     = eps_SiO2;
  mfdata.eps_strip     = eps_Si;
  mfdata.eps_ambient   = 1.0;  // ambient medium is vacuum
  meep_geom::material_type my_material
   = meep_geom::make_user_material(material_func, (void *)&mfdata);
  bool use_anisotropic_averaging=true;
  double tol  = DEFAULT_SUBPIXEL_TOL;
  int maxeval = DEFAULT_SUBPIXEL_MAXEVAL;
  bool ensure_periodicity = false;
  bool verbose            = false;
  geometric_object_list g={0,0};
  meep_geom::set_materials_from_geometry(&the_structure, g,
                                         use_anisotropic_averaging,
                                         tol, maxeval, ensure_periodicity,
                                         verbose, my_material);

  fields f(&the_structure);
  f.output_hdf5(Dielectric,f.total_volume(),0,false,false,filebase);

  /***************************************************************/
  /* add source                                                  */
  /***************************************************************/
  double fcen = 0.41;      // center frequency
  double df   = 0.5*fcen;  // bandwidth
  int nfreq   = 1;         // number of frequency points
  gaussian_src_time gsrc(fcen, df);

  double LXP=LX-dpml;
  double LYP=LY-dpml;
  double LZP=LZ-dpml;
  volume fvA( vec(-0.5*LX, -LYP, -LZP), vec(-0.5*LX, +LYP, +LZP) );
  volume fvB( vec(+0.5*LX, -LYP, -LZP), vec(+0.5*LX, +LYP, +LZP) );
  volume fvC( vec(-LXP,0,-LZP),     vec(LXP,0,LZP) );

  direction dA = f.normal_direction(fvA);
  direction dB = f.normal_direction(fvB);

  if (point_source)
   { 
     double z_source = -LZ + h_substrate + h_oxide + 0.5*h_strip;
     f.add_point_source(Ex, gsrc, vec(-0.5*LX, 0.0, z_source));
   }
  else
   {
     bool match_frequency = true;
     int parity = 0; // NO_PARITY
     double resolution=f.a;
     double eigensolver_tol=1.0e-4;
     vec kpoint=k_guess((void *)&wA, fcen, band_num);
     f.add_eigenmode_source(Dielectric, gsrc, dA, fvA, fvA, band_num,
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
  /* timestep                                                    */
  /***************************************************************/
  double NextFileTime=f.round_time();
  while( f.round_time() < (f.last_source_time() + 5.0/df) )
   { f.step();
     if (frame_interval>0.0 && f.round_time()>NextFileTime)
      { NextFileTime += frame_interval;
        f.output_hdf5(Ex,f.total_volume(),0,false,false,filebase);
        f.output_hdf5(Hx,f.total_volume(),0,false,false,filebase);
        f.output_hdf5(Sx,f.total_volume(),0,false,false,filebase);
      };
   };
 
  /***************************************************************/
  /* write output files ******************************************/
  /***************************************************************/
  if (plot_flux)
   { char filename[100];
     snprintf(filename,100,"%s_fluxA",filebase);
     f.output_flux_fields(&fluxA, fvA, filename);
     snprintf(filename,100,"%s_fluxB",filebase);
     f.output_flux_fields(&fluxB, fvB, filename);
     snprintf(filename,100,"%s_fluxC",filebase);
     f.output_flux_fields(&fluxC, fvC, filename);
   };

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  if (plot_modes)
   { 
     void **vedata = (void **)malloc(num_bands * sizeof(void *));
     double *vgrp  = new double[num_bands];
     char filename[100];
     for(int nb=0; nb<num_bands; nb++)
      { int band_num=nb+1;
        vedata[nb]=f.get_eigenmode(fcen, dB, fvB, fvB,
                                   band_num, k_guess((void *)&wB,
                                   fcen,band_num),
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
     double vol=4.0*LYP*LZP; // 2-dimensional "volume" of flux plane
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
   f.get_eigenmode_coefficients(&fluxB, dB, fvB, bands, k_guess, (void *)&wB);

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
         double atot=0.0;
         for(int nbb=0; nbb<num_bands; nbb++)
          for(int pm=0; pm<2; pm++)
           atot += norm( coeffs[2*nbb*num_freqs + 2*nf + pm] );
   
         cdouble aP = coeffs[2*nb*num_freqs + 2*nf + 0];
         cdouble aM = coeffs[2*nb*num_freqs + 2*nf + 1];
         printf("%2i  %2i  (+)  %e {%+e,%+e} (%e %%)\n",nf,nb,abs(aP),real(aP),imag(aP),100.0*norm(aP)/atot);
         printf("%2i  %2i  (-)  %e {%+e,%+e} (%e %%)\n",nf,nb,abs(aM),real(aM),imag(aM),100.0*norm(aM)/atot);
         fprintf(ff,"%2i  %2i   %e {%+e,%+e} (%e %%)  %e  {%+e, %+e} (%e %%)\n",nf,nb,
                     abs(aP),real(aP),imag(aP),norm(aP)/atot, abs(aM),real(aM),imag(aM),norm(aM)/atot);
      };
     fclose(ff);
   };
   
  return 0;
}
