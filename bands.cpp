/* Copyright (C) 2003 Massachusetts Institute of Technology  
%
%  This program is free software; you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation; either version 2, or (at your option)
%  any later version.
%
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%
%  You should have received a copy of the GNU General Public License
%  along with this program; if not, write to the Free Software Foundation,
%  Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "dactyl.h"
#include "dactyl_internals.h"
#include "harminv.h"

#define BAND(b,r,t) ((b)[(r)+(t)*nr])

bandsdata::bandsdata() {
  verbosity = 0;
  maxbands = -1;
  tstart = nr = z = 0;
  tend = -1;
  hr = hp = hz = er = ep = ez = NULL;
  P = NULL;
}

bandsdata::~bandsdata() {
  delete[] hr;
  delete[] hp;
  delete[] hz;
  delete[] er;
  delete[] ep;
  delete[] ez;
  delete[] P;
}

int src::find_last_source(int sofar) {
  if (peaktime + cutoff > sofar) sofar = (int)peaktime + cutoff;
  if (next == NULL) return sofar;
  return next->find_last_source(sofar);
}

int fields::find_last_source() {
  int last_source = 0;
  if (e_sources != NULL)
    last_source = e_sources->find_last_source();
  if (h_sources != NULL)
    last_source = max(last_source, h_sources->find_last_source());
  return last_source;  
}

void fields::prepare_for_bands(int z, double endtime, double fmax,
                               double qmin, double frac_pow_min) {
  int last_source = find_last_source();
  last_source = max(last_source, t + phasein_time);
  if (fmax == 0) fmax = preferred_fmax;
  else preferred_fmax = fmax;
  if (!bands) bands = new bandsdata;
  bands->tstart = last_source+1;
  if (bands->tstart < t) bands->tstart = t;
  bands->tend = t + (int)(endtime*a/c) - 1;

  if (z >= nz) {
    printf("Specify a lower z for your band structure! (%d > %d)\n",
           z, nz);
    exit(1);
  }
  bands->fpmin = frac_pow_min;
  bands->z = z;
  bands->nr = nr;

  // Set fmin properly...
  double epsmax = 1;
  for (int r=0;r<nr;r++) {
    for (z=0;z<nz;z++) {
      if (MA(ma->eps,r,z) > epsmax) epsmax = MA(ma->eps,r,z);
    }
  }
  const double cutoff_freq = 1.84*c/(2*pi)/nr/sqrt(epsmax);
  bands->fmin = sqrt(cutoff_freq*cutoff_freq + k*k*c*c/epsmax);
  bands->fmin = cutoff_freq/(c*inva);
  bands->qmin = qmin;
  // Set fmax and determine how many timesteps to skip over...
  bands->fmax = fmax;
  {
    // for when there are too many data points...
    double decayconst = bands->fmax*(c*inva)/qmin*8.0;
    double smalltime = 1./(decayconst + bands->fmax*(c*inva));
    bands->scale_factor = (int)(0.06*smalltime);
    if (bands->scale_factor < 1) bands->scale_factor = 1;
    if (verbosity) printf("scale_factor is %d (%lg,%lg)\n",
                          bands->scale_factor, bands->fmax*(c*inva), decayconst);
  }

  if (bands->tend <= bands->tstart) {
    printf("Oi, we don't have any time to take a fourier transform!\n");
    printf("FT start is %d and end is %d\n", bands->tstart, bands->tend);
    exit(1);
  }
  bands->ntime = (1+(bands->tend-bands->tstart)/bands->scale_factor);
  bands->a = a;
  bands->inva = inva;
  bands->hr = new cmplx[nr*bands->ntime];
  bands->hp = new cmplx[nr*bands->ntime];
  bands->hz = new cmplx[nr*bands->ntime];
  bands->er = new cmplx[nr*bands->ntime];
  bands->ep = new cmplx[nr*bands->ntime];
  bands->ez = new cmplx[nr*bands->ntime];
  bands->P = new cmplx[bands->ntime];
  for (int i=0;i<bands->ntime;i++) bands->P[i] = 0.0;
  if (bands->ez == NULL) {
    printf("Unable to allocate bandstructure array!\n");
    exit(1);
  }
  bands->verbosity = verbosity;
}

void fields::record_bands() {
  if (t > bands->tend || t < bands->tstart) return;
  if (t % bands->scale_factor != 0) return;
  int thet = (t-bands->tstart)/bands->scale_factor;
  if (thet >= bands->ntime) return;
  for (int r=0;r<nr;r++) {
    BAND(bands->hr,r,thet) =
      cmplx(RE(hr,r,bands->z), IM(hr,r,bands->z));
    BAND(bands->hp,r,thet) =
      cmplx(RE(hp,r,bands->z), IM(hp,r,bands->z));
    BAND(bands->hz,r,thet) =
      cmplx(RE(hz,r,bands->z), IM(hz,r,bands->z));
    BAND(bands->er,r,thet) =
      cmplx(RE(er,r,bands->z), IM(er,r,bands->z));
    BAND(bands->ep,r,thet) =
      cmplx(RE(ep,r,bands->z), IM(ep,r,bands->z));
    BAND(bands->ez,r,thet) =
      cmplx(RE(ez,r,bands->z), IM(ez,r,bands->z));
  }
  if (pol) {
    int zmean = nz/2;
    int rmean = nr/2;
    // The 1.137 below is just a pathetic attempt to avoid the situation
    // where some symmetry makes Pp and Pz equal and opposite.  The better
    // solution would be to store the two of them separately.
    bands->P[thet] = cmplx(RE(pol->Pp, rmean, zmean)+1.137*RE(pol->Pz, rmean, zmean),
                           IM(pol->Pp, rmean, zmean)+1.137*IM(pol->Pz, rmean, zmean));
  }
}

#define HARMOUT(o,r,n,f) ((o)[(r)+(n)*nr+(f)*nr*maxbands])

complex<double> fields::get_band(int nn, int maxbands) {
  complex<double> *fad = get_the_bands(maxbands);
  complex<double> thef = fad[nn-1];
  delete[] fad;
  return thef;
}

extern "C" {
  extern void zgesvd_(char *, char *, int *m, int *n, 
                     cmplx *A, int*lda, double *S, cmplx *U, int*ldu,
                     cmplx *VT,int*ldvt,
                     cmplx *WORK, int *lwork,double *RWORK, int*);
}

void bandsdata::get_fields(cmplx *eigen, cmplx *fad,
                           int nbands, int n) {
  const int dofourier = 1;
  // eigen has dimensions of nr*maxbands*6
  // n here is the total time
  double unitconvert = (2*pi)*c*scale_factor*inva;

  // Now let's zero out any decays that seem unreasonably small...
  for (int i = 0; i < nbands; i++) {
    if (abs(unitconvert*imag(fad[i])*n) < 1.0) {
      printf("I think that band with freq %lg and decay %lg has no decay.\n",
             real(fad[i]), imag(fad[i]));
      fad[i] = real(fad[i]);
    } else {
      printf("It seems the band with freq %lg and decay %lg does decay.\n",
             real(fad[i]), imag(fad[i]));
    }
  }
  // First, we want to take a fourier transform of each mode that has zero
  // decay, and subtract off its power.
  int numfit = 0;
  if (dofourier) {
    for (int r=0;r<nr;r++) {
      for (int whichf=0;whichf<6;whichf++) {
        cmplx *bdata;
        switch (whichf) {
        case 0: bdata = er; break;
        case 1: bdata = ep; break;
        case 2: bdata = ez; break;
        case 3: bdata = hr; break;
        case 4: bdata = hp; break;
        case 5: bdata = hz; break;
        }
        {
          cmplx mean = 0;
          for (int t=0;t<n;t++) {
            mean += BAND(bdata,r,t);
          }
          for (int t=0;t<n;t++) BAND(bdata,r,t) -= (1.0/n)*mean;
        }
        for (int j=0;j<nbands;j++) {
          if (imag(fad[j]) == 0.0) {
            cmplx posfreq = 0;
            cmplx negfreq = 0;
            double energy = 0;
            for (int t=0;t<n;t++) {
              complex<double> phase = exp(cmplx(0.0,unitconvert*real(fad[j])*t));
              posfreq += BAND(bdata,r,t)*phase;
              negfreq += BAND(bdata,r,t)*conj(phase);
              energy += abs(BAND(bdata,r,t))*abs(BAND(bdata,r,t));
            }
            double oldenergy = energy/n;
            posfreq *= 1.0/n;
            negfreq *= 1.0/n;
            HARMOUT(eigen,r,j,whichf) = posfreq;
            for (int t=0;t<n;t++) {
              complex<double> phase = exp(cmplx(0.0,unitconvert*real(fad[j])*t));
              BAND(bdata,r,t) -= posfreq*conj(phase);
              BAND(bdata,r,t) -= negfreq*phase;
              energy += abs(BAND(bdata,r,t))*abs(BAND(bdata,r,t));
            }
          } else numfit++;
        }
      }
    }
    numfit /= 6*nr;
  } else {
    numfit = nbands;
  }
  if (verbosity) printf("Looking at %d bands using least squares...\n", numfit);
  if (numfit < 2) {
    printf("It looks like we can't do the least squares fit, since we\n");
    printf("only have %d frequencies to fit.\n", numfit);
  } else {
    int num_param = numfit;
    cmplx *A = new cmplx[num_param*n];

    int ifit = 0;
    for (int i=0;i<nbands;i++) {
      if (imag(fad[i]) != 0.0 || !dofourier) {
        for (int t=0;t<n;t++) {
          A[ifit*n+t] = exp(-unitconvert*fad[i]*t);
        }
        ifit++;
      }
    }
    int one = 1, zero=0, worksize = num_param*n, dummy;
    double *S = new double[num_param];
    cmplx *VT = new cmplx[num_param*num_param];
    cmplx *V = new cmplx[num_param*num_param];
    cmplx *work = new cmplx[worksize];
    double *rwork = new double[5*num_param];
    zgesvd_("O","A",&n,&num_param,A,&n,S,NULL,&one,VT,&num_param,
            work,&worksize,rwork,&dummy);
    // Unhermitian conjugate V
    for (int j=0;j<num_param;j++) {
      if (verbosity) printf("S[%d] value is %lg.\n", j, S[j]);
      for (int i=0;i<num_param;i++) {
        V[i*num_param+j] = conj( VT[j*num_param+i]);
      }
    }
    for (int r=0;r<nr;r++) {
      for (int whichf=0;whichf<6;whichf++) {
        cmplx *bdata;
        switch (whichf) {
        case 0: bdata = er; break;
        case 1: bdata = ep; break;
        case 2: bdata = ez; break;
        case 3: bdata = hr; break;
        case 4: bdata = hp; break;
        case 5: bdata = hz; break;
        }
        for (int j=0;j<nbands;j++) {
          cmplx possofar = 0, negsofar = 0;
          if (imag(fad[j]) != 0.0 || !dofourier) {
            for (int i=0;i<numfit;i++) {
              cmplx possum=0, negsum=0;
              if (S[i] > 1e-6) {
                for (int t=0;t<n;t++) {
                  possum += A[i*n+t]*BAND(bdata,r,t);
                  //negsum += A[(i+numfit)*n+t]*bandhere;
                }
                possofar += possum*V[i*numfit+j]/S[i];
                //negsofar += negsum*V[(i+numfit)*numfit+j]/S[(i+numfit)];
              }
            }
            if (abs(negsofar) > abs(possofar))
              HARMOUT(eigen,r,j,whichf) = negsofar;
            else
              HARMOUT(eigen,r,j,whichf) = possofar;
          }
        }
      }
    }
    delete[] VT;
    delete[] V;
    delete[] S;
    delete[] rwork;
    delete[] work;
    delete[] A;
  }
}

void fields::output_bands_and_modes(FILE *o, const char *name, int maxbands) {
  out_bands(o, name, maxbands, 1);
}

void fields::grace_bands(grace *g, int maxbands) {
  double *approx_power = new double[maxbands];
  //complex<double> *fad = get_the_bands(maxbands, approx_power);
  complex<double> *fad = clever_cluster_bands(maxbands, approx_power);

  int num_found = 0;
  for (int i=0;i<maxbands;i++) if (fad[i] != 0) num_found = i+1;

  for (int i = 0; i < num_found; ++i) {
    g->output_out_of_order(i, k, fabs(real(fad[i])), fabs(imag(fad[i])),
                           approx_power[i]);
  }
  delete[] fad;
  delete[] approx_power;
}

void fields::output_bands(FILE *o, const char *name, int maxbands) {
  out_bands(o, name, maxbands, 0);
}

void fields::out_bands(FILE *o, const char *name, int maxbands, int and_modes) {
  double *approx_power = new double[maxbands];
  complex<double> *fad = get_the_bands(maxbands, approx_power);

  cmplx *eigen = new cmplx[nr*maxbands*6];
  if (!eigen) {
    printf("Error allocating...\n");
    exit(1);
  }

  for (int r=0;r<nr;r++) {
    for (int whichf = 0; whichf < 6; whichf++) {
      for (int n=0;n<maxbands;n++) {
        HARMOUT(eigen,r,n,whichf) = 0;
      }
    }
  }
  int num_found = 0;
  for (int i=0;i<maxbands;i++) if (fad[i] != 0) num_found = i+1;
  if (and_modes) bands->get_fields(eigen,fad,num_found,bands->ntime);

  for (int i = 0; i < num_found; ++i) {
    // k m index freq decay Q
    fprintf(o, "%s %lg %d %d %lg %lg %lg %lg\n", name,
            k, m, i, fabs(real(fad[i])), imag(fad[i]),
            -fabs(real(fad[i])) / (2 * imag(fad[i])), approx_power[i]);
    if (and_modes) {
      for (int r=0;r<nr;r++) {
        fprintf(o, "%s-fields %lg %d %d %lg", name, k, m, i, r*inva);
        for (int whichf = 0; whichf < 6; whichf++) {
          fprintf(o, " %lg %lg", real(HARMOUT(eigen,r,i,whichf)),
                  imag(HARMOUT(eigen,r,i,whichf)));
        }
        fprintf(o, "\n");
      }
    }
  }
  delete[] approx_power;
  delete[] fad;
}

int bandsdata::look_for_more_bands(complex<double> *simple_data,
                                   double *reff, double *refd,
                                   complex<double> *refa,
                                   complex<double> *refdata,
                                   int numref) {
  const double total_time = (tend-tstart)*c/a;
  const double deltaf = 2.0/total_time;
  if (numref == 0) { // Have no reference bands so far...
    numref = get_freqs(simple_data, ntime, refa, reff, refd);
    for (int n=0;n<numref;n++) {
      if (verbosity > 1) printf("Here's a mode (%10lg,%10lg) (%10lg,%10lg) -- %d\n",
                                reff[n], refd[n], real(refa[n]),imag(refa[n]), n);
      for (int t=0;t<ntime;t++)
        refdata[t+n*ntime] = simple_data[t];
    }
  } else if (numref < 0) { // Don't try to match these new bands...
    int tempmaxbands = maxbands;
    if (tempmaxbands+numref > 0) {
      maxbands = maxbands+numref;
      int moreref = get_freqs(simple_data, ntime,
                              refa+abs(numref), reff+abs(numref), refd+abs(numref));
      maxbands = tempmaxbands;
      numref = moreref - numref;
      for (int n=numref-moreref;n<numref;n++) {
        if (verbosity > 1)
          printf("HERE's a mode (%10lg,%10lg) (%10lg,%10lg) -- %d\n",
                 reff[n], refd[n], real(refa[n]),imag(refa[n]), n);
        for (int t=0;t<ntime;t++)
          refdata[t+n*ntime] = simple_data[t];
      }
    } else if (verbosity) {
      printf("Yikes, running out of room in bandsdata!\n");
    }
  } else {
    double *tf = new double[maxbands];
    double *heref = new double[maxbands];
    double *td = new double[maxbands];
    double *hered = new double[maxbands];
    cmplx *ta = new cmplx[maxbands];
    cmplx *herea = new cmplx[maxbands];
    for (int n=0;n<numref;n++) {
      int num_match = get_both_freqs(refdata+n*ntime,simple_data,ntime,
                                     ta, herea, tf, td);
      if (verbosity && !num_match) printf("I don't see any matches for ref %d here!\n",n);
      int best_match=-1;
      double err_best = 1e300;
      for (int i=0;i<num_match;i++) {
        double err = (abs(tf[i]-reff[n])+0.25*abs(td[i]-refd[n]));
        if (err < err_best) {
          best_match = i;
          err_best = err;
        }
      }
      //printf("Setting amp to %lg (vs %lg)\n",
      //       abs(herea[best_match]), abs(ta[best_match]));
      const double err_limit = deltaf;
      if (err_best > err_limit && err_best < 1e299) {
        if (verbosity > 2 && abs(herea[best_match]) != 0.0)
          printf("Missed a match at %d between %lg and %lg (%lg vs %lg) -- %lg vs %lg\n",
                 n, tf[best_match], reff[n], abs(herea[best_match]), abs(refa[n]),
                 err_best, err_limit);
        /*printf("OOOOOOOOO\n");
        printf("---------\n");
        if (err_best > 0.02) {
          printf("Didn't find a nice frequency! (%lg)\n", err_best);
        } else {
          printf("Found a nice frequency! (%lg)\n", err_best);
        }
        printf("Ref %d: %10lg %10lg\t(%10lg ,%10lg)\n",
               n, reff[n], refd[n], refa[n]);
        for (int i=0;i<num_match;i++) {
          printf("%5d: %10lg %10lg\t(%10lg ,%10lg)\n",
                 i, tf[i], td[i], ta[i]);
        } 
        printf("---------\n");
        printf("OOOOOOOOO\n");*/
      } else if (err_best < err_limit) {
        if (verbosity > 2 && abs(herea[best_match]) != 0.0)
          printf("Found a match at %d between %lg and %lg (%lg vs %lg) -- %lg vs %lg\n",
                 n, tf[best_match], reff[n], abs(herea[best_match]), abs(refa[n]),
                 err_best, err_limit);
        if (abs(herea[best_match]) > abs(refa[n])) { // Change reference...
          if (verbosity > 1) printf("Changing reference %d... (%lg vs %lg)\n",
                                    n, err_best, err_limit);
          if (verbosity > 1) printf("Freq goes from (%lg,%lg) to (%lg,%lg).\n",
                                    reff[n], refd[n], tf[best_match], td[best_match]);
          //printf("best_err is %lg\n", err_best);
          //printf("amp (%lg,%lg) (%lg,%lg)\n", 
          //       real(refa[n]),imag(refa[n]),
          //       real(ta[best_match]), imag(ta[best_match]));
          reff[n] = tf[best_match];
          refd[n] = td[best_match];
          refa[n] = herea[best_match];
          for (int t=0;t<ntime;t++)
            refdata[t+n*ntime] = simple_data[t];
        }
      }
    }
    int num_here = get_freqs(simple_data, ntime, herea, heref, hered);
    if (num_here > numref || 1) {
      if (verbosity > 1) printf("Deltaf here is %lg\n", deltaf);
      // It looks like we see a new mode at this point...
      int *refnum = new int[maxbands];
      for (int i=0;i<num_here;i++) refnum[i] = -1;
      for (int n=0;n<numref;n++) {
        int best_match=-1;
        double err_best = 1e300;
        for (int i=0;i<num_here;i++) {
          double err = abs(heref[i]-reff[n])+0.1*abs(hered[i]-refd[n]);
          //printf("heref[%d] %lg vs reff[%d] %lg gives %lg\n",
          //       i, heref[i], n, reff[n], err);
          if (err < err_best && refnum[i] == -1) {
            best_match = i;
            err_best = err;
          }
        }
        const double err_limit = deltaf;
        if (err_best < err_limit) {
          refnum[best_match] = n;
          if (verbosity > 1)
            printf("Matched %d: %10lg Got a best err of %8lg on an f of %lg %d (%lg)\n",
                   n, heref[best_match], err_best, reff[n], n, abs(herea[best_match]));
        } else if (err_best < 1e299) {
          if (verbosity > 1)
            printf("Missed %d:  %10lg Got a best err of %8lg on an f of %lg %d (%lg)\n",
                   n, heref[best_match], err_best, reff[n], n, abs(herea[best_match]));
        }
      }
      for (int i=0;i<num_here;i++) {
        if (refnum[i] == -1) { // New mode!!! Change reference...
          reff[numref] = heref[i];
          refd[numref] = hered[i];
          if (verbosity > 1)
            printf("Found one more mode (was i == %d)! (%10lg,%10lg) -- %d\n",
                   i, heref[i], refd[numref], numref);
          refa[numref] = herea[i];
          for (int t=0;t<ntime;t++) 
            refdata[t+numref*ntime] = simple_data[t];
          numref++;
          if (numref > maxbands-2) numref = maxbands-2;
        }
      }
      delete[] refnum;
    }
    delete[] ta;
    delete[] tf;
    delete[] td;
    delete[] herea;
    delete[] heref;
    delete[] hered;
  }
  return numref;
}

static inline int get_closest(double f[], int fmax) {
  double deltamin = 1e300;
  for (int i=0;i<fmax-1;i++) {
    if (f[i+1]-f[i] < deltamin) deltamin = f[i+1]-f[i];
  }
  for (int i=0;i<fmax-1;i++) {
    if (f[i+1]-f[i] == deltamin) return i;
  }
  return 0;
}

static inline int go_higher(double f[], int fmax, int lo, int hi) {
  if (lo == 0) return 1;
  if (hi == fmax-1) return 0;
  if (f[lo]-f[lo-1] > f[hi+1]-f[hi]) return 1;
  else return 0;
}

static inline int am_done(double f[], int fmax, int lo, int hi) {
  double wid = f[hi]-f[lo] + 0.001;
  return f[lo]-f[lo-1] > wid && f[hi+1]-f[hi] > wid;
}

static void get_cluster(double f[], int fmax, int maxsize, double maxwid,
                        int *out_lo, int *out_hi) {
  int lo = get_closest(f,fmax);
  int hi = lo+1;
  int minsize = maxsize/2+1;
  if (minsize < 3) minsize = 3;
  for (int i=0;i<minsize-2;i++) {
    if (go_higher(f,fmax,lo,hi)) hi++;
    else lo--;
  }
  while (hi + 1 - lo < maxsize) {
    if (am_done(f,fmax,lo,hi)) break;
    if (go_higher(f,fmax,lo,hi)) {
      if (f[hi+1]-f[lo] > maxwid) break;
      hi++;
    } else {
      if (f[hi]-f[lo-1] > maxwid) break;
      lo--;
    }
  }
  *out_lo = lo;
  *out_hi = hi;
}

int fields::cluster_some_bands_cleverly(double *tf, double *td, complex<double> *ta,
                                        int num_freqs, int fields_considered, int maxbands,
                                        complex<double> *fad, double *approx_power) {
  const double total_time = (bands->tend-bands->tstart)*c/a;
  const double deltaf = 1.0/total_time;
  int freqs_so_far = num_freqs;
  printf("About to sort by frequency... (%d frequencies)\n", freqs_so_far);
  // Sort by frequency...
  for (int i = 1; i < freqs_so_far; i++) {
    for (int j=i; j>0;j--) {
      if (tf[j]<tf[j-1]) {
        double t1 = tf[j], t2 = td[j];
        tf[j] = tf[j-1];
        td[j] = td[j-1];
        tf[j-1] = t1;
        td[j-1] = t2;
        complex<double> temp = ta[j];
        ta[j] = ta[j-1];
        ta[j-1] = temp;
      }
    }
  }
  printf("Looking for clusters...\n");
  int num_found = 0;
  double totwid = 0.001;
  while (freqs_so_far >= fields_considered/2 + 1) {
    int hi, lo;
    get_cluster(tf,freqs_so_far,fields_considered,deltaf,&lo,&hi);
    int mid = lo + (hi-lo)/2;
    if (tf[hi]-tf[lo] < deltaf) {
      printf("Got a cluster from %lg to %lg (%d freqs)\n", tf[lo], tf[hi], 1+hi-lo);
      fad[num_found] = complex<double>(tf[mid],td[mid]);
      if (approx_power) {
        approx_power[num_found] = 0;
        for (int i=lo;i<=hi;i++) {
          if (abs(ta[i])*abs(ta[i]) > approx_power[num_found]) {
            approx_power[num_found] = abs(ta[i])*abs(ta[i]);
          }
        }
      }
      totwid += tf[hi]-tf[lo];
      num_found++;
      if (num_found >= maxbands) num_found--;
    } else {
      printf("Rejected a cluster from %lg to %lg (%d freqs)\n", tf[lo], tf[hi], 1+hi-lo);
      if (verbosity > 1) printf("width is %g vs %g\n", tf[hi] - tf[lo], deltaf);
      lo = get_closest(tf,freqs_so_far);
      hi = lo+1;
      if (verbosity > 1) printf("dropping %g and %g\n", tf[hi], tf[lo]);
    }
    freqs_so_far -= 1 + hi - lo;
    for (int i=lo;i<freqs_so_far;i++) {
      tf[i] = tf[i + 1 + hi - lo];
      td[i] = td[i + 1 + hi - lo];
      ta[i] = ta[i + 1 + hi - lo];
    }
  }
  for (int i=0;i<freqs_so_far;i++) {
    if (verbosity > 1) printf("Have a leftover freq: %g\n", tf[i]);
  }
  return num_found;
}

complex<double> *fields::clever_cluster_bands(int maxbands, double *approx_power) {
  const double total_time = (bands->tend-bands->tstart)*c/a;
  const double deltaf = 1.0/total_time;
  bands->maxbands = maxbands;
  const int max_harminvs = 120;
  const int max_freqs = max_harminvs*maxbands;
  double *tf = new double[max_freqs];
  double *td = new double[max_freqs];
  cmplx *ta = new cmplx[max_freqs];
  const int ntime = bands->ntime;
  cmplx *simple_data = new cmplx[ntime];
  if (!simple_data || !ta) {
    printf("Error allocating...\n");
    exit(1);
  }
  /*
  freqs_here = bands->get_freqs(bands->P, ntime,
                               ta+freqs_so_far,tf+freqs_so_far,td+freqs_so_far);
  if (freqs_here) {
    fields_considered++;
    freqs_so_far += freqs_here;
  }
  */
  int num_found = 0;
  complex<double> *fad = new complex<double>[maxbands];
  for (int i=0;i<maxbands;i++) fad[i] = 0.0;
  if (k==0.0 && nz == 1) {
    // Treat special case when we're looking at pure TE/TM
    int freqs_so_far = 0;
    int fields_considered = 0;

    for (int r=0;r<nr;r+=1+(int)(bands->scale_factor/c*1.99)) {
      for (int t=0;t<ntime;t++) simple_data[t] = BAND(bands->er,r,t);
      int freqs_here = bands->get_freqs(simple_data, ntime,
                                        ta+freqs_so_far,tf+freqs_so_far,td+freqs_so_far);
      if (freqs_here) {
        fields_considered++;
        freqs_so_far += freqs_here;
      }
      for (int t=0;t<ntime;t++) simple_data[t] = BAND(bands->ep,r,t);
      freqs_here = bands->get_freqs(simple_data, ntime,
                                    ta+freqs_so_far,tf+freqs_so_far,td+freqs_so_far);
      if (freqs_here) {
        fields_considered++;
        freqs_so_far += freqs_here;
      }
      for (int t=0;t<ntime;t++) simple_data[t] = BAND(bands->hz,r,t);
      freqs_here = bands->get_freqs(simple_data, ntime,
                                    ta+freqs_so_far,tf+freqs_so_far,td+freqs_so_far);
      if (freqs_here) {
        fields_considered++;
        freqs_so_far += freqs_here;
      }
    }
    num_found = cluster_some_bands_cleverly(tf, td, ta, freqs_so_far, fields_considered,
                                            maxbands, fad, approx_power);
    freqs_so_far = 0;
    fields_considered = 0;
    for (int r=0;r<nr;r+=1+(int)(bands->scale_factor/c*1.99)) {
      for (int t=0;t<ntime;t++) simple_data[t] = BAND(bands->hr,r,t);
      int freqs_here = bands->get_freqs(simple_data, ntime,
                                        ta+freqs_so_far,tf+freqs_so_far,td+freqs_so_far);
      if (freqs_here) {
        fields_considered++;
        freqs_so_far += freqs_here;
      }
      for (int t=0;t<ntime;t++) simple_data[t] = BAND(bands->hp,r,t);
      freqs_here = bands->get_freqs(simple_data, ntime,
                                    ta+freqs_so_far,tf+freqs_so_far,td+freqs_so_far);
      if (freqs_here) {
        fields_considered++;
        freqs_so_far += freqs_here;
      }
      for (int t=0;t<ntime;t++) simple_data[t] = BAND(bands->ez,r,t);
      freqs_here = bands->get_freqs(simple_data, ntime,
                                    ta+freqs_so_far,tf+freqs_so_far,td+freqs_so_far);
      if (freqs_here) {
        fields_considered++;
        freqs_so_far += freqs_here;
      }
    }
    num_found += cluster_some_bands_cleverly(tf, td, ta, freqs_so_far, fields_considered,
                                             maxbands-num_found,
                                             fad+num_found, approx_power+num_found);
  } else {
    int freqs_so_far = 0;
    int fields_considered = 0;

    for (int r=0;r<nr;r+=1+(int)(bands->scale_factor/c*1.99)) {
      cmplx *bdata;
      for (int whichf = 0; whichf < 6; whichf++) {
        if (verbosity>1) printf("Looking at r == %lg, field %d\n", r*inva, whichf);
        switch (whichf) {
        case 0: bdata = bands->er; break;
        case 1: bdata = bands->ep; break;
        case 2: bdata = bands->ez; break;
        case 3: bdata = bands->hr; break;
        case 4: bdata = bands->hp; break;
        case 5: bdata = bands->hz; break;
        }
        for (int t=0;t<ntime;t++) simple_data[t] = BAND(bdata,r,t);
        int freqs_here = bands->get_freqs(simple_data, ntime,
                                          ta+freqs_so_far,tf+freqs_so_far,td+freqs_so_far);
        if (freqs_here) {
          fields_considered++;
          freqs_so_far += freqs_here;
        }
        if (freqs_so_far + maxbands > max_freqs) break;
      }
    }
    num_found = cluster_some_bands_cleverly(tf, td, ta, freqs_so_far, fields_considered,
                                            maxbands, fad, approx_power);
  }
    
  delete[] ta;
  delete[] tf;
  delete[] td;
  // Sorting by frequency again...
  for (int i=0;i<num_found;i++) {
    for (int j=i; j>0;j--) {
      if (real(fad[j])<real(fad[j-1])) {
        complex<double> t1 = fad[j];
        fad[j] = fad[j-1];
        fad[j-1] = t1;
        double temp = approx_power[j];
        approx_power[j] = approx_power[j-1];
        approx_power[j-1] = temp;
      }
    }
  }
  return fad;
}

complex<double> *fields::get_the_bands(int maxbands, double *approx_power) {
  bands->maxbands = maxbands;
  double *tf = new double[maxbands];
  double *td = new double[maxbands];
  cmplx *ta = new cmplx[maxbands];
  double *heref = new double[maxbands];
  double *hered = new double[maxbands];
  cmplx *herea = new cmplx[maxbands];
  const int ntime = bands->ntime;

  cmplx *eigen = new cmplx[nr*maxbands*6];
  cmplx *simple_data = new cmplx[ntime];
  if (!eigen || !simple_data) {
    printf("Error allocating...\n");
    exit(1);
  }

  for (int r=0;r<nr;r++) {
    for (int whichf = 0; whichf < 6; whichf++) {
      for (int n=0;n<maxbands;n++) {
        HARMOUT(eigen,r,n,whichf) = 0;
      }
    }
  }
  int *refnum = new int[maxbands];
  int *refr = new int[maxbands], *refw = new int[maxbands], numref = 0;
  double *reff = new double[maxbands], *refd = new double[maxbands];
  cmplx *refa = new complex<double>[maxbands];
  cmplx *refdata = new complex<double>[maxbands*ntime];
  if ((k==0.0 || m==0) && nz==1) { // We have pure TE and TM modes!
    for (int t=0;t<ntime;t++)
      simple_data[t] = BAND(bands->ez,nr/2,t);
    numref = bands->look_for_more_bands(simple_data, reff, refd, refa, refdata, numref);
    for (int t=0;t<ntime;t++)
      simple_data[t] = BAND(bands->hz,nr/2,t);
    numref = bands->look_for_more_bands(simple_data, reff, refd, refa, refdata,-numref);
  }
  numref = bands->look_for_more_bands(bands->P, reff, refd, refa, refdata, numref);
  if (numref && verbosity) printf("I found %d bands in the polarization...\n", numref);
  for (int r=0;r<nr;r+=1+(int)(bands->scale_factor/c*3.99)) {
    cmplx *bdata;
    for (int whichf = 0; whichf < 6; whichf++) {
      if (verbosity>1) printf("Looking at r == %lg, field %d\n", r*inva, whichf);
      switch (whichf) {
      case 0: bdata = bands->er; break;
      case 1: bdata = bands->ep; break;
      case 2: bdata = bands->ez; break;
      case 3: bdata = bands->hr; break;
      case 4: bdata = bands->hp; break;
      case 5: bdata = bands->hz; break;
      }
      for (int t=0;t<ntime;t++) {
        simple_data[t] = BAND(bdata,r,t);
      }
      numref = bands->look_for_more_bands(simple_data, reff, refd, refa, refdata, numref);
    }
  }
  delete[] refdata;
  delete[] tf;
  delete[] td;
  delete[] ta;
  delete[] herea;
  delete[] heref;
  delete[] hered;

  // Sort by frequency...
  for (int i = 1; i < numref; i++) {
    for (int j=i; j>0;j--) {
      if (reff[j]<reff[j-1]) {
        double t1 = reff[j], t2 = refd[j];
        reff[j] = reff[j-1];
        refd[j] = refd[j-1];
        reff[j-1] = t1;
        refd[j-1] = t2;
        complex<double> temp = refa[j];
        refa[j] = refa[j-1];
        refa[j-1] = temp;
      }
    }
  }

  // Now get rid of any really low power solutions...
  for (int i=0;i<numref;i++)
    refa[i] *= sqrt((refd[i]<0)?exp(refd[i]*2*pi*ntime/a*c):1);
  double powmax = 0.0;
  for (int i=0;i<numref;i++)
    powmax = max(powmax,abs(refa[i]*refa[i]));
  double powmin = powmax*bands->fpmin;
  for (int i=numref-1;i>=0;i--) {
    if (abs(refa[i]*refa[i]) < powmin) {
      numref--;
      if (verbosity > 2)
        printf("Trashing a low power solution with freq %lg %lg (%lg vs %lg)\n",
               reff[i], refd[i], abs(refa[i]*refa[i]), powmin);
      for (int j=i;j<numref;j++) {
        reff[j] = reff[j+1];
        refd[j] = refd[j+1];
        refa[j] = refa[j+1];
      }
    }
  }

  complex<double> *fad = new complex<double>[maxbands];
  for (int i=0;i<maxbands;i++) fad[i] = 0.0;
  for (int i=0;i<numref;i++) fad[i] = complex<double>(reff[i],refd[i]);
  if (approx_power) for (int i=0;i<numref;i++) approx_power[i] = abs(refa[i])*abs(refa[i]);

  delete[] reff;
  delete[] refd;
  delete[] eigen;
  delete[] refa;
  return fad;
}

int bandsdata::get_both_freqs(cmplx *data1, cmplx *data2, int n,
                              cmplx *amps1, cmplx *amps2, 
                              double *freqs, double *decays) {
  const double total_time = (tend-tstart)*c/a;
  const double deltaf = 2.0/total_time;
  double phi = 0.5;
  int numfound = 0;
  double mag1 = 0, mag2 = 0;
  for (int i=0;i<n;i++)
    mag1 += norm(data1[i]); // norm(a) is actually sqr(abs(a))
  for (int i=0;i<n;i++)
    mag2 += norm(data2[i]); // norm(a) is actually sqr(abs(a))
  int times_through = 0;
  do {
    complex<double> shift = polar(2.0,phi);
    complex<double> unshift = polar(0.5,-phi);
    if (verbosity > 3 && phi != 0.5) printf("CHANGING PHI! (1.0,%lg)\n",phi);
    cmplx *plus = new complex<double>[n];
    cmplx *minus = new complex<double>[n];
    cmplx *Ap = new complex<double>[n];
    cmplx *Am = new complex<double>[n];
    double *fp = new double[maxbands];
    double *fm = new double[maxbands];
    double *dp = new double[maxbands];
    double *dm = new double[maxbands];
    int plusboring = 1;
    int minusboring = 1;
    for (int i=0;i<n;i++) {
      plus[i] = data1[i]+shift*data2[i];
      minus[i] = data1[i]-shift*data2[i];
      if (plus[i] != 0) plusboring = 0;
      if (minus[i] != 0) minusboring = 0;
    }
    if (verbosity > 2 && (plusboring || minusboring))
      printf("This is a boring point...\n");
    if (!plusboring && !minusboring) {
      int numplus = get_freqs(plus, n, Ap, fp, dp);
      int numminus = get_freqs(minus, n, Am, fm, dm);
      if (verbosity > 3) printf("We got %d plus and %d minus...\n", numplus, numminus);
      if (numplus == numminus) {
        // Looks like we agree on the number of bands...
        numfound = numplus;
        //printf("Mags: %10lg/%10lg\n", mag1, mag2);
        //for (int i=0;i<numfound;i++) {
        //  printf("%10lg/%10lg %10lg/%10lg\t(%11lg,%11lg)/(%11lg,%11lg)\n",
        //         fp[i], fm[i], dp[i], dm[i],
        //         0.5*(Ap[i]+Am[i]), 0.5*(Ap[i]-Am[i])*unshift);
        //}
        for (int i=0;i<numfound;i++) {
          if (0.5*(fp[i]-fm[i]) > 0.01*deltaf) {
            numplus--; // Try another phi...
            if (verbosity > 3)
              printf("We've got some weird frequencies: %lg and %lg\n",fp[i],fm[i]);
          }
          freqs[i]  = (abs(Ap[i])>abs(Am[i]))?fp[i]:fm[i];//0.5*(fp[i]+fm[i]);
          decays[i] = (abs(Ap[i])>abs(Am[i]))?dp[i]:dm[i];//0.5*(dp[i]+dm[i]);
          amps1[i] = 0.5*(Ap[i]+Am[i]);
          amps2[i] = 0.5*(Ap[i]-Am[i])*unshift;
        }
      }
    }
    delete[] plus;
    delete[] minus;
    delete[] Ap;
    delete[] Am;
    delete[] fp;
    delete[] fm;
    delete[] dp;
    delete[] dm;
    phi += 0.1;
    times_through++;
    if (times_through > 50) return 0;
  } while (numfound == 0);
  return numfound;
}

int bandsdata::get_freqs(cmplx *data, int n,
                         cmplx *amps, double *freq_re, double *freq_im) {
  
  int num = do_harminv(data, n, scale_factor, a, fmin, fmax, maxbands, amps, 
		       freq_re, freq_im);
  // First deal with any negative frequency solutions.
  const double total_time = n*scale_factor*c/a;
  for (int i=0;i<num-1;i++) {
    if (freq_re[i] < 0) {
      for (int j=i+1;j<num;j++) {
        if (abs(freq_re[j]+freq_re[i]) < 2.0/total_time) {
          if (verbosity > 2 && freq_re[i] != 0.0) {
            printf("Got a plus/minus freq match at %lg\n",freq_re[j]); 
            printf("Total time: %lg and delta freq limit %lg\n",
                   total_time, 2.0/total_time);
          }
          freq_re[i] = -0.0; // It will get cleaned up later...
        }
      }
      freq_re[i] = -freq_re[i];
      if (verbosity > 2 && freq_re[i] != 0.0)
        printf("Flipping sign of a negative freq:  %lg %lg\n", freq_re[i], freq_im[i]);
    }
  }
  // Now sort the silly solutions again...
  for (int i=0;i<num-1;i++) { // This is a really bad sort algorithm...
    for (int j=i+1;j<num;j++) {
      if (freq_re[i] > freq_re[j]) {
        double t = freq_re[i];
        freq_re[i] = freq_re[j];
        freq_re[j] = t;
        t = freq_im[i];
        freq_im[i] = freq_im[j];
        freq_im[j] = t;
        complex<double> tc = amps[i];
        amps[i] = amps[j];
        amps[j] = tc;
      }
    }
  }
  // Now get rid of any spurious low frequency solutions...
  int orignum = num;
  for (int i=0;i<orignum;i++) {
    if (freq_re[0] < fmin*.9) {
      if (verbosity > 2 && freq_re[0] != 0.0) {
        printf("Trashing a spurious low frequency solution with freq %lg %lg\n",
               freq_re[0], freq_im[0]);
        //printf("For your info, fmin is %lg\n", fmin);
      }
      for (int j=0;j<num-1;j++) {
        freq_re[j]=freq_re[j+1];
        freq_im[j]=freq_im[j+1];
        amps[j]=amps[j+1];
      }
      num--;
    }
  }
  // Now get rid of any spurious transient solutions...
  for (int i=num-1;i>=0;i--) {
    double qminhere = 1.0/(1.0/qmin + 0.25/(freq_re[i]*total_time));
    double qhere = 0.5*fabs(fabs(freq_re[i])/freq_im[i]);
    if (qhere < qminhere) {
      num--;
      if (verbosity > 2) {
        printf("Trashing a spurious low Q solution with freq %lg %lg (%lg vs %lg)\n",
               freq_re[i], freq_im[i], qhere, qminhere);
      }
      for (int j=i;j<num;j++) {
        freq_re[j] = freq_re[j+1];
        freq_im[j] = freq_im[j+1];
        amps[j] = amps[j+1];
      }
    }
  }
  return num;
}

int do_harminv(cmplx *data, int n, int sampling_rate, double a, 
	       double fmin, double fmax, int maxbands,
	       cmplx *amps, double *freq_re, double *freq_im, double *errors) {
  // data is a size n array.

  // check for all zeros in input
  {
    int all_zeros = 1;
    for (int i=0; i<n; i++)
      if (data[i] != 0) all_zeros = 0;
    if (all_zeros)
      return 0;
  }

  harminv_data hd = 
    harminv_data_create(n, data, fmin*sampling_rate*c/a, fmax*sampling_rate*c/a, maxbands);

  int prev_nf, cur_nf;
  harminv_solve(hd);
  prev_nf = cur_nf = harminv_get_num_freqs(hd);

  /* keep re-solving as long as spurious solutions are eliminated */
  do {
    prev_nf = cur_nf;
    harminv_solve_again(hd);
    cur_nf = harminv_get_num_freqs(hd);
  } while (cur_nf < prev_nf);
  if (cur_nf > prev_nf)
    fprintf(stderr,
            "harminv: warning, number of solutions increased from %d to %d!\n",
            prev_nf, cur_nf);
  
  cmplx *tmpamps = harminv_compute_amplitudes(hd);
  double *tmperrors = harminv_compute_frequency_errors(hd);

  freq_re[0] = a*harminv_get_freq(hd, 0)/c/sampling_rate;
  freq_im[0] = -1/(2*pi)*a*harminv_get_decay(hd, 0)/c/sampling_rate;
  for (int i = 1; i < harminv_get_num_freqs(hd); ++i) {
    freq_re[i] = a*harminv_get_freq(hd, i)/c/sampling_rate;
    freq_im[i] = -1/(2*pi)*a*harminv_get_decay(hd, i)/c/sampling_rate;
    for (int j=i; j>0;j--) {
      if (freq_re[j]<freq_re[j-1]) {
        double t1 = freq_re[j], t2 = freq_im[j], e = tmperrors[j];
        cmplx a = tmpamps[j];
        tmpamps[j] = tmpamps[j-1];
        tmperrors[j] = tmperrors[j-1];
        freq_re[j] = freq_re[j-1];
        freq_im[j] = freq_im[j-1];
        freq_re[j-1] = t1;
        freq_im[j-1] = t2;
        tmpamps[j-1] = a;
        tmperrors[j-1] = e;
      }
    }
  }
  int num = harminv_get_num_freqs(hd);
  for (int i = 0; i < num; ++i) {
    amps[i] = tmpamps[i];
    if (errors)
      errors[i] = tmperrors[i];
  }
  free(tmpamps);
  free(tmperrors);
  harminv_data_destroy(hd);
  return num;
}


