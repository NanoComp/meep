#include "meep.hpp"
#include "meep_internals.hpp"
#include "config.h"

#define DPR double * restrict
#define RPR realnum * restrict

/* These macros get into the guts of the LOOP_OVER_VOL loops to efficiently
   construct the index k into a PML sigma array.  Basically, k needs to
   increment by 2 for each increment of one of LOOP's for-loops, starting
   at the appropriate corner of the volume, and these macros define the
   relevant strides etc. for each loop.  If sigsize_dsig <= 1, however,
   our index k should always be zero (sigma array is trivial length-1). 
   KSTRIDE_DEF defines the relevant strides etc. and goes outside the 
   LOOP, wheras KDEF defines the k index and goes inside the LOOP. */
#define SZ0(dsig, expr) (sigsize_ ## dsig > 1 ? (expr) : 0)
#define KSTRIDE_DEF(dsig, k, corner)					  \
     const int k##0 = SZ0(dsig, corner.in_direction(dsig)		  \
                              - v.little_corner().in_direction(dsig));	  \
     const int s##k##1 = SZ0(dsig, v.yucky_direction(0) == dsig ? 2 : 0); \
     const int s##k##2 = SZ0(dsig, v.yucky_direction(1) == dsig ? 2 : 0); \
     const int s##k##3 = SZ0(dsig, v.yucky_direction(2) == dsig ? 2 : 0)
#define KDEF(k,dsig) const int k = ((k##0 + s##k##1*loop_i1) + s##k##2*loop_i2) + s##k##3*loop_i3
#define DEF_k KDEF(k,dsig)

namespace meep {

#define SWAP(t,a,b) { t xxxx = a; a = b; b = xxxx; }

/* update step for df/dt = curl g,
   i.e. f += dt curl g = dt/dx (dg1 - dg2)
   where dgk = gk[i] - gk[i+sk].

   g = (g1,g2), where g1 or g2 may be NULL.  Note that dt/dx and/or s1
   and s2 may be negative to flip signs of derivatives.

   PML: sig[k] = sigma[k]*dt/2, siginv[k] = 1 / (1 + sigma[k]*dt/2).
   Here, k is the index in the dsig direction.  if dsig ==
   NO_DIRECTION, then PML is not used.  (dsig is the sigma direction.)

   if non-NULL, then cnd is an array of conductivity values, changing
   the underlying PDE to:
       df/dt = curl g - cnd f
   which is updated as:
       f = [ dt * curl g + (1 - dt cnd/2) f ] / (1 + dt cnd/2)
   cndinv should be an array of 1 / (1 + dt cnd/2).  In the case
   of PML, cndinv should contain 1 / (1 + dt (cnd + sigma)/2).

   fcnd is an auxiliary field used ONLY when we simultaneously have
   PML and conductivity, in which case fcnd solves 
       dfcnd/dt = curl g - cnd*fcnd
   and f satisfies 
       df/dt = dfcnd/dt - sigma*f.
*/
void step_curl(RPR f, component c, const RPR g1, const RPR g2,
	       int s1, int s2, // strides for g1/g2 shift
	       const volume &v, double dtdx,
	       direction dsig, const DPR sig, const DPR siginv,
	       double dt, 
	       const RPR cnd, const RPR cndinv, RPR fcnd)
{
  if (!g1) { // swap g1 and g2
    SWAP(const RPR, g1, g2);
    SWAP(int, s1, s2);
    dtdx = -dtdx; // need to flip derivative sign
  }
  if (dsig == NO_DIRECTION) { // no PML
    if (cnd) {
      double dt2 = dt * 0.5;
      if (g2) {
	LOOP_OVER_VOL_OWNED0(v, c, i)
	  f[i] = ((1 - dt2 * cnd[i]) * f[i] - 
		  dtdx * (g1[i+s1] - g1[i] + g2[i] - g2[i+s2])) * cndinv[i];
      }
      else {
	LOOP_OVER_VOL_OWNED0(v, c, i)
	  f[i] = ((1 - dt2 * cnd[i]) * f[i] 
		  - dtdx * (g1[i+s1] - g1[i])) * cndinv[i];
      }
    }
    else { // no conductivity
      if (g2) {
	LOOP_OVER_VOL_OWNED0(v, c, i)
	  f[i] -= dtdx * (g1[i+s1] - g1[i] + g2[i] - g2[i+s2]);
      }
      else {
	LOOP_OVER_VOL_OWNED0(v, c, i)
	  f[i] -= dtdx * (g1[i+s1] - g1[i]);
      }
    }
  }
  else { /* PML */
    const int sigsize_dsig=2; KSTRIDE_DEF(dsig, k, v.little_owned_corner0(c));
    if (cnd) {
      double dt2 = dt * 0.5;
      if (g2) {
	LOOP_OVER_VOL_OWNED0(v, c, i) {
	  DEF_k;
	  realnum fcnd_prev = fcnd[i];
	  fcnd[i] = ((1 - dt2 * cnd[i]) * fcnd[i] - 
		     dtdx * (g1[i+s1] - g1[i] + g2[i] - g2[i+s2])) * cndinv[i];
	  f[i] = ((1 - sig[k]) * f[i] + (fcnd[i] - fcnd_prev)) * siginv[k];
	}
      }
      else {
	LOOP_OVER_VOL_OWNED0(v, c, i) {
	  DEF_k;
	  realnum fcnd_prev = fcnd[i];
	  fcnd[i] = ((1 - dt2 * cnd[i]) * fcnd[i] - 
		     dtdx * (g1[i+s1] - g1[i])) * cndinv[i];
	  f[i] = ((1 - sig[k]) * f[i] + (fcnd[i] - fcnd_prev)) * siginv[k];
	}
      }
    }
    else { // no conductivity (other than PML conductivity)
      if (g2) {
	LOOP_OVER_VOL_OWNED0(v, c, i) {
	  DEF_k;
	  f[i] = ((1 - sig[k]) * f[i] -
		  dtdx * (g1[i+s1] - g1[i] + g2[i] - g2[i+s2])) * siginv[k];
	}
      }
      else {
	LOOP_OVER_VOL_OWNED0(v, c, i) {
	  DEF_k;
	  f[i] = ((1 - sig[k]) * f[i] - dtdx * (g1[i+s1] - g1[i])) * siginv[k];
	}
      }
    }
  }
}

/* field-update equation f += betadt * g (plus variants for conductivity 
   and/or PML).  This is used in 2d calculations to add an exp(i beta z)
   time dependence, which gives an additional i \beta \hat{z} \times
   cross-product in the curl equations. */
void step_beta(RPR f, component c, const RPR g,
	       const volume &v, double betadt,
	       direction dsig, const DPR siginv,
	       const RPR cndinv, RPR fcnd)
{
  if (!g) return;
  if (dsig != NO_DIRECTION) { // PML
    const int sigsize_dsig=2;
    KSTRIDE_DEF(dsig, k, v.little_owned_corner0(c));
    if (cndinv) { // conductivity + PML
      LOOP_OVER_VOL_OWNED0(v, c, i) {
	DEF_k;
	double dfcnd = betadt * g[i] * cndinv[i];
	fcnd[i] += dfcnd;
	f[i] += dfcnd * siginv[k];
      }
    }
    else { // PML only
      LOOP_OVER_VOL_OWNED0(v, c, i) {
	DEF_k;
	f[i] += betadt * g[i] * siginv[k];
      }
    }
  }
  else { // no PML
    if (cndinv) { // conductivity, no PML
      LOOP_OVER_VOL_OWNED0(v, c, i)
	f[i] += betadt * g[i] * cndinv[i];
    }
    else { // no conductivity or PML
      LOOP_OVER_VOL_OWNED0(v, c, i)
	f[i] += betadt * g[i];
    }
  }
}

/* Given Dsqr = |D|^2 and Di = component of D, compute the factor f so
   that Ei = chi1inv * f * Di.   In principle, this would involve solving
   a cubic equation, but instead we use a Pade approximant that is 
   accurate to several orders.  This is inaccurate if the nonlinear
   index change is large, of course, but in that case the chi2/chi3
   power-series expansion isn't accurate anyway, so the cubic isn't
   physical there either. */
inline double calc_nonlinear_u(const double Dsqr, 
			       const double Di,
			       const double chi1inv,
			       const double chi2, const double chi3) {
  double c2 = Di*chi2*(chi1inv*chi1inv);
  double c3 = Dsqr*chi3*(chi1inv*chi1inv*chi1inv);
  return (1 + c2 + 2*c3)/(1 + 2*c2 + 3*c3);
}

/* Update E from D using epsilon and PML, *or* update H from B using
   mu and PML.

   To be generic, here we set f = u * g (for the non-PML), where u may
   be a tensor, and we also have a nonlinear susceptibility chi.
   Here, g = (g,g1,g2) where g1 and g2 are the off-diagonal
   components, if any (g2 may be NULL).

   Note that for the common case of mu = 1 and no PML, we don't even
   call this routine. 

   Here, sig = sigma[k]*dt/2, and siginv[k] = 1 / (1 + sig[k]), and
   sig2 is the other sigma array.  gb etc. are the backups of g from
   the previous time step. */
void step_update_EDHB(RPR f, component fc, const volume &v, 
		      const RPR g, const RPR g1, const RPR g2,
		      const RPR gb, const RPR g1b, const RPR g2b,
		      const RPR u, const RPR u1, const RPR u2,
		      int s, int s1, int s2,
		      const RPR chi2, const RPR chi3,
		      direction dsig, const DPR sig, const DPR siginv,
		      direction dsigg, const DPR sigg,
		      direction dsig1, const DPR sig1,
		      direction dsig1inv, const DPR sig1inv,
		      direction dsig2, const DPR sig2,
		      direction dsig2inv, const DPR sig2inv,
		      int sigsize_dsig,int sigsize_dsigg,int sigsize_dsig1)
{
  if (!f) return;
  int sigsize_dsig1inv = sigsize_dsigg;
  int sigsize_dsig2 = sigsize_dsig;
  int sigsize_dsig2inv = sigsize_dsig1;
  
  if ((!g1 && g2) || (g1 && g2 && !u1 && u2)) { /* swap g1 and g2 */
    SWAP(const RPR, g1, g2);
    SWAP(const RPR, g1b, g2b);
    SWAP(const RPR, u1, u2);
    SWAP(int, s1, s2);
    SWAP(const DPR, sig1, sig2);
    SWAP(const DPR, sig1inv, sig2inv);
    SWAP(direction, dsig1, dsig2);
    SWAP(direction, dsig1inv, dsig2inv);
    SWAP(int, sigsize_dsig1, sigsize_dsig2);
    SWAP(int, sigsize_dsig1inv, sigsize_dsig2inv);
  }

  // stable averaging of offdiagonal components
#define OFFDIAG(u,g,sx) (0.25 * ((g[i]+g[i-sx])*u[i] \
		   	       + (g[i+s]+g[(i+s)-sx])*u[i+s]))
  
  if (sigsize_dsig <= 1 && sigsize_dsigg <= 1 && 
      (!u1 || (sigsize_dsig1 <= 1 && sigsize_dsig1inv <= 1))) { // no PML
    if (u1 && u2) { // 3x3 off-diagonal u
      if (chi3) {
	LOOP_OVER_VOL_OWNED(v, fc, i) {
	  double g1s = g1[i]+g1[i+s]+g1[i-s1]+g1[i+(s-s1)];
	  double g2s = g2[i]+g2[i+s]+g2[i-s2]+g2[i+(s-s2)];
	  double gs = g[i]; double us = u[i];
	  f[i] = (gs * us + OFFDIAG(u1,g1,s1) + OFFDIAG(u2,g2,s2))
	    * calc_nonlinear_u(gs * gs + 0.0625 * (g1s*g1s + g2s*g2s),
			       gs, us, chi2[i], chi3[i]);	  
	}
      }
      else {
	LOOP_OVER_VOL_OWNED(v, fc, i) {
	  double gs = g[i]; double us = u[i];
	  f[i] = (gs * us + OFFDIAG(u1,g1,s1) + OFFDIAG(u2,g2,s2));
	}
      }
    }
    else if (u1) { // 2x2 off-diagonal u
      if (chi3) {
	LOOP_OVER_VOL_OWNED(v, fc, i) {
	  double g1s = g1[i]+g1[i+s]+g1[i-s1]+g1[i+(s-s1)];
	  double gs = g[i]; double us = u[i];
	  f[i] = (gs * us + OFFDIAG(u1,g1,s1))
	    * calc_nonlinear_u(gs * gs + 0.0625 * (g1s*g1s),
			       gs, us, chi2[i], chi3[i]);
	}
      }
      else {
	LOOP_OVER_VOL_OWNED(v, fc, i) {
	  double gs = g[i]; double us = u[i];
	  f[i] = (gs * us + OFFDIAG(u1,g1,s1));
	}
      }
    }
    else if (u2) { // 2x2 off-diagonal u
      abort("bug - didn't swap off-diagonal terms!?");
    }
    else { // diagonal u
      if (chi3) {
	if (g1 && g2) {
	  LOOP_OVER_VOL_OWNED(v, fc, i) {
	    double g1s = g1[i]+g1[i+s]+g1[i-s1]+g1[i+(s-s1)];
	    double g2s = g2[i]+g2[i+s]+g2[i-s2]+g2[i+(s-s2)];
	    double gs = g[i]; double us = u[i];
	    f[i] = (gs*us)*calc_nonlinear_u(gs*gs+0.0625*(g1s*g1s+g2s*g2s),
					gs, us, chi2[i], chi3[i]);
	  }
	}
	else if (g1) {
	  LOOP_OVER_VOL_OWNED(v, fc, i) {
	    double g1s = g1[i]+g1[i+s]+g1[i-s1]+g1[i+(s-s1)];
	    double gs = g[i]; double us = u[i];
	    f[i] = (gs*us)*calc_nonlinear_u(gs*gs + 0.0625*(g1s*g1s),
						gs, us, chi2[i], chi3[i]);
	  }
	}
	else if (g2) {
	  abort("bug - didn't swap off-diagonal terms!?");
	}
	else {
	  LOOP_OVER_VOL_OWNED(v, fc, i) {
	    double gs = g[i]; double us = u[i];
	    f[i] = (gs*us)*calc_nonlinear_u(gs*gs, gs,us, chi2[i],chi3[i]);
	  }
	}
      }
      else if (u) {
	LOOP_OVER_VOL_OWNED(v, fc, i) {
	  double gs = g[i]; double us = u[i];
	  f[i] = (gs * us);
	}
      }
      else
	LOOP_OVER_VOL_OWNED(v, fc, i) f[i] = g[i];
    }
  }
  else { // PML
    
    // strides etc. for updating each sig[] index inside the LOOP:
    KSTRIDE_DEF(dsig, k, v.little_owned_corner(fc));
    KSTRIDE_DEF(dsigg, kg, v.little_owned_corner(fc));
    KSTRIDE_DEF(dsig1, k1, v.little_owned_corner(fc));
    KSTRIDE_DEF(dsig2, k2, v.little_owned_corner(fc));
    KSTRIDE_DEF(dsig1inv, k1inv, v.little_owned_corner(fc));
    KSTRIDE_DEF(dsig2inv, k2inv, v.little_owned_corner(fc));

   // the following definitions are used over and over

   // indices into sigma arrays:
#  define DEF_kx(x) KDEF(k ## x, dsig ## x)

   // fields
#  define DEF_gs  double gs0 = g[i]; double gs = ((1+sigg[kg])*siginv[k])*gs0
#  define DEF_gbs double gbs = ((1-sigg[kg])*siginv[k]) * gb[i]
#  define DEF_gss  double gs0 = g[i]; double gss = (1+sigg[kg]) * gs0
#  define DEF_gbss double gbss = (1-sigg[kg]) * gb[i]
#  define DEF_us  double us = u[i]
#  define DEF_g1s0 double g1s0 = g1[i]+g1[i+s]+g1[i-s1]+g1[i+(s-s1)];
#  define DEF_g2s0 double g2s0 = g2[i]+g2[i+s]+g2[i-s2]+g2[i+(s-s2)];
#  define SIG1 ((1+sig1[k1]) * sig1inv[k1inv])
#  define SIG2 ((1+sig2[k2]) * sig2inv[k2inv])
#  define SIG1b ((1-sig1[k1]) * sig1inv[k1inv])
#  define SIG2b ((1-sig2[k2]) * sig2inv[k2inv])

   if (u1 && u2) { // 3x3 off-diagonal u
     if (chi3) {
       LOOP_OVER_VOL_OWNED(v, fc, i) {
	  DEF_kx(1); DEF_kx(1inv); DEF_g1s0;
	  DEF_kx(2); DEF_kx(2inv); DEF_g2s0;
	  DEF_k; DEF_kx(g); DEF_gs; DEF_gbs;
	  DEF_us;
	  f[i] = ((1-sig[k])*siginv[k]) * f[i] + 
	    ((gs-gbs) * us
	     + SIG1*OFFDIAG(u1,g1,s1) - SIG1b*OFFDIAG(u1,g1b,s1)
	     + SIG2*OFFDIAG(u2,g2,s2) - SIG1b*OFFDIAG(u2,g2b,s2))
	    * calc_nonlinear_u(gs0 * gs0 + 0.0625 * (g1s0*g1s0 + g2s0*g2s0),
			       gs0, us, chi2[i], chi3[i]);
	}
      } else {
	LOOP_OVER_VOL_OWNED(v, fc, i) {
	  DEF_kx(1); DEF_kx(1inv);
          DEF_kx(2); DEF_kx(2inv);
          DEF_k; DEF_kx(g); DEF_gs; DEF_gbs;
          DEF_us;
	  f[i] = ((1-sig[k])*siginv[k]) * f[i] + 
	    ((gs-gbs) * us
	     + SIG1*OFFDIAG(u1,g1,s1) - SIG1b*OFFDIAG(u1,g1b,s1)
	     + SIG2*OFFDIAG(u2,g2,s2) - SIG1b*OFFDIAG(u2,g2b,s2));
	}
      }
    } else if (u1) { // 2x2 off-diagonal u
      if (chi3) {
	LOOP_OVER_VOL_OWNED(v, fc, i) {
	  DEF_kx(1); DEF_kx(1inv); DEF_g1s0;
	  DEF_k; DEF_kx(g); DEF_gs; DEF_gbs;
          DEF_us;
	  f[i] = ((1-sig[k])*siginv[k]) * f[i] + 
	    ((gs-gbs) * us
	     + SIG1*OFFDIAG(u1,g1,s1) - SIG1b*OFFDIAG(u1,g1b,s1))
	    * calc_nonlinear_u(gs0 * gs0 + 0.0625 * (g1s0*g1s0),
			       gs0, us, chi2[i], chi3[i]);
	}
      } else {
	LOOP_OVER_VOL_OWNED(v, fc, i) {
	  DEF_kx(1); DEF_kx(1inv);
          DEF_k; DEF_kx(g); DEF_gs; DEF_gbs;
          DEF_us;
	  f[i] = ((1-sig[k])*siginv[k]) * f[i] + 
	    ((gs-gbs) * us
	     + SIG1*OFFDIAG(u1,g1,s1) - SIG1b*OFFDIAG(u1,g1b,s1));
	}	
      }
    } else if (u2) { // 2x2 off-diagonal u
      abort("bug - didn't swap off-diagonal terms!?");
    } else if (u) { // diagonal u
      if (chi3) {
	if (g1 && g2) {
	  LOOP_OVER_VOL_OWNED(v, fc, i) {
	    DEF_g1s0; DEF_g2s0;
	    DEF_k; DEF_kx(g); DEF_gss; DEF_gbss;
	    DEF_us;
	    f[i] = siginv[k] * ((1-sig[k])*f[i] + (gss-gbss)*us *
	      calc_nonlinear_u(gs0 * gs0 + 0.0625 * (g1s0*g1s0 + g2s0*g2s0),
			       gs0, us, chi2[i], chi3[i]));
	  }
	} else if (g1) {
	  LOOP_OVER_VOL_OWNED(v, fc, i) {
	    DEF_g1s0;
            DEF_k; DEF_kx(g); DEF_gss; DEF_gbss;
            DEF_us;
	    f[i] = siginv[k] * ((1-sig[k])*f[i] + (gss-gbss)*us *
     	      calc_nonlinear_u(gs0 * gs0 + 0.0625 * (g1s0*g1s0),
     			       gs0, us, chi2[i], chi3[i]));
	  }
	} else if (g2) {
	  abort("bug - didn't swap off-diagonal terms!?");
	} else {
	  LOOP_OVER_VOL_OWNED(v, fc, i) {
	    DEF_k; DEF_kx(g); DEF_gss; DEF_gbss;
            DEF_us;
	    f[i] = siginv[k] * ((1-sig[k])*f[i] + (gss-gbss)*us *
     	      calc_nonlinear_u(gs0 * gs0, gs0, us, chi2[i], chi3[i]));
	  }
	}
      } else { //linear, diagonal u
	LOOP_OVER_VOL_OWNED(v, fc, i) {
	  DEF_k; DEF_kx(g); DEF_gss; DEF_gbss;
	  DEF_us;
	  f[i] = siginv[k] * ((1-sig[k])*f[i] + (gss-gbss) * us);
	}
      }
    }
    else { // NULL u array, corresponding to u = 1 everywhere
      if (chi3) abort("bug - should not have chi3 without chi1");
      // since this case is so common, do a few special cases:
      if (sigsize_dsig > 1 && sigsize_dsigg > 1)
	LOOP_OVER_VOL_OWNED(v, fc, i) {
	  DEF_k; DEF_kx(g);
	  DEF_gss; DEF_gbss;
	  f[i] = siginv[k] * ((1-sig[k])*f[i] + (gss-gbss));
	}
      else if (sigsize_dsig > 1 && sigsize_dsigg <= 1)
	LOOP_OVER_VOL_OWNED(v, fc, i) {
	  DEF_k;
	  f[i] = siginv[k] * ((1-sig[k])*f[i] + (g[i]-gb[i]));
	}
      else if (sigsize_dsig <= 1 && sigsize_dsigg > 1)
	LOOP_OVER_VOL_OWNED(v, fc, i) {
	  DEF_kx(g);
	  f[i] = f[i] + ((1+sigg[kg])*g[i]-(1-sigg[kg])*gb[i]);
	}
      else abort("bug - non-PML case in PML-only code");
    }
  }
}

} // namespace meep
