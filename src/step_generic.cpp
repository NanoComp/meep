#include "meep.hpp"

namespace meep {

/* update step for df/dt = curl g,
   i.e. f += dt curl g = dt/dx (dg1 - dg2)
   where dgk = g[i+sk] - g[i].

   g = (g1,g2), where g2 may be NULL.  Note that dt/dx and/or s1 and s2 may
   be negative to flip signs of derivatives. 

   PML: sig[k] = sigma[k]*dt/2, siginv[k] = 1 / (1 + sigma[k]*dt/2).
   Here, k is the dsig direction.  if dsig == NO_DIRECTION, then PML
   is not used.  (dsig is the sigma direction.)

   if non-NULL, then cnd is an array of conductivity values, changing
   the underlying PDE to:
       df/dt = curl g - cnd f
   which is updated as:
       f = [ dt * curl g + (1 - dt cnd/2) f ] / (1 + dt cnd/2)
   cndinv should be an array of 1 / (1 + dt cnd/2).  In the case
   of PML, cndinv should contain 1 / (1 + dt (cnd + sigma)/2).
*/

void step_curl(double *f, component c, const double *g1, const double *g2,
	       int s1, int s2, // strides for g1/g2 shift
	       const volume &v, double dtdx,
	       direction dsig, const double *sig, const double *siginv,
	       double dt, const double *cnd, const double *cndinv)
{
  if (dsig == NO_DIRECTION) { // no PML
    if (cnd) {
      double dt2 = dt * 0.5;
      if (g2) {
	LOOP_OVER_VOL_OWNED(v, c, i)
	  f[i] = ((1 - dt2 * cnd[i]) * f[i] - 
		  dtdx * (g1[i+s1] - g1[i] + g2[i] - g2[i+s2])) * cndinv[i];
      }
      else {
	LOOP_OVER_VOL_OWNED(v, c, i)
	  f[i] = ((1 - dt2 * cnd[i]) * f[i] 
		  - dtdx * (g1[i+s1] - g1[i])) * cndinv[i];
      }
    }
    else { // no conductivity
      if (g2) {
	LOOP_OVER_VOL_OWNED(v, c, i)
	  f[i] -= dtdx * (g1[i+s1] - g1[i] + g2[i] - g2[i+s2]);
      }
      else {
	LOOP_OVER_VOL_OWNED(v, c, i)
	  f[i] -= dtdx * (g1[i+s1] - g1[i]);
      }
    }
  }
  else { /* PML */
    int k0 = v.little_corner().in_direction(dsig);
    if (cnd) {
      double dt2 = dt * 0.5;
      if (g2) {
	LOOP_OVER_VOL_OWNED(v, c, i) {
	  IVEC_LOOP_ILOC(v, iloc);
	  int k = iloc.in_direction(dsig) - k0;
	  f[i] = ((1 - dt2 * cnd[i] - sig[k]) * f[i] - 
		  dtdx * (g1[i+s1] - g1[i] + g2[i] - g2[i+s2])) * cndinv[i];
	}
      }
      else {
	LOOP_OVER_VOL_OWNED(v, c, i) {
	  IVEC_LOOP_ILOC(v, iloc);
	  int k = iloc.in_direction(dsig) - k0;
	  f[i] = ((1 - dt2 * cnd[i] - sig[k]) * f[i] 
		  - dtdx * (g1[i+s1] - g1[i])) * cndinv[i];
	}
      }
    }
    else { // no conductivity (other than PML conductivity)
      if (g2) {
	LOOP_OVER_VOL_OWNED(v, c, i) {
	  IVEC_LOOP_ILOC(v, iloc);
	  int k = iloc.in_direction(dsig) - k0;
	  f[i] = ((1 - sig[k]) * f[i] -
		  dtdx * (g1[i+s1] - g1[i] + g2[i] - g2[i+s2])) * siginv[k];
	}
      }
      else {
	LOOP_OVER_VOL_OWNED(v, c, i) {
	  IVEC_LOOP_ILOC(v, iloc);
	  int k = iloc.in_direction(dsig) - k0;
	  f[i] = ((1 - sig[k]) * f[i] - dtdx * (g1[i+s1] - g1[i])) * siginv[k];
	}
      }
    }
  }
}

/* Given Dsqr = |D|^2 and Di = component of D, compute the factor f so
   that Ei = inveps * f * Di.   In principle, this would involve solving
   a cubic equation, but instead we use a Pade approximant that is 
   accurate to several orders.  This is inaccurate if the nonlinear
   index change is large, of course, but in that case the chi2/chi3
   power-series expansion isn't accurate anyway, so the cubic isn't
   physical there either. */
inline double calc_nonlinear_u(const double Dsqr, 
			       const double Di,
			       const double inveps,
			       const double chi2, const double chi3) {
  double c2 = Di*chi2*(inveps*inveps);
  double c3 = Dsqr*chi3*(inveps*inveps*inveps);
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
  void step_update_EDHB(double *f, component fc, const volume &v, 
		     const double *g, const double *g1, const double *g2,
		     const double *gb, const double *g1b, const double *g2b,
		     const double *u, const double *u1, const double *u2,
		     int s, int s1, int s2,
		     const double *chi2, const double *chi3,
		     direction dsig, const double *sig, const double *siginv,
		     direction dsigg, const double *sigg,
		     direction dsig1, const double *sig1,
		     direction dsig1inv, const double *sig1inv,
		     direction dsig2, const double *sig2,
		     direction dsig2inv, const double *sig2inv,
		     int sigsize_dsig,int sigsize_dsigg,int sigsize_dsig1)
{
  if (!f) return;
  int sigsize_dsig1inv = sigsize_dsigg;
  int sigsize_dsig2 = sigsize_dsig;
  int sigsize_dsig2inv = sigsize_dsig1;
  
# define SWAP(t,a,b) { t xxxx = a; a = b; b = xxxx; }
  if ((!g1 && g2) || (g1 && g2 && !u1 && u2)) { /* swap g1 and g2 */
    SWAP(const double*, g1, g2);
    SWAP(const double*, g1b, g2b);
    SWAP(const double*, u1, u2);
    SWAP(int, s1, s2);
    SWAP(const double*, sig1, sig2);
    SWAP(const double*, sig1inv, sig2inv);
    SWAP(direction, dsig1, dsig2);
    SWAP(direction, dsig1inv, dsig2inv);
    SWAP(int, sigsize_dsig1, sigsize_dsig2);
    SWAP(int, sigsize_dsig1inv, sigsize_dsig2inv);
  }
  
  if (sigsize_dsig <= 1 && sigsize_dsigg <= 1 && sigsize_dsig1 <= 1) { // no PML
    if (u1 && u2) { // 3x3 off-diagonal u
      if (chi3) {
	LOOP_OVER_VOL_OWNED(v, fc, i) {
	  double g1s = g1[i]+g1[i+s]+g1[i-s1]+g1[i+(s-s1)];
	  double g2s = g2[i]+g2[i+s]+g2[i-s2]+g2[i+(s-s2)];
	  double gs = g[i]; double us = u[i];
	  f[i] = (gs * us + 0.25 * (u1[i]*g1s + u2[i]*g2s)) *
	    calc_nonlinear_u(gs * gs + 0.0625 * (g1s*g1s + g2s*g2s),
			     gs, us, chi2[i], chi3[i]);
	}
      }
      else {
	LOOP_OVER_VOL_OWNED(v, fc, i) {
	  double g1s = g1[i]+g1[i+s]+g1[i-s1]+g1[i+(s-s1)];
	  double g2s = g2[i]+g2[i+s]+g2[i-s2]+g2[i+(s-s2)];
	  double gs = g[i]; double us = u[i];
	  f[i] = (gs * us + 0.25 * (u1[i]*g1s + u2[i]*g2s));
	}
      }
    }
    else if (u1) { // 2x2 off-diagonal u
      if (chi3) {
	LOOP_OVER_VOL_OWNED(v, fc, i) {
	  double g1s = g1[i]+g1[i+s]+g1[i-s1]+g1[i+(s-s1)];
	  double gs = g[i]; double us = u[i];
	  f[i] = (gs * us + 0.25 * (u1[i]*g1s)) *
	    calc_nonlinear_u(gs * gs + 0.0625 * (g1s*g1s),
			     gs, us, chi2[i], chi3[i]);
	}
      }
      else {
	LOOP_OVER_VOL_OWNED(v, fc, i) {
	  double g1s = g1[i]+g1[i+s]+g1[i-s1]+g1[i+(s-s1)];
	  double gs = g[i]; double us = u[i];
	  f[i] = (gs * us + 0.25 * (u1[i]*g1s));
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
      else {
	LOOP_OVER_VOL_OWNED(v, fc, i) {
	  double gs = g[i]; double us = u[i];
	  f[i] = (gs * us);
	}
      }
    }
  }
  else { // PML
    
    // offset so that we can index by iloc
    /* The sigsize_dsig>1 check seems to be redudant, since there
       is an identical check when k0 (etc.) is used.  However,
       for some reason we don't understand, removing that check
       here causes valgrind to complain about some uninitialized
       pointer value when running bragg_transmission.dac ... since
       it is harmless, we keep it here to shut valgrind up
       (and on the chance that valgrind is identifying a real bug
       that we can't figure out). */
   int k0 = (sigsize_dsig > 1)?
     v.little_corner().in_direction(dsig):0;
   int kg0 = (sigsize_dsigg > 1)?
     v.little_corner().in_direction(dsigg):0;
   int k10 = (sigsize_dsig1 > 1)?
     v.little_corner().in_direction(dsig1):0;
   int k1inv0 = (sigsize_dsigg > 1)?
     v.little_corner().in_direction(dsig1inv):0;
   int k20 = (sigsize_dsig > 1)?
     v.little_corner().in_direction(dsig2):0;
   int k2inv0 = (sigsize_dsig1 > 1)?
     v.little_corner().in_direction(dsig2inv):0;

   // the following definitions are used over and over

   // indices into sigma arrays:
#  define KDEF(k,dsig,k0) int k = (sigsize_ ## dsig > 1) * \
                                  (iloc.in_direction(dsig) - k0)
#  define DEF_k KDEF(k,dsig,k0)
#  define DEF_kx(x) KDEF(k ## x, dsig ## x, k ## x ## 0)

   // fields
#  define DEF_gs  double gs = ((1+sigg[kg])*siginv[k]) * g[i]
#  define DEF_gbs double gbs = ((1-sigg[kg])*siginv[k]) * gb[i]
#  define DEF_us  double us = u[i]
#  define DEF_g1s double g1s = ((1+sig1[k1]) * sig1inv[k1inv]) * \
                               (g1[i]+g1[i+s]+g1[i-s1]+g1[i+(s-s1)])
#  define DEF_g2s double g2s = ((1+sig2[k2]) * sig2inv[k2inv]) * \
                               (g2[i]+g2[i+s]+g2[i-s2]+g2[i+(s-s2)])
#  define DEF_g1bs double g1bs = ((1-sig1[k1]) * sig1inv[k1inv]) * \
                                 (g1b[i]+g1b[i+s]+g1b[i-s1]+g1b[i+(s-s1)])
#  define DEF_g2bs double g2bs = ((1-sig2[k2]) * sig2inv[k2inv]) * \
                                 (g2b[i]+g2b[i+s]+g2b[i-s2]+g2b[i+(s-s2)])


   if (u1 && u2) { // 3x3 off-diagonal u
     if (chi3) {
       LOOP_OVER_VOL_OWNED(v, fc, i) {
	  IVEC_LOOP_ILOC(v, iloc);
	  DEF_kx(1); DEF_kx(1inv); DEF_g1s; DEF_g1bs;
	  DEF_kx(2); DEF_kx(2inv); DEF_g2s; DEF_g2bs;
	  DEF_k; DEF_kx(g); DEF_gs; DEF_gbs;
	  DEF_us;
	  f[i] = ((1-sig[k])*siginv[k]) * f[i] + 
	    ((gs-gbs) * us + 0.25 * (u1[i]*(g1s-g1bs) + u2[i]*(g2s-g2bs))) *
	    calc_nonlinear_u(gs * gs + 0.0625 * (g1s*g1s + g2s*g2s),
			     gs, us, chi2[i], chi3[i]);
	}
      } else {
	LOOP_OVER_VOL_OWNED(v, fc, i) {
	  IVEC_LOOP_ILOC(v, iloc);
	  DEF_kx(1); DEF_kx(1inv); DEF_g1s; DEF_g1bs;
          DEF_kx(2); DEF_kx(2inv); DEF_g2s; DEF_g2bs;
          DEF_k; DEF_kx(g); DEF_gs; DEF_gbs;
          DEF_us;
	  f[i] = ((1-sig[k])*siginv[k]) * f[i] + 
	    ((gs-gbs) * us + 0.25 * (u1[i]*(g1s-g1bs) + u2[i]*(g2s-g2bs)));
	}
      }
    } else if (u1) { // 2x2 off-diagonal u
      if (chi3) {
	LOOP_OVER_VOL_OWNED(v, fc, i) {
	  IVEC_LOOP_ILOC(v, iloc);
	  DEF_kx(1); DEF_kx(1inv); DEF_g1s; DEF_g1bs;
	  DEF_k; DEF_kx(g); DEF_gs; DEF_gbs;
          DEF_us;
	  f[i] = ((1-sig[k])*siginv[k]) * f[i] + 
	    ((gs - gbs) * us + 0.25 * (u1[i]*(g1s-g1bs))) *
	    calc_nonlinear_u(gs * gs + 0.0625 * (g1s*g1s),
			     gs, us, chi2[i], chi3[i]);
	}
      } else {
	LOOP_OVER_VOL_OWNED(v, fc, i) {
	  IVEC_LOOP_ILOC(v, iloc);
	  DEF_kx(1); DEF_kx(1inv); DEF_g1s; DEF_g1bs;
          DEF_k; DEF_kx(g); DEF_gs; DEF_gbs;
          DEF_us;
	  f[i] = ((1-sig[k])*siginv[k]) * f[i] + 
	    ((gs - gbs) * us + 0.25 * (u1[i]*(g1s-g1bs)));
	}	
      }
    } else if (u2) { // 2x2 off-diagonal u
      abort("bug - didn't swap off-diagonal terms!?");
    } else { // diagonal u
      if (chi3) {
	if (g1 && g2) {
	  LOOP_OVER_VOL_OWNED(v, fc, i) {
	    IVEC_LOOP_ILOC(v, iloc);
	    DEF_kx(1); DEF_kx(1inv); DEF_g1s;
	    DEF_kx(2); DEF_kx(2inv); DEF_g2s;
	    DEF_k; DEF_kx(g); DEF_gs; DEF_gbs;
	    DEF_us;
	    f[i] = ((1-sig[k])*siginv[k])*f[i] + (gs-gbs)*us *
	      calc_nonlinear_u(gs * gs + 0.0625 * (g1s*g1s + g2s*g2s),
			       gs, us, chi2[i], chi3[i]);
	  }
	} else if (g1) {
	  LOOP_OVER_VOL_OWNED(v, fc, i) {
	    IVEC_LOOP_ILOC(v, iloc);
	    DEF_kx(1); DEF_kx(1inv); DEF_g1s;
            DEF_k; DEF_kx(g); DEF_gs; DEF_gbs;
            DEF_us;
	    f[i] = ((1-sig[k])*siginv[k])*f[i] + (gs-gbs)*us *
     	      calc_nonlinear_u(gs * gs + 0.0625 * (g1s*g1s),
     			       gs, us, chi2[i], chi3[i]);
	  }
	} else if (g2) {
	  abort("bug - didn't swap off-diagonal terms!?");
	} else {
	  LOOP_OVER_VOL_OWNED(v, fc, i) {
	    IVEC_LOOP_ILOC(v, iloc);
	    DEF_k; DEF_kx(g); DEF_gs; DEF_gbs;
            DEF_us;
	    f[i] = ((1-sig[k])*siginv[k])*f[i] + (gs-gbs)*us *
     	      calc_nonlinear_u(gs * gs, gs, us, chi2[i], chi3[i]);
	  }
	}
      } else { //linear, diagonal u
	LOOP_OVER_VOL_OWNED(v, fc, i) {
	  IVEC_LOOP_ILOC(v, iloc);
	  DEF_k; DEF_kx(g); DEF_gs; DEF_gbs;
	  DEF_us;
	  f[i] = ((1-sig[k])*siginv[k])*f[i] + (gs-gbs) * us;
	}
      }
    }
  }
}

} // namespace meep
