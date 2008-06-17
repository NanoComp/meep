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
*/

void step_curl(double *f, component c, const double *g1, const double *g2,
	       int s1, int s2, // strides for g1/g2 shift
	       const volume &v, double dtdx,
	       direction dsig, const double *sig, const double *siginv)
{
  if (dsig == NO_DIRECTION) { // no PML
    if (g2) {
      LOOP_OVER_VOL_OWNED(v, c, i)
	f[i] -= dtdx * (g1[i+s1] - g1[i] + g2[i] - g2[i+s2]);
    }
    else {
      LOOP_OVER_VOL_OWNED(v, c, i)
	f[i] -= dtdx * (g1[i+s1] - g1[i]);
    }
  }
  else { /* PML */
    int k0 = v.little_corner().in_direction(dsig);
    if (g2) {
      LOOP_OVER_VOL_OWNED(v, c, i) {
	IVEC_LOOP_ILOC(v, iloc);
	int k = iloc.in_direction(dsig) - k0;
	f[i] = ((1 - sig[k]) * siginv[k]) * f[i] -
	  (siginv[k] * dtdx) * (g1[i+s1] - g1[i] + g2[i] - g2[i+s2]);
      }
    }
    else {
      LOOP_OVER_VOL_OWNED(v, c, i) {
	IVEC_LOOP_ILOC(v, iloc);
	int k = iloc.in_direction(dsig) - k0;
	f[i] = ((1 - sig[k]) * siginv[k]) * f[i] -
	  (siginv[k] * dtdx) * (g1[i+s1] - g1[i]);
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
    else if (u1 && !u2) { // 2x2 off-diagonal u
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
    else if (!u1 && u2) { // 2x2 off-diagonal u
      if (chi3) {
	LOOP_OVER_VOL_OWNED(v, fc, i) {
	  double g2s = g2[i]+g2[i+s]+g2[i-s2]+g2[i+(s-s2)];
	  double gs = g[i]; double us = u[i];
	  f[i] = (gs * us + 0.25 * (u2[i]*g2s)) *
	    calc_nonlinear_u(gs * gs + 0.0625 * (g2s*g2s),
			     gs, us, chi2[i], chi3[i]);
	}
      }
      else {
	LOOP_OVER_VOL_OWNED(v, fc, i) {
	  double g2s = g2[i]+g2[i+s]+g2[i-s2]+g2[i+(s-s2)];
	  double gs = g[i]; double us = u[i];
	  f[i] = (gs * us + 0.25 * (u2[i]*g2s));
	}
      }
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
	  LOOP_OVER_VOL_OWNED(v, fc, i) {
	    double g2s = g2[i]+g2[i+s]+g2[i-s2]+g2[i+(s-s2)];
	    double gs = g[i]; double us = u[i];
	    f[i] = (gs*us)*calc_nonlinear_u(gs*gs + 0.0625*(g2s*g2s),
						gs, us, chi2[i], chi3[i]);
	  }
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

   if (u1 && u2) { // 3x3 off-diagonal u
     if (chi3) {
       LOOP_OVER_VOL_OWNED(v, fc, i) {
	  IVEC_LOOP_ILOC(v, iloc);
	  int k = iloc.in_direction(dsig) - k0;
	  int k1 = iloc.in_direction(dsig1) - k10;
	  int k1inv = iloc.in_direction(dsig1inv) - k1inv0;
	  double g1s = ((1+sig1[k1]) * sig1inv[k1inv]) *
	    (g1[i]+g1[i+s]+g1[i-s1]+g1[i+(s-s1)]);
	  double g1bs = ((1-sig1[k1]) * sig1inv[k1]) *
	    (g1b[i]+g1b[i+s]+g1b[i-s1]+g1b[i+(s-s1)]);
	  int k2 = iloc.in_direction(dsig2);
	  int k2inv = iloc.in_direction(dsig2inv);
	  double g2s = ((1+sig2[k2]) * sig2inv[k2inv]) *
	    (g2[i]+g2[i+s]+g2[i-s2]+g2[i+(s-s2)]);
	  double g2bs = ((1-sig2[k2]) * sig2inv[k2]) *
	    (g2b[i]+g2b[i+s]+g2b[i-s2]+g2b[i+(s-s2)]);
	  int kg = iloc.in_direction(dsigg);
	  double gs = ((1+sigg[kg])*siginv[k]) * g[i];
	  double gbs = ((1-sigg[kg])*siginv[k]) * gb[i];
	  double us = u[i];
	  f[i] = ((1-sig[k])*siginv[k]) * f[i] + 
	    ((gs-gbs) * us + 0.25 * (u1[i]*(g1s-g1bs) + u2[i]*(g2s-g2bs))) *
	    calc_nonlinear_u(gs * gs + 0.0625 * (g1s*g1s + g2s*g2s),
			     gs, us, chi2[i], chi3[i]);
	}
      } else {
	LOOP_OVER_VOL_OWNED(v, fc, i) {
	  IVEC_LOOP_ILOC(v, iloc);
	  int k = iloc.in_direction(dsig) - k0;
	  int k1 = iloc.in_direction(dsig1) - k10;
	  int k1inv = iloc.in_direction(dsig1inv) - k1inv0;
	  double g1s = ((1+sig1[k1]) * sig1inv[k1inv]) *
	    (g1[i]+g1[i+s]+g1[i-s1]+g1[i+(s-s1)]);
	  double g1bs = ((1-sig1[k1]) * sig1inv[k1]) *
	    (g1b[i]+g1b[i+s]+g1b[i-s1]+g1b[i+(s-s1)]);
	  int k2 = iloc.in_direction(dsig2) - k20;
	  int k2inv = iloc.in_direction(dsig2inv) - k2inv0;
	  double g2s = ((1+sig2[k2]) * sig2inv[k2inv]) *
	    (g2[i]+g2[i+s]+g2[i-s2]+g2[i+(s-s2)]);
	  double g2bs = ((1-sig2[k2]) * sig2inv[k2]) *
	    (g2b[i]+g2b[i+s]+g2b[i-s2]+g2b[i+(s-s2)]);
	  int kg = iloc.in_direction(dsigg) - kg0;
	  double gs = ((1+sigg[kg])*siginv[k]) * g[i];
	  double gbs = ((1-sigg[kg])*siginv[k]) * gb[i];
	  double us = u[i];
	  f[i] = ((1-sig[k])*siginv[k]) * f[i] + 
	    ((gs-gbs) * us + 0.25 * (u1[i]*(g1s-g1bs) + u2[i]*(g2s-g2bs)));
	}
      }
    } else if (u1 && !u2) { // 2x2 off-diagonal u
      if (chi3) {
	LOOP_OVER_VOL_OWNED(v, fc, i) {
	  IVEC_LOOP_ILOC(v, iloc);
	  int k = iloc.in_direction(dsig) - k0;
	  int k1 = iloc.in_direction(dsig1) - k10;
	  int k1inv = iloc.in_direction(dsig1inv) - k1inv0;
	  double g1s = ((1+sig1[k1]) * sig1inv[k1inv]) *
	    (g1[i]+g1[i+s]+g1[i-s1]+g1[i+(s-s1)]);
	  double g1bs = ((1-sig1[k1]) * sig1inv[k1]) *
	    (g1b[i]+g1b[i+s]+g1b[i-s1]+g1b[i+(s-s1)]);
	  int kg = iloc.in_direction(dsigg) - kg0;
	  double gs = ((1+sigg[kg])*siginv[k]) * g[i];
	  double gbs = ((1-sigg[kg])*siginv[k]) * gb[i];
	  double us = u[i];
	  f[i] = ((1-sig[k])*siginv[k]) * f[i] + 
	    ((gs - gbs) * us + 0.25 * (u1[i]*(g1s-g1bs))) *
	    calc_nonlinear_u(gs * gs + 0.0625 * (g1s*g1s),
			     gs, us, chi2[i], chi3[i]);
	}
      } else {
	LOOP_OVER_VOL_OWNED(v, fc, i) {
	  IVEC_LOOP_ILOC(v, iloc);
	  int k = iloc.in_direction(dsig) - k0;
	  int k1 = iloc.in_direction(dsig1) - k10;
	  int k1inv = iloc.in_direction(dsig1inv) - k1inv0;
	  double g1s = ((1+sig1[k1]) * sig1inv[k1inv]) *
	    (g1[i]+g1[i+s]+g1[i-s1]+g1[i+(s-s1)]);
	  double g1bs = ((1-sig1[k1]) * sig1inv[k1]) *
	    (g1b[i]+g1b[i+s]+g1b[i-s1]+g1b[i+(s-s1)]);
	  int kg = iloc.in_direction(dsigg) - kg0;
	  double gs = ((1+sigg[kg])*siginv[k]) * g[i];
	  double gbs = ((1-sigg[kg])*siginv[k]) * gb[i];
	  double us = u[i];
	  f[i] = ((1-sig[k])*siginv[k]) * f[i] + 
	    ((gs - gbs) * us + 0.25 * (u1[i]*(g1s-g1bs)));
	}	
      }
    } else if (!u1 && u2) { // 2x2 off-diagonal u
      if (chi3) {
	LOOP_OVER_VOL_OWNED(v, fc, i) {
	  IVEC_LOOP_ILOC(v, iloc);
	  int k = iloc.in_direction(dsig) - k0;
	  int k2 = iloc.in_direction(dsig2) - k20;
	  int k2inv = iloc.in_direction(dsig2inv) - k2inv0;
	  double g2s = ((1+sig2[k2]) * sig2inv[k2inv]) *
	    (g2[i]+g2[i+s]+g2[i-s2]+g2[i+(s-s2)]);
	  double g2bs = ((1-sig2[k2]) * sig2inv[k2]) *
	    (g2b[i]+g2b[i+s]+g2b[i-s2]+g2b[i+(s-s2)]);
	  int kg = iloc.in_direction(dsigg) - kg0;
	  double gs = ((1+sigg[kg])*siginv[k]) * g[i];
	  double gbs = ((1-sigg[kg])*siginv[k]) * gb[i];
	  double us = u[i];
	  f[i] = ((1-sig[k])*siginv[k]) * f[i] + 
	    ((gs - gbs) * us + 0.25 * (u2[i]*(g2s-g2bs))) *
	    calc_nonlinear_u(gs * gs + 0.0625 * (g2s*g2s),
			     gs, us, chi2[i], chi3[i]);
	}
      } else {
	LOOP_OVER_VOL_OWNED(v, fc, i) {
	  IVEC_LOOP_ILOC(v, iloc);
	  int k = iloc.in_direction(dsig) - k0;
	  int k2 = iloc.in_direction(dsig1) - k20;
	  int k2inv = iloc.in_direction(dsig2inv) - k2inv0;
	  double g2s = ((1+sig2[k2]) * sig2inv[k2inv]) *
	    (g2[i]+g2[i+s]+g2[i-s2]+g2[i+(s-s2)]);
	  double g2bs = ((1-sig2[k2]) * sig2inv[k2]) *
	    (g2b[i]+g2b[i+s]+g2b[i-s2]+g2b[i+(s-s2)]);
	  int kg = iloc.in_direction(dsigg) - kg0;
	  double gs = ((1+sigg[kg])*siginv[k]) * g[i];
	  double gbs = ((1-sigg[kg])*siginv[k]) * gb[i];
	  double us = u[i];
	  f[i] = ((1-sig[k])*siginv[k]) * f[i] + 
	    ((gs - gbs) * us + 0.25 * (u2[i]*(g2s-g2bs)));
	}	
      }
    } else { // diagonal u
      if (chi3) {
	if (g1 && g2) {
	  LOOP_OVER_VOL_OWNED(v, fc, i) {
	    IVEC_LOOP_ILOC(v, iloc);
	    int k = iloc.in_direction(dsig) - k0;
	    int k1 = iloc.in_direction(dsig1) - k10;
	    int k1inv = iloc.in_direction(dsig1inv) - k1inv0;
	    double g1s = ((1+sig1[k1]) * sig1inv[k1inv]) *
	      (g1[i]+g1[i+s]+g1[i-s1]+g1[i+(s-s1)]);
	    int k2 = iloc.in_direction(dsig2) - k20;
	    int k2inv = iloc.in_direction(dsig2inv) - k2inv0;
	    double g2s = ((1+sig2[k2]) * sig2inv[k2inv]) *
	      (g2[i]+g2[i+s]+g2[i-s2]+g2[i+(s-s2)]);
	    int kg = iloc.in_direction(dsigg) - kg0;
	    double gs = ((1+sigg[kg])*siginv[k]) * g[i];
	    double gbs = ((1-sigg[kg])*siginv[k]) * gb[i];
	    double us = u[i];
	    f[i] = ((1-sig[k])*siginv[k])*f[i] + (gs-gbs)*us *
	      calc_nonlinear_u(gs * gs + 0.0625 * (g1s*g1s + g2s*g2s),
			       gs, us, chi2[i], chi3[i]);
	  }
	} else if (g1 && !g2) {
	  LOOP_OVER_VOL_OWNED(v, fc, i) {
	    IVEC_LOOP_ILOC(v, iloc);
	    int k = iloc.in_direction(dsig) - k0;
	    int k1 = iloc.in_direction(dsig1) - k10;
	    int k1inv = iloc.in_direction(dsig1inv) - k1inv0;
	    double g1s = ((1+sig1[k1]) * sig1inv[k1inv]) *
	      (g1[i]+g1[i+s]+g1[i-s1]+g1[i+(s-s1)]);
	    int kg = iloc.in_direction(dsigg) - kg0;
	    double gs = ((1+sigg[kg])*siginv[k]) * g[i];
	    double gbs = ((1-sigg[kg])*siginv[k]) * gb[i];
	    double us = u[i];
	    f[i] = ((1-sig[k])*siginv[k])*f[i] + (gs-gbs)*us *
     	      calc_nonlinear_u(gs * gs + 0.0625 * (g1s*g1s),
     			       gs, us, chi2[i], chi3[i]);
	  }
	} else if (!g1 && g2) {
	  LOOP_OVER_VOL_OWNED(v, fc, i) {
	    IVEC_LOOP_ILOC(v, iloc);
	    int k = iloc.in_direction(dsig) - k0;
	    int k2 = iloc.in_direction(dsig2) - k20;
	    int k2inv = iloc.in_direction(dsig2inv) - k2inv0;
	    double g2s = ((1+sig2[k2]) * sig2inv[k2inv]) *
	      (g2[i]+g2[i+s]+g2[i-s2]+g2[i+(s-s2)]);
	    int kg = iloc.in_direction(dsigg) - kg0;
	    double gs = ((1+sigg[kg])*siginv[k]) * g[i];
	    double gbs = ((1-sigg[kg])*siginv[k]) * gb[i];
	    double us = u[i];
	    f[i] = ((1-sig[k])*siginv[k])*f[i] + (gs-gbs)*us *
     	      calc_nonlinear_u(gs * gs + 0.0625 * (g2s*g2s),
     			       gs, us, chi2[i], chi3[i]);
	  }
	} else {
	  LOOP_OVER_VOL_OWNED(v, fc, i) {
	    IVEC_LOOP_ILOC(v, iloc);
	    int k = (sigsize_dsig>1)*(iloc.in_direction(dsig) - k0);
	    int kg = (sigsize_dsigg>1)*(iloc.in_direction(dsigg) - kg0);
	    double gs = ((1+sigg[kg])*siginv[k]) * g[i];
	    double gbs = ((1-sigg[kg])*siginv[k]) * gb[i];
	    double us = u[i];
	    f[i] = ((1-sig[k])*siginv[k])*f[i] + (gs-gbs)*us *
     	      calc_nonlinear_u(gs * gs, gs, us, chi2[i], chi3[i]);
	  }
	}
      } else { //linear, diagonal u
	LOOP_OVER_VOL_OWNED(v, fc, i) {
	  IVEC_LOOP_ILOC(v, iloc);
	  int k = (sigsize_dsig>1)*(iloc.in_direction(dsig) - k0);
	  int kg = (sigsize_dsigg>1)*(iloc.in_direction(dsigg) - kg0); 
	  double gs = ((1+sigg[kg])*siginv[k]) * g[i];
	  double gbs = ((1-sigg[kg])*siginv[k]) * gb[i];
	  double us = u[i];
	  f[i] = ((1-sig[k])*siginv[k])*f[i] + (gs-gbs) * us;
	}
      }
    }
  }
}

} // namespace meep
