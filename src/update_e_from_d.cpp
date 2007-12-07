/* Copyright (C) 2005-2007 Massachusetts Institute of Technology
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

#include "meep.hpp"
#include "meep_internals.hpp"

namespace meep {

/* Given Dsqr = |D|^2 and Di = component of D, compute the factor f so
   that Ei = inveps * f * Di.   In principle, this would involve solving
   a cubic equation, but instead we use a Pade approximant that is 
   accurate to several orders.  This is inaccurate if the nonlinear
   index change is large, of course, but in that case the chi2/chi3
   power-series expansion isn't accurate anyway, so the cubic isn't
   physical there either. */
inline double calc_nonlinear_inveps(const double Dsqr, 
				    const double Di,
				    const double inveps,
				    const double chi2, const double chi3) {
  double c2 = Di*chi2*(inveps*inveps);
  double c3 = Dsqr*chi3*(inveps*(inveps*inveps));
  return (1 + c2 + 2*c3)/(1 + 2*c2 + 3*c3);
}

void fields::update_e_from_d() {
  for (int i=0;i<num_chunks;i++)
    if (chunks[i]->is_mine()) {
      src_vol *save_e_sources = chunks[i]->e_sources;
      if (disable_sources) chunks[i]->e_sources = NULL; // temporary
      chunks[i]->update_e_from_d();
      chunks[i]->e_sources = save_e_sources;
    }
}

void fields_chunk::update_e_from_d() {
  bool have_int_sources = false;
  for (src_vol *sv = e_sources; sv; sv = sv->next)
    if (sv->t->is_integrated) {
      have_int_sources = true;
      break;
    }

  FOR_ELECTRIC_COMPONENTS(ec) DOCMP 
    if (!d_minus_p[ec][cmp] && f[ec][cmp] && (pol || have_int_sources)) {
      d_minus_p[ec][cmp] = new double[v.ntot()];
      have_d_minus_p = true;
    }

  const int ntot = s->v.ntot();

  //////////////////////////////////////////////////////////////////////////
  // First, initialize d_minus_p to D - P, if necessary

  if (have_d_minus_p) {
    if (pol) {
      FOR_E_AND_D(ec,dc) if (f[ec][0]) {
	for (polarization *np=pol,*op=olpol; np; np=np->next,op=op->next) {
	  if (is_real) for (int i = 0; i < ntot; ++i) {
	    np->energy[ec][i] = op->energy[ec][i] +
	      (0.5)*(np->P[ec][0][i] - op->P[ec][0][i])
              * f[ec][0][i];
	  }
	  else for (int i = 0; i < ntot; ++i) {
            np->energy[ec][i] = op->energy[ec][i] +
              (0.5)*(np->P[ec][0][i] - op->P[ec][0][i])
              * f[ec][0][i] +
              (0.5)*(np->P[ec][1][i] - op->P[ec][1][i])
              * f[ec][1][i];
          }
	}
	DOCMP {
	  for (int i=0;i<ntot;i++) {
	    double sum = f[dc][cmp][i];
            for (polarization *p = pol; p; p = p->next) {
              sum -= p->P[ec][cmp][i];
            }	  
            d_minus_p[ec][cmp][i] = sum;
	  }
	}
      }
    }
    else {
      FOR_E_AND_D(ec,dc) if (f[ec][0]) DOCMP
	memcpy(d_minus_p[ec][cmp], f[dc][cmp], ntot * sizeof(double));
    }
  }

  //////////////////////////////////////////////////////////////////////////
  // Next, subtract time-integrated sources (i.e. polarizations, not currents)

  if (have_d_minus_p) {
    for (src_vol *sv = e_sources; sv; sv = sv->next) {  
      if (sv->t->is_integrated && f[sv->c][0]) {
	for (int j = 0; j < sv->npts; ++j) { 
	  const complex<double> A = sv->dipole(j);
	  DOCMP {
	    d_minus_p[sv->c][cmp][sv->index[j]] -= 
	      (cmp) ? imag(A) :  real(A);
	  }
	}
      }
    }
  }

  //////////////////////////////////////////////////////////////////////////
  // Finally, compute E = inveps * D
  
  double *dmp[NUM_FIELD_COMPONENTS][2];
  if (have_d_minus_p) {
    FOR_ELECTRIC_COMPONENTS(ec) DOCMP2 dmp[ec][cmp] = d_minus_p[ec][cmp];
  } else {
    FOR_E_AND_D(ec,dc) DOCMP2 dmp[ec][cmp] = f[dc][cmp];
  }

  FOR_E_AND_D(ec,dc) if (f[ec][0]) {
    const int d_ec = component_direction(ec);
    const int s_ec = stride_any_direction[d_ec];
    const direction d_1 = direction(v.dim != Dcyl 
				    ? (d_ec+1)%3 : ((d_ec-2)+1)%3+2);
    const component ec_1 = direction_component(ec,d_1);
    const int s_1 = stride_any_direction[d_1];
    const direction d_2 = direction(v.dim != Dcyl ? (d_ec+2)%3 : d_ec%3+2);
    const component ec_2 = direction_component(ec,d_2);
    const int s_2 = stride_any_direction[d_2];

    const double *ieps = s->inveps[ec][d_ec];
    const double *ieps1 = dmp[ec_1][0] ? s->inveps[ec][d_1] : 0;
    const double *ieps2 = dmp[ec_2][0] ? s->inveps[ec][d_2] : 0;

    if (ieps1 && ieps2) { // have 3x3 off-diagonal inveps
      if (s->chi3[ec]) { // nonlinear
	const double *chi2 = s->chi2[ec];
	const double *chi3 = s->chi3[ec];
	DOCMP {
	  double *efield = f[ec][cmp];
	  const double *dfield = dmp[ec][cmp];
	  const double *df1 = dmp[ec_1][cmp];
	  const double *df2 = dmp[ec_2][cmp];
	  LOOP_OVER_VOL_OWNED(v, ec, i) {
	    double df1s = df1[i]+df1[i+s_ec]+df1[i-s_1]+df1[i+(s_ec-s_1)];
	    double df2s = df2[i]+df2[i+s_ec]+df2[i-s_2]+df2[i+(s_ec-s_2)];
	    double df = dfield[i]; double iep = ieps[i];
	    efield[i] = (df * iep + 0.25 * (ieps1[i]*df1s + ieps2[i]*df2s)) *
	      calc_nonlinear_inveps(df * df + 0.0625 * (df1s*df1s + df2s*df2s),
				    df, iep, chi2[i], chi3[i]);
	  }
	}
      }
      else { // linear, 3x3 off-diagonal inveps
	DOCMP {
	  double *efield = f[ec][cmp];
	  const double *dfield = dmp[ec][cmp];
	  const double *df1 = dmp[ec_1][cmp];
	  const double *df2 = dmp[ec_2][cmp];
	  LOOP_OVER_VOL_OWNED(v, ec, i)
	    if (ieps1[i] * ieps2[i] != 0)
	      efield[i] = dfield[i] * ieps[i] +
		0.25 * (ieps1[i] * (df1[i] + df1[i+s_ec]
				    + df1[i-s_1] + df1[i+(s_ec-s_1)]) +
			ieps2[i] * (df2[i] + df2[i+s_ec]
				    + df2[i-s_2] + df2[i+(s_ec-s_2)]));
	    else
	      efield[i] = dfield[i] * ieps[i];
	}
      }
    }
    else if (ieps1 || ieps2) { // 2x2 off-diagonal inveps
      int s_o = ieps1 ? s_1 : s_2;
      const double *iepso = ieps1 ? ieps1 : ieps2;
      if (s->chi3[ec]) { // nonlinear
	const double *chi2 = s->chi2[ec];
	const double *chi3 = s->chi3[ec];
	DOCMP {
	  double *efield = f[ec][cmp];
	  const double *dfield = dmp[ec][cmp];
	  const double *dfo = dmp[ieps1 ? ec_1 : ec_2][cmp];
	  LOOP_OVER_VOL_OWNED(v, ec, i) {
	    double dfos = dfo[i]+dfo[i+s_ec]+dfo[i-s_o]+dfo[i+(s_ec-s_o)];
	    double df = dfield[i]; double iep = ieps[i];
	    efield[i] = (df * iep + 0.25 * (iepso[i]*dfos)) *
	      calc_nonlinear_inveps(df * df + 0.0625*dfos*dfos, 
				    df, iep, chi2[i], chi3[i]);
	  }
	}
      }
      else { // linear, 2x2 off-diagonal inveps
	DOCMP {
	  double *efield = f[ec][cmp];
	  const double *dfield = dmp[ec][cmp];
	  const double *dfo = dmp[ieps1 ? ec_1 : ec_2][cmp];
	  LOOP_OVER_VOL_OWNED(v, ec, i)
	    if (iepso[i] != 0)
	      efield[i] = dfield[i] * ieps[i] +
		0.25 * (iepso[i] * (dfo[i] + dfo[i+s_ec]
				    + dfo[i-s_o] + dfo[i+(s_ec-s_o)]));
	    else
	      efield[i] = dfield[i] * ieps[i];

	}
      }
    }
    else { // inveps is diagonal
      if (s->chi3[ec]) { // nonlinear
	const double *chi2 = s->chi2[ec];
	const double *chi3 = s->chi3[ec];
	if (dmp[ec_1][0] && dmp[ec_2][0]) {
	  DOCMP {
	    double *efield = f[ec][cmp];
	    const double *dfield = dmp[ec][cmp];
	    const double *df1 = dmp[ec_1][cmp];
	    const double *df2 = dmp[ec_2][cmp];
	    LOOP_OVER_VOL_OWNED(v, ec, i) {
	      double df1s = df1[i]+df1[i+s_ec]+df1[i-s_1]+df1[i+(s_ec-s_1)];
	      double df2s = df2[i]+df2[i+s_ec]+df2[i-s_2]+df2[i+(s_ec-s_2)];
	      double df = dfield[i]; double iep = ieps[i];
	      efield[i] = df * iep *
		calc_nonlinear_inveps(df * df +
				      0.0625 * (df1s*df1s + df2s*df2s),
				      df, iep, chi2[i], chi3[i]);
	    }
	  }
	}
	else if (dmp[ec_1][0] || dmp[ec_2][0]) {
	  int s_o = dmp[ec_1][0] ? s_1 : s_2;
	  DOCMP {
	    double *efield = f[ec][cmp];
	    const double *dfield = dmp[ec][cmp];
	    const double *dfo = dmp[dmp[ec_1][0] ? ec_1 : ec_2][cmp];
	    LOOP_OVER_VOL_OWNED(v, ec, i) {
	      double dfos = dfo[i]+dfo[i+s_ec]+dfo[i-s_o]+dfo[i+(s_ec-s_o)];
	      double df = dfield[i]; double iep = ieps[i];
	      efield[i] = df * iep *
		calc_nonlinear_inveps(df * df + 0.0625 * dfos*dfos,
				      df, iep, chi2[i], chi3[i]);
	    }
	  }
	}
	else {
	  DOCMP {
	    double *efield = f[ec][cmp];
	    const double *dfield = dmp[ec][cmp];
	    for (int i = 0; i < ntot; ++i) {
	      double df = dfield[i]; double iep = ieps[i];
	      efield[i] = df * iep * 
		calc_nonlinear_inveps(df * df, df, iep, chi2[i], chi3[i]);
	    }
	  }
	}
      }
      else { // linear, diagonal inveps
	DOCMP {
	  double *efield = f[ec][cmp];
	  const double *dfield = dmp[ec][cmp];
	  for (int i = 0; i < ntot; ++i)
	    efield[i] = dfield[i] * ieps[i];
	}
      }
    }
  }

  /* Do annoying special cases for r=0 in cylindrical coords.  Note
     that this only really matters for field output; the Ez and Ep
     components at r=0 don't usually affect the fields elsewhere
     because of the form of Maxwell's equations in cylindrical coords. */
  // (FIXME: handle Kerr case?).
  if (v.dim == Dcyl && v.origin_r() == 0.0)
    DOCMP FOR_E_AND_D(ec,dc) if (f[ec][cmp] && ec != Er) {
      const int yee_idx = v.yee_index(ec);
      const int d_ec = component_direction(ec);
      const int sR = stride_any_direction[R];
      const double *D = have_d_minus_p ? d_minus_p[ec][cmp] : f[dc][cmp];
      for (int iZ=0; iZ<num_any_direction[Z]; iZ++) {
	const int i = yee_idx + iZ - sR;
	f[ec][cmp][i] = s->inveps[ec][d_ec][i] * D[i];
      }
    }

}

} // namespace meep
