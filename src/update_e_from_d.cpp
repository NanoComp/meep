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
	  if (np->energy[ec] && op->energy[ec]) {
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

  DOCMP FOR_E_AND_D(ec,dc) if (f[ec][cmp]) {
    const int d_ec = component_direction(ec);
    const int s_ec = stride_any_direction[d_ec];
    const direction d_1 = direction((d_ec+1)%3);
    const component ec_1 = direction_component(ec,d_1);
    const int s_1 = stride_any_direction[d_1];
    const direction d_2 = direction((d_ec+2)%3);
    const component ec_2 = direction_component(ec,d_2);
    const int s_2 = stride_any_direction[d_2];

    component dc_1 = direction_component(dc,d_1);
    component dc_2 = direction_component(dc,d_2);

    direction dsig = (direction)((d_ec+2)%3);
    direction dsigg = (direction)(d_ec);
    direction dsig1 = (direction)((d_ec+1)%3);
    direction dsig1inv = (direction)(d_ec);
    direction dsig2 = (direction)((d_ec+2)%3);
    direction dsig2inv = (direction)((d_ec+1)%3);

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
