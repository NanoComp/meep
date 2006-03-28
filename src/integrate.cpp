/* Copyright (C) 2006 Massachusetts Institute of Technology.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

#include "meep.hpp"
#include "meep_internals.hpp"

/* generic integration and related routines, based fields::loop_in_chunk */

namespace meep {

struct integrate_data {
  int num_fields;
  const component *components;
  component *cS;
  complex<double> *ph;
  complex<double> *fields;
  int *offsets;
  int ninveps;
  component inveps_cs[3];
  direction inveps_ds[3];
  int inveps_offsets[6];
  complex<long double> sum;
  double maxabs;
  field_function integrand;
  void *integrand_data_;
};

static void integrate_chunkloop(fields_chunk *fc, int ichunk, component cgrid,
				ivec is, ivec ie,
				vec s0, vec s1, vec e0, vec e1,
				double dV0, double dV1,
				ivec shift, complex<double> shift_phase,
				const symmetry &S, int sn,
				void *data_)
{
  integrate_data *data = (integrate_data *) data_;
  int *off = data->offsets;
  component *cS = data->cS;
  complex<double> *fields = data->fields, *ph = data->ph;
  complex<long double> sum = 0.0;
  double maxabs = 0;
  const component *iecs = data->inveps_cs;
  const direction *ieds = data->inveps_ds;
  int *ieos = data->inveps_offsets;

  for (int i = 0; i < data->num_fields; ++i) {
    cS[i] = S.transform(data->components[i], -sn);
    if (cS[i] == Dielectric)
      ph[i] = 1.0;
    else {
      if (cgrid == Dielectric)
	fc->v.yee2diel_offsets(cS[i], off[2*i], off[2*i+1]);
      ph[i] = shift_phase * S.phase_shift(cS[i], sn);
    }
  }
  for (int k = 0; k < data->ninveps; ++k)
    fc->v.yee2diel_offsets(iecs[k], ieos[2*k], ieos[2*k+1]);

  vec rshift(shift * (0.5*fc->v.inva));
  LOOP_OVER_IVECS(fc->v, is, ie, idx) {
    IVEC_LOOP_LOC(fc->v, loc);
    loc = S.transform(loc, sn) + rshift;

    for (int i = 0; i < data->num_fields; ++i) {
      if (cS[i] == Dielectric) {
	double tr = 0.0;
	for (int k = 0; k < data->ninveps; ++k)
	  tr += (fc->s->inveps[iecs[k]][ieds[k]][idx]
		 + fc->s->inveps[iecs[k]][ieds[k]][idx+ieos[2*k]]
		 + fc->s->inveps[iecs[k]][ieds[k]][idx+ieos[1+2*k]]
		 + fc->s->inveps[iecs[k]][ieds[k]][idx+ieos[2*k]+ieos[1+2*k]]);
	fields[i] = (4 * data->ninveps) / tr;
      }
      else {
	double f[2];
	for (int k = 0; k < 2; ++k)
	  if (fc->f[cS[i]][k])
	    f[k] = 0.25 * (fc->f[cS[i]][k][idx]
			   + fc->f[cS[i]][k][idx+off[2*i]]
			   + fc->f[cS[i]][k][idx+off[2*i+1]]
			   + fc->f[cS[i]][k][idx+off[2*i]+off[2*i+1]]);
	  else
	    f[k] = 0;
	fields[i] = complex<double>(f[0], f[1]) * ph[i];
      }
    }

    complex<double> integrand = 
      data->integrand(fields, loc, data->integrand_data_);
    maxabs = max(maxabs, abs(integrand));
    sum += integrand * IVEC_LOOP_WEIGHT(s0, s1, e0, e1, dV0 + dV1 * loop_i2);
  }

  data->maxabs = max(data->maxabs, maxabs);
  data->sum += sum;
}

complex<double> fields::integrate(int num_fields, const component *components,
				  field_function integrand,
				  const geometric_volume &where,
				  void *integrand_data_,
				  double *maxabs)
{
  // check if components are all on the same grid:
  bool same_grid = true;
  for (int i = 1; i < num_fields; ++i)
    if (v.iyee_shift(components[i]) != v.iyee_shift(components[0])) {
      same_grid = false;
      break;
    }

  component cgrid = Dielectric;
  if (same_grid && num_fields > 0)
    cgrid = components[0];

  integrate_data data;
  data.num_fields = num_fields;
  data.components = components;
  data.cS = new component[num_fields];
  data.ph = new complex<double>[num_fields];
  data.fields = new complex<double>[num_fields];
  data.sum = 0;
  data.maxabs = 0;
  data.integrand = integrand;
  data.integrand_data_ = integrand_data_;

  /* compute inverse-epsilon directions for computing Dielectric fields */
  data.ninveps = 0;
  bool needs_dielectric = false;
  for (int i = 0; i < num_fields; ++i)
    if (components[i] == Dielectric) { needs_dielectric = true; break; }
  if (needs_dielectric) 
    FOR_ELECTRIC_COMPONENTS(c) if (v.has_field(c)) {
      if (data.ninveps == 3) abort("more than 3 field components??");
      data.inveps_cs[data.ninveps] = c;
      data.inveps_ds[data.ninveps] = component_direction(c);
      ++data.ninveps;
    }
  
  data.offsets = new int[2 * num_fields];
  for (int i = 0; i < 2 * num_fields; ++i)
    data.offsets[i] = 0;

  loop_in_chunks(integrate_chunkloop, (void *) &data, where, cgrid);

  delete[] data.offsets;
  delete[] data.fields;
  delete[] data.ph;
  delete[] data.cS;

  if (maxabs)
    *maxabs = max_to_all(data.maxabs);
  data.sum = sum_to_all(data.sum);

  return complex<double>(real(data.sum), imag(data.sum));
}

typedef struct { 
  field_rfunction integrand; void *integrand_data; 
} rfun_wrap_data;
static complex<double> rfun_wrap(const complex<double> *fields,
				 const vec &loc, void *data_) {
  rfun_wrap_data *data = (rfun_wrap_data *) data_;
  return data->integrand(fields, loc, data->integrand_data);
}

double fields::integrate(int num_fields, const component *components,
			 field_rfunction integrand,
			 const geometric_volume &where,
			 void *integrand_data_,
			 double *maxabs)
{
  rfun_wrap_data data;
  data.integrand = integrand;
  data.integrand_data = integrand_data_;
  return real(integrate(num_fields, components, rfun_wrap,
			where, &data, maxabs));
}

double fields::max_abs(int num_fields, const component *components,
		       field_function integrand,
		       const geometric_volume &where,
		       void *integrand_data_)
{
  double maxabs;
  integrate(num_fields, components, integrand, where, integrand_data_,
	    &maxabs);
  return maxabs;
}

double fields::max_abs(int num_fields, const component *components,
		       field_rfunction integrand,
		       const geometric_volume &where,
		       void *integrand_data_)
{
  rfun_wrap_data data;
  data.integrand = integrand;
  data.integrand_data = integrand_data_;
  return max_abs(num_fields, components, rfun_wrap, where, &data);
}

static complex<double> return_the_field(const complex<double> *fields,
					const vec &loc,
					void *integrand_data_)
{
  (void) integrand_data_; (void) loc; // unused
  return fields[0];
}

double fields::max_abs(int c, const geometric_volume &where)
{
  if (is_derived(c))
    return max_abs(derived_component(c), where);
  else
    return max_abs(component(c), where);
}

double fields::max_abs(component c, const geometric_volume &where)
{
  if (is_derived(int(c)))
    return max_abs(derived_component(c), where);
  return max_abs(1, &c, return_the_field, where, 0);
}

double fields::max_abs(derived_component c, const geometric_volume &where)
{
  if (!is_derived(int(c)))
    return max_abs(component(c), where);
  int nfields;
  component cs[6];
  field_rfunction fun = derived_component_func(c, v, nfields, cs);
  return max_abs(nfields, cs, fun, where, &nfields);
}

} // namespace meep
