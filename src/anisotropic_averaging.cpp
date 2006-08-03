#include <math.h>

#include "meep_internals.hpp"

/* This file contains routines to compute the "average" or "effective"
   dielectric constant for a pixel, using an anisotropic averaging
   procedure described in an upcoming paper (similar to the one in
   MPB). */

namespace meep {

////////////////////////////////////////////////////////////////////////////

#include "sphere-quad.h"

static vec sphere_pt(const vec &cent, double R, int n, double &weight) {
     switch (cent.dim) {
	 case D1:
	 {
	      weight = sphere_quad[0][n][3];
	      vec pt(sphere_quad[0][n][2]);
	      return cent + pt * R;
	 }
	 case D2:
	 {
	      weight = sphere_quad[1][n][3];
	      vec pt(sphere_quad[1][n][0], sphere_quad[1][n][1]);
	      return cent + pt * R;
	 }
	 case D3:
	 {
	      weight = sphere_quad[2][n][3];
	      vec pt(sphere_quad[2][n][0], sphere_quad[2][n][1],
		     sphere_quad[2][n][2]);
	      return cent + pt * R;
	 }
	 case Dcyl:
	 {
	      weight = sphere_quad[1][n][3];
	      return cent 
		+ veccyl(sphere_quad[1][n][0], sphere_quad[1][n][1]) * R;
	 }
         default:
	   abort("unknown dimensions in sphere_pt\n");
     }
}

////////////////////////////////////////////////////////////////////////////

vec material_function::normal_vector(const geometric_volume &gv) {
  // detect 2D line averages
  if (gv.get_min_corner().x()==gv.get_max_corner().x() || gv.get_min_corner().y()==gv.get_max_corner().y())
    return zero_vec(gv.dim);
  vec p(gv.center());
  double R = gv.diameter();  
  vec gradient = zero_vec(p.dim);
  for (int i = 0; i < num_sphere_quad[number_of_directions(gv.dim) - 1]; ++i) {
    double weight;
    vec pt = sphere_pt(p, R, i, weight);
    gradient += (pt - p) * (weight * eps(pt));
  }
  return gradient;
}

/* default: simple numerical integration of surfaces/cubes, relative
   tolerance 'tol'.   This is superseded by the routines in the libctl
   interface, which either use a semi-analytical average or can
   use a proper adaptive cubature. */
void material_function::meaneps(double &meps, double &minveps, vec &gradient, const geometric_volume &gv, double tol, int maxeval) {
  gradient = normal_vector(gv);
  if (sqrt(gradient & gradient) < 1e-8) {
    meps = eps(gv.center());
    minveps = 1/meps;
    return;
  }
  vec d = gv.get_max_corner() - gv.get_min_corner();
  int ms = 10; 
  double old_meps=0, old_minveps=0;
  int iter = 0;
  meps=1; minveps=1;
  switch(gv.dim) {
  case D3:
    while ((fabs(meps - old_meps) > tol*fabs(old_meps)) && (fabs(minveps - old_minveps) > tol*fabs(old_minveps))) {
      old_meps=meps; old_minveps=minveps;
      meps = minveps = 0;
      for (int k=0; k < ms; k++)
	for (int j=0; j < ms; j++)
	  for (int i=0; i < ms; i++) {
	    double ep = eps(gv.get_min_corner() + vec(i*d.x()/ms, j*d.y()/ms, k*d.z()/ms));
	    meps += ep; minveps += 1/ep;
	  }
      meps /= ms*ms*ms;
      minveps /= ms*ms*ms;
      ms *= 2;
      if (maxeval && (iter += ms*ms*ms) >= maxeval) return;
    }
    break;
  case D2:
    while ((fabs(meps-old_meps) > tol*old_meps) && (fabs(minveps-old_minveps) > tol*old_minveps)) {
      old_meps=meps; old_minveps=minveps;
      meps = minveps = 0;
      for (int j=0; j < ms; j++)
	for (int i=0; i < ms; i++) {
	  double ep = eps(gv.get_min_corner() + vec(i*d.x()/ms, j*d.y()/ms));
	  meps += ep; minveps += 1/ep;
	}
      meps /= ms*ms;
      minveps /= ms*ms;
      ms *= 2; 
      if (maxeval && (iter += ms*ms) >= maxeval) return;
    }
    break;
  case Dcyl:
    while ((fabs(meps-old_meps) > tol*old_meps) && (fabs(minveps-old_minveps) > tol*old_minveps)) {
      old_meps=meps; old_minveps=minveps;
      meps = minveps = 0;
      for (int j=0; j < ms; j++)
	for (int i=0; i < ms; i++) {
	  double ep = eps(gv.get_min_corner() + veccyl(i*d.r()/ms, j*d.z()/ms));
	  meps += ep; minveps += 1/ep;
	}
      meps /= ms*ms;
      minveps /= ms*ms;
      ms *= 2; 
      if (maxeval && (iter += ms*ms) >= maxeval) return;
    }
    break;
  case D1:
    while ((fabs(meps-old_meps) > tol*old_meps) && (fabs(minveps-old_minveps) > tol*old_minveps)) {
      old_meps=meps; old_minveps=minveps;
      meps = minveps = 0;
      for (int i=0; i < ms; i++) {
	double ep = eps(gv.get_min_corner() + vec(i*d.z()/ms));
	meps += ep; minveps += 1/ep;
      }
      meps /= ms;
      minveps /= ms;
      ms *= 2; 
      if (maxeval && (iter += ms*ms) >= maxeval) return;
    }
    break;
  }
}

////////////////////////////////////////////////////////////////////////////

static void anisoaverage(material_function &epsilon, const geometric_volume dV,
			 component ec, double invepsrow[3],
			 double tol, int maxeval) {    
  double meps = 0, minveps = 0;
  vec norm(dV.dim);

  epsilon.meaneps(meps,minveps,norm,dV,tol,maxeval); 

  double n[3] = {0,0,0};
  LOOP_OVER_DIRECTIONS(norm.dim, k) n[k%3] = norm.in_direction(k);
  if (abs(norm) < 1e-8) { /* couldn't find normal: just use meps */
    minveps = 1/meps;
    n[0] = 1; n[1] = 0; n[2] = 0;
  }
  else {
    double nabsinv = 1/abs(norm);
    for (int i=0; i<3; ++i) n[i] *= nabsinv;
  }

  /* get rownum'th row of effective tensor
          P * minveps + (I-P) * 1/meps = P * (minveps-1/meps) + I * 1/meps
     where I is the identity and P is the projection matrix
     P_{ij} = n[i] * n[j]. */
  int rownum = component_direction(ec) % 3;
  for (int i=0; i<3; ++i) 
    invepsrow[i] = n[rownum] * n[i] * (minveps - 1/meps);
  invepsrow[rownum] += 1/meps;
}

void structure_chunk::set_epsilon(material_function &epsilon,
				  bool use_anisotropic_averaging,
				  double tol, int maxeval) {
  if (!is_mine()) return;

  epsilon.set_volume(v.pad().surroundings());

  if (!eps) eps = new double[v.ntot()];
  LOOP_OVER_VOL(v, v.eps_component(), i) {
    IVEC_LOOP_LOC(v, here);
    eps[i] = epsilon.eps(here);
  }
  
  if (!use_anisotropic_averaging) {
    FOR_ELECTRIC_COMPONENTS(c)
      if (v.has_field(c)) {
#if 1 // legacy method: very simplistic averaging (TODO: delete this?)
	bool have_other_direction = false;
	vec dxa = zero_vec(v.dim);
	vec dxb = zero_vec(v.dim);
	direction c_d = component_direction(c);
	LOOP_OVER_DIRECTIONS(v.dim,da)
	  if (da != c_d) {
	    dxa.set_direction(da,0.5/a);
	    LOOP_OVER_DIRECTIONS(v.dim,db)
	      if (db != c_d && db != da) {
		dxb.set_direction(db,0.5/a);
		have_other_direction = true;
	      }
	    break;
	  }
	if (!inveps[c][c_d]) inveps[c][c_d] = new double[v.ntot()];
	LOOP_OVER_VOL(v, c, i) {
	  IVEC_LOOP_LOC(v, here);
	  if (!have_other_direction)
	    inveps[c][c_d][i] =
	      2.0/(epsilon.eps(here + dxa) + epsilon.eps(here - dxa));
	  else
	    inveps[c][c_d][i] = 4.0/(epsilon.eps(here + dxa + dxb) +
				     epsilon.eps(here + dxa - dxb) +
				     epsilon.eps(here - dxa + dxb) +
				     epsilon.eps(here - dxa - dxb));
	}
#else // really no averaging at all
	direction c_d = component_direction(c);
        if (!inveps[c][c_d]) inveps[c][c_d] = new double[v.ntot()]; 
        LOOP_OVER_VOL(v, c, i) {
          IVEC_LOOP_LOC(v, here);
          inveps[c][c_d][i] = 1/epsilon.eps(here);
	}
#endif
      }
  } else {
    const double smoothing_diameter = 1.0; // FIXME: make this user-changable?

    // may take a long time in 3d, so prepare to print status messages
    int npixels = 0, ipixel = 0;
    FOR_ELECTRIC_COMPONENTS(c) if (v.has_field(c)) {
      int loop_npixels = 0;
      LOOP_OVER_VOL(v, c, i) {
	loop_npixels = loop_n1 * loop_n2 * loop_n3;
	goto breakout;
      }
    breakout: npixels += loop_npixels;
    }
    double last_output_time = wall_time();

    FOR_ELECTRIC_COMPONENTS(c)
      if (v.has_field(c)) {
	FOR_ELECTRIC_COMPONENTS(c2) if (v.has_field(c2)) {
	  direction d = component_direction(c2);
	  if (!inveps[c][d]) inveps[c][d] = new double[v.ntot()];
	  if (!inveps[c][d]) abort("Memory allocation error.\n");
	}
	direction d0 = X, d1 = Y, d2 = Z;
	if (v.dim == Dcyl) { d0 = R; d1 = P; }
	LOOP_OVER_VOL(v, c, i) {
	  double invepsrow[3];
	  IVEC_LOOP_ILOC(v, here);
	  anisoaverage(epsilon, v.dV(here, smoothing_diameter), 
		       c, invepsrow, tol, maxeval);
	  if (inveps[c][d0]) inveps[c][d0][i] = invepsrow[0];
	  if (inveps[c][d1]) inveps[c][d1][i] = invepsrow[1];
	  if (inveps[c][d2]) inveps[c][d2][i] = invepsrow[2];

	  if (!quiet && (ipixel+1) % 1000 == 0
	      && wall_time() > last_output_time + MIN_OUTPUT_TIME) {
	    master_printf("subpixel-averaging is %g%% done, %g s remaining\n", 
			  ipixel * 100.0 / npixels,
			  (npixels - ipixel) *
			  (wall_time() - last_output_time) / ipixel);
	    last_output_time = wall_time();
	  }
	  ++ipixel;
	}
      }
  }

  update_pml_arrays(); // PML stuff depends on epsilon
}

} // namespace meep
