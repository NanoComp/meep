#include <math.h>

#include "meep_internals.hpp"

/* This file contains routines to compute the "average" or "effective"
   dielectric constant for a pixel, using an anisotropic averaging
   procedure described in an upcoming paper (similar to the one in
   MPB). */

using namespace std;

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

vec material_function::normal_vector(field_type ft, const volume &v)
{
  vec gradient(zero_vec(v.dim));
  vec p(v.center());
  double R = v.diameter();

  for (int i = 0; i < num_sphere_quad[number_of_directions(v.dim)-1]; ++i) {
    double weight;
    vec pt = sphere_pt(p, R, i, weight);
    gradient += (pt - p) * (weight * chi1p1(ft,pt));
  }
  return gradient;
}

/* default: simple numerical integration of surfaces/cubes, relative
   tolerance 'tol'.   This is superseded by the routines in the libctl
   interface, which either use a semi-analytical average or can
   use a proper adaptive cubature. */
void material_function::eff_chi1inv_row(component c, double chi1inv_row[3],
					const volume &v,
					double tol, int maxeval) {
  field_type ft = type(c);
  if (!maxeval) {
  trivial:
    chi1inv_row[0] = chi1inv_row[1] = chi1inv_row[2] = 0.0;
    chi1inv_row[component_direction(c) % 3] = 1/chi1p1(ft,v.center());
    return;
  }

  vec gradient(normal_vector(ft, v));
  if (abs(gradient) < 1e-8) goto trivial;

  double meps=1, minveps=1;
  vec d = v.get_max_corner() - v.get_min_corner();
  int ms = 10; 
  double old_meps=0, old_minveps=0;
  int iter = 0;
  switch(v.dim) {
  case D3:
    while ((fabs(meps - old_meps) > tol*fabs(old_meps)) && (fabs(minveps - old_minveps) > tol*fabs(old_minveps))) {
      old_meps=meps; old_minveps=minveps;
      meps = minveps = 0;
      for (int k=0; k < ms; k++)
	for (int j=0; j < ms; j++)
	  for (int i=0; i < ms; i++) {
	    double ep = chi1p1(ft,v.get_min_corner() + vec(i*d.x()/ms, j*d.y()/ms, k*d.z()/ms));
	    if (ep < 0) goto trivial;
	    meps += ep; minveps += 1/ep;
	  }
      meps /= ms*ms*ms;
      minveps /= ms*ms*ms;
      ms *= 2;
      if (maxeval && (iter += ms*ms*ms) >= maxeval) goto done;
    }
    break;
  case D2:
    while ((fabs(meps-old_meps) > tol*old_meps) && (fabs(minveps-old_minveps) > tol*old_minveps)) {
      old_meps=meps; old_minveps=minveps;
      meps = minveps = 0;
      for (int j=0; j < ms; j++)
	for (int i=0; i < ms; i++) {
	  double ep = chi1p1(ft,v.get_min_corner() + vec(i*d.x()/ms, j*d.y()/ms));
	  if (ep < 0) goto trivial;
	  meps += ep; minveps += 1/ep;
	}
      meps /= ms*ms;
      minveps /= ms*ms;
      ms *= 2; 
      if (maxeval && (iter += ms*ms) >= maxeval) goto done;
    }
    break;
  case Dcyl:
    while ((fabs(meps-old_meps) > tol*old_meps) && (fabs(minveps-old_minveps) > tol*old_minveps)) {
      old_meps=meps; old_minveps=minveps;
      meps = minveps = 0;
      double sumvol = 0;
      for (int j=0; j < ms; j++)
	for (int i=0; i < ms; i++) {
	  double r = v.get_min_corner().r() + i*d.r()/ms;
	  double ep = chi1p1(ft,v.get_min_corner() + veccyl(i*d.r()/ms, j*d.z()/ms));
	  if (ep < 0) goto trivial;
	  sumvol += r;
	  meps += ep * r; minveps += r/ep;
	}
      meps /= sumvol;
      minveps /= sumvol;
      ms *= 2; 
      if (maxeval && (iter += ms*ms) >= maxeval) goto done;
    }
    break;
  case D1:
    while ((fabs(meps-old_meps) > tol*old_meps) && (fabs(minveps-old_minveps) > tol*old_minveps)) {
      old_meps=meps; old_minveps=minveps;
      meps = minveps = 0;
      for (int i=0; i < ms; i++) {
	double ep = chi1p1(ft,v.get_min_corner() + vec(i*d.z()/ms));
	if (ep < 0) {
	  meps = chi1p1(ft,v.center());
	  minveps = 1/meps;
	  goto done;
	}
	meps += ep; minveps += 1/ep;
      }
      meps /= ms;
      minveps /= ms;
      ms *= 2; 
      if (maxeval && (iter += ms*ms) >= maxeval) goto done;
    }
    break;
  }

 done:
  {
    double n[3] = {0,0,0};
    double nabsinv = 1.0/abs(gradient);
    LOOP_OVER_DIRECTIONS(gradient.dim, k)
      n[k%3] = gradient.in_direction(k) * nabsinv;
    
    /* get rownum'th row of effective tensor
       P * minveps + (I-P) * 1/meps = P * (minveps-1/meps) + I * 1/meps
       where I is the identity and P is the projection matrix
       P_{ij} = n[i] * n[j]. */
    int rownum = component_direction(c) % 3;
    for (int i=0; i<3; ++i) 
      chi1inv_row[i] = n[rownum] * n[i] * (minveps - 1/meps);
    chi1inv_row[rownum] += 1/meps;
  }
}

void structure_chunk::set_chi1inv(component c,
				  material_function &medium,
				  bool use_anisotropic_averaging,
				  double tol, int maxeval) {
  if (!is_mine() || !gv.has_field(c)) return;
  field_type ft = type(c);
  if (ft != E_stuff && ft != H_stuff) abort("only E or H can have chi");
  medium.set_volume(gv.pad().surroundings());

  if (!use_anisotropic_averaging) maxeval = 0;

  const double smoothing_diameter = 1.0; // FIXME: make user-changable?
      
  // may take a long time in 3d, so prepare to print status messages
  int npixels = 0, ipixel = 0;
  int loop_npixels = 0;
  LOOP_OVER_VOL(gv, c, i) {
    loop_npixels = loop_n1 * loop_n2 * loop_n3;
    goto breakout; // hack to use loop-size computation from LOOP_OVER_VOL
  }
 breakout: npixels += loop_npixels;
  double last_output_time = wall_time();

  FOR_FT_COMPONENTS(ft,c2) if (gv.has_field(c2)) {
    direction d = component_direction(c2);
    if (!chi1inv[c][d]) chi1inv[c][d] = new realnum[gv.ntot()];
    if (!chi1inv[c][d]) abort("Memory allocation error.\n");
  }
  direction dc = component_direction(c);
  direction d0 = X, d1 = Y, d2 = Z;
  if (gv.dim == Dcyl) { d0 = R; d1 = P; }
  int idiag = component_index(c);
  bool trivial[3] = {true,true,true};
  double trivial_val[3] = {0,0,0};
  trivial_val[idiag] = 1.0;
  ivec shift1(unit_ivec(gv.dim,component_direction(c))
	      * (ft == E_stuff ? 1 : -1));
  LOOP_OVER_VOL(gv, c, i) {
    double chi1invrow[3], chi1invrow_offdiag[3];
    IVEC_LOOP_ILOC(gv, here);
    medium.eff_chi1inv_row(c, chi1invrow,
			   gv.dV(here,smoothing_diameter), tol,maxeval);
    medium.eff_chi1inv_row(c, chi1invrow_offdiag,
			   gv.dV(here-shift1,smoothing_diameter), tol,maxeval);
    if (chi1inv[c][d0]) {
      chi1inv[c][d0][i] = (d0 == dc) ? chi1invrow[0] : chi1invrow_offdiag[0];
      trivial[0] = trivial[0] && (chi1inv[c][d0][i] == trivial_val[0]);
    }
    if (chi1inv[c][d1]) {
      chi1inv[c][d1][i] = (d1 == dc) ? chi1invrow[1] : chi1invrow_offdiag[1];
      trivial[1] = trivial[1] && (chi1inv[c][d1][i] == trivial_val[1]);
    }
    if (chi1inv[c][d2]) {
      chi1inv[c][d2][i] = (d2 == dc) ? chi1invrow[2] : chi1invrow_offdiag[2];
      trivial[2] = trivial[2] && (chi1inv[c][d2][i] == trivial_val[2]);
    }
    
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
  direction ds[3]; ds[0] = d0; ds[1] = d1; ds[2] = d2;
  for (int i = 0; i < 3; ++i) {
    trivial_chi1inv[c][ds[i]] = trivial[i];
    if (i != idiag && trivial[i]) { // deallocate trivial offdiag
      delete[] chi1inv[c][ds[i]]; 
      chi1inv[c][ds[i]] = 0; 
    }
  }
  // only deallocate trivial diag if entire tensor is trivial
  if (trivial[0] && trivial[1] && trivial[2]) {
    delete[] chi1inv[c][dc];
    chi1inv[c][dc] = 0;
  }
  medium.unset_volume();
}


void structure_chunk::add_susceptibility(material_function &sigma,
                                         field_type ft,
                                         const susceptibility &sus)
{
    if (ft != E_stuff && ft != H_stuff)
        abort("susceptibilities must be for E or H fields");

    sigma.set_volume(gv.pad().surroundings());

    susceptibility *newsus = sus.clone();
    newsus->next = NULL;
    newsus->ntot = gv.ntot();
    // get rid of previously allocated sigma, normally not the case here:
    FOR_COMPONENTS(c) FOR_DIRECTIONS(d) if (newsus->sigma[c][d]) {
        delete[] newsus->sigma[c][d];
        newsus->sigma[c][d] = NULL;
        newsus->trivial_sigma[c][d] = true;
    }

    // if we own this chunk, set up the sigma array(s):
    if (is_mine()) FOR_FT_COMPONENTS(ft,c) if (gv.has_field(c)) {
        FOR_FT_COMPONENTS(ft,c2) if (gv.has_field(c2)) {
            direction d = component_direction(c2);
            if (!newsus->sigma[c][d]) newsus->sigma[c][d] = new realnum[gv.ntot()];
            if (!newsus->sigma[c][d]) abort("Memory allocation error.\n");
        }
        bool trivial[3] = {true, true, true};
        direction dc = component_direction(c);
        direction d0 = X, d1 = Y, d2 = Z;
        if (gv.dim == Dcyl) { d0 = R; d1 = P; }
        int idiag = component_index(c);
        realnum *s0 = newsus->sigma[c][d0];
        realnum *s1 = newsus->sigma[c][d1];
        realnum *s2 = newsus->sigma[c][d2];
        vec shift1(gv[unit_ivec(gv.dim,component_direction(c))
        * (ft == E_stuff ? 1 : -1)]);
        LOOP_OVER_VOL(gv, c, i) {
            double sigrow[3], sigrow_offdiag[3];
            IVEC_LOOP_LOC(gv, here);
            sigma.sigma_row(c, sigrow, here);
            sigma.sigma_row(c, sigrow_offdiag, here - shift1);
            sigrow[(idiag+1) % 3] = sigrow_offdiag[(idiag+1) % 3];
            sigrow[(idiag+2) % 3] = sigrow_offdiag[(idiag+2) % 3];
            if (s0 && (s0[i] = sigrow[0]) != 0.) trivial[0] = false;
            if (s1 && (s1[i] = sigrow[1]) != 0.) trivial[1] = false;
            if (s2 && (s2[i] = sigrow[2]) != 0.) trivial[2] = false;
        }

        direction ds[3]; ds[0] = d0; ds[1] = d1; ds[2] = d2;
        for (int i = 0; i < 3; ++i) {
            newsus->trivial_sigma[c][ds[i]] = trivial[i];
            if (i != idiag && trivial[i]) { // deallocate trivial offdiag
                delete[] newsus->sigma[c][ds[i]];
                newsus->sigma[c][ds[i]] = 0;
            }
        }
        // only deallocate trivial diag if entire tensor is trivial
        if (trivial[0] && trivial[1] && trivial[2]) {
            delete[] newsus->sigma[c][dc];
            newsus->sigma[c][dc] = 0;
        }
    }

    // finally, add to the beginning of the chiP list:
    newsus->next = chiP[ft];
    chiP[ft] = newsus;

    sigma.unset_volume();
}


////////////////////////////////////////////////////////////
//--- Begin: defs for polygon-based material functions ---//
////////////////////////////////////////////////////////////

// Returns the absolute area of a polygon.
// For self-intersecting polygons, this will return the difference between clockwise
// parts and anti-clockwise parts, not the real area.
// Last point must equal first point!
static double get_polygon_area(const double* const* matrix, std::size_t dim1, std::size_t dim2) {
    if (dim1 < 3)
        return 0;

    if (dim2 != 2) {
        abort("get_polygon_area: Array with polygon points should have a shape X,2. The second dimension is not 2.");
    }

    //2A = \sum_{i=0}^{n-1}( (x_{i-1} - x_{i+1}) * y_i ) # Gauss's area formula (aka shoelace algorithm)
    double sum = (matrix[dim1 - 1][0] - matrix[1][0]) * matrix[0][1]; // (i = 0)
    for (std::size_t i = 1; i < dim1 - 1; ++i)
        sum += (matrix[i - 1][0] - matrix[i + 1][0]) * matrix[i][1];
    sum += (matrix[dim1 - 2][0] - matrix[0][0]) * matrix[dim1 - 1][1]; // (i = n-1)
    return fabs(0.5 * sum);
}

simple_polygon::simple_polygon(const double (* const points)[2], std::size_t num_points) {
    if (num_points > 0) {
        //make sure last point equals first point. If not, append copy of first point to end:
        if (points[0][0] != points[num_points - 1][0] || points[0][1] != points[num_points - 1][1])
            _number_of_points = num_points + 1;
        else
            _number_of_points = num_points;

        // copy matrix:
        _polygon_points = new double*[_number_of_points];
        for (std::size_t i = 0; i < num_points; i++) {
            _polygon_points[i] = new double[2];
            _polygon_points[i][0] = points[i][0];
            _polygon_points[i][1] = points[i][1];
        }

        //if last point does not equal first point -> append copy of first point to end:
        if (_number_of_points == num_points + 1) {
            _polygon_points[num_points] = new double[2];
            _polygon_points[num_points][0] = points[0][0];
            _polygon_points[num_points][1] = points[0][1];
        }
    }
    else {
        _number_of_points = 0;
        _polygon_points = 0;
    }

}

simple_polygon::simple_polygon(const double* matrix, std::size_t dim1, std::size_t dim2) {
    if (dim1 > 0) {
        if (dim2 != 2) {
            abort("Array with polygon points should have a shape X,2. The second dimension is not 2.");
        }

        //make sure last point equals first point. If not, append copy of first point to end:
        if (matrix[0] != matrix[2 * (dim1 - 1)] || matrix[1] != matrix[2 * dim1 - 1])
            this->_number_of_points = dim1 + 1;
        else
            this->_number_of_points = dim1;

        // copy matrix:
        _polygon_points = new double*[_number_of_points];
        for (std::size_t i = 0; i < dim1; i++) {
            _polygon_points[i] = new double[2];
            _polygon_points[i][0] = matrix[2 * i];
            _polygon_points[i][1] = matrix[2 * i + 1];
        }

        //if last point does not equal first point -> append copy of first point to end:
        if (this->_number_of_points == dim1 + 1) {
            _polygon_points[dim1] = new double[2];
            _polygon_points[dim1][0] = matrix[0];
            _polygon_points[dim1][1] = matrix[1];
        }
    }
    else {
        _number_of_points = 0;
        _polygon_points = 0;
    }
};

simple_polygon::simple_polygon(const simple_polygon &pol) {
    this->_number_of_points = pol._number_of_points;
    // copy matrix:
    if (_number_of_points > 0) {
        _polygon_points = new double*[_number_of_points];
        for (std::size_t i = 0; i < _number_of_points; i++) {
            _polygon_points[i] = new double[2];
            _polygon_points[i][0] = pol._polygon_points[i][0];
            _polygon_points[i][1] = pol._polygon_points[i][1];
        }
    }
    else {
        _polygon_points = 0;
    }
}

void swap(simple_polygon& first, simple_polygon& second)
{
    using std::swap;
    swap(first._number_of_points, second._number_of_points);
    swap(first._polygon_points, second._polygon_points);
}

simple_polygon& simple_polygon::operator=(simple_polygon other)
{
    swap(*this, other);
    return *this;
}

simple_polygon::~simple_polygon() {
    for (std::size_t i = 0; i < _number_of_points; i++) {
        delete [] _polygon_points[i];
    }
    delete [] _polygon_points;
};

// Clip the polygon by a rectangular boundary:
void simple_polygon::clip_polygon(
    double left, double right, double bottom, double top) {
    // It's easier to do this in two steps: First clip the polygon by left and
    // right only, then clip the resulting polygon by bottom and top.

    // Note: This whole method was only created so we can compare the area of
    // the original (apart from this clipping), unsplitted polygons with the
    // total area of the splitted polygons, which must be equal, in
    // polygonal_material::update_splitting(). The splitting also
    // works fine if we don't clip the polygons to the computational cell
    // before, only it will then output a warning that the areas don't match.
    // This probably was not worth the effort, but anyway, now the method is
    // already written and can be used. ;)

    if (_number_of_points == 0)
        return;

    // Go through points, keep points inside boundaries,
    // add intersecting points when transitioning from outside->inside or
    // vice versa.

    // First step, clip by left and right:
    vector<double*> poly1;
    double *intp_pt, *pt, *prev_pt = _polygon_points[0];
    double slope;
    // -1 for left, 0 for inside, +1 for right of border:
    int pos, prev_pos = int(prev_pt[0] > right) - int(prev_pt[0] < left);
    for (std::size_t i = 0; i < _number_of_points; i++) {
        pt = _polygon_points[i];
        // -1 for left, 0 for inside, +1 for right of border:
        pos = int(pt[0] > right) - int(pt[0] < left);
        if (prev_pos != 0 && prev_pos != pos) {
            // The polygon entered or crossed the inner area.
            // Calc intersection point at the entry to inner area:
            slope = (pt[1] - prev_pt[1]) / (pt[0] - prev_pt[0]);
            intp_pt = new double[2];
            intp_pt[0] = prev_pos == -1 ? left : right;
            intp_pt[1] = prev_pt[1] + (intp_pt[0] - prev_pt[0]) * slope;
            // add intersection point:
            poly1.push_back(intp_pt);
        }
        if (pos == 0) {
            // add pt (which is inside):
            poly1.push_back(pt);
            // This point is now owned by poly1:
            _polygon_points[i] = 0;
        }
        else if (prev_pos != pos) {
            // The polygon left or crossed the inner area.
            // Calc intersection point at the exit from inner area:
            slope = (pt[1] - prev_pt[1]) / (pt[0] - prev_pt[0]);
            intp_pt = new double[2];
            intp_pt[0] = pos == -1 ? left : right;
            intp_pt[1] = prev_pt[1] + (intp_pt[0] - prev_pt[0]) * slope;
            // add intersection point:
            poly1.push_back(intp_pt);
        }
        prev_pt = pt;
        prev_pos = pos;
    }
    // End of first step; Add first point if needed:
    if (poly1.size() > 1 && (
        poly1.back()[0] != poly1.front()[0] ||
        poly1.back()[1] != poly1.front()[1]))
    {
        pt = new double[2];
        pt[0] = poly1.front()[0];
        pt[1] = poly1.front()[1];
        poly1.push_back(pt);
    }

    if (poly1.empty()){
        for (std::size_t i = 0; i < _number_of_points; i++) {
            delete[] _polygon_points[i];
        }
        delete[] _polygon_points;
        _number_of_points = 0;
        _polygon_points = 0;
        return;
    }

    // Second step, clip by top and bottom:
    vector<double*> poly2;
    prev_pt = poly1[0];
    // -1 for below, 0 for inside, +1 for above border:
    prev_pos = int(prev_pt[1] > top) - int(prev_pt[1] < bottom);
    for (std::size_t i = 0; i < poly1.size(); i++) {
        pt = poly1[i];
        // -1 for below, 0 for inside, +1 for above border:
        pos = int(pt[1] > top) - int(pt[1] < bottom);
        if (prev_pos != 0 && pos != prev_pos) {
            // The polygon entered or crossed the inner area.
            // Calc intersection point at the entry to inner area:
            slope = (pt[0] - prev_pt[0]) / (pt[1] - prev_pt[1]);
            intp_pt = new double[2];
            intp_pt[1] = prev_pos == -1 ? bottom : top;
            intp_pt[0] = prev_pt[0] + (intp_pt[1] - prev_pt[1]) * slope;
            // add intersection point:
            poly2.push_back(intp_pt);
        }
        if (pos == 0) {
            // add pt (which is inside):
            poly2.push_back(pt);
            // This point is now owned by poly2:
            poly1[i] = 0;
        }
        else if (pos != prev_pos) {
            // The polygon left or crossed the inner area.
            // Calc intersection point at the exit from inner area:
            slope = (pt[0] - prev_pt[0]) / (pt[1] - prev_pt[1]);
            intp_pt = new double[2];
            intp_pt[1] = pos == -1 ? bottom : top;
            intp_pt[0] = prev_pt[0] + (intp_pt[1] - prev_pt[1]) * slope;
            // add intersection point:
            poly2.push_back(intp_pt);
        }
        prev_pt = pt;
        prev_pos = pos;
    }

    // End of second step; Add first point if needed:
    if (poly2.size() > 1 && (
        poly2.back()[0] != poly2.front()[0] ||
        poly2.back()[1] != poly2.front()[1]))
    {
        pt = new double[2];
        pt[0] = poly2.front()[0];
        pt[1] = poly2.front()[1];
        poly2.push_back(pt);
    }

    // finally set the polygon's points to the clipped ones in poly2:
    for (std::size_t i = 0; i < _number_of_points; i++) {
        delete[] _polygon_points[i];
    }
    delete[] _polygon_points;
    _number_of_points = poly2.size();
    if (_number_of_points > 0) {
        _polygon_points = new double*[_number_of_points];
        for (std::size_t i = 0; i < _number_of_points; i++) {
            _polygon_points[i] = poly2[i];
        }
    }
    else
        _polygon_points = NULL;
    // Now all points in poly2 belong to _polygons.

    // clean up:
    for (std::size_t i = 0; i < poly1.size(); i++) {
        delete[] poly1[i];
    }
    poly1.clear();
    poly2.clear();
}

double simple_polygon::get_area() {
    return get_polygon_area(_polygon_points, _number_of_points, 2);
}

polygon::polygon(const polygon &pol) : simple_polygon(pol) {
    for(vector<simple_polygon*>::const_iterator p = pol._inner_polygons.begin();
        p != pol._inner_polygons.end(); ++p)
    {
        _inner_polygons.push_back(new simple_polygon(**p));
    }
}

void swap(polygon& first, polygon& second)
{
    using std::swap;
    swap(static_cast<simple_polygon&>(first), static_cast<simple_polygon&>(second));
    swap(first._inner_polygons, second._inner_polygons);
}

polygon& polygon::operator=(polygon other)
{
    swap(*this, other);
    return *this;
}

polygon::~polygon() {
    //_polygon_points from parent class are deleted in parent class destructor.
    for(vector<simple_polygon*>::iterator p = _inner_polygons.begin(); p != _inner_polygons.end(); ++p) {
        delete (*p);
    }
};

// TODO? all add_inner_polygon: check if inner polygon is real inner polygon of parent:
// first point is inside parent polygon AND inner polygon's edges does not cross parent polygon's edges
void polygon::add_inner_polygon(const double (* const points)[2], std::size_t num_points) {
    _inner_polygons.push_back(new simple_polygon(points, num_points));
}

void polygon::add_inner_polygon(const double* matrix, std::size_t dim1, std::size_t dim2) {
    _inner_polygons.push_back(new simple_polygon(matrix, dim1, dim2));
}

void polygon::add_inner_polygon(const simple_polygon &pol) {
    _inner_polygons.push_back(new simple_polygon(pol));
}

// Clip the polygon by a rectangular boundary:
void polygon::clip_polygon(
    double left, double right, double bottom, double top) {
    for(vector<simple_polygon*>::iterator inner_p = _inner_polygons.begin();
        inner_p != _inner_polygons.end(); ++inner_p)
    {
        (*inner_p)->clip_polygon(left, right, bottom, top);
    }
    simple_polygon::clip_polygon(left, right, bottom, top);
}

double polygon::get_area() {
    double parea = simple_polygon::get_area();
    for(vector<simple_polygon*>::iterator inner_p = _inner_polygons.begin();
        inner_p != _inner_polygons.end(); ++inner_p)
    {
        //all inner polygons must be real inner polygons (i.e. completely inside parent polygon)
        parea -= (*inner_p)->get_area();
    }
    return parea;
}

material_stack::material_stack(
    const double* layer_thicknesses, const double* layer_epsilons,
    const int number_of_layers)
{
    _number_of_layers = number_of_layers;
    // see if we can remove or join some layers:
    for (int i = 0; i < number_of_layers; ++i )
    {
        if (
            layer_thicknesses[i] <= 0 ||
            (i > 0 && layer_epsilons[i] == layer_epsilons[i - 1]))
            _number_of_layers--;
    }

    if (_number_of_layers == 0)
        abort("Please specify at least one layer in material stack.");

    _high_z = new double[_number_of_layers];
    _epsilon = new double[_number_of_layers];
    double current_z = 0;
    int j = 0;
    for (int i = 0; i < number_of_layers; ++i) {
        if (layer_thicknesses[i] <= 0) {
            // don't include empty layer:
            continue;
        }
        if (j > 0 && layer_epsilons[j-1] == layer_epsilons[i]) {
            // join layers:
            current_z += layer_thicknesses[i];
            _high_z[j-1] = current_z;
            continue;
        }
        current_z += layer_thicknesses[i];
        _high_z[j] = current_z;
        _epsilon[j] = layer_epsilons[i];
        j++;
    }
};

material_stack::material_stack(const material_stack& stack) :
        _number_of_layers(stack._number_of_layers)
    {
    if (_number_of_layers == 0)
        abort("Please specify at least one layer in material stack.");

    // make copy of arrays:
    _high_z = new double[_number_of_layers];
    _epsilon = new double[_number_of_layers];
    for (int i = 0; i < _number_of_layers; ++i) {
        _high_z[i] = stack._high_z[i];
        _epsilon[i] = stack._epsilon[i];
    }
}

material_stack::~material_stack(){
    delete[] _high_z;
    delete[] _epsilon;
};

double material_stack::chi1p1(const double &z) const
{
    // if z <= 0, return epsilon of first layer
    for (int i = 0; i < _number_of_layers - 1; ++i) {
        if (z <= _high_z[i])
            return _epsilon[i];
    }
    // if z > total thickness, return epsilon of last layer:
    return _epsilon[_number_of_layers - 1];
}

double material_stack::normal(const double &lower_z, const double &higher_z) const
{
    double lower_eps, higher_eps;
    for (int i = 0; i < _number_of_layers - 1; ++i) {
        if (lower_z <= _high_z[i]) {
            lower_eps = _epsilon[i];
            if (higher_z <= _high_z[i]) {
                // both z values are in same layer:
                return 0.0;
            }
            for (int j = i + 1; j < _number_of_layers - 1; ++j) {
                if (higher_z <= _high_z[j]) {
                    higher_eps = _epsilon[j];
                    return 0.5 * (lower_eps - higher_eps);
                }
            }
            return 0.5 * (lower_eps - _epsilon[_number_of_layers - 1]);
        }
    }
    // both lower_z and higher_z are in topmost layer or above:
    return 0.0;

}

double material_stack::mean_eps(const double &lower_z, const double &higher_z) const
{
    double distance = higher_z - lower_z;
    double result = 0;
    for (int i = 0; i < _number_of_layers - 1; ++i) {
        if (lower_z <= _high_z[i]) {
            if (higher_z <= _high_z[i]) {
                // both z values are in same layer:
                return _epsilon[i];
            }
            result += (_high_z[i] - lower_z) * _epsilon[i];
            for (int j = i + 1; j < _number_of_layers - 1; ++j) {
                if (higher_z <= _high_z[j]) {
                    result += (higher_z - _high_z[j - 1]) * _epsilon[j];
                    return result / distance;
                }
                else {
                    // Add the thin layer in between. This should normally not
                    // occur, only with bad resolution and very thin layers.
                    result += (_high_z[j] - _high_z[j - 1]) * _epsilon[j];
                }
            }
            // higher_z is in topmost layer or above:
            result += (higher_z - _high_z[_number_of_layers - 2]) * _epsilon[_number_of_layers - 1];
            return result / distance;
        }
    }
    // both lower_z and higher_z are in topmost layer or above:
    return _epsilon[_number_of_layers - 1];
}

double material_stack::mean_inveps(const double &lower_z, const double &higher_z) const
{
    double distance = higher_z - lower_z;
    double result = 0;
    for (int i = 0; i < _number_of_layers - 1; ++i) {
        if (lower_z <= _high_z[i]) {
            if (higher_z <= _high_z[i]) {
                // both z values are in same layer:
                return 1.0 / _epsilon[i];
            }
            result += (_high_z[i] - lower_z) / _epsilon[i];
            for (int j = i + 1; j < _number_of_layers - 1; ++j) {
                if (higher_z <= _high_z[j]) {
                    result += (higher_z - _high_z[j - 1]) / _epsilon[j];
                    return result / distance;
                }
                else {
                    // Add the thin layer in between. This should normally not
                    // occur, only with bad resolution or very thin layers.
                    result += (_high_z[j] - _high_z[j - 1]) / _epsilon[j];
                }
            }
            // higher_z is in topmost layer or above:
            result += (higher_z - _high_z[_number_of_layers - 2]) / _epsilon[_number_of_layers - 1];
            return result / distance;
        }
    }
    // both lower_z and higher_z are in topmost layer or above:
    return 1.0 / _epsilon[_number_of_layers - 1];
}

polygonal_material::polygonal_material(
    const material_stack * mat_stack, double block_size, int nx, int ny, int nz)
{
    init_polymat(mat_stack, 0, block_size, nx, ny, nz);
}

polygonal_material::polygonal_material(
    double epsilon, double block_size, int nx, int ny)
{
    init_polymat(0, epsilon, block_size, nx, ny, 0);
}

void polygonal_material::init_polymat(
    const material_stack * mat_stack,
    double epsilon, double block_size, int nx, int ny, int nz)
{
    if (mat_stack) {
        _material_stack = new material_stack(*mat_stack);
    }
    else {
        _material_stack = 0;
    }
    _epsilon = epsilon;
    _block_size = block_size;
    _nx = nx;
    _ny = ny;
    _nz = nz;
    splitting_done = false;
    _areas = new double*[nx];
    for (int i = 0; i < nx; ++i){
        _areas[i] = new double[ny];
        for (int j = 0; j < ny; ++j) {
            _areas[i][j] = 0;
        }
    }
};

polygonal_material::~polygonal_material() {
    //delete polygons list:
    for (size_t i = 0; i < _polygons.size(); ++i) {
        //delete polygon object
        delete _polygons[i];
    }
    // remove pointers to deleted polygon objects:
    _polygons.clear();

    delete _material_stack;

    // delete _areas:
    for (int i = 0; i < _nx; ++i)
        delete[] _areas[i];
    delete[] _areas;
};

void polygonal_material::add_polygon(const polygon& pol)
{
    if (splitting_done)
        abort("Cannot add polygons after splitting has been done.");
    // make a copy that belongs to this class:
    polygon *cpol = new polygon(pol);
    // clip the polygon to the computational domain:
    cpol->clip_polygon(0, _nx * _block_size, 0, _ny * _block_size);
    if (cpol->get_number_of_points() > 2)
        _polygons.push_back(cpol);
    else {
        delete cpol;
    }
}

bool polygonal_material::is_trivial(const ivec &small_ivec, const ivec &big_ivec)
{
    if (!splitting_done)
        update_splitting();
    int i2, j2;
    double val = -2;
    for (int i = small_ivec.x(); i < big_ivec.x(); ++i)
    {
        // If pixel is at outer edge, overlapping the undefined area outside
        // the grid volume, wrap the structure to the opposite side. If a PML
        // boundary is used, there shouldn't be a structure close to the
        // boundary anyways, and if periodic boundary conditions are used, we
        // would want to wrap the structure:
        i2 = i;
        if (i2 < 0)
            i2 += _nx;
        else if (i2 >= _nx)
            i2 -= _nx;
        for (int j = small_ivec.y(); j < big_ivec.y(); ++j)
        {
            j2 = j;
            if (j2 < 0)
                j2 += _ny;
            else if (j2 >= _ny)
                j2 -= _ny;
            if (_areas[i2][j2] != 0 && _areas[i2][j2] != 1)
                return false;
            else if (val == -2)
                val = _areas[i2][j2];
            else if (val != _areas[i2][j2])
                return false;
        }
    }
    return true;
}

vec polygonal_material::normal_vector(const ivec &center)
{
    if (!splitting_done)
        update_splitting();

    // These 4 pointers will point to the 4x4 blocks (=2x2 pixels) that
    // will be used as parameter in get_2D_gradient:
    double* area_cols[4];
    int xi, yi;
    vec gradient;
    double inner_area;

    if (center.y() < 2 || center.y() + 1 >= _ny) {
        /**
            * If pixel is at outer edge, overlapping the undefined area
            * outside the grid volume, wrap the structure to the opposite side.
            * If a PML boundary is used, there shouldn't be a structure close
            * to the boundary anyways, and if periodic boundary conditions are
            * used, we would want to wrap the structure.
            */
        for (int i = 0; i < 4; ++i)
        {
            xi = center.x() - 2 + i;
            if (xi < 0)
                xi += _nx;
            else if (xi >= _nx)
                xi -= _nx;
            /**
                * If the 4x4 area wraps around the grid volume, we cannot
                * point to the existing memory location, but need to create
                * new arrays: */
            area_cols[i] = new double[4];
            for (int j = 0; j < 4; ++j)
            {
                yi = center.y() - 2 + j;
                if (yi < 0)
                    yi += _ny;
                else if (yi >= _ny)
                    yi -= _ny;
                area_cols[i][j] = _areas[xi][yi];
            }
        }
        inner_area = (
            area_cols[1][1] + area_cols[1][2] +
            area_cols[2][1] + area_cols[2][2]);
        gradient = get_2D_gradient(area_cols);
        for (int i = 0; i < 4; ++i)
        {
            delete[] area_cols[i];
        }
    }
    else
    {
        for (int i = 0; i < 4; ++i) {
            xi = center.x() - 2 + i;
            if (xi < 0)
                xi += _nx;
            else if (xi >= _nx)
                xi -= _nx;
            area_cols[i] = &_areas[xi][center.y() - 2];
        }
        inner_area = (
            area_cols[1][1] + area_cols[1][2] +
            area_cols[2][1] + area_cols[2][2]);
        gradient = get_2D_gradient(area_cols);
    }

    // Normalize gradient such that it can be added to gradient of
    // any neighboring polygonal_material in the same pixel. That is,
    // the returned gradient must have same length than it would have
    // if it was calculated with the Lebedev quadrature scheme, or at
    // least approximately.
    double meps = is_3D() ? _material_stack->mean_eps(
        (center.z() - 1) * _block_size,
        (center.z() + 1) * _block_size) : _epsilon;
    double len = (meps - 1) * sqrt(fabs((4 - inner_area) * inner_area));

    if (is_3D()){
        double grad_z = _material_stack->normal(
            (center.z() - 1) * _block_size, (center.z() + 1) * _block_size);
        if (abs(grad_z) < 1e-8)
        {
            // There is no interface in z direction.
            // The (len / 8)-factor is to make the result comparable with the
            // Lebedev scheme implementation for gradienQuad below.
            return gradient * len / 8;
        }
        else
        {
            // There is an interface in z direction
            if (abs(gradient) < 1e-8) {
                // there is ONLY an interface in z direction:
                return vec(0.0, 0.0, grad_z);
            }
            else
            {
                /**
                 * There's an interface in z and transversal direction. Such a
                 * corner cannot be solved exactly, it also isn't clear how to
                 * calculate the anisotropic averaging here. The best solution
                 * is to do the same than in material_function::normal_vector:
                 * That is, integrate chi*r over the unit sphere using the
                 * Lebedev quadrature scheme:
                 */
                return get_gradient_from_sphere_quad(center, gradient);
            }
        }
    }
    return gradient * len;
}

double polygonal_material::get_area(int i, int j) {
    if (!splitting_done)
        update_splitting();
    while (i < 0)
        i += _nx;
    while (i >= _nx)
        i -= _nx;
    while (j < 0)
        j += _ny;
    while (j >= _ny)
        j -= _ny;
    return _areas[i][j];
}

//     void polygonal_material::message_truncate(double* pt) {
//
//      I have deactivated this output, since we now clip the polygon to the
//      computational cell in add_polygon.
//         if (pt[0] + 1e-8 < 0 || pt[0] - 1e-8 > _nx * _block_size || //hardcoded epsilon: yuck!
//             pt[1] + 1e-8 < 0 || pt[1] - 1e-8 > _ny * _block_size)
//             // if point is exactly on upper border, this will also be called,
//             // but point will not really be truncated, so check bounds again.
//             master_printf(
//                 "point (%.8f, %.8f) outside computational grid (%f, %f) - polygon will be truncated.\n",
//                 pt[0], pt[1], _nx * _block_size, _ny * _block_size);
//     }

void polygonal_material::update_splitting()
{
    /**
     * Splits all _polygons into blocks with size _blocksize x _blocksize.
     *
     * This is a quick and efficient algorithm, which can be used for both
     * simple and for selfintersecting polygons, both convex and concave.
     *
     * It might produce ghost lines in the case of concave polygons (i.e.
     * concave polygons will never be split into more than two polygons by one
     * splitting border, with ghost lines connecting otherwise seperate
     * polygons). The additional area created by these ghost lines is zero, so
     * this should not be an issue, as the eps function is solely calculated
     * from areas. */


    // First, split polygons into columns, then afterwards split these columns
    // also into rows. If it is NOT done in this way, blocks completely
    // enclosed by a single polygon will be left out.

    // The column-wise split polygons will be temporary saved in col_arr,
    // the row-wise split col_arr elements will be temporary saved in row_arr:
    // (xxx_arr item = vector of points = polygon)
    vector<double*>* col_arr = new vector<double*>[_nx];
    vector<double*>* row_arr = new vector<double*>[_ny];

    polygon* pol;
    simple_polygon* spol;
    double area = 0; // area calculated directly from given polygons
    double area2 = 0; // area calculated from split polygons, for comparison
    // go through each polygon separately:
    for (size_t i = 0; i < _polygons.size(); ++i) {
        pol = _polygons[i];
        if (pol->get_number_of_points() > 2) {
            area += pol->get_area();
            // split polygon into blocks, save area sizes in _areas:
            area2 += split_polygon(pol, false, col_arr, row_arr);
            // repeat for inner polygons:
            for (std::size_t i = 0; i < pol->get_number_inner_polygons(); ++i){
                spol = pol->get_inner_polygon(i);
                if (spol->get_number_of_points() > 2) {
                    // split polygon into blocks, update area sizes in _areas:
                    area2 -= split_polygon(spol, true, col_arr, row_arr);
                }
            }
        }
        delete pol;
        _polygons[i] = 0;
    }
    _polygons.clear();

    // Make sure the arr items are empty:
    for (int i = 0; i < _nx; ++i) {
        if (!col_arr[i].empty()) {
            // Not empty yet. Delete leftover points to prevent memory leaks:
            for (size_t j = 0; j < col_arr[i].size(); ++j) {
                delete[] col_arr[i][j];
            }
            col_arr[i].clear();
        }
    }
    delete[] col_arr;
    for (int i = 0; i < _ny; ++i) {
        if (!row_arr[i].empty()) {
            // Not empty yet. Delete leftover points to prevent memory leaks:
            for (size_t j = 0; j < row_arr[i].size(); ++j) {
                delete[] row_arr[i][j];
            }
            row_arr[i].clear();
        }
    }
    delete[] row_arr;
    // area is real area, but must be relative area occupied per block
    // like area2, so normalize with block area:
    double block_area = _block_size * _block_size;
    // hardcoded epsilon, but it's okay, since the areas are normalized:
    if (fabs(area / block_area - area2) > 1e-6) {
        master_printf("input polygons area: %f\n", area);
        master_printf("area after splitting: %f\n", area2 * block_area);
        master_printf("WARNING, total polygon area differs after splitting, which could lead to errors in defined structure.\n");
        master_printf("Please provide only simple, non self-intersecting polygons\n");
    }

    // Try to find overlapping polygons and print a warning:
    bool overlap_message_printed = false;
    for (int i = 0; i < _nx; ++i){
        for (int j = 0; j < _ny; ++j) {
            if ((_areas[i][j] < 0) || (_areas[i][j] > 1)) {
                if (!overlap_message_printed) {
                    // print the message only once (per material)
                    master_printf("WARNING: There are overlapping polygons. Structure will be inexact.\n");
                    overlap_message_printed = true;
                }
                if (_areas[i][j] < 0)
                    _areas[i][j] = 0.0;
                else
                    _areas[i][j] = 1.0;
            }
        }
    }
    splitting_done = true;
};

double polygonal_material::split_polygon(
    simple_polygon* pol, bool inner_polygon,
    vector<double*>* col_arr, vector<double*>* row_arr)
{
    std::size_t pol_size;
    double pol_area = 0;
    double split_area = 0;

    // split polygon into columns and save in col_arr:
    split_in_one_dimension(
        col_arr, _nx,
        &pol->get_polygon_point(0), pol->get_number_of_points(), 0);

    // now, go through columns (col_arr) separately and split each
    // polygon into rows:
    for (int c = 0; c < _nx; ++c) {
        if (!col_arr[c].empty()){
            // split col_arr[c] in rows and save in row_arr:
            split_in_one_dimension(
                row_arr, _ny,
                col_arr[c].data(), col_arr[c].size(), 1);

            // go through row_arr, get and save areas,
            // delete split polygon points:
            for (int j = 0; j < _ny; ++j){
                pol_size = row_arr[j].size();
                if (pol_size > 3) {
                    // get area:
                    pol_area = get_polygon_area(row_arr[j].data(), pol_size, 2);
                    if (inner_polygon)
                        _areas[c][j] -= pol_area;
                    else
                        _areas[c][j] += pol_area;
                    split_area += pol_area;

                    // delete polygon points:
                    for (std::size_t k = 0; k < pol_size; ++k){
                        delete[] row_arr[j].at(k);
                    }
                    row_arr[j].clear();
                }
            }
        }
    }
    return split_area;
}

void polygonal_material::split_in_one_dimension(
    std::vector< double* >* arr, int arr_size,
    double** matrix, int matrix_len, int axis){
    // The split polygons will be added to arr:
    // (arr item = polygon = vector of points)
    // arr_size is the length of arr (the number of cols/rows to split the polygon into).

    double *pt, *pt2, *ept, *ept2;
    int col, col2;
    double u;
    int other_axis = axis == 0 ? 1 : 0;

    // Make sure the arr items are empty:
    for (int i = 0; i < arr_size; ++i) {
        if (!arr[i].empty()) {
            // Not empty yet. Delete leftover points to prevent memory leaks:
            for (size_t j = 0; j < arr[i].size(); ++j) {
                delete[] arr[i][j];
            }
            arr[i].clear();
        }
    }

    if (matrix_len == 0)
        return;

    // Note: In the following, I will only talk about columns, but think
    // of them as rows if axis=1.

    // Note 2: The resulting split polygons will be scaled and shifted along
    // axis to only have values between 0 and 1. The reason for this is that
    // at the end after both diretions have been split, trivial polygons which
    // fill a whole block will have an area of exactly one and we don't need
    // to care about floating point rounding errors which we would have
    // otherwise.

    // Go through all points. Check inside which column
    // each point is, and add the point to the polygon of
    // this column in arr. If column changes from one
    // point to next, calculate edge point (= intersection of line
    // connecting the two points and border between columns) and add
    // this point to polygons in both columns in arr, in this order:
    // arr[n]: (other points, previous point, edge point);
    // arr[n+1]: (other points, edge point, current point).
    // Algorithm only works if polygon's first point equals last point.
    pt = matrix[0]; //get first point
    // Calculate inside which column pt is:
    col = static_cast<int>(floor(pt[axis] / _block_size));
    if (col >= 0 && col < arr_size) {
        // Add point to polygon in column:
        arr[col].push_back(pt);
        // This point now belongs to arr[col], it should not be
        // deleted when matrix' items are deleted:
        matrix[0] = 0;
    }
//     else
//         message_truncate(pt);

    //go through remaining points:
    for (int i = 1; i < matrix_len; ++i)
    {
        pt2 = matrix[i]; //get next point
        //calculate inside which column pt2 is:
        col2 = static_cast<int>(floor(pt2[axis] / _block_size));
        if (col == col2){
            // pt2 is in same column than previous point;
            // add point to polygon in column :
            if (col2 >= 0 && col2 < arr_size) {
                arr[col2].push_back(pt2);
                // This point now belongs to arr[col2], it should not be
                // deleted when matrix' items are deleted:
                matrix[i] = 0;
            }
//             else
//                 message_truncate(pt2);
        }
        else {
            // pt2 is in different column than previous point ->
            // add edge point to polygons in both columns
            // (and to additional columns in between, if neccessary);

            // calculate slope of line:
            u = (pt2[other_axis] - pt[other_axis]) / (pt2[axis] - pt[axis]);

            if (col2 > col){ // line is going from left to right
                for (int c = col; c < col2; ++c){
                    // calculate edge point:
                    ept = new double[2];
                    ept[axis] = _block_size * (c + 1); //right border
                    ept[other_axis] = pt[other_axis] + (ept[axis] - pt[axis]) * u;
                    ept[axis] = 1; // right border, scaled and shifted to first block
                    if (c >= 0 && c < arr_size) {
                        // add edge_point to right side of current c:
                        arr[c].push_back(ept);
                        if (c + 1 != arr_size) {
                            // add copy of edge_point to left side of next c:
                            ept2 = new double[2];
                            ept2[axis] = 0; //left border
                            ept2[other_axis] = ept[other_axis];
                            arr[c + 1].push_back(ept2);
                        }
                    }
                    else if (c == -1) {
                        // add edge_point to left side of leftmost c:
                        ept[axis] = 0;
                        arr[c + 1].push_back(ept);
                    }
                    else {
                        // both this and next column are outside computational volume
                        delete[] ept;
                    }
                }
            }
            else { // line is going from right to left
                for (int c = col; c > col2; --c){
                    // calculate edge point:
                    ept = new double[2];
                    ept[axis] = _block_size * c; //left border
                    ept[other_axis] = pt[other_axis] + (ept[axis] - pt[axis]) * u;
                    ept[axis] = 0; //left border, scaled and shifted to first block
                    if (c >= 0 && c < arr_size) {
                        // add edge_point to left side of current c:
                        arr[c].push_back(ept);
                        // if previous point was exactly on border, it will not differ
                        if (c != 0) {
                            // add copy of edge_point to right side of next c:
                            ept2 = new double[2];
                            ept2[axis] = 1; // right border
                            ept2[other_axis] = ept[other_axis];
                            arr[c - 1].push_back(ept2);
                        }
                    }
                    else if (c == arr_size) {
                        // add edge_point to right side of rightmost c:
                        ept[axis] = 1;
                        arr[c - 1].push_back(ept);
                    }
                    else {
                        // both this and next column are outside computational volume
                        delete[] ept;
                    }
                }
            }
            //add point2 after last edge point:
            if (col2 >= 0 && col2 < arr_size) {
                // pt2 will be scaled and shifted later, at the end of the
                // next iteration, but add pointer to it now:
                arr[col2].push_back(pt2);
                // This point now belongs to arr[col2], it should not be
                // deleted when matrix' items are deleted:
                matrix[i] = 0;
            }
//             else
//                 message_truncate(pt2);
        }
        // Scale and shift point to first block:
        pt[axis] = pt[axis] / _block_size - col;
        pt = pt2; // update current point
        col = col2; // update current column
    }
    // Scale and shift point to first block:
    pt[axis] = pt[axis] / _block_size - col;

    // Finalize all new split polygons in arr, i.e. they must all be closed:
    for (int i = 0; i < arr_size; ++i) {
        if (arr[i].size() > 0) {
            pt = arr[i].front();
            pt2 = arr[i].back();
            if (pt[0] != pt2[0] || pt[1] != pt2[1]) {
                pt2 = new double[2];
                pt2[0] = pt[0];
                pt2[1] = pt[1];
                arr[i].push_back(pt2);
            }
        }
    }
}

vec polygonal_material::rightanglerotate_2Dvec(
    vec v, int num, bool mirror_x_before, bool mirror_y_before) const
{
    while (num < 0)
        num += 4;
    num %= 4;
    double x = v.x() * (mirror_x_before ? -1 : 1);
    double y = v.y() * (mirror_y_before ? -1 : 1);
    switch (num)
    {
        case 0 : v.set_direction(X, x); v.set_direction(Y, y); break;
        case 1 : v.set_direction(X, -y); v.set_direction(Y, x); break;
        case 2 : v.set_direction(X, -x); v.set_direction(Y, -y); break;
        case 3 : v.set_direction(X, y); v.set_direction(Y, -x); break;
        default : abort("strange error in rightanglerotate_2Dvec");
    }
    return v;
};

vec polygonal_material::calc_grad_case1(const double A1, const double A2) const
/*
 * Calculate the normal of the line intersecting the two 1x1 squares. The line
 * crosses the squares as shown in the diagram, with the constraints below:
 * (Note: the aspect ratio is not right, think of the blocks of being square)
 *
 *     -|--------------|-
 *      |              |   A1: area of upper triangle
 *   y1 |\             |   A2: area of lower trapezoid
 *      |.\            |   y1: height of upper triangle
 *      |..\           |   y2 + y1: height of combined triangle A1 + A2
 *      |A1.\          |   x1: width of upper triangle
 *     -|----\---------|-  x2: width of combined triangle A1 + A2
 *      |...x1\        |
 *      |......\       |   constraints:
 *   y2 |.......\      |   A1 > 0             A2 > 0
 *      |..A2....\     |   0 < x1 <= 0.5      0 < x2 <= 1
 *      |.........\    |   0 < y1 <= 1        y2 == 1
 *     -|----------\---|-
 *                x2
 *
 * The gradient (unnormalized) is: vec(1, x2 - x1)
 *
 *   (I): A1 = 0.5 * x1 * y1      <=>   (I'): y1 = 2*A1 / x1
 *  (II): A2 = 0.5 * (x1 + x2)    <=>  (II'): x2 = 2*A2 - x1
 * (III): x2 / x1 = (1 + y1) / y1 <=> (III'): x1 + x1*y1 = x2*y1
 */
{
    // Substitute y1 and x2 in (III') with (I') and (II'), respectively, to
    // get an equation for x1:
    double x1 = 2.0 * (sqrt(A1 * (A1 + A2)) - A1);
    // and get x2 from (II):
    double x2 = 2.0 * A2 - x1;
    double dy = x2 - x1;
    double a = sqrt(1.0 + dy * dy);
    return is_3D() ? vec(1.0 / a, dy / a, 0.0) : vec(1.0 / a, dy / a);
};

vec polygonal_material::calc_grad_case2(const double A1, const double A2) const
/*
 * Calculate the normal of the line intersecting the two 1x1 squares. The line
 * crosses the squares as shown in the diagram, with the constraints below:
 * (Note: the aspect ratio is not right, think of the blocks of being square)
 *
 *     -|---------|-       A1: area of upper triangle
 *      |         |        A2: area of lower trapezoid
 *      |\        |        y1: height of upper triangle
 *      |.\       |        y2: height of the lower unfilled triangle (the
 *   y1 |..\      |            triangle with area (1.0-A2) )
 *      |A1.\     |        x1: width of upper triangle
 *     -|----\----|-       x2: width of lower trapezoid (at base)
 *      |..x1.\   |
 *      |......\  |y2      constraints:
 *    1 |..A2...\ |        A1 > 0            A2 > 0.5
 *      |........\|        0 < x1 < 1        x2 == 1
 *      |.........|        0 < y1 <= 1       0 < y2 <= 1
 *     -|---------|-
 *          x2
 *
 * The gradient (unnormalized) is: vec(y1 + y2, 1)
 *
 *   (I):      A1 = 0.5 * x1 * y1       <=>   (I'): x1 = 2*A1 / y1
 *  (II):  1 - A2 = 0.5 * y2 * (1 - x1) <=>  (II'): y2 = 2 * (1-A2) / (1-x1)
 * (III): y1 / y2 = x1 / (1 - x1)       <=> (III'): y1 * (1-x1) = y2 * x1
 */
{
    double y1, y2, dx, a;
    if (A1 != 0) {
        // Substitute y2 in (III') with (II'), then substitute x1 with (I') to
        // get an equation for y1:
        y1 = 2.0 * (A1 + sqrt(A1 * (1.0 - A2)));
        // Now substitute x1 in (III') with (I'):
        y2 = y1 * (0.5 * y1 / A1 - 1.0);
    }
    else {
        // In get_2D_gradient() it is already made sure that this can never
        // be called with A1 == 0, but this is just to make double sure we
        // never get a division by zero.
        y1 = 0.0;
        y2 = 2.0 * (1.0 - A2);
    }

    dx = y1 + y2;
    a = sqrt(1.0 + dx * dx);
    return is_3D() ? vec(dx / a, 1.0 / a, 0.0) : vec(dx / a, 1.0 / a);
};

vec polygonal_material::calc_grad_case3(const double A1, const double A2) const
/*
 * Calculate the normal of the line intersecting the two 1x1 squares. The line
 * crosses the squares as shown in the diagram, with the constraints below:
 * (Note: the aspect ratio is not right, think of the blocks of being square)
 *
 *     -|--\-------------|-  A1: area of upper trapezoid
 *      |.x0\            |   A2: area of lower trapezoid
 *      |....\           |   y1: height of upper trapezoid
 *   y1 |.....\          |   y2: height of lower trapezoid
 *      |..A1..\         |   x0: width of upper trapezoid (at top)
 *      |.......\        |   x1: width of upper trapezoid (at base)
 *     -|--------\-------|-  x2: width of lower trapezoid (at base)
 *      |.......x1\      |
 *      |........,,\     |   constraints:
 *   y2 |....A2....,\    |   A1 > 0         A2 > 0
 *      |............\   |   0 < x0 <= 1    0 < x1 <= 1    0 < x2 <= 1
 *      |.............\  |   y1 == 1        y2 == 1
 *     -|--------------\-|-
 *                    x2
 *
 * The gradient (unnormalized) is: vec(1, x2 - x1)
 *
 *   (I): A2 = 0.5 * (x1 + x2)          <=>   (I'): x2 = 2 * A2 - x1
 *  (II): A1 + A2 = (x0 + x2) = 2 * x1  <=>  (II'): x1 = (A1 + A2) / 2
 *
 *          (I')             (II')
 *  x2 - x1 ==== 2*A2 - 2*x1 ==== 2*A2 - A1 - A2 = A2 - A1
 */
{
    double dy = A2 - A1;
    double a = sqrt(1.0 + dy * dy);
    return is_3D() ? vec(1.0 / a, dy / a, 0.0) : vec(1.0 / a, dy / a);
};

vec polygonal_material::get_2D_gradient(double** areas4x4) const
/*
 * These are the three fundamentally different cases how a straight line might
 * cross two blocks. Every other case can be transformed to one of these by
 * (90-degree-)rotating and/or flipping the problem and transforming the
 * resulting gradient back afterwards. (Note: The aspect ratios in these
 * scetches are not the same and not correct. Think of all blocks being exactly
 * square). The crucial part here is on which side the line crosses a block's
 * outer boundary:
 *
 *   case 1:                case 2:           case 3:
 *  -|--------------|-     -|---------|-     -|--\-------------|-
 *   |              |       |         |       |...\            |
 *   |\             |       |\        |       |....\           |
 *   |.\            |       |.\       |       |.....\          |
 *   |..\           |       |..\      |       |..A1..\         |
 *   |A1.\          |       |A1.\     |       |.......\        |
 *  -|----\---------|-     -|----\----|-     -|--------\-------|-
 *   |.....\        |       |.....\   |       |.........\      |
 *   |......\       |       |......\  |       |........,,\     |
 *   |..A2...\      |       |..A2...\ |       |....A2....,\    |
 *   |........\     |       |........\|       |............\   |
 *   |.........\    |       |.........|       |.............\  |
 *  -|----------\---|-     -|---------|-     -|--------------\-|-
 *
 * Given only the areas A1 and A2, the exact normal vector on the surface can
 * be calculated with calc_grad_case1(A1, A2) to calc_grad_case3(A1, A2),
 * depending on the case.
 *
 * In this method, we will start with the following setup of interesting
 * squares amongst the 4x4 input matrix:
 *
 *  ^ y
 *  |
 *  |
 *   ----> x
 *
 *  |-------------|-------------|-------------|-------------|
 *  |             |             |             |             |
 *  |             |             |             |             |
 *  |             | above_left  | above_right |             |
 *  |             |             |             |             |
 *  |             |             |             |             |
 *  |-------------|-------------|-------------|-------------|
 *  |             |             |             |             |
 *  |             |             |             |             |
 *  |             |     A12     |     A22     |  right_top  |
 *  |             |             |             |             |
 *  |             |             |             |             |
 *  |-------------|-------------|-------------|-------------|
 *  |             |             |             |             |
 *  |             |             |             |             |
 *  | left_bottom |     A11     |     A21     |             |
 *  |             |             |             |             |
 *  |             |             |             |             |
 *  |-------------|-------------|-------------|-------------|
 *  |             |             |             |             |
 *  |             |             |             |             |
 *  |             | below_left  | below_right |             |
 *  |             |             |             |             |
 *  |             |             |             |             |
 *  |-------------|-------------|-------------|-------------|
 *
 * Then, by checking which of the inner squares (A11 - A22) have an interface
 * (they have one when their area is > 0 and < 1) and also checking the area of
 * their neighboring squares (that's where we need the outer squares), we can
 * find out which of the 3 cases and which 2 squares to use for the actual
 * calculation of the normalized gradient.
 *
 * Since there are many possibilities, first check all possible cases
 * that are unique for one orientation, then rotate the whole problem by 90
 * degrees and try again (rotate the resulting vector back 90 degrees
 * accordingly if it could be calculated).
 *
 */
{
    double A11, A12, A21, A22, above_left, below_left, left_bottom;
    double above_right, below_right, right_top;

    for (int rotate = 0; rotate < 4; rotate++)
    {
        switch (rotate)
        {
            case 0:
                A11 = areas4x4[1][1];
                A12 = areas4x4[1][2];
                A21 = areas4x4[2][1];
                A22 = areas4x4[2][2];
                above_left = areas4x4[1][3];
                below_left = areas4x4[1][0];
                left_bottom = areas4x4[0][1];
                above_right = areas4x4[2][3];
                below_right = areas4x4[2][0];
                right_top = areas4x4[3][2];
                break;
            case 1: // rotate blocks in anticlockwise-direction:
                A11 = areas4x4[2][1];
                A12 = areas4x4[1][1];
                A21 = areas4x4[2][2];
                A22 = areas4x4[1][2];
                above_left = areas4x4[0][1];
                below_left = areas4x4[3][1];
                left_bottom = areas4x4[2][0];
                above_right = areas4x4[0][2];
                below_right = areas4x4[3][2];
                right_top = areas4x4[1][3];
                break;
            case 2: // rotate blocks in anticlockwise-direction:
                A11 = areas4x4[2][2];
                A12 = areas4x4[2][1];
                A21 = areas4x4[1][2];
                A22 = areas4x4[1][1];
                above_left = areas4x4[2][0];
                below_left = areas4x4[2][3];
                left_bottom = areas4x4[3][2];
                above_right = areas4x4[1][0];
                below_right = areas4x4[1][3];
                right_top = areas4x4[0][1];
                break;
            default: // rotate == 3: // rotate blocks in anticlockwise-direction:
                A11 = areas4x4[1][2];
                A12 = areas4x4[2][2];
                A21 = areas4x4[1][1];
                A22 = areas4x4[2][1];
                above_left = areas4x4[3][2];
                below_left = areas4x4[0][2];
                left_bottom = areas4x4[1][3];
                above_right = areas4x4[3][1];
                below_right = areas4x4[0][1];
                right_top = areas4x4[2][0];
                break;
        }

        // check case:
        if ((A11 == 0 && A12 == 0 && A21 == 0 && A22 == 0) ||
            (A11 == 1 && A12 == 1 && A21 == 1 && A22 == 1)) {
            // no border in inner 4 blocks
            return vec(ndim(D2 + int(is_3D())), 0.0);
        }
        if ((A11 == 0 && A22 == 0 && A12 != 0 && A21 != 0) ||
            (A12 == 0 && A21 == 0 && A11 != 0 && A22 != 0) ||
            (A11 != 1 && A22 != 1 && A12 == 1 && A21 == 1) ||
            (A12 != 1 && A21 != 1 && A11 == 1 && A22 == 1)) {
            // These cases can only mean there are multiple separate polygons
            // inside the 2x2 center block; i.e. the resolution is way too bad.
            return vec(ndim(D2 + int(is_3D())), 0.0);
        }
        if (A12 == 0 && A21 == 0 && A22 == 0) // only A11 != 0
        {
            if (left_bottom < below_left) {
                return rightanglerotate_2Dvec(
                    calc_grad_case2(A11, below_left),
                    rotate
                );
            }
            else if (left_bottom == 0) {
                // which means below_left and left_bottom are zero,
                // and we have an isolated polygon only inside A11:
                return vec(ndim(D2 + int(is_3D())), 0.0);
            }
            else {
                return rightanglerotate_2Dvec(
                    calc_grad_case2(A11, left_bottom),
                    rotate - 1, true, false
                );
            }
        }
        if (A11 != 0 && A12 != 0)
        {
            if (A21 == 0 && A22 == 0)
            {
                if ((above_left != 0 && below_left != 0) ||
                    (above_left == 0 && below_left == 0)) {
                    // First condition: solution is exact if no corner;
                    // Second one means there is a corner somewhere,
                    // then this is just approx. solution
                    return rightanglerotate_2Dvec(
                        calc_grad_case3(A12, A11), rotate);
                }
                else if (above_left == 0) { // below_left != 0
                    return rightanglerotate_2Dvec(
                        calc_grad_case1(A12, A11),
                        rotate
                    );
                }
                else { // above_left != 0; below_left == 0
                    return rightanglerotate_2Dvec(
                        calc_grad_case1(A11, A12),
                        rotate, false, true
                    );
                }
            }
            else if (A21 != 0 && A22 == 0) // only A22 == 0
            {
                if (A12 < A21) {
                    return rightanglerotate_2Dvec(
                        calc_grad_case2(A12, A11),
                        rotate
                    );
                }
                else {
                    return rightanglerotate_2Dvec(
                        calc_grad_case2(A21, A11),
                        rotate - 1, true, false
                    );
                }
            }
            else if (A21 != 0 && A22 != 0) // all != 0
            {
                if (A11 == 1 && A12 != 1 && A21 != 1 && A22 != 1)
                {
                    if (A12 < A21) {
                        return rightanglerotate_2Dvec(
                            calc_grad_case2(A22, A21),
                            rotate
                        );
                    }
                    else {
                        return rightanglerotate_2Dvec(
                            calc_grad_case2(A22, A12),
                            rotate - 1, true, false
                        );
                    }
                }
                else if (A11 == 1 && A12 == 1 && A21 != 1 && A22 != 1)
                {
                    if ((above_right != 1 && below_right != 1) ||
                        (above_right == 1 && below_right == 1)) {
                        // First condition: solution is exact if no corner;
                        // Last one means there is a corner somewhere,
                        // then this is just approx. solution
                        return rightanglerotate_2Dvec(
                            calc_grad_case3(A22, A21), rotate);
                    }
                    else if (above_right != 1) { // below_right == 1
                        return rightanglerotate_2Dvec(
                            calc_grad_case1(1.0 - A21, 1.0 - A22),
                            rotate
                        );
                    }
                    else { // above_right == 1; below_right != 1
                        return rightanglerotate_2Dvec(
                            calc_grad_case1(1.0 - A22, 1.0 - A21),
                            rotate, false, true
                        );
                    }
                }
                else if (A11 == 1 && A12 == 1 && A21 == 1 && A22 != 1)
                {
                    if (above_right < right_top) {
                        return rightanglerotate_2Dvec(
                            //calc_grad_case2(above_right, A22),
                            calc_grad_case2(1.0 - A22, 1.0 - above_right),
                            rotate
                        );
                    }
                    else if (right_top == 1) {
                        // which means right_top and above_right are both one,
                        // and we have an isolated hole only inside A22:
                        return vec(ndim(D2 + int(is_3D())), 0.0);
                    }
                    else {
                        return rightanglerotate_2Dvec(
                            calc_grad_case2(1.0 - A22, 1.0 - right_top),
                            rotate - 1, true, false
                        );
                    }
                }
            } // (A21 != 0 && A22 != 0) // all != 0
        } // (A11 != 0 && A12 != 0)
    } // for (int rotate = 0; rotate < 4; rotate++)

    // If this point is reached, it can only mean there is a corner somewhere
    // or even multiple polygons in one block. To fix this, the resolution must
    // be increased. For now, an approximate solution based on the Sobel
    // operator (expanded for 4x4 pixel, to be symmetric around the interesting
    // center point of v) for edge detection is calculated:

    // Define the two kernels of the 'expanded' Sobel operator:
    int G_avg_weights[4] = {1, 2, 2, 1};
    int G_diff[4] = {+1, 0, 0, -1};
    double Gx = 0;
    double Gy = 0;
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j){
            Gx += G_avg_weights[j] * G_diff[i] * areas4x4[i][j];
            Gy += G_avg_weights[i] * G_diff[j] * areas4x4[i][j];
        }
        double a = sqrt(Gx * Gx + Gy * Gy);
    if (a == 0.0)
        return is_3D() ? vec(0.0, 0.0, 0.0) : vec(0.0, 0.0);
    else
        return is_3D() ? vec(Gx / a, Gy / a, 0.0) : vec(Gx / a, Gy / a);
};

vec polygonal_material::get_gradient_from_sphere_quad(
    const ivec &pixel_center, const vec &gradient2D)
{
    /**
     * We want to do the same than in material_function::normal_vector:
     * That is, integrate chi1*r over the unit sphere using the
     * Lebedev quadrature scheme.
     *
     * Before we can do that, we need to know chi at arbitrary points inside
     * the pixel: From the 2D-normal and the area occupied by a polygon in the
     * pixel (=2x2 ivec-pixels), a straight interface can be (re)constructed
     * and with this (and the stack in 3D) we can quickly determine chi at any
     * point in the pixel.
     *
     * Note: I could just have kept the original polygons before splitting
     * (after splitting, they will take up too much memory) and use a
     * winding number search to find out whether a specific point is inside a
     * polygon, but that takes longer and I need the areas and gradients
     * anyways.
     */
    if (!splitting_done)
        update_splitting();

    double gx = gradient2D.x();
    double gy = gradient2D.y();
    int rot_num = 0;
    int i2, j2;
    double slope = 1;
    double y0 = 0;

    // sum total area in pixel:
    double total_area = 0;
    for (int i = pixel_center.x() - 1; i <= pixel_center.x(); ++i)
    {
        i2 = i;
        if (i2 < 0) i2 += _nx;
        else if (i2 >= _nx) i2 -= _nx;
        for (int j = pixel_center.y() - 1; j <= pixel_center.y(); ++j)
        {
            j2 = j;
            if (j2 < 0) j2 += _ny;
            else if (j2 >= _ny) j2 -= _ny;
            total_area += _areas[i2][j2];
        }
    }

    // If more than half the area is filled, invert and rotate the whole
    // problem by 180 degrees (then it's easier to solve later on):
    bool invert = total_area > 2.0;
    if (invert) {
        total_area = 4.0 - total_area;
        rot_num = 2;
    }

    bool trivial3D = false;
    /**
     * Check into which (rough) direction the normal points, then rotate
     * everything (the slope of the interface and later the point of interest)
     * such that we can work with a normal pointing roughly upwards
     * (slope between -1 and 1):
     */
    if ((gy > 0) && (fabs(gx) <= gy)) {
        // normal points roughly upwards
        rot_num += 0;
        slope = -gx / gy;
    }
    else if ((gx > 0) && (fabs(gy) <= gx)) {
        // normal points roughly to the right
        rot_num += 1;
        slope = gy / gx;
    }
    else if ((gy < 0) && (fabs(gx) <= -gy)) {
        // normal points roughly downwards
        rot_num += 2;
        slope = -gx / gy;
    }
    else if ((gx < 0) && (fabs(gy) <= -gx)) {
        // normal points roughly to the left
        rot_num += 3;
        slope = gy / gx;
    }
    else {
        // There is no polygon interface here: ((gx == 0) && (gy == 0))
        if (!is_3D()) {
            return zero_vec(gradient2D.dim);
        }
        else if (!invert) {
            // total_area <= 2.0, i.e. we are (mostly) outside any polygons:
            return zero_vec(gradient2D.dim);
        }
        else {
            // we are completely inside the material
            trivial3D = true;
            rot_num = 0;
        }
    }

    if (!trivial3D) {
        // The interface is a line with the equation y = slope*x + y0. Find y0:
        if (slope < 0 && total_area <= -2 * slope) {
            // total_area is a triangle at the lower left edge of the pixel
            y0 = slope - 1 + sqrt(-2 * slope * total_area);
        }
        else if (slope > 0 && total_area <= 2 * slope) {
            // total_area is a triangle at the lower right edge of the pixel
            y0 = -slope - 1 + sqrt(2 * slope * total_area);
        }
        else {
            // total_area is a trapezoid at the bottom of the pixel
            y0 = total_area * 0.5 - 1;
        }
    }

    // Increased radius so the pixel corners are also included:
    const double radius = sqrt(gradient2D.dim + 1);
    vec zerovec(zero_vec(gradient2D.dim));
    vec result(zerovec);
    for (int i = 0;
         i < num_sphere_quad[number_of_directions(gradient2D.dim) - 1];
         ++i)
        {
            double weight, chi1p1;
            vec pt = sphere_pt(zerovec, radius, i, weight);

            bool on_edge = fabs(pt.y() - slope * pt.x() + y0) < 0.001;

            bool inside = (
                on_edge ||
                trivial3D ||
                (!invert && pt.y() <= slope * pt.x() + y0) ||
                (invert && pt.y() > slope * pt.x() + y0));

            if (inside)
            {
                if (is_3D()) {
                    chi1p1 = _material_stack->chi1p1(
                        (pt.z() + pixel_center.z()) * _block_size);
                }
                else
                   chi1p1 = _epsilon;
            }
            else {
                chi1p1 = 1;
            }

            /**
             * Note1: chi1 = epsilon-1 instead of epsilon is used so that all
             * gradients returned from different polygonal_materials in same
             * pixel can simply be added.
             *
             * Note2: Contrary to the implementation in
             * material_function::normal_vector, here we _substract_
             * the individual weighted points. The reason is that
             * the gradients in polygonal_material point away from
             * the polygons, and therefore also away from areas of
             * higher eps (background is always air), whereas the
             * gradients in material_function::normal_vector point
             * towards areas of higher eps. While the sign doesn't
             * matter at the end in eff_chi1inv_row, we need to stay
             * consistent within one material_function.
             */
            if (on_edge)
                result -= pt * weight * (chi1p1 - 1) * 0.5;
            else
                result -= pt * weight * (chi1p1 - 1);
        }

    // Counter the rotation applied to the slope:
    result = rightanglerotate_2Dvec(result, -rot_num) / radius;

    // Now, we actually only needed the z component, since the input gradient2D
    // is supposed to be more accurate than the lateral part of the result. So
    // we fix the result with the better lateral vector:
    double a = gx * gx + gy * gy;
    double b = result.x() * result.x() + result.y() * result.y();
    double c = sqrt(b/a);
    result.set_direction(X, gx * c);
    result.set_direction(Y, gy * c);
    return result;
}

material_function_for_polygons::material_function_for_polygons(
    const grid_volume &thegv)
{
    if ((thegv.dim != D2) && (thegv.dim != D3))
        abort("material_function_for_polygons only available for 2D and 3D.\n");
    // Number of blocks to split polygons into:
    // (One block should have halve the size of one pixel to allow
    // calculating averaged eps for different components (half pixel
    // apart due to yee lattice). The smooting diameter should stay
    // 1 pixel (Meep default) or be an integer multiple thereof,
    // in which case material_function_for_polygons::eff_chi1inv_row()
    // must be changed though.)
    nx = thegv.nx() * 2;
    ny = thegv.ny() * 2;
    nz = thegv.nz() * 2;

    block_size = thegv.inva * 0.5;
};

material_function_for_polygons::~material_function_for_polygons() {
    for (vector<polygonal_material*>::const_iterator it = _polymats.begin();
            it != _polymats.end();
    ++it)
        {
            delete (*it); //delete polygonal_material object
        }
        _polymats.clear(); //delete all pointers
}

unsigned int material_function_for_polygons::add_material_stack(
    double* layer_thicknesses, int dim1, double* layer_epsilons, int dim2)
{
    if (!is_3D()) {
        abort("Specifying material stack is only available in 3D.\n%s",
            "For 2D, please specify epsilon in add_polygon().\n");
    }
    if (dim1 != dim2)
        abort("layer_thicknesses and layer_epsilons must have same length, "
              "i. e. the number of layers in the material stack.");
    material_stack mat_stack(layer_thicknesses, layer_epsilons, dim1);
    _polymats.push_back(
        new polygonal_material(
            &mat_stack, block_size, nx, ny, nz));
    return _polymats.size() - 1;
};

// TODO? all add_polygon: check if polygon overlaps other already added polygons
void material_function_for_polygons::add_polygon(const polygon &pol, unsigned int mat_ID){
    if (!is_3D()) {
        abort("Specifying material stack is only available in 3D.\n%s",
            "For 2D, please specify epsilon in add_polygon().\n");
    }
    if (mat_ID >= _polymats.size())
        abort("Wrong mat_ID specified in add_polygon. Please add material stacks before adding polygons.");
    _polymats[mat_ID]->add_polygon(pol);
};

void material_function_for_polygons::add_polygon(const polygon &pol, double eps){
    if (is_3D()) {
        abort("Specifying polygon with epsilon value is only available in 2D.\n%s",
            "For 3D, please specify material_stack ID instead.\n");
    }
    // try to use id of previously added eps:
    int id = add_epsilon(eps);
    _polymats[id]->add_polygon(pol);
};

vec material_function_for_polygons::normal_vector(const ivec &center)
{
    vec result = vec(is_3D()? D3 : D2, 0.0);
    for (
        vector<polygonal_material*>::const_iterator it = _polymats.begin();
        it != _polymats.end();
        ++it)
    {
        if (!(*it)->is_trivial(
            center - one_ivec(center.dim),
            center + one_ivec(center.dim)))
        {
            /**
             * The gradients returned from polygonal_material::normal_vector
             * are normalized such that neighboring materials' gradients in
             * same pixel can simply be added:
             */
            result += (*it)->normal_vector(center);
        }
        else
        if (is_3D()) {
            /**
             * If the material is trivial in pixel, i.e. the pixel is completely
             * empty or completely filled, the transversal gradient (gradient
             * projected onto x-y plane) is (0, 0).
             * In 3D, we must find out the gradient's z-component, which is
             * only non-zero if there is an interface in z, which in turn can
             * only be if the pixel is filled: */
            if ((*it)->get_area(center) == 1) {
                double grad_z = (*it)->get_material_stack()->normal(
                    (center.z() - 1) * block_size,
                    (center.z() + 1) * block_size);
                if (abs(grad_z) > 1e-8)
                {
                    /**
                     * There is an interface in z direction. We can directly
                     * return it, since the structure won't make sense if
                     * there are any other interfaces in the same pixel:
                     */
                    return vec(0, 0, grad_z);
                }
            }
        }
    }
    return result;
}

double material_function_for_polygons::mean_eps(const ivec &center){
    // start with complete area (made up of 4 subpixel), multiplied by 1 (air):
    double result = 4.0;
    for (
        vector<polygonal_material*>::const_iterator it = _polymats.begin();
        it != _polymats.end();
        ++it)
    {
        // epsilon for this block (or stack of blocks in 3D, averaged along z):
        double eps = (
            is_3D() ?
            (*it)->get_material_stack()->mean_eps
            (
                (center.z() - 1) * block_size,
                (center.z() + 1) * block_size
            ) :
            (*it)->get_material_epsilon()
        );

        // add up all 4 areas (in xy):
        for (int i = center.x() - 1; i <  center.x() + 1; ++i)
            for (int j = center.y() - 1; j < center.y() + 1; ++j) {
                // add epsilon, weighted by area,
                // minus 1 (= air: already added in beginning):
                result += (*it)->get_area(i, j) * (eps - 1.0);
            }
    }
    return result / 4.0;
}

double material_function_for_polygons::mean_inveps(const ivec &center){
    // start with complete area (made up of 4 subpixel), multiplied by 1 (air):
    double result = 4.0;
    for (
        vector<polygonal_material*>::const_iterator it = _polymats.begin();
        it != _polymats.end();
        ++it)
    {
        // inverse epsilon for this block (or stack of blocks in 3D, averaged along z):
        double inveps = (
            is_3D() ?
            (*it)->get_material_stack()->mean_inveps
            (
                (center.z() - 1) * block_size,
                (center.z() + 1) * block_size
            ) :
            1.0 / (*it)->get_material_epsilon()
        );

        // add up all 4 areas (in xy):
        for (int i = center.x() - 1; i <  center.x() + 1; ++i)
            for (int j = center.y() - 1; j < center.y() + 1; ++j) {
                // add 1/epsilon, weighted by area,
                // minus 1 (= air: already added in beginning):
                result += (*it)->get_area(i, j) * (inveps - 1.0);
            }
    }
    return result / 4.0;
}

ivec material_function_for_polygons::round_to_block_edge(const vec &v) {
    ivec iv(v.dim, 0.0);
    LOOP_OVER_DIRECTIONS(v.dim, d)
    {
        iv.set_direction(
            d, my_round(v.in_direction(d) / block_size));
    }
    return iv;
}

unsigned int material_function_for_polygons::add_epsilon(double eps) {
    if (!is_3D()){ //2D
        for (std::size_t i = 0; i < _polymats.size(); ++i){
            if (_polymats[i]->get_material_epsilon() == eps) {
                //epsilon value already added before. Just return ID.
                return i;
            }
        }
        // add new polygon_for_single_material with epsilon value and return
        // new ID:
        _polymats.push_back(new polygonal_material(
            eps, block_size, nx, ny));
        return _polymats.size() - 1;
    }
    else { //3D
        abort("Specifying polygon with epsilon value is only available in 2D.\n%s",
              "For 3D, please specify material_stack ID instead.\n");
    }
};

#define UNUSED(x) (void) x // silence compiler warnings
void material_function_for_polygons::eff_chi1inv_row(
    component c, double chi1inv_row[3],
    const volume &v, double tol, int maxeval)
{
    UNUSED(tol);

    ivec center_ivec(round_to_block_edge(v.center()));
    if (!maxeval) {
        trivial:
        chi1inv_row[0] = chi1inv_row[1] = chi1inv_row[2] = 0.0;
        chi1inv_row[component_direction(c) % 3] = 1.0 / mean_eps(center_ivec);
        return;
    }

    vec gradient(normal_vector(center_ivec));
    if (abs(gradient) < 1e-8) {
        goto trivial;
    }

    double meps = mean_eps(center_ivec);
    double minveps = mean_inveps(center_ivec);

    double n[3] = {0,0,0};
    double nabsinv = 1.0/abs(gradient);
    LOOP_OVER_DIRECTIONS(gradient.dim, k)
    n[k%3] = gradient.in_direction(k) * nabsinv;

    /* get rownum'th row of effective tensor
    *  P * minveps + (I-P) * 1/meps = P * (minveps-1/meps) + I * 1/meps
    *  where I is the identity and P is the projection matrix
    *  P_{ij} = n[i] * n[j]. */
    int rownum = component_direction(c) % 3;
    for (int i=0; i<3; ++i)
        chi1inv_row[i] = n[rownum] * n[i] * (minveps - 1/meps);
    chi1inv_row[rownum] += 1/meps;
}

////////////////////////////////////////////////////////////
//---- End: defs for polygon-based material functions ----//
////////////////////////////////////////////////////////////



} // namespace meep
