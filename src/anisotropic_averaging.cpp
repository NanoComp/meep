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


////////////////////////////////////////////////////////////
//--- Begin: defs for polygon-based material functions ---//
////////////////////////////////////////////////////////////

// polygon defined by a numpy array with dimensions N*2 (containing the points) and defining an area with a certain epsilon
simple_polygon::simple_polygon(double* matrix, int dimX, int dimY) {
    if (dimY != 2) {
        abort("Array with polygon points should have a shape X,2. The second dimension is not 2.");
    }

    //make sure last point equals first point. If not, append copy of first point to end:
    if (dimX > 0 && (matrix[0] != matrix[dimY * (dimX - 1)] || matrix[1] != matrix[dimY * dimX - 1]))
        this->_number_of_points = dimX + 1;
    else
        this->_number_of_points = dimX;

    // copy matrix:
    _polygon_points = new double*[_number_of_points];
    for (int i = 0; i < dimX; i++) {
        _polygon_points[i] = new double[2];
        _polygon_points[i][0] = matrix[dimY * i];
        _polygon_points[i][1] = matrix[dimY * i + 1];
    }

    //if last point does not equal first point -> append copy of first point to end:
    if (this->_number_of_points == dimX + 1) {
        _polygon_points[dimX] = new double[2];
        _polygon_points[dimX][0] = matrix[0];
        _polygon_points[dimX][1] = matrix[1];
    }
};

simple_polygon::simple_polygon(double** matrix, int dimX) {
    // if this constructor is used, care must be taken that matrix is not deleted
    // externally. The memory pointed to by matrix will be deleted with ~polygon()
    // First point of polygon must equal last point.

    //make sure last point equals first point:
    if (dimX > 0 && (matrix[0][0] != matrix[dimX - 1][0] || matrix[0][1] != matrix[dimX - 1][1]))
        abort("Error creating polygon: First point should equal last point.");

    this->_number_of_points = dimX;
    _polygon_points = matrix;
    //the caller should set matrix = 0;
}

simple_polygon::~simple_polygon() {
    for (int i = 0; i < _number_of_points; i++) {
        delete [] _polygon_points[i];
    }
    delete [] _polygon_points;
};

double simple_polygon::get_area() {
    //Returns the absolute area of polygon.
    //For self-intersecting polygons, this will return the difference between clockwise
    //parts and anti-clockwise parts, not the real area.
    int n = _number_of_points;
    if (n < 3)
        return 0;
    //2A = \sum_{i=0}^{n-1}( (x_{i-1} - x_{i+1}) * y_i ) # Gauss's area formula (aka shoelace algorithm)
    double sum = (_polygon_points[n - 1][0] - _polygon_points[1][0]) * _polygon_points[0][1]; // (i = 0)
    for (int i = 1; i < n - 1; ++i)
        sum += (_polygon_points[i - 1][0] - _polygon_points[i + 1][0]) * _polygon_points[i][1];
    sum += (_polygon_points[n - 2][0] - _polygon_points[0][0]) * _polygon_points[n - 1][1]; // (i = n-1)
    return fabs(0.5 * sum);
}

polygon::~polygon() {
    for(vector<simple_polygon*>::iterator p = _inner_polygons.begin(); p != _inner_polygons.end(); ++p) {
        delete (*p);
    }
};

void polygon::add_inner_polygon(double* matrix, int dimX, int dimY) {
    // TODO? check if inner polygon is real inner polygon of parent:
    // first point is inside parent polygon AND inner polygon's edges does not cross parent polygon's edges
    _inner_polygons.push_back(new simple_polygon(matrix, dimX, dimY));
}

double polygon::get_area() {
    double parea = simple_polygon::get_area();
    for(vector<simple_polygon*>::iterator inner_p = _inner_polygons.begin(); inner_p != _inner_polygons.end(); ++inner_p) {
        parea -= (*inner_p)->get_area(); //all inner polygons must be real inner polygons (i.e. completely inside parent polygon)
    }
    return parea;
}

material_stack::material_stack(const double* material_heights,
                               const double* epsilon_values, const int number_of_layers)
{
    init_material_stack(material_heights, epsilon_values, number_of_layers);
};

material_stack::material_stack(const double* material_heights, const int mh_dim,
                               const double* epsilon_values, const int ev_dim)
{
    if (mh_dim != ev_dim)
        abort("material_heights and epsilon_values must have same length, \
        i. e. the number of layers in the material stack.");
    init_material_stack(material_heights, epsilon_values, mh_dim);
}

void material_stack::init_material_stack(const double* material_heights,
                                         const double* epsilon_values, const int number_of_layers)
{
    if (number_of_layers == 0)
        abort("Please specify at least one layer in material stack.");

    // make copy of arrays:
    // It is really important that the user supplies the right length of the two
    // arrays with number_of_layers, else some undefined memory outside of the
    // arrays could be accessed. Unfortunately, this can't be checked here.
    this->_material_heights = new double[number_of_layers];
    this->_epsilon_values = new double[number_of_layers];
    for (int i = 0; i < number_of_layers; ++i) {
        this->_material_heights[i] = material_heights[i];
        this->_epsilon_values[i] = epsilon_values[i];
    }
    _number_of_layers = number_of_layers;
}

material_stack::~material_stack(){
    delete[] _material_heights;
    delete[] _epsilon_values;
};

bool material_stack::interface_inside_block(
    double center_z, double half_block_size,
    double &lower_layer_eps, double &upper_layer_eps,
    double &percentage_upper_layer) const
    {
        if (_number_of_layers == 1) {
            // no interface
            lower_layer_eps = _epsilon_values[0];
            upper_layer_eps = 0.0;
            percentage_upper_layer = 0.0;
            return false;
        }

        double lower = center_z - half_block_size;
        double higher = center_z + half_block_size;
        double height;
        double z = 0.0;
        for (int i = 0; i < _number_of_layers - 1; i++)
        {
            z += _material_heights[i];
            if (lower < z) {
                lower_layer_eps = _epsilon_values[i];
                if (higher <= z) {
                    // no interface
                    upper_layer_eps = 0.0;
                    percentage_upper_layer = 0.0;
                    return false;
                }
                else {
                    height = higher - z;
                    if (i < _number_of_layers - 2 && // upmost layer may be thinner than resolution
                        height > _material_heights[i + 1])
                        abort("Layer specified in material_stack too thin or resolution too poor.");
                    upper_layer_eps = _epsilon_values[i + 1];
                    percentage_upper_layer = height / half_block_size / 2.0;
                    return true;
                }
            }
        }
        // lower > z, i.e. block is in highest layer or above:
        lower_layer_eps = _epsilon_values[_number_of_layers - 1];
        upper_layer_eps = 0.0;
        percentage_upper_layer = 0.0;
        return false;
    }

    polygons_for_single_material::polygons_for_single_material(
        material_stack * mat_stack, double block_size, int nx, int ny, int nz)
    {
        init_pfsm(mat_stack, 0, block_size, nx, ny, nz);
    }

    polygons_for_single_material::polygons_for_single_material(
        double epsilon, double block_size, int nx, int ny)
    {
        init_pfsm(0, epsilon, block_size, nx, ny, 0);
    }

    void polygons_for_single_material::init_pfsm(material_stack * mat_stack,
                                                    double epsilon, double block_size, int nx, int ny, int nz)
    {
        _material_stack = mat_stack;
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

    polygons_for_single_material::~polygons_for_single_material() {
        //delete polygons list:
        while (!_polygons.empty()){
            delete _polygons.front(); //delete polygon object
            _polygons.pop(); //delete pointer to polygon object
        }

        delete _material_stack;

        // delete _areas:
        for (int i = 0; i < _nx; ++i)
            delete[] _areas[i];
        delete[] _areas;
    };

    polygon* polygons_for_single_material::add_polygon(polygon* pol) {
        if (splitting_done)
            abort("Cannot add polygons after splitting has been done.");
        _polygons.push(pol);
        return pol;
    };

    bool polygons_for_single_material::is_trivial(const ivec &small_ivec, const ivec &big_ivec)
    {
        if (!splitting_done)
            update_splitting();
        int i2, j2;
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
            }
        }
        return true;
    }

    vec polygons_for_single_material::normal_vector(const ivec &center)
    {
        if (!splitting_done)
            update_splitting();
        // Initialize the 4x4 blocks (=2x2 pixels) that will be used as parameter in
        // calc_gradient_4x4:
        double** area_cols = new double*[4];
        bool created_area_cols = false;
        int xi, yi;
        vec result;
        if (center.y() < 2 || center.y() + 1 >= _ny) {
            // If pixel is at outer edge, overlapping the undefined area outside
            // the grid volume, wrap the structure to the opposite side. If a PML
            // boundary is used, there shouldn't be a structure close to the
            // boundary anyways, and if periodic boundary conditions are used, we
            // would want to wrap the structure.
            for (int i = 0; i < 4; ++i)
            {
                xi = center.x() - 2 + i;
                if (xi < 0)
                    xi += _nx;
                else if (xi >= _nx)
                    xi -= _nx;
                area_cols[i] = new double[4];
                created_area_cols = true;
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


            result = calc_gradient_4x4(area_cols);
            if (center == ivec(216, 0, 0)) {
                master_printf("singmat normvec: %f, %f, %f; area cols, wrapping: \n", result.x(), result.y(), result.z() );
                for(int i = 0; i < 4; ++i)
                    master_printf("%f, %f, %f, %f\n", area_cols[i][0], area_cols[i][1], area_cols[i][2], area_cols[i][3]);
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
            result = calc_gradient_4x4(area_cols);
        }
        if (is_3D()){
            double lower_layer_eps;
            double upper_layer_eps;
            double percentage_upper_layer;
            double z;
            if (_material_stack->interface_inside_block(
                center.z() * _block_size, _block_size,
                                                        lower_layer_eps, upper_layer_eps,
                                                        percentage_upper_layer))
            {
                z = upper_layer_eps < lower_layer_eps ? 1.0     : -1.0;
                if (abs(result) < 1e-8)
                    // interface only in z direction
                    return vec(0.0, 0.0, z);
                else
                {
                    // Interface in z and transversal direction;
                    // Need to return approximate vector for corner.
                    // This is the only time the sign of the normal vector is
                    // important, because the transversal vector (from polygons)
                    // and the z vector (from material stack) must be combined,
                    // and the resulting vector will point in a different direction
                    // if transversal and z vector have opposite sign than when
                    // both have same sign.
                    // This is difficult, because a neigboring polygon (on same
                    // block) with different material function (= different
                    // polygons_for_single_material) could have same epsilon on
                    // different height, e.g. going lower and higher up than the
                    // currently regarded point, which would mean the transversal
                    // part of this polygon is actually pointing inward.
                    // The only possible way to solve this is on a higher level,
                    // so this problem is moved to
                    // material_function_for_polygons::normal_vector().
                    // For now, just return the combined vector ignoring possible
                    // neighbors, with vector pointing from high to low eps:


                    // height of material with higher eps:
                    double h = upper_layer_eps < lower_layer_eps ?
                    1 - percentage_upper_layer :
                    percentage_upper_layer;
                    // square for better comparability with area:
                    h *= h;
                    // area of current 2x2 block (= 1 pixel):
                    double area = area_cols[1][1] + area_cols[1][2] +
                    area_cols[2][1] + area_cols[2][2];
                    // scale transversal vector:
                    result = result * (h * (1.0 - area));
                    // add scaled z component:
                    result.set_direction(Z, area * (1.0 - h) * z);
                    // normalize:
                    if (abs(result) != 0.0)
                        result = result / abs(result);
                }
            }
        }
        // delete area_cols[] if arrays have been created earlier:
        if (created_area_cols)
            for (int i = 0; i < 4; ++i)
                delete[] area_cols[i];
            delete[] area_cols;
        return result;
    }

    double polygons_for_single_material::get_area(int i, int j) {
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

    void polygons_for_single_material::message_truncate(double* pt) {
        if (pt[0] < 0 || pt[0] - 1e-8 > _nx * _block_size || //TODO hardcoded epsilon: yuck!
            pt[1] < 0 || pt[1] - 1e-8 > _ny * _block_size)
            // if point is exactly on upper border, this will also be called,
            // but point will not really be truncated, so check bounds again.
            master_printf(//"~");
        "point (%f, %f) outside computational grid (%f, %f) - polygon will be truncated.\n",
                            pt[0], pt[1], _nx * _block_size, _ny * _block_size);
    }

    void polygons_for_single_material::split_in_one_dimension(queue<double*>* arr, int arr_size, queue<double*>* qpol, int axis){
        if (qpol->empty())
            return;
        double *f = qpol->front();
        double *b = qpol->back();
        // check if last point equals first point
        if (f[0] != b[0] || f[1] != b[1]){
            // add copy of first point to end:
            b = new double[2];
            b[0] = f[0];
            b[1] = f[1];
            qpol->push(b);
        }
        std::size_t pol_size = qpol->size();
        double** matrix = new double*[pol_size];
        for (std::size_t i = 0; i < pol_size; ++i) {
            matrix[i] = qpol->front();
            qpol->pop();
        }
        simple_polygon* pol = new simple_polygon(matrix, pol_size);
        matrix = 0; //the memory pointed to by matrix is now owned by pol;
        split_in_one_dimension(arr, arr_size, pol, axis);
        delete pol;
    }

    void polygons_for_single_material::split_in_one_dimension(queue<double*>* arr, int arr_size, simple_polygon* pol, int axis){
        // The split polygons will be added to arr:
        // (arr item = polygon = queue of points)
        // arr_size is the length of arr (the number of cols/rows to split the polygon into).
        // pol will only be split along axis (0 or 1).
        // The array arr must be initialized before call, but the items (the queues)
        // must be empty.

        double *pt, *pt2, *ept, *ept2;
        int col, col2;
        double u;
        int other_axis = axis == 0 ? 1 : 0;

        /* Go through all points. Check inside which column
            *                                                e ach point is, and add the point to the polygon o*f
            *                                                this column in arr. If column changes from one
            *                                                point to next, calculate edge point (= intersection of line
            *                                                connecting the two points and border between columns) and add
            *                                                this point to polygons in both columns in arr:
            *                                                arr[n]: (other points, previous point, edge point);
            *                                                arr[n+1]: (other points, edge point, current point)
            *                                                Algorithm only works if polygon's first point equals last point,
            *                                                but this is always true for simple_polygon class.*/

        pt = pol->get_polygon_point(0); //get first point
        // note: There can't be any empty polygons (with zero points)
        // at this point, since they have been sorted out in
        // add_polygon already.
        col = static_cast<int>(floor(pt[axis] / _block_size)); //TODO: special handling if point is exactly on border:
        // such a point belongs to current col
        // compare to col = ceil(pt) - 1
        // don't add point to neighbouring column - it will happen automatically
        // if next point is inside neighbouring column.
        //calculate inside which column pt is
        if (col >= 0 && col < arr_size) {
            //add point to last polygon in column:
            arr[col].push(pt);
            // steal point (=double* with length 2) from pol,
            // so it will not be deleted together with pol:
            pol->take_ownership_of_point(0);
        }
        else
            message_truncate(pt);

        //go through remaining points:
        for (int i = 1; i < pol->get_number_of_points(); ++i){
            pt2 = pol->get_polygon_point(i); //get next point
            //calculate inside which column pt2 is:
            col2 = static_cast<int>(floor(pt2[axis] / _block_size));//TODO: special handling if point is exactly on border:
            // such a point belongs to current col
            // compare to col = ceil(pt) - 1
            // don't add point to neighbouring column - it will happen automatically
            // if next point is inside neighbouring column.
            if (col == col2){
                // pt2 is in same column than previous point;
                // add point to last polygon in column :
                if (col2 >= 0 && col2 < arr_size) {
                    arr[col2].push(pt2);
                    // steal point (=double* with length 2) from pol,
                    // so it will not be deleted together with pol:
                    pol->take_ownership_of_point(i);
                }
                else
                    message_truncate(pt2);
                pt = pt2;
            }
            else {
                // pt2 is in different column than previous point ->
                // add edge point to last polygon in both columns
                // (and to additional columns in between, if neccessary);

                // calculate slope of line:
                u = (pt2[other_axis] - pt[other_axis]) / (pt2[axis] - pt[axis]);

                if (col2 > col){ // line is going from left to right
                    for (int c = col; c < col2; ++c){
                        // calculate edge point:
                        ept = new double[2];
                        ept[axis] = _block_size * (c + 1); //right border
                        ept[other_axis] = pt[other_axis] + (ept[axis] - pt[axis]) * u;
                        if (c >= 0 && c < arr_size) {
                            // add edge_point to right side of current c:
                            arr[c].push(ept); // TODO: only add if ept differs from previous point
                            // if previous point was exactly on border, it will not differ
                            if (c + 1 != arr_size) {
                                // add copy of edge_point to left side of next c:
                                ept2 = new double[2];
                                ept2[0] = ept[0];
                                ept2[1] = ept[1];
                                arr[c + 1].push(ept2);
                            }
                        }
                        else if (c == -1)
                            // add edge_point to left side of leftmost c:
                            arr[c + 1].push(ept);
                        else
                            // both this and next column are outside computational volume
                            delete[] ept;
                    }
                }
                else { // line is going from right to left
                    for (int c = col; c > col2; --c){
                        // calculate edge point:
                        ept = new double[2];
                        ept[axis] = _block_size * c; //left border
                        ept[other_axis] = pt[other_axis] + (ept[axis] - pt[axis]) * u;
                        if (c >= 0 && c < arr_size) {
                            // add edge_point to left side of current c:
                            arr[c].push(ept);  // TODO: only add if ept differs from previous point
                            // if previous point was exactly on border, it will not differ
                            if (c != 0) {
                                // add copy of edge_point to right side of next c:
                                ept2 = new double[2];
                                ept2[0] = ept[0];
                                ept2[1] = ept[1];
                                arr[c - 1].push(ept2);
                            }
                        }
                        else if (c == arr_size)
                            // add edge_point to right side of rightmost c:
                            arr[c - 1].push(ept);
                        else
                            // both this and next column are outside computational volume
                            delete[] ept;
                    }
                }
                //add point2 after last edge point:
                if (col2 >= 0 && col2 < arr_size) {
                    arr[col2].push(pt2);
                    // steal point (=double* with length 2) from pol,
                    // so it will not be deleted together with pol:
                    pol->take_ownership_of_point(i);
                }
                else
                    message_truncate(pt2);
                pt = pt2; // update current point
                col = col2; // update current column
            }
        }
        // TODO?: remove ghost lines & points: if polpoints enter and leave column again on the same side, start new polygon in that column.
        // at end, go through first added point, if in same column than last added point - polygons should be merged.
    }

    double polygons_for_single_material::split_polygon(simple_polygon* pol,
                                                        bool inner_polygon, queue<double*>* col_arr, queue<double*>* row_arr)
    {
        std::size_t pol_size;
        double** matrix;
        simple_polygon *spol;
        double pol_area = 0;
        double split_area = 0;

        // split polygon into columns and save in col_arr:
        split_in_one_dimension(col_arr, _nx, pol, 0);
        // now, go through columns (col_arr) separately and split each
        // polygon into rows:
        for (int c = 0; c < _nx; ++c) {
            if (!col_arr[c].empty()){
                // split col_arr[c] in rows and save in row_arr:
                split_in_one_dimension(row_arr, _ny, &col_arr[c], 1);
                // after the call to split_in_one_dimension, all polygons in
                // col_arr[c] are empty.

                // go through row_arr, make polygons,
                // get and save areas, delete polygons again:
                for (int j = 0; j < _ny; ++j){
                    pol_size = row_arr[j].size();
                    if (pol_size > 0) { //don't add empty polygons
                        double *f = row_arr[j].front();
                        double *b = row_arr[j].back();
                        // check if last point equals first point
                        if (f[0] != b[0] || f[1] != b[1]){
                            // add copy of first point to end:
                            b = new double[2];
                            b[0] = f[0];
                            b[1] = f[1];
                            pol_size ++;
                        }
                        else
                            b = 0;
                        matrix = new double*[pol_size];
                        for (std::size_t k = 0; k < pol_size - int(b!=0); ++k){
                            //get first point in queue (double[2] point):
                            matrix[k] = row_arr[j].front();
                            //pop this point:
                            row_arr[j].pop();
                        }
                        if (b)
                            matrix[pol_size - 1] = b;

                        // now, row_arr should be empty
                        spol = new simple_polygon(matrix, pol_size);
                        //the memory pointed to by matrix is now owned by spol:
                        matrix = 0;
                        pol_area = spol->get_area();
                        if (inner_polygon)
                            _areas[c][j] -= pol_area;
                        else
                            _areas[c][j] += pol_area;
                        split_area += pol_area;
                        delete spol;
                    }
                }
            }
        }
        return split_area;
    }

    void polygons_for_single_material::update_splitting() {
        /*Splits all _polygons in blocks with size _blocksize x _blocksize.
         * 
         *  This is a quick and efficient algorithm, which can be us*ed for both simple
         *  and for selfintersecting polygons, both convex and concave.
         *
         *  It might produce ghost lines in the case of concave polygons (i.e. concave
         *  polygons will never be split in more than two polygons by one splitting border,
         *  with ghost lines connecting otherwise seperate polygons). The additional area
         *  created by these ghost lines is zero, so this will not be an issue for most
         *  applications.*/


        // first, split polygons into columns, then afterwards split these columns also in the rows.
        // If it is not done in this way, blocks completely enclosed by a single polygon will be left out.

        // the column-wise split polygons will be temporary saved in col_arr,
        // the row-wise split col_arr elements will be temporary saved in split_arr:
        // (column item = polygon = queue of points)
        queue<double*>* col_arr = new queue<double*>[_nx];
        queue<double*>* row_arr = new queue<double*>[_ny];

        polygon* pol;
        simple_polygon* spol;
        double area = 0; // area calculated directly from given polygons
        double area2 = 0; // area calculated from split polygons, for comparison
        // go through each polygon separately:
        while (!_polygons.empty()){
            pol = _polygons.front();
            area += pol->get_area();
            // split polygon into blocks, save area sizes in _areas:
            area2 += split_polygon(pol, false, col_arr, row_arr);
            //col_arr and row_arr are empty after split_polygon()
            // repeat for inner polygons:
            for (int i = 0; i < pol->get_number_inner_polygons(); ++i){
                spol = pol->get_inner_polygon(i);
                // split polygon into blocks, save area sizes in _areas:
                area2 -= split_polygon(spol, true, col_arr, row_arr);
                //col_arr and row_arr are empty after split_polygon()
            }
            delete pol;
            _polygons.pop();
        }

        delete[] col_arr; //delete array, free memory
        delete[] row_arr; //delete array, free memory

        if (fabs(area - area2) > 1e-6) { //TODO: hardcoded epsilon: yuck!
            master_printf("input polygons area: %f\n", area);
            master_printf("area after splitting: %f\n", area2);
            master_printf("WARNING, total polygon area differs after splitting, which could lead to errors in defined structure.\n");
            master_printf("Please provide only simple, non self-intersecting polygons\n");
        }

        bool overlap_message_printed = false;
        // areas must be relative area occupied per block, so normalize with
        // block area:
        double block_area = _block_size * _block_size;
        for (int i = 0; i < _nx; ++i){
            for (int j = 0; j < _ny; ++j) {
                _areas[i][j] /= block_area;
                // try to round to trivial cases:
                if (fabs(_areas[i][j]) < 1E-12)
                    _areas[i][j] = 0.0;
                else if (fabs(1.0 - _areas[i][j]) < 1E-12)
                    _areas[i][j] = 1.0;
                // Catch overlapping polygons. This will not catch all overlapping polygons, but at least the user gets a message:
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

    vec polygons_for_single_material::rightanglerotate_2Dvec(
        vec v, int num, bool mirror_x_before, bool mirror_y_before) const
    {
        while (num < 0) // modulo works strange for negative numbers
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

    vec polygons_for_single_material::calc_grad_case1(const double A1, const double A2) const
    {
        double x1 = 2.0 * (sqrt(A1 * (A1 + A2)) - A1);
        double x2 = 2.0 * A2 - x1;
        double dy = x2 - x1;
        double a = sqrt(1.0 + dy * dy);
        return is_3D() ? vec(1.0 / a, dy / a, 0.0) : vec(1.0 / a, dy / a);
    };

    vec polygons_for_single_material::calc_grad_case2(const double A1, const double A2) const
    {
        double y1, y2, dx, a;
        if (A1 == 0.0) {
            y1 = 0.0;
            y2 = 2.0 * (1.0 - A2);
        }
        else {
            y1 = 2.0 * (A1 + sqrt(A1 * (1.0 - A2)));
            y2 = y1 * (0.5 * y1 / A1 - 1.0);
        }
        dx = y1 + y2;
        a = sqrt(1.0 + dx * dx);
        return is_3D() ? vec(dx / a, 1.0 / a, 0.0) : vec(dx / a, 1.0 / a);
    };

    vec polygons_for_single_material::calc_grad_case3(const double A1, const double A2) const
    {
        double dy = A2 - A1;
        double a = sqrt(1.0 + dy * dy);
        return is_3D() ? vec(1.0 / a, dy / a, 0.0) : vec(1.0 / a, dy / a);
    };

    vec polygons_for_single_material::calc_gradient_4x4(double** areas) const
    {
        double A11, A12, A21, A22, above_left, below_left, left_bottom;
        double above_right, below_right, right_top;
        // First, find out which rough position the interface has relative to the
        // squares, then we know which of the three cases to calculate and which two
        // squares to use for the calculation.
        // Since there are many possibilities, first check all possible positions
        // that are unique for one orientation. If the rough position could not be
        // found, rotate the whole problem by 90 degrees and try again, rotate the
        // resulting vector back 90 degrees accordingly if it could be calculated,
        // otherwise rotate more, and so forth:
        for (int rotate = 0; rotate < 4; rotate++)
        {
            switch (rotate)
            {
                case 0:
                    A11 = areas[1][1];
                    A12 = areas[1][2];
                    A21 = areas[2][1];
                    A22 = areas[2][2];
                    above_left = areas[1][3];
                    below_left = areas[1][0];
                    left_bottom = areas[0][1];
                    above_right = areas[2][3];
                    below_right = areas[2][0];
                    right_top = areas[3][2];
                    break;
                case 1: // rotate blocks in anticlockwise-direction:
                    A11 = areas[2][1];
                    A12 = areas[1][1];
                    A21 = areas[2][2];
                    A22 = areas[1][2];
                    above_left = areas[0][1];
                    below_left = areas[3][1];
                    left_bottom = areas[2][0];
                    above_right = areas[0][2];
                    below_right = areas[3][2];
                    right_top = areas[1][3];
                    break;
                case 2: // rotate blocks in anticlockwise-direction:
                    A11 = areas[2][2];
                    A12 = areas[2][1];
                    A21 = areas[1][2];
                    A22 = areas[1][1];
                    above_left = areas[2][0];
                    below_left = areas[2][3];
                    left_bottom = areas[3][2];
                    above_right = areas[1][0];
                    below_right = areas[1][3];
                    right_top = areas[0][1];
                    break;
                default: // rotate == 3: // rotate blocks in anticlockwise-direction:
                    A11 = areas[1][2];
                    A12 = areas[2][2];
                    A21 = areas[1][1];
                    A22 = areas[2][1];
                    above_left = areas[3][2];
                    below_left = areas[0][2];
                    left_bottom = areas[1][3];
                    above_right = areas[3][1];
                    below_right = areas[0][1];
                    right_top = areas[2][0];
                    break;
            }
            // check case:
            if ((A11 == 0 && A12 == 0 && A21 == 0 && A22 == 0) ||
                (A11 == 1 && A12 == 1 && A21 == 1 && A22 == 1))
                // no border in center block
                return vec(D2 + int(is_3D()), 0.0);
            if ((A11 == 0 && A22 == 0 && A12 != 0 && A21 != 0) ||
                (A12 == 0 && A21 == 0 && A11 != 0 && A22 != 0))
                // two separate polygons inside block;
                // the resolution should be increased
                return vec(D2 + int(is_3D()), 0.0);
            if (A12 == 0 && A21 == 0 && A22 == 0) // A11 != 0
            {
                if (left_bottom < below_left)
                    return rightanglerotate_2Dvec(
                        calc_grad_case2(A11, below_left),
                                                rotate);
                    else
                        return rightanglerotate_2Dvec(
                            calc_grad_case2(A11, left_bottom),
                                                    rotate - 1, true, false);
            }
            if (A11 != 0 && A12 != 0)
            {
                if (A21 == 0 && A22 == 0)
                {
                    if ((above_left != 0 && below_left != 0) ||
                        (above_left == 0 && below_left == 0))
                        // last one means there is a corner somewhere;
                        // then calculate approx. easy solution
                        return rightanglerotate_2Dvec(
                            calc_grad_case3(A12, A11), rotate);
                        else if (above_left == 0) // below_left != 0
                            return rightanglerotate_2Dvec(calc_grad_case1(A12, A11),
                                                        rotate);
                            else // above_left != 0; below_left == 0
                                return rightanglerotate_2Dvec(calc_grad_case1(A11, A12),
                                                            rotate, false, true);
                }
                else if (A21 != 0 && A22 == 0)
                {
                    if (A12 < A21)
                        return rightanglerotate_2Dvec(
                            calc_grad_case2(A12, A11),
                                                    rotate);
                        else
                            return rightanglerotate_2Dvec(
                                calc_grad_case2(A21, A11),
                                                        rotate - 1, true, false);
                }
                else if (A21 != 0 && A22 != 0) // all != 0
                {
                    if (A11 == 1 && A12 != 1 && A21 != 1 && A22 != 1)
                    {
                        if (A12 < A21)
                            return rightanglerotate_2Dvec(
                                calc_grad_case2(A22, A21),
                                                        rotate);
                            else
                                return rightanglerotate_2Dvec(
                                    calc_grad_case2(A22, A12),
                                                            rotate - 1, true, false);
                    }
                    else if (A11 == 1 && A12 == 1 && A21 != 1 && A22 != 1)
                    {
                        if ((above_right != 1 && below_right != 1) ||
                            (above_right == 1 && below_right == 1))
                            // last one means there is a corner somewhere;
                            // then calculate approx. easy solution
                            return rightanglerotate_2Dvec(
                                calc_grad_case3(A22, A21), rotate);
                            else if (above_right != 1) // below_right == 1
                                return rightanglerotate_2Dvec(
                                    calc_grad_case1(1.0 - A21, 1.0 - A22),
                                                            rotate);
                                else // above_right == 1; below_right != 1
                                    return rightanglerotate_2Dvec(
                                        calc_grad_case1(1.0 - A22, 1.0 - A21),
                                                                rotate, false, true);
                    }
                    else if (A11 == 1 && A12 == 1 && A21 == 1 && A22 != 1)
                    {
                        if (above_right < right_top)
                            return rightanglerotate_2Dvec(
                                calc_grad_case2(above_right, A22),
                                                        rotate);
                            else
                                return rightanglerotate_2Dvec(
                                    calc_grad_case2(right_top, A22),
                                                            rotate - 1, true, false);
                    }
                }
            }
        }
        // if this is reached, it can only mean there is a corner somewhere or even
        // multiple polygons in one block. To fix this, the resolution must be
        // increased. For now, an approximate solution based on the Sobel operator
        // (expanded for 4x4 pixel, to be symmetric around the interesting center
        // point of v) for edge detection is calculated:


        // Define the two kernels of the 'expanded' Sobel operator:
        int G_avg_weights[4] = {1, 2, 2, 1};
        int G_diff[4] = {+1, 0, 0, -1};
        double Gx = 0;
        double Gy = 0;
        for (int i = 0; i < 4; ++i)
            for (int j = 0; j < 4; ++j){
                Gx += G_avg_weights[j] * G_diff[i] * areas[i][j];
                Gy += G_avg_weights[i] * G_diff[j] * areas[i][j];
            }
            double a = sqrt(Gx * Gx + Gy * Gy);
        if (a == 0.0)
            return is_3D() ? vec(0.0, 0.0, 0.0) : vec(0.0, 0.0);
        else
            return is_3D() ? vec(Gx / a, Gy / a, 0.0) : vec(Gx / a, Gy / a);
    };

    bool polygons_for_single_material::index_in_bounds(int x, int y){
        return (x >= 0 && x < _nx && y >= 0 && y < _ny);
    }

    material_function_for_polygons::material_function_for_polygons(
        const grid_volume &thegv)
    {
        if ((thegv.dim != D2) && (thegv.dim != D3))
            abort("material_function_for_polygons only available for 2D and 3D.\n");
        nx = thegv.nx() * 2;
        ny = thegv.ny() * 2; // Number of blocks to split polygons into. One block should have halve the size of one pixel
        // to allow calculating averaged eps for different components (half pixel apart due to yee lattice).
        // The smooting diameter should stay 1 pixel (Meep default) or be an integer multiple thereof,
        // in which case material_function_for_polygons::eff_chi1inv_row() must be changed.
        nz = thegv.nz() * 2;

        block_size = thegv.inva * 0.5;
    };

    material_function_for_polygons::~material_function_for_polygons() {
        for (vector<polygons_for_single_material*>::const_iterator it = _polygons.begin();
            it != _polygons.end();
        ++it)
            {
                delete (*it); //delete polygons_for_single_material object
            }
            _polygons.clear(); //delete all pointers
    }

    //Adds a material stack to the material function, for use in 3D case. Returns mat_ID, which must
    //be used for add_polygon(polygon, mat_ID)
    unsigned int material_function_for_polygons::add_material_stack(
        double* material_heights, int mh_dim, double* epsilon_values, int ev_dim)
    {
        if (!is_3D()) {
            abort("Specifying material stack is only available in 3D.\n%s",
                "For 2D, please specify epsilon in add_polygon().\n");
        }
        material_stack * mat_stack = new material_stack(material_heights, mh_dim,
                                                        epsilon_values, ev_dim);
        _polygons.push_back(new polygons_for_single_material(mat_stack, block_size,
                                                            nx, ny, nz));
        //_material_stacks.push_back(new material_stack(material_heights,
        //              refractive_indices, number_of_layers));
        //_polygons.push_back(new queue<polygon*>());
        return _polygons.size() - 1;
    };

    // Adds a polygon with predefined stack (->mat_ID) to the simulation.
    // Care must be taken that none of the polygons overlap. Only for 3D simulation.
    polygon& material_function_for_polygons::add_polygon(polygon* pol, unsigned int mat_ID){
        // TODO? check if polygon overlaps other already added polygon:
        // first point is outside other polygon AND polygon's edges does not cross other polygon's edges
        if (!is_3D()) {
            abort("Specifying material stack is only available in 3D.\n%s",
                "For 2D, please specify epsilon in add_polygon().\n");
        }
        if (mat_ID >= _polygons.size())
            abort("Wrong mat_ID specified in add_polygon. Please add material stacks before adding polygons.");
        polygon* result = 0;
        if (pol->get_area() > 0){
            result = _polygons[mat_ID]->add_polygon(pol);
        }
        else {
            master_printf("polygon with zero area not added.\n");
        }
        return *result;
    };

    // Adds a polygon with predefined stack (->mat_ID) to the simulation.
    // Care must be taken that none of the polygons overlap. Only for 3D simulation.
    polygon& material_function_for_polygons::add_polygon(double* matrix, int dimX, int dimY, unsigned int mat_ID) {
        if (!is_3D()) {
            abort("Specifying material stack is only available in 3D.\n%s",
                "For 2D, please specify epsilon in add_polygon().\n");
        }
        return add_polygon(new polygon(matrix, dimX, dimY), mat_ID);
    };

    // Adds a polygon with epsilon value eps to the simulation.
    // Care must be taken that none of the polygons overlap.
    polygon& material_function_for_polygons::add_polygon(polygon* pol, double eps){
        // TODO? check if polygon overlaps other already added polygon:
        // first point is outside other polygon AND polygon's edges does not cross other polygon's edges
        if (is_3D()) {
            abort("Specifying polygon with epsilon value is only available in 2D.\n%s",
                "For 3D, please specify material_stack ID instead.\n");
        }
        polygon* result = 0;
        if (pol->get_area() > 0){
            // only adds eps if neccessary, else use id
            // of previously added eps:
            int id = add_epsilon(eps);
            result = _polygons[id]->add_polygon(pol);
        }
        else {
            master_printf("polygon with zero area not added.\n");
        }
        return *result;
    };

    // Adds a polygon with epsilon value eps to the simulation.
    // Care must be taken that none of the polygons overlap.
    polygon& material_function_for_polygons::add_polygon(double* matrix, int dimX, int dimY, double eps) {
        if (is_3D()) {
            abort("Specifying polygon with epsilon value is only available in 2D.\n%s",
                "For 3D, please specify material_stack ID instead.\n");
        }
        return add_polygon(new polygon(matrix, dimX, dimY), eps);
    };

    vec material_function_for_polygons::normal_vector(const ivec &center)
    {
        int debugn = 216; //108;
        for (vector<polygons_for_single_material*>::const_iterator it = _polygons.begin();
            it != _polygons.end();
        ++it)
            {
                if ((*it)->is_trivial(center - one_ivec(center.dim),
                    center + one_ivec(center.dim)))
                {
                    if (center == ivec(debugn, 0, 0))
                        master_printf("pol for single mat is trivial\n");
                    continue;
                }
                else {
                    if (center == ivec(debugn, 0, 0))
                        master_printf("pol for single mat is not trivial\n");
                    // Return normal vector of first material that's non-trivial in v.
                    // If a second material is non-trivial, it's (absolute) normal
                    // vector should be the same. If that's not the case, the polygons
                    // have been defined badly by the user, or the user should increase
                    // the resolution.
                    // TODO: This changes slightly in 3D: The sign of the normal vector's
                    // xy components could change for example in following situation:
                    // An interface dividing the pixel into two halves left and right,
                    // both sides with different material-stacks. The left side has
                    // an interface in z direction (high-eps below, 1 above), the right
                    // side has no interface, just the whole pixel height with high-eps.
                    // The resulting surface vector should point to high z and to the
                    // left, but when only the left side is regarded, it will point to
                    // high z and to the right. If only the right side is regarded, it
                    // will only point to the left.
                    // Possible solution: calculate normal_vector for all non-trivial
                    // polygons, then logically combine the results.
                    // Since it is not exactly clear how to handle corners and edges
                    // (also not in original MEEP-1.2.1), I will postpone the
                    // implementation until later.
                    return (*it)->normal_vector(center);
                }
            }
            if (!is_3D())
                // no interface of polygon in volume
                return vec(D2, 0.0);
            else
                // return vector pointing in z direction if on interface of material
                // stack whose polygon area is one (i.e. not zero):
                for (vector<polygons_for_single_material*>::const_iterator it = _polygons.begin();
                    it != _polygons.end();
                ++it)
                    {
                        // polygons must all be trivial (see above), so just check at one
                        // point whether area == 1:
                        if ((*it)->get_area(center) == 1) {
                            if (center == ivec(debugn, 0, 0))
                                master_printf("found area = 1\n");
                            double lower_layer_eps;
                            double upper_layer_eps;
                            double percentage_upper_layer;
                            if ((*it)->get_material_stack()->interface_inside_block(
                                center.z() * block_size, block_size,
                                                                                    lower_layer_eps, upper_layer_eps,
                                                                                    percentage_upper_layer))
                            {
                                if (center == ivec(debugn, 0, 0))
                                    master_printf("interface in z\n");
                                // interface in z only; return vector pointing to lower eps:
                                return vec(0.0, 0.0,
                                            upper_layer_eps < lower_layer_eps ? 1.0 : -1.0);
                            }
                            else {
                                if (center == ivec(debugn, 0, 0))
                                    master_printf("no interface in z\n");
                                // no interface, also not in z
                                return vec(D3, 0.0);
                            }
                        }
                    }
                    // no interface in block (all areas are trivial and no areas are one):
                    if (center == ivec(debugn, 0, 0))
                        master_printf("all blocks zero, i.e no interface\n");
                    return vec(D3, 0.0);
    }

    double material_function_for_polygons::mean_eps(const ivec &center){
        // start with complete area/volume, multiplied by 1 (air):
        double result = is_3D() ? 8.0 : 4.0;
        for (vector<polygons_for_single_material*>::const_iterator it = _polygons.begin();
            it != _polygons.end();
        ++it)
            {
                if (center == ivec(216, 0, 0)) {
                    master_printf("meaneps areas: \n" );
                    for (int i = center.x() - 1; i <  center.x() + 1; ++i)
                        master_printf("%f, %f\n", (*it)->get_area(i, center.y() - 1), (*it)->get_area(i, center.y()));
                }
                // add up all 4 areas (in xy):
                for (int i = center.x() - 1; i <  center.x() + 1; ++i)
                    for (int j = center.y() - 1; j < center.y() + 1; ++j)
                        if (!is_3D()) {
                            // add epsilon, weighted by area,
                            // minus 1 (=air - already added in beginning)
                            result += (*it)->get_area(i, j) * ((*it)->get_material_epsilon() - 1.0);
                        }
                        else // 3D
                            for (int k = center.z() - 1; k < center.z() + 1; ++k){
                                double lower_layer_eps;
                                double upper_layer_eps;
                                double percentage_upper_layer;
                                if ((*it)->get_material_stack()->interface_inside_block(
                                    k * block_size, block_size / 2.0,
                                    lower_layer_eps, upper_layer_eps,
                                    percentage_upper_layer))
                                {
                                    // add epsilon, weighted by volume,
                                    // minus 1 (=air - already added in beginning)
                                    result += (*it)->get_area(i, j) *
                                    ((1.0 - percentage_upper_layer) * (lower_layer_eps - 1.0)
                                    + percentage_upper_layer * (upper_layer_eps - 1.0));
                                }
                                else {
                                    result += (*it)->get_area(i, j) * (lower_layer_eps - 1.0);
                                    if (center == ivec(216, 0, 0))
                                        master_printf("meps result so far: %f (area %i, %i: %f; eps: %f)\n", result, i, j, (*it)->get_area(i, j), lower_layer_eps);
                                }
                            }
            }
            return (is_3D() ? result / 8.0 : result / 4.0);
        }

        double material_function_for_polygons::mean_inveps(const ivec &center){
            // start with complete area/volume, multiplied by 1 (air):
            double result = is_3D() ? 8.0 : 4.0;
            for (vector<polygons_for_single_material*>::const_iterator it = _polygons.begin();
                it != _polygons.end();
            ++it)
                {
                    // add up all 4 areas (in xy):
                    for (int i = center.x() - 1; i <  center.x() + 1; ++i)
                        for (int j = center.y() - 1; j < center.y() + 1; ++j)
                            if (!is_3D())
                                // add 1/epsilon, weighted by area,
                                // minus 1 (=air - already added in beginning)
                                result += (*it)->get_area(i, j) * (1.0 / (*it)->get_material_epsilon() - 1.0);
                            else // 3D
                                for (int k = center.z() - 1; k < center.z() + 1; ++k){
                                    double lower_layer_eps;
                                    double upper_layer_eps;
                                    double percentage_upper_layer;
                                    if ((*it)->get_material_stack()->interface_inside_block(
                                        k * block_size, block_size / 2.0,
                                        lower_layer_eps, upper_layer_eps,
                                        percentage_upper_layer))
                                    {
                                        // add 1/epsilon, weighted by volume,
                                        // minus 1 (=air - already added in beginning)
                                        result += (*it)->get_area(i, j) *
                                        ((1.0 - percentage_upper_layer) * (1.0 / lower_layer_eps - 1.0)
                                        + percentage_upper_layer * (1.0 / upper_layer_eps - 1.0));
                                    }
                                    else
                                        result += (*it)->get_area(i, j) * (1.0 / lower_layer_eps - 1.0);
                                }
                }
                return is_3D() ? result / 8.0 : result / 4.0;
        }

        ivec material_function_for_polygons::round_to_block_edge(const vec &v) {
            ivec iv(v.dim, 0.0);
            LOOP_OVER_DIRECTIONS(v.dim, d)
            {
                iv.set_direction(d,
                                my_round(v.in_direction(d) / block_size));
            }
            return iv;
        }

        unsigned int material_function_for_polygons::add_epsilon(double eps) {
            if (!is_3D()){ //2D
                for (std::size_t i = 0; i < _polygons.size(); ++i){
                    if (_polygons[i]->get_material_epsilon() == eps) {
                        //epsilon value already added before. Just return ID.
                        return i;
                    }
                }
                // add new polygon_for_single_material with epsilon value and return
                // new ID:
                _polygons.push_back(new polygons_for_single_material(eps, block_size,
                                                                    nx, ny));
                return _polygons.size() - 1;
            }
            else { //3D
                //TODO to be implemented: create a simple material stack with just one eps an no (inf?) thickness
                return 0;
            }
        };

        #define UNUSED(x) (void) x // silence compiler warnings
        void material_function_for_polygons::eff_chi1inv_row(component c,
                                                            double chi1inv_row[3],
                                                            const volume &v,
                                                            double tol, int maxeval)
        {
            UNUSED(tol);
            ivec center_ivec(round_to_block_edge(v.center()));
            if (!maxeval) {
                trivial:
                chi1inv_row[0] = chi1inv_row[1] = chi1inv_row[2] = 0.0;
                chi1inv_row[component_direction(c) % 3] = 1.0 / mean_eps(center_ivec);
                return;
            }

            if (v.center() == vec(108.0/20, 0.0, 0.0)) {
                master_printf("ivec: %i, %i, %i", center_ivec.x(), center_ivec.y(), center_ivec.z());
                master_printf("normalx: %f\n", normal_vector(center_ivec).x());
            }

            vec gradient(normal_vector(center_ivec));
            if (abs(gradient) < 1e-8)
                goto trivial;

            double meps = mean_eps(center_ivec);
            double minveps = mean_inveps(center_ivec);

            double n[3] = {0,0,0};
            //      double nabsinv = 1.0/abs(gradient);
            LOOP_OVER_DIRECTIONS(gradient.dim, k)
            n[k%3] = gradient.in_direction(k);// * nabsinv;

            /* get rownum'th row of effective tensor
            *                                                                                                        P * minveps + (I-P) * 1/meps = P * (minveps-1/mep*s) + I * 1/meps
            *                                                                                                        where I is the identity and P is the projection matrix
            *                                                                                                        P_{ij} = n[i] * n[j]. */
            int rownum = component_direction(c) % 3;
            for (int i=0; i<3; ++i)
                chi1inv_row[i] = n[rownum] * n[i] * (minveps - 1/meps);
            chi1inv_row[rownum] += 1/meps;
            if (v.center() == vec(108.0/20, 0.0, 0.0)) {
                master_printf("gradient: %f, %f, %f - ", gradient.x(), gradient.y(), gradient.z());
                master_printf("meps: %f, minveps: %f, nvec: (%f, %f, %f)\n", meps, minveps, n[0], n[1], n[2]);
            }
        }

////////////////////////////////////////////////////////////
//---- End: defs for polygon-based material functions ----//
////////////////////////////////////////////////////////////

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



} // namespace meep
