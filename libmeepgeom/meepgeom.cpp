/* Copyright (C) 2005-2015 Massachusetts Institute of Technology
%
%  This program is free software; you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation; either version 2, or (at your option)
%  any later version.
%
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.  %
%  You should have received a copy of the GNU General Public License
%  along with this program; if not, write to the Free Software Foundation,
%  Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
*/

#include <vector>
#include "meepgeom.hpp"

namespace meep_geom {

#define master_printf meep::master_printf

/***************************************************************/
/* global variables for default material                       */
/***************************************************************/
material_data vacuum_material_data;
material_type vacuum = &vacuum_material_data;

void check_offdiag(medium_struct *m) {
    if (m->epsilon_offdiag.x.im != 0 ||
        m->epsilon_offdiag.y.im != 0 ||
        m->epsilon_offdiag.z.im != 0 ||
        m->mu_offdiag.x.im != 0 ||
        m->mu_offdiag.y.im != 0 ||
        m->mu_offdiag.z.im != 0) {

        meep::abort("Found non-zero imaginary part of epsilon or mu offdiag.\n");
    }
}

bool susceptibility_equal(const susceptibility &s1, const susceptibility &s2)
{
    return (vector3_equal(s1.sigma_diag, s2.sigma_diag) &&
            vector3_equal(s1.sigma_offdiag, s2.sigma_offdiag) &&
            s1.frequency == s2.frequency && s1.gamma == s2.gamma &&
            s1.noise_amp == s2.noise_amp && s1.drude == s2.drude &&
            s1.is_file == s2.is_file);
}

bool susceptibility_list_equal(const susceptibility_list &s1, const susceptibility_list &s2)
{
    if (s1.num_items != s2.num_items) return false;
    for (int i = 0; i < s1.num_items; ++i)
      if (!susceptibility_equal(s1.items[i], s2.items[i]))
        return false;
    return true;
}

bool medium_struct_equal(const medium_struct *m1, const medium_struct *m2)
{
  return (vector3_equal(m1->epsilon_diag, m2->epsilon_diag) &&
          cvector3_equal(m1->epsilon_offdiag, m2->epsilon_offdiag) &&
          vector3_equal(m1->mu_diag, m2->mu_diag) &&
          cvector3_equal(m1->mu_offdiag, m2->mu_offdiag) &&
          vector3_equal(m1->E_chi2_diag, m2->E_chi2_diag) &&
          vector3_equal(m1->E_chi3_diag, m2->E_chi3_diag) &&
          vector3_equal(m1->H_chi2_diag, m2->H_chi2_diag) &&
          vector3_equal(m1->D_conductivity_diag, m2->D_conductivity_diag) &&
          vector3_equal(m1->B_conductivity_diag, m2->B_conductivity_diag) &&
          susceptibility_list_equal(m1->E_susceptibilities, m2->E_susceptibilities) &&
          susceptibility_list_equal(m1->H_susceptibilities, m2->H_susceptibilities));
}

// garbage collection for susceptibility_list structures.
// Assumes that the 'items' field, if non-empty, was allocated using new[];
// this is automatically the case for python code but is not checked
// for c++ code and will yield runtime errors if a user's user_material_func
// uses e.g. malloc() instead.
static void susceptibility_list_gc(susceptibility_list *sl)
{ if ( !sl || !(sl->num_items) ) return;
  delete[] sl->items;
  sl->items = NULL;
  sl->num_items=0;
}

// garbage collection for material structures: called to deallocate memory
// allocated for susceptibilities in user-defined materials.
void material_gc(material_type m)
 { if ( !m || m->which_subclass!=material_data::MATERIAL_USER ) return;
   susceptibility_list_gc( &(m->medium.E_susceptibilities) );
   susceptibility_list_gc( &(m->medium.H_susceptibilities) );
 }

bool material_type_equal(const material_type m1, const material_type m2)
{
    if (m1 == m2) return true;
    if (m1->which_subclass != m2->which_subclass) return false;
    switch (m1->which_subclass) {
      case material_data::MATERIAL_FILE: case material_data::PERFECT_METAL: return true;
      case material_data::MATERIAL_USER: return m1->user_func == m2->user_func && m1->user_data == m2->user_data;
      case material_data::MEDIUM: return medium_struct_equal( &(m1->medium), &(m2->medium));
      default: return false;
    }
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
typedef struct {
  double m00, m01, m02,
              m11, m12,
                   m22;
} symmetric_matrix;

/* rotate A by a unitary (real) rotation matrix R:
      RAR = transpose(R) * A * R
*/
void sym_matrix_rotate(symmetric_matrix *RAR,
			       const symmetric_matrix *A_,
			       const double R[3][3])
{
     int i,j;
     double A[3][3], AR[3][3];
     A[0][0] = A_->m00;
     A[1][1] = A_->m11;
     A[2][2] = A_->m22;
     A[0][1] = A[1][0] = A_->m01;
     A[0][2] = A[2][0] = A_->m02;
     A[1][2] = A[2][1] = A_->m12;
     for (i = 0; i < 3; ++i) for (j = 0; j < 3; ++j)
	  AR[i][j] = A[i][0]*R[0][j] + A[i][1]*R[1][j] + A[i][2]*R[2][j];
     for (i = 0; i < 3; ++i) for (j = i; j < 3; ++j)
	  A[i][j] = R[0][i]*AR[0][j] + R[1][i]*AR[1][j] + R[2][i]*AR[2][j];
     RAR->m00 = A[0][0];
     RAR->m11 = A[1][1];
     RAR->m22 = A[2][2];
     RAR->m01 = A[0][1];
     RAR->m02 = A[0][2];
     RAR->m12 = A[1][2];
}

/* Set Vinv = inverse of V, where both V and Vinv are real-symmetric matrices.*/
void sym_matrix_invert(symmetric_matrix *Vinv,
			       const symmetric_matrix *V)
{
     double m00 = V->m00, m11 = V->m11, m22 = V->m22;
     double m01 = V->m01, m02 = V->m02, m12 = V->m12;

     if (m01 == 0.0 && m02 == 0.0 && m12 == 0.0) {
	  /* optimize common case of a diagonal matrix: */
	  Vinv->m00 = 1.0 / m00;
	  Vinv->m11 = 1.0 / m11;
	  Vinv->m22 = 1.0 / m22;
	  Vinv->m01 = Vinv->m02 = Vinv->m12 = 0.0;
     }
     else {
	  double detinv;

	  /* compute the determinant: */
	  detinv = m00*m11*m22 - m02*m11*m02 + 2.0 * m01*m12*m02 -
	       m01*m01*m22 - m12*m12*m00;

	  if (detinv == 0.0) meep::abort( "singular 3x3 matrix");

	  detinv = 1.0/detinv;

	  Vinv->m00 = detinv * (m11*m22 - m12*m12);
	  Vinv->m11 = detinv * (m00*m22 - m02*m02);
	  Vinv->m22 = detinv * (m11*m00 - m01*m01);

	  Vinv->m02 = detinv * (m01*m12 - m11*m02);
	  Vinv->m01 = detinv * (m12*m02 - m01*m22);
	  Vinv->m12 = detinv * (m01*m02 - m00*m12);
     }
}

/* Returns whether or not V is positive-definite. */
int sym_matrix_positive_definite(symmetric_matrix *V)
{
     double det2, det3;
     double m00 = V->m00, m11 = V->m11, m22 = V->m22;

#if defined(WITH_HERMITIAN_EPSILON)
     scalar_complex m01 = V->m01, m02 = V->m02, m12 = V->m12;

     det2 = m00*m11 - CSCALAR_NORMSQR(m01);
     det3 = det2*m22 - m11*CSCALAR_NORMSQR(m02) - CSCALAR_NORMSQR(m12)*m00 +
	  2.0 * ((m01.re * m12.re - m01.im * m12.im) * m02.re +
		 (m01.re * m12.im + m01.im * m12.re) * m02.im);
#else /* real matrix */
     double m01 = V->m01, m02 = V->m02, m12 = V->m12;

     det2 = m00*m11 - m01*m01;
     det3 = det2*m22 - m02*m11*m02 + 2.0 * m01*m12*m02 - m12*m12*m00;
#endif /* real matrix */

     return (m00 > 0.0 && det2 > 0.0 && det3 > 0.0);
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
static meep::ndim dim = meep::D3;
void set_dimensions(int dims)
{
  if (dims == CYLINDRICAL) {
    dimensions = 2;
    dim = meep::Dcyl;
  }
  else {
    dimensions = dims;
    dim = meep::ndim(dims - 1);
  }
}

vector3 vec_to_vector3(const meep::vec &pt)
{
  vector3 v3;

  switch (pt.dim) {
  case meep::D1:
    v3.x = 0;
    v3.y = 0;
    v3.z = pt.z();
    break;
  case meep::D2:
    v3.x = pt.x();
    v3.y = pt.y();
    v3.z = 0;
    break;
  case meep::D3:
    v3.x = pt.x();
    v3.y = pt.y();
    v3.z = pt.z();
    break;
  case meep::Dcyl:
    v3.x = pt.r();
    v3.y = 0;
    v3.z = pt.z();
    break;
  }
  return v3;
}

meep::vec vector3_to_vec(const vector3 v3)
{
  switch (dim) {
  case meep::D1:
    return meep::vec(v3.z);
  case meep::D2:
    return meep::vec(v3.x, v3.y);
  case meep::D3:
    return meep::vec(v3.x, v3.y, v3.z);
  case meep::Dcyl:
    return meep::veccyl(v3.x, v3.z);
  default:
    meep::abort("unknown dimensionality in vector3_to_vec");
  }
}

static geom_box gv2box(const meep::volume &v)
{
  geom_box box;
  box.low = vec_to_vector3(v.get_min_corner());
  box.high = vec_to_vector3(v.get_max_corner());
  return box;
}

bool is_variable(material_type mt)
{
  return (mt->which_subclass == material_data::MATERIAL_USER);
}
bool is_variable(void* md) { return is_variable((material_type) md); }

bool is_file(material_type md)
{
  return (md->which_subclass == material_data::MATERIAL_FILE);
}
bool is_file(void* md) { return is_file((material_type) md); }

bool is_medium(material_type md, medium_struct **m)
{
  if (md->which_subclass == material_data::MEDIUM)
   { *m = &(md->medium);
     return true;
   };
  return false;
}

bool is_medium(void* md, medium_struct **m)
 { return is_medium((material_type) md, m); }

bool is_metal(meep::field_type ft, const material_type *material) {
  material_data *md = *material;
  if (ft == meep::E_stuff)
    switch (md->which_subclass) {
    case material_data::MEDIUM:
      return (md->medium.epsilon_diag.x < 0 ||
              md->medium.epsilon_diag.y < 0 ||
              md->medium.epsilon_diag.z < 0);
    case material_data::PERFECT_METAL:
      return true;
    default:
      meep::abort("unknown material type");
  }
  else
    switch (md->which_subclass) {
    case material_data::MEDIUM:
      return (md->medium.mu_diag.x < 0 ||
              md->medium.mu_diag.y < 0 ||
              md->medium.mu_diag.z < 0);
    case material_data::PERFECT_METAL:
      return false; // is an electric conductor, but not a magnetic conductor
    default:
      meep::abort("unknown material type");
  }
}

/* Linearly interpolate a given point in a 3d grid of data.  The point
   coordinates should be in the range [0,1], or at the very least [-1,2]
   ... anything outside [0,1] is *mirror* reflected into [0,1] */
meep::realnum linear_interpolate(
		     meep::realnum rx, meep::realnum ry, meep::realnum rz,
		     meep::realnum *data, int nx, int ny, int nz, int stride)
{
     int x, y, z, x2, y2, z2;
     meep::realnum dx, dy, dz;

     /* mirror boundary conditions for r just beyond the boundary */
     if (rx < 0.0) rx = -rx; else if (rx > 1.0) rx = 1.0 - rx;
     if (ry < 0.0) ry = -ry; else if (ry > 1.0) ry = 1.0 - ry;
     if (rz < 0.0) rz = -rz; else if (rz > 1.0) rz = 1.0 - rz;

     /* get the point corresponding to r in the epsilon array grid: */
     x = rx * nx; if (x == nx) --x;
     y = ry * ny; if (y == ny) --y;
     z = rz * nz; if (z == nz) --z;

     /* get the difference between (x,y,z) and the actual point
        ... we shift by 0.5 to center the data points in the pixels */
     dx = rx * nx - x - 0.5;
     dy = ry * ny - y - 0.5;
     dz = rz * nz - z - 0.5;

     /* get the other closest point in the grid, with mirror boundaries: */
     x2 = (dx >= 0.0 ? x + 1 : x - 1);
     if (x2 < 0) x2++; else if (x2 == nx) x2--;
     y2 = (dy >= 0.0 ? y + 1 : y - 1);
     if (y2 < 0) y2++; else if (y2 == ny) y2--;
     z2 = (dz >= 0.0 ? z + 1 : z - 1);
     if (z2 < 0) z2++; else if (z2 == nz) z2--;

     /* take abs(d{xyz}) to get weights for {xyz} and {xyz}2: */
     dx = fabs(dx);
     dy = fabs(dy);
     dz = fabs(dz);

     /* define a macro to give us data(x,y,z) on the grid,
	in row-major order (the order used by HDF5): */
#define D(x,y,z) (data[(((x)*ny + (y))*nz + (z)) * stride])

     return(((D(x,y,z)*(1.0-dx) + D(x2,y,z)*dx) * (1.0-dy) +
	     (D(x,y2,z)*(1.0-dx) + D(x2,y2,z)*dx) * dy) * (1.0-dz) +
	    ((D(x,y,z2)*(1.0-dx) + D(x2,y,z2)*dx) * (1.0-dy) +
	     (D(x,y2,z2)*(1.0-dx) + D(x2,y2,z2)*dx) * dy) * dz);

#undef D
}

// return material of the point p from the file (assumed already read)
void epsilon_file_material(material_data *md, vector3 p)
{
  default_material = (void*) md;

  if (md->which_subclass != material_data::MATERIAL_FILE)
    meep::abort("epsilon-input-file only works with a type=file default-material");

  if (!(md->epsilon_data)) return;
  medium_struct *mm=&(md->medium);
  double rx = geometry_lattice.size.x == 0
    ? 0 : 0.5 + (p.x-geometry_center.x) / geometry_lattice.size.x;
  double ry = geometry_lattice.size.y == 0
    ? 0 : 0.5 + (p.y-geometry_center.y) / geometry_lattice.size.y;
  double rz = geometry_lattice.size.z == 0
    ? 0 : 0.5 + (p.z-geometry_center.z) / geometry_lattice.size.z;
  mm->epsilon_diag.x = mm->epsilon_diag.y = mm->epsilon_diag.z =
    linear_interpolate(rx, ry, rz, md->epsilon_data,
		       md->epsilon_dims[0],
                       md->epsilon_dims[1],
                       md->epsilon_dims[2], 1);
  mm->epsilon_offdiag.x.re = mm->epsilon_offdiag.y.re = mm->epsilon_offdiag.z.re = 0;
}

struct pol {
  susceptibility user_s;
  struct pol *next;
};

// structure to hold a conductivity profile (for scalar absorbing layers)
struct cond_profile {
    double L; // thickness
    int N; // number of points prof[n] from 0..N corresponding to 0..L
    double *prof; // (NULL if none)
};


class geom_epsilon : public meep::material_function {
  geometric_object_list geometry;
  geom_box_tree geometry_tree;
  geom_box_tree restricted_tree;

  cond_profile cond[5][2]; // [direction][side]

public:
  geom_epsilon(geometric_object_list g, material_type_list mlist,
	       const meep::volume &v);
  virtual ~geom_epsilon();

  virtual void set_cond_profile(meep::direction, meep::boundary_side,
                                double L, double dx,
                                double (*prof)(int,double*,void*), void*,
                                double R);

  virtual void set_volume(const meep::volume &v);
  virtual void unset_volume(void);

  bool has_chi(meep::component c, int p);
  virtual bool has_chi3(meep::component c);
  virtual bool has_chi2(meep::component c);

  double chi(meep::component c, const meep::vec &r, int p);
  virtual double chi3(meep::component c, const meep::vec &r);
  virtual double chi2(meep::component c, const meep::vec &r);

  virtual bool has_mu();

  virtual bool has_conductivity(meep::component c);
  virtual double conductivity(meep::component c, const meep::vec &r);

  virtual double chi1p1(meep::field_type ft, const meep::vec &r);
  virtual void eff_chi1inv_row(meep::component c, double chi1inv_row[3],
                               const meep::volume &v,
                               double tol, int maxeval);

  void eff_chi1inv_matrix(meep::component c, symmetric_matrix *chi1inv_matrix,
                          const meep::volume &v, double tol, int maxeval);

  void fallback_chi1inv_row(meep::component c, double chi1inv_row[3],
			      const meep::volume &v,
			      double tol, int maxeval);

  virtual void sigma_row(meep::component c, double sigrow[3],
			 const meep::vec &r);
  void add_susceptibilities(meep::structure *s);
  void add_susceptibilities(meep::field_type ft, meep::structure *s);

  static bool verbose;

private:
  void get_material_pt(material_type &material, const meep::vec &r);

  material_type_list extra_materials;
  pol *current_pol;
};


/***********************************************************************/
bool geom_epsilon::verbose=false;

geom_epsilon::geom_epsilon(geometric_object_list g,
                           material_type_list mlist,
			   const meep::volume &v)
{
  geometry = g; // don't bother making a copy, only used in one place
  extra_materials = mlist;
  current_pol = NULL;

  FOR_DIRECTIONS(d) FOR_SIDES(b) cond[d][b].prof = NULL;

  if (meep::am_master()) {
    for (int i = 0; i < geometry.num_items; ++i) {

      display_geometric_object_info(5, geometry.items[i]);

      medium_struct *mm;
      if ( is_medium(geometry.items[i].material, &mm) ) {
        check_offdiag(mm);
        printf("%*sdielectric constant epsilon diagonal "
               "= (%g,%g,%g)\n", 5 + 5, "",
               mm->epsilon_diag.x,
               mm->epsilon_diag.y,
               mm->epsilon_diag.z
              );
      }
    }
  }

  geom_fix_objects0(geometry);
  geom_box box = gv2box(v);
  geometry_tree = create_geom_box_tree0(geometry, box);
  if (verbose && meep::am_master()) {
    printf("Geometric-object bounding-box tree:\n");
    display_geom_box_tree(5, geometry_tree);

    int tree_depth, tree_nobjects;
    geom_box_tree_stats(geometry_tree, &tree_depth, &tree_nobjects);
    master_printf("Geometric object tree has depth %d "
		  "and %d object nodes (vs. %d actual objects)\n",
		  tree_depth, tree_nobjects, geometry.num_items);
  }

  restricted_tree = geometry_tree;
}

geom_epsilon::~geom_epsilon()
{
  unset_volume();
  destroy_geom_box_tree(geometry_tree);
  FOR_DIRECTIONS(d) FOR_SIDES(b) if (cond[d][b].prof)
      delete[] cond[d][b].prof;
}

void geom_epsilon::set_cond_profile(meep::direction dir,
                                    meep::boundary_side side,
                                    double L, double dx,
                                    double (*P)(int,double*,void*), void *data,
                                    double R)
{
    if (cond[dir][side].prof) delete[] cond[dir][side].prof;

    int N = int(L / dx + 0.5);
    cond[dir][side].L = L;
    cond[dir][side].N = N;
    double *prof = cond[dir][side].prof = new double[N+1];

    double umin = 0, umax = 1, esterr;
    int errflag;
    double prof_int = adaptive_integration(P, &umin,&umax, 1, data, 1e-9, 1e-4,
                                           50000, &esterr, &errflag);

    double prefac = (-log(R)) / (4*L*prof_int);
    for (int i = 0; i <= N; ++i) {
        double u = double(i)/N;
        prof[i] = prefac * P(1, &u, data);
    }
}

void geom_epsilon::unset_volume(void)
{
  if (restricted_tree != geometry_tree) {
    destroy_geom_box_tree(restricted_tree);
    restricted_tree = geometry_tree;
  }
}

void geom_epsilon::set_volume(const meep::volume &v)
{
  unset_volume();

  geom_box box = gv2box(v);
  restricted_tree = create_geom_box_tree0(geometry, box);
}

static void material_epsmu(meep::field_type ft, material_type material,
		    symmetric_matrix *epsmu, symmetric_matrix *epsmu_inv) {

  material_data *md=material;
  if (ft == meep::E_stuff)
    switch (md->which_subclass) {

    case material_data::MEDIUM:
    case material_data::MATERIAL_FILE:
    case material_data::MATERIAL_USER:
      epsmu->m00 = md->medium.epsilon_diag.x;
      epsmu->m11 = md->medium.epsilon_diag.y;
      epsmu->m22 = md->medium.epsilon_diag.z;
      epsmu->m01 = md->medium.epsilon_offdiag.x.re;
      epsmu->m02 = md->medium.epsilon_offdiag.y.re;
      epsmu->m12 = md->medium.epsilon_offdiag.z.re;
      sym_matrix_invert(epsmu_inv,epsmu);
      break;

    case material_data::PERFECT_METAL:
      epsmu->m00 = -meep::infinity;
      epsmu->m11 = -meep::infinity;
      epsmu->m22 = -meep::infinity;
      epsmu_inv->m00 = -0.0;
      epsmu_inv->m11 = -0.0;
      epsmu_inv->m22 = -0.0;
      epsmu->m01 = epsmu->m02 = epsmu->m12 = 0.0;
      epsmu_inv->m01 = epsmu_inv->m02 = epsmu_inv->m12 = 0.0;
      break;

    default:
      meep::abort("unknown material type");
  }
  else
    switch (md->which_subclass) {
    case material_data::MEDIUM:
    case material_data::MATERIAL_FILE:
    case material_data::MATERIAL_USER:
      epsmu->m00 = md->medium.mu_diag.x;
      epsmu->m11 = md->medium.mu_diag.y;
      epsmu->m22 = md->medium.mu_diag.z;
      epsmu->m01 = md->medium.mu_offdiag.x.re;
      epsmu->m02 = md->medium.mu_offdiag.y.re;
      epsmu->m12 = md->medium.mu_offdiag.z.re;
      sym_matrix_invert(epsmu_inv,epsmu);
      break;

    case material_data::PERFECT_METAL:
      epsmu->m00 = 1.0;
      epsmu->m11 = 1.0;
      epsmu->m22 = 1.0;
      epsmu_inv->m00 = 1.0;
      epsmu_inv->m11 = 1.0;
      epsmu_inv->m22 = 1.0;
      epsmu->m01 = epsmu->m02 = epsmu->m12 = 0.0;
      epsmu_inv->m01 = epsmu_inv->m02 = epsmu_inv->m12 = 0.0;
      break;

    default:
      meep::abort("unknown material type");
  }
}

// the goal of this routine is to fill in the 'medium' field
// within the material structure as appropriate for the
// material properties at r.
void geom_epsilon::get_material_pt(material_type &material, const meep::vec &r)
{
  vector3 p = vec_to_vector3(r);
  boolean inobject;
  material = (material_type)
   material_of_unshifted_point_in_tree_inobject(p, restricted_tree,&inobject);
  material_data *md = material;

  switch(md->which_subclass)
   {
     // material read from file: interpolate to get properties at r
     case material_data::MATERIAL_FILE:
      if (md->epsilon_data)
       epsilon_file_material(md, p);
      else
       material = (material_type) default_material;
      return;

     // material specified by user-supplied function: call user
     // function to get properties at r.
     // Note that we initialize the medium to vacuum, so that
     // the user's function only needs to fill in whatever is
     // different from vacuum.
     case material_data::MATERIAL_USER:
      md->medium = medium_struct();
      md->user_func(p, md->user_data, &(md->medium));
      check_offdiag(&md->medium);
      return;

     // position-independent material or metal: there is nothing to do
     case material_data::MEDIUM:
     case material_data::PERFECT_METAL:
      return;

     default:
      meep::abort("unknown material type");
   };
}

// returns trace of the tensor diagonal
double geom_epsilon::chi1p1(meep::field_type ft, const meep::vec &r)
{
  symmetric_matrix chi1p1, chi1p1_inv;

#ifdef DEBUG
  vector3 p = vec_to_vector3(r);
  if (p.x < restricted_tree->b.low.x ||
      p.y < restricted_tree->b.low.y ||
      p.z < restricted_tree->b.low.z ||
      p.x > restricted_tree->b.high.x ||
      p.y > restricted_tree->b.high.y ||
      p.z > restricted_tree->b.high.z)
    meep::abort("invalid point (%g,%g,%g)\n", p.x,p.y,p.z);
#endif

  material_type material;
  get_material_pt(material, r);
  material_epsmu(ft, material, &chi1p1, &chi1p1_inv);
  material_gc(material);

  return (chi1p1.m00 + chi1p1.m11 + chi1p1.m22)/3;
}

/* Find frontmost object in v, along with the constant material behind it.
   Returns false if material behind the object is not constant.

   Requires moderately horrifying logic to figure things out properly,
   stolen from MPB. */
static bool get_front_object(const meep::volume &v,
			     geom_box_tree geometry_tree,
			     vector3 &pcenter,
			     const geometric_object **o_front,
			     vector3 &shiftby_front,
			     material_type &mat_front,
			     material_type &mat_behind) {
  vector3 p;
  const geometric_object *o1 = 0, *o2 = 0;
  vector3 shiftby1 = {0,0,0}, shiftby2 = {0,0,0};
  geom_box pixel;
  material_type mat1, mat2;
  int id1 = -1, id2 = -1;
  const int num_neighbors[3] = { 3, 5, 9 };
  const int neighbors[3][9][3] = {
    { {0,0,0}, {0,0,-1}, {0,0,1},
      {0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0} },
    { {0,0,0},
      {-1,-1,0}, {1,1,0}, {-1,1,0}, {1,-1,0},
      {0,0,0},{0,0,0},{0,0,0},{0,0,0} },
    { {0,0,0},
      {1,1,1},{1,1,-1},{1,-1,1},{1,-1,-1},
      {-1,1,1},{-1,1,-1},{-1,-1,1},{-1,-1,-1} }
  };
  pixel = gv2box(v);
  pcenter = p = vec_to_vector3(v.center());
  double d1, d2, d3;
  d1 = (pixel.high.x - pixel.low.x) * 0.5;
  d2 = (pixel.high.y - pixel.low.y) * 0.5;
  d3 = (pixel.high.z - pixel.low.z) * 0.5;
  for (int i = 0; i < num_neighbors[dimensions - 1]; ++i) {
    const geometric_object *o;
    material_type mat;
    vector3 q, shiftby;
    int id;
    q.x = p.x + neighbors[dimensions - 1][i][0] * d1;
    q.y = p.y + neighbors[dimensions - 1][i][1] * d2;
    q.z = p.z + neighbors[dimensions - 1][i][2] * d3;
    o = object_of_point_in_tree(q, geometry_tree, &shiftby, &id);
    if ((id == id1 && vector3_equal(shiftby, shiftby1)) ||
	(id == id2 && vector3_equal(shiftby, shiftby2)))
      continue;

    mat = (material_type) default_material;
    if (o)
     { material_data *md = (material_data *)o->material;
       if (md->which_subclass != material_data::MATERIAL_FILE)
        mat = md;
     }
    if (id1 == -1) {
      o1 = o;
      shiftby1 = shiftby;
      id1 = id;
      mat1 = mat;
    }
    else if (id2 == -1 || ( (id >= id1 && id >= id2) &&
			    (id1 == id2
			    || material_type_equal(mat1,mat2)))) {
      o2 = o;
      shiftby2 = shiftby;
      id2 = id;
      mat2 = mat;
    }
    else if (!(id1 < id2 &&
	       (id1 == id || material_type_equal(mat1,mat))) &&
	     !(id2 < id1 &&
	       (id2 == id || material_type_equal(mat2,mat))))
      return false;
  }

  // CHECK(id1 > -1, "bug in object_of_point_in_tree?");
  if (id2 == -1) { /* only one nearby object/material */
    id2 = id1;
    o2 = o1;
    mat2 = mat1;
    shiftby2 = shiftby1;
  }


  if (    (o1 && is_variable(o1->material))
       || (o2 && is_variable(o2->material))
       || ( (is_variable(default_material) || is_file(default_material))
            && (    !o1 || is_file(o1->material)
                 || !o2 || is_file(o2->material)
               )
          )
     )
    return false;

  if (id1 >= id2) {
    *o_front = o1;
    shiftby_front = shiftby1;
    mat_front = mat1;
    if (id1 == id2) mat_behind = mat1; else mat_behind = mat2;
  }
  if (id2 > id1) {
    *o_front = o2;
    shiftby_front = shiftby2;
    mat_front = mat2;
    mat_behind = mat1;
  }
  return true;
}

void geom_epsilon::eff_chi1inv_row(meep::component c, double chi1inv_row[3],
                                   const meep::volume &v, double tol, int maxeval) {
  symmetric_matrix meps_inv;
  eff_chi1inv_matrix(c, &meps_inv, v, tol, maxeval);

  switch (component_direction(c)) {
    case meep::X:
    case meep::R:
      chi1inv_row[0] = meps_inv.m00;
      chi1inv_row[1] = meps_inv.m01;
      chi1inv_row[2] = meps_inv.m02;
      break;
    case meep::Y:
    case meep::P:
      chi1inv_row[0] = meps_inv.m01;
      chi1inv_row[1] = meps_inv.m11;
      chi1inv_row[2] = meps_inv.m12;
      break;
    case meep::Z:
      chi1inv_row[0] = meps_inv.m02;
      chi1inv_row[1] = meps_inv.m12;
      chi1inv_row[2] = meps_inv.m22;
      break;
    case meep::NO_DIRECTION:
      chi1inv_row[0] = chi1inv_row[1] = chi1inv_row[2] = 0;
      break;
  }
}

void geom_epsilon::eff_chi1inv_matrix(meep::component c, symmetric_matrix *chi1inv_matrix,
                                      const meep::volume &v, double tol, int maxeval) {
  const geometric_object *o;
  material_type mat, mat_behind;
  symmetric_matrix meps;
  vector3 p, shiftby, normal;

  if (maxeval == 0 || !get_front_object(v, geometry_tree, p, &o, shiftby, mat, mat_behind)) {
  noavg:
    get_material_pt(mat, v.center());
  trivial:
    material_epsmu(meep::type(c), mat, &meps, chi1inv_matrix);
    material_gc(mat);
    return;
  }

  // FIXME: reimplement support for fallback integration, without
  //        messing up anisotropic support
  //  if (!get_front_object(v, geometry_tree,
  //                        p, &o, shiftby, mat, mat_behind)) {
  //     fallback_chi1inv_row(c, chi1inv_row, v, tol, maxeval);
  //     return;
  //  }

  /* check for trivial case of only one object/material */
  if (material_type_equal(mat, mat_behind)) goto trivial;

  // it doesn't make sense to average metals (electric or magnetic)
  if (is_metal(meep::type(c), &mat) || is_metal(meep::type(c), &mat_behind))
    goto noavg;

  normal = unit_vector3(normal_to_fixed_object(vector3_minus(p, shiftby), *o));
  if (normal.x == 0 && normal.y == 0 && normal.z == 0)
    goto noavg; // couldn't get normal vector for this point, punt
  geom_box pixel = gv2box(v);
  pixel.low = vector3_minus(pixel.low, shiftby);
  pixel.high = vector3_minus(pixel.high, shiftby);

  double fill = box_overlap_with_object(pixel, *o, tol, maxeval);

  material_epsmu(meep::type(c), mat, &meps, chi1inv_matrix);
  symmetric_matrix eps2, epsinv2;
  symmetric_matrix eps1, delta;
  double Rot[3][3];
  material_epsmu(meep::type(c), mat_behind, &eps2, &epsinv2);
  eps1 = meps;

  Rot[0][0] = normal.x;
  Rot[1][0] = normal.y;
  Rot[2][0] = normal.z;
  if (fabs(normal.x) > 1e-2 || fabs(normal.y) > 1e-2) {
    Rot[0][2] = normal.y;
    Rot[1][2] = -normal.x;
    Rot[2][2] = 0;
  }
  else { /* n is not parallel to z direction, use (x x n) instead */
    Rot[0][2] = 0;
    Rot[1][2] = -normal.z;
    Rot[2][2] = normal.y;
  }
  { /* normalize second column */
    double s = Rot[0][2]*Rot[0][2]+Rot[1][2]*Rot[1][2]+Rot[2][2]*Rot[2][2];
    s = 1.0 / sqrt(s);
    Rot[0][2] *= s;
    Rot[1][2] *= s;
    Rot[2][2] *= s;
  }
  /* 1st column is 2nd column x 0th column */
  Rot[0][1] = Rot[1][2] * Rot[2][0] - Rot[2][2] * Rot[1][0];
  Rot[1][1] = Rot[2][2] * Rot[0][0] - Rot[0][2] * Rot[2][0];
  Rot[2][1] = Rot[0][2] * Rot[1][0] - Rot[1][2] * Rot[0][0];

  /* rotate epsilon tensors to surface parallel/perpendicular axes */
  sym_matrix_rotate(&eps1, &eps1, Rot);
  sym_matrix_rotate(&eps2, &eps2, Rot);

#define AVG (fill * (EXPR(eps1)) + (1-fill) * (EXPR(eps2)))
#define SQR(x) ((x) * (x))

#define EXPR(eps) (-1 / eps.m00)
  delta.m00 = AVG;
#undef EXPR
#define EXPR(eps) (eps.m11 - SQR(eps.m01) / eps.m00)
  delta.m11 = AVG;
#undef EXPR
#define EXPR(eps) (eps.m22 - SQR(eps.m02) / eps.m00)
  delta.m22 = AVG;
#undef EXPR

#define EXPR(eps) (eps.m01 / eps.m00)
  delta.m01 = AVG;
#undef EXPR
#define EXPR(eps) (eps.m02 / eps.m00)
  delta.m02 = AVG;
#undef EXPR
#define EXPR(eps) (eps.m12 - eps.m02 * eps.m01 / eps.m00)
  delta.m12 = AVG;
#undef EXPR

  meps.m00 = -1/delta.m00;
  meps.m11 = delta.m11 - SQR(delta.m01) / delta.m00;
  meps.m22 = delta.m22 - SQR(delta.m02) / delta.m00;
  meps.m01 = -delta.m01/delta.m00;
  meps.m02 = -delta.m02/delta.m00;
  meps.m12 = delta.m12 - (delta.m02 * delta.m01) / delta.m00;

#undef SQR

#define SWAP(a,b) { double xxx = a; a = b; b = xxx; }
  /* invert rotation matrix = transpose */
  SWAP(Rot[0][1], Rot[1][0]);
  SWAP(Rot[0][2], Rot[2][0]);
  SWAP(Rot[2][1], Rot[1][2]);
  sym_matrix_rotate(&meps, &meps, Rot); /* rotate back */
#undef SWAP

#ifdef DEBUG
  if(!sym_matrix_positive_definite(&meps))
    meep::abort("negative mean epsilon from Kottke algorithm");
#endif

  sym_matrix_invert(chi1inv_matrix, &meps);
}

static int eps_ever_negative = 0;
static meep::field_type func_ft = meep::E_stuff;

#ifdef CTL_HAS_COMPLEX_INTEGRATION
static cnumber ceps_func(int n, number *x, void *geomeps_)
{
  geom_epsilon *geomeps = (geom_epsilon *) geomeps_;
  vector3 p = {0,0,0};
  p.x = x[0]; p.y = n > 1 ? x[1] : 0; p.z = n > 2 ? x[2] : 0;
  double s = 1;
  if (dim == meep::Dcyl) { double py = p.y; p.y = p.z; p.z = py; s = p.x; }
  cnumber ret;
  double ep = geomeps->chi1p1(func_ft, vector3_to_vec(p));
  if (ep < 0) eps_ever_negative = 1;
  ret.re = ep * s;
  ret.im = s / ep;
  return ret;
}
#else
static number eps_func(int n, number *x, void *geomeps_)
{
  geom_epsilon *geomeps = (geom_epsilon *) geomeps_;
  vector3 p = {0,0,0};
  double s = 1;
  p.x = x[0]; p.y = n > 1 ? x[1] : 0; p.z = n > 2 ? x[2] : 0;
  if (dim == meep::Dcyl) { double py = p.y; p.y = p.z; p.z = py; s = p.x; }
  double ep = geomeps->chi1p1(func_ft, vector3_to_vec(p));
  if (ep < 0) eps_ever_negative = 1;
  return ep * s;
}
static number inveps_func(int n, number *x, void *geomeps_)
{
  geom_epsilon *geomeps = (geom_epsilon *) geomeps_;
  vector3 p = {0,0,0};
  double s = 1;
  p.x = x[0]; p.y = n > 1 ? x[1] : 0; p.z = n > 2 ? x[2] : 0;
  if (dim == meep::Dcyl) { double py = p.y; p.y = p.z; p.z = py; s = p.x; }
  double ep = geomeps->chi1p1(func_ft, vector3_to_vec(p));
  if (ep < 0) eps_ever_negative = 1;
  return s / ep;
}
#endif

// fallback meaneps using libctl's adaptive cubature routine
void geom_epsilon::fallback_chi1inv_row(meep::component c,
					double chi1inv_row[3],
					const meep::volume &v,
					double tol, int maxeval)
{

  symmetric_matrix chi1p1, chi1p1_inv;
  material_type material;
  get_material_pt(material, v.center());
  material_epsmu(meep::type(c), material, &chi1p1, &chi1p1_inv);
  material_gc(material);
  if (chi1p1.m01 != 0 || chi1p1.m02 != 0 || chi1p1.m12 != 0
      || chi1p1.m00 != chi1p1.m11 || chi1p1.m11 != chi1p1.m22 ||
      chi1p1.m00 != chi1p1.m22) {
    int rownum = meep::component_direction(c) % 3;
    if (rownum == 0) {
      chi1inv_row[0] = chi1p1.m00;
      chi1inv_row[1] = chi1p1.m01;
      chi1inv_row[2] = chi1p1.m02;
    }
    else if (rownum == 1) {
      chi1inv_row[0] = chi1p1.m01;
      chi1inv_row[1] = chi1p1.m11;
      chi1inv_row[2] = chi1p1.m12;
    }
    else {
      chi1inv_row[0] = chi1p1.m02;
      chi1inv_row[1] = chi1p1.m12;
      chi1inv_row[2] = chi1p1.m22;
    }
    return;
      }

  number esterr;
  integer errflag, n;
  number xmin[3], xmax[3];
  vector3 gvmin, gvmax;
  gvmin = vec_to_vector3(v.get_min_corner());
  gvmax = vec_to_vector3(v.get_max_corner());
  xmin[0] = gvmin.x; xmax[0] = gvmax.x;
  if (dim == meep::Dcyl) {
    xmin[1] = gvmin.z; xmin[2] = gvmin.y; xmax[1] = gvmax.z; xmax[2] = gvmax.y;
  }
  else{
    xmin[1] = gvmin.y; xmin[2] = gvmin.z; xmax[1] = gvmax.y; xmax[2] = gvmax.z;
  }
  if (xmin[2] == xmax[2])
    n = xmin[1] == xmax[1] ? 1 : 2;
  else
    n = 3;
  double vol = 1;
  for (int i = 0; i < n; ++i) vol *= xmax[i] - xmin[i];
  if (dim == meep::Dcyl) vol *= (xmin[0] + xmax[0]) * 0.5;
  eps_ever_negative = 0;
  func_ft = meep::type(c);
  double meps, minveps;
#ifdef CTL_HAS_COMPLEX_INTEGRATION
  cnumber ret = cadaptive_integration(ceps_func, xmin, xmax, n, (void*) this,
				      0, tol, maxeval, &esterr, &errflag);
  meps = ret.re / vol;
  minveps = ret.im / vol;
#else
  meps = adaptive_integration(eps_func, xmin, xmax, n, (void*) this,
			      0, tol, maxeval, &esterr, &errflag) / vol;
  minveps = adaptive_integration(inveps_func, xmin, xmax, n, (void*) this,
				 0, tol, maxeval, &esterr, &errflag) / vol;
#endif
  if (eps_ever_negative) // averaging negative eps causes instability
    minveps = 1.0 / (meps = eps(v.center()));

  {
    meep::vec gradient(normal_vector(meep::type(c), v));
    double n[3] = {0,0,0};
    double nabsinv = 1.0/meep::abs(gradient);
    LOOP_OVER_DIRECTIONS(gradient.dim, k)
      n[k%3] = gradient.in_direction(k) * nabsinv;
    int rownum = meep::component_direction(c) % 3;
    for (int i=0; i<3; ++i)
      chi1inv_row[i] = n[rownum] * n[i] * (minveps - 1/meps);
    chi1inv_row[rownum] += 1/meps;
  }
}

static double get_chi3(meep::component c, const medium_struct *m) {
  switch (c) {
  case meep::Er: case meep::Ex: return m->E_chi3_diag.x;
  case meep::Ep: case meep::Ey: return m->E_chi3_diag.y;
  case meep::Ez: return m->E_chi3_diag.z;
  case meep::Hr: case meep::Hx: return m->H_chi3_diag.x;
  case meep::Hp: case meep::Hy: return m->H_chi3_diag.y;
  case meep::Hz: return m->H_chi3_diag.z;
  default: return 0;
  }
}

static double get_chi2(meep::component c, const medium_struct *m) {
  switch (c) {
  case meep::Er: case meep::Ex: return m->E_chi2_diag.x;
  case meep::Ep: case meep::Ey: return m->E_chi2_diag.y;
  case meep::Ez: return m->E_chi2_diag.z;
  case meep::Hr: case meep::Hx: return m->H_chi2_diag.x;
  case meep::Hp: case meep::Hy: return m->H_chi2_diag.y;
  case meep::Hz: return m->H_chi2_diag.z;
  default: return 0;
  }
}

static double get_chi(meep::component c, const medium_struct *m, int p)
 { return ((p==2) ? get_chi2(c,m) : get_chi3(c,m)); }

// the following routines consolidate the former has_chi2, has_chi3 and chi2, chi3
// p=2,3 for chi2,chi3
bool geom_epsilon::has_chi(meep::component c, int p)
{
  medium_struct *mm;

  for (int i = 0; i < geometry.num_items; ++i)
   if (is_medium(geometry.items[i].material, &mm) )
    if (get_chi(c, mm, p)!=0)
     return true;

  for (int i = 0; i < extra_materials.num_items; ++i)
   if (is_medium(extra_materials.items[i], &mm ))
    if (get_chi(c, mm, p )!=0)
     return true;

  return (    is_medium(default_material, &mm)
           && get_chi(c, mm, p) != 0
         );
}

bool geom_epsilon::has_chi3(meep::component c)
{ return has_chi(c, 3); }

bool geom_epsilon::has_chi2(meep::component c)
{ return has_chi(c, 2); }

double geom_epsilon::chi(meep::component c, const meep::vec &r, int p) {
  material_type material;
  get_material_pt(material, r);

  double chi_val;

  material_data *md = material;
  switch (md->which_subclass)
   {
     case material_data::MEDIUM:
     case material_data::MATERIAL_USER:
       chi_val = get_chi(c, &(md->medium), p );
       break;

     default:
       chi_val = 0;
   };

  material_gc(material);
  return chi_val;
}

double geom_epsilon::chi3(meep::component c, const meep::vec &r)
 { return chi(c, r, 3); }

double geom_epsilon::chi2(meep::component c, const meep::vec &r)
 { return chi(c, r, 2); }

static bool mu_not_1(material_type m)
{
  medium_struct *mm;
  return (    is_medium(m, &mm)
           && (     mm->mu_diag.x!=1
                ||  mm->mu_diag.y!=1
                ||  mm->mu_diag.z!=1
                ||  mm->mu_offdiag.x.re!=0
                ||  mm->mu_offdiag.y.re!=0
                ||  mm->mu_offdiag.z.re!=0
              )
        );
}
static bool mu_not_1(void *m) { return mu_not_1((material_type) m); }

bool geom_epsilon::has_mu()
{
  for (int i = 0; i < geometry.num_items; ++i)
   if (mu_not_1(geometry.items[i].material))
    return true;

  for (int i = 0; i < extra_materials.num_items; ++i)
   if (mu_not_1(extra_materials.items[i]))
    return true;

  return (mu_not_1(default_material));
}

/* a global scalar conductivity to add to all materials; this
   is mostly for the convenience of Casimir calculations where
   the global conductivity corresponds to a rotation to
   complex frequencies */
static double global_D_conductivity = 0, global_B_conductivity = 0;

static double get_cnd(meep::component c, const medium_struct *m) {
  switch (c) {
  case meep::Dr: case meep::Dx: return m->D_conductivity_diag.x + global_D_conductivity;
  case meep::Dp: case meep::Dy: return m->D_conductivity_diag.y + global_D_conductivity;
  case meep::Dz: return m->D_conductivity_diag.z + global_D_conductivity;
  case meep::Br: case meep::Bx: return m->B_conductivity_diag.x + global_B_conductivity;
  case meep::Bp: case meep::By: return m->B_conductivity_diag.y + global_B_conductivity;
  case meep::Bz: return m->B_conductivity_diag.z + global_B_conductivity;
  default: return 0;
  }
}

bool geom_epsilon::has_conductivity(meep::component c)
{
  medium_struct *mm;

  FOR_DIRECTIONS(d) FOR_SIDES(b) if (cond[d][b].prof) return true;

  for (int i = 0; i < geometry.num_items; ++i)
   if (    is_medium(geometry.items[i].material,&mm)
        && get_cnd(c,mm)
      ) return true;

  for (int i = 0; i < extra_materials.num_items; ++i)
   if (    is_medium(extra_materials.items[i],&mm)
        && get_cnd(c,mm)
      ) return true;

  return (    is_medium(default_material, &mm)
	   && get_cnd(c, mm) !=0
         );
}

static meep::vec geometry_edge; // geometry_lattice.size / 2
double geom_epsilon::conductivity(meep::component c, const meep::vec &r) {
  material_type material;
  get_material_pt(material, r);

  double cond_val;
  material_data *md = material;
  switch (md->which_subclass) {
  case material_data::MEDIUM:
    cond_val = get_cnd(c, &(md->medium));
    break;
  default:
    cond_val = 0;
  }
  material_gc(material);

  // if the user specified scalar absorbing layers, add their conductivities
  // to cond_val (isotropically, for both magnetic and electric conductivity).
  LOOP_OVER_DIRECTIONS(r.dim, d) {
      double x = r.in_direction(d);
      double edge = geometry_edge.in_direction(d) - cond[d][meep::High].L;
      if (cond[d][meep::High].prof && x >= edge) {
          int N = cond[d][meep::High].N;
          double ui = N * (x-edge) / cond[d][meep::High].L;
          int i = int(ui);
          if (i >= N)
              cond_val += cond[d][meep::High].prof[N];
          else {
              double di = ui - i;
              cond_val += cond[d][meep::High].prof[i] * (1-di)
                        + cond[d][meep::High].prof[i+1] * di;
          }
      }
      edge = cond[d][meep::Low].L - geometry_edge.in_direction(d);
      if (cond[d][meep::Low].prof && x <= edge) {
          int N = cond[d][meep::Low].N;
          double ui = N * (edge-x) / cond[d][meep::Low].L;
          int i = int(ui);
          if (i >= N)
              cond_val += cond[d][meep::Low].prof[N];
          else {
              double di = ui - i;
              cond_val += cond[d][meep::Low].prof[i] * (1-di)
                        + cond[d][meep::Low].prof[i+1] * di;
          }
      }
  }

  return cond_val;
}

/* like susceptibility_equal in ctl-io.cpp, but ignores sigma and id
   (must be updated manually, re-copying from ctl-io.cpp), if we
   add new susceptibility subclasses) */
static bool susceptibility_equiv(const susceptibility *o0,
				 const susceptibility *o)
{
  if (o0->frequency != o->frequency) return 0;
  if (o0->gamma     != o->gamma)     return 0;
  if (o0->noise_amp != o->noise_amp) return 0;
  if (o0->drude     != o->drude)     return 0;
  if (o0->is_file   != o->is_file)   return 0;

  return 1;
}

void geom_epsilon::sigma_row(meep::component c, double sigrow[3],
			     const meep::vec &r) {

  vector3 p = vec_to_vector3(r);

  boolean inobject;
  material_type mat =
      (material_type)material_of_unshifted_point_in_tree_inobject(p, restricted_tree, &inobject);

  if (mat->which_subclass == material_data::MATERIAL_USER) {
      mat->medium = medium_struct();
      mat->user_func(p, mat->user_data, &(mat->medium));
      check_offdiag(&mat->medium);
  }

  sigrow[0] = sigrow[1] = sigrow[2] = 0.0;

  if (mat->which_subclass == material_data::MATERIAL_USER ||
      mat->which_subclass == material_data::MEDIUM) {

    susceptibility_list slist =
      type(c) == meep::E_stuff ? mat->medium.E_susceptibilities : mat->medium.H_susceptibilities;
    for (int j = 0; j < slist.num_items; ++j) {
      if (susceptibility_equiv(&slist.items[j], &current_pol->user_s)) {
        int ic = meep::component_index(c);
        switch (ic) { // which row of the sigma tensor to return
        case 0:
          sigrow[0] = slist.items[j].sigma_diag.x;
          sigrow[1] = slist.items[j].sigma_offdiag.x;
          sigrow[2] = slist.items[j].sigma_offdiag.y;
          break;
        case 1:
          sigrow[0] = slist.items[j].sigma_offdiag.x;
          sigrow[1] = slist.items[j].sigma_diag.y;
          sigrow[2] = slist.items[j].sigma_offdiag.z;
          break;
        default: // case 2:
          sigrow[0] = slist.items[j].sigma_offdiag.y;
          sigrow[1] = slist.items[j].sigma_offdiag.z;
          sigrow[2] = slist.items[j].sigma_diag.z;
          break;
        }
        break;
      }
    }
  }
  material_gc(mat);
}

// add a polarization to the list if it is not already there
static pol *add_pol(pol *pols, const susceptibility *user_s)
{
  struct pol *p = pols;
  while (p && !susceptibility_equiv(user_s, &p->user_s)) p = p->next;
  if (!p) {
    p = new pol;
    p->user_s = *user_s;
    p->next = pols;
    pols = p;
  }
  return pols;
}

static pol *add_pols(pol *pols, const susceptibility_list slist) {
  for (int j = 0; j < slist.num_items; ++j)
    pols = add_pol(pols, &slist.items[j]);
  return pols;
}

void geom_epsilon::add_susceptibilities(meep::structure *s) {
  add_susceptibilities(meep::E_stuff, s);
  add_susceptibilities(meep::H_stuff, s);
}

void geom_epsilon::add_susceptibilities(meep::field_type ft,
					meep::structure *s) {
  pol *pols = 0;
  medium_struct *mm;

  // construct a list of the unique susceptibilities in the geometry:
  for (int i = 0; i < geometry.num_items; ++i)
   if ( is_medium( geometry.items[i].material, &mm) )
    pols = add_pols(pols, ft == meep::E_stuff
		          ? mm->E_susceptibilities
		          : mm->H_susceptibilities
                   );

  for (int i = 0; i < extra_materials.num_items; ++i)
   if ( is_medium(extra_materials.items[i], &mm ) )
    pols = add_pols(pols, ft == meep::E_stuff
		          ? mm->E_susceptibilities
		          : mm->H_susceptibilities
                   );

  if ( is_medium(default_material, &mm ) )
    pols = add_pols(pols, ft == meep::E_stuff
	          	  ? mm->E_susceptibilities
	  	          : mm->H_susceptibilities
                   );

  for (struct pol *p = pols; p; p = p->next) {

    susceptibility *ss=&(p->user_s);
    if (ss->is_file)
     meep::abort("unknown susceptibility");
    bool noisy = (ss->noise_amp!=0.0);

    meep::susceptibility *sus =
     noisy ? new meep::noisy_lorentzian_susceptibility(ss->noise_amp,
                                                       ss->frequency,
                                                       ss->gamma,
                                                       ss->drude
                                                      )
           : new meep::lorentzian_susceptibility(ss->frequency,
                                                 ss->gamma,
                                                 ss->drude
                                                );

    master_printf("%s%s susceptibility: frequency=%g, gamma=%g",
                  noisy     ? "noisy " : "",
                  ss->drude ? "drude"  : "lorentzian",
                  ss->frequency, ss->gamma);
    if(noisy) master_printf(" amp=%g ",ss->noise_amp);
    master_printf("\n");

    current_pol = p;
    if(sus) {
      s->add_susceptibility(*this, ft, *sus);
      delete sus;
    }
  }
  current_pol = NULL;

  while (pols) {
    struct pol *p = pols;
    pols = pols->next;
    delete p;
  }
}

typedef struct pml_profile_thunk
{
 meep::pml_profile_func func;
 void *func_data;
} pml_profile_thunk;

double pml_profile_wrapper(int dim, double *u, void *user_data)
{ (void )dim; // unused
  pml_profile_thunk *mythunk = (pml_profile_thunk *)user_data;
  return mythunk->func(u[0], mythunk->func_data);
}

/***************************************************************/
/* mechanism for allowing users to specify non-PML absorbing   */
/* layers.                                                     */
/***************************************************************/
absorber_list create_absorber_list()
{
  absorber_list alist = new absorber_list_type;
  return alist;
}

void destroy_absorber_list(absorber_list alist)
{ delete alist;
}

void add_absorbing_layer(absorber_list alist,
                         double thickness, int direction, int side,
                         double R_asymptotic, double mean_stretch,
                         meep::pml_profile_func func, void *func_data)
{
  absorber myabsorber;
  myabsorber.thickness=thickness;
  myabsorber.direction=direction;
  myabsorber.side=side;
  myabsorber.R_asymptotic=R_asymptotic;
  myabsorber.mean_stretch=mean_stretch;
  myabsorber.pml_profile=func;
  myabsorber.pml_profile_data=func_data;

  if (alist==0)
   meep::abort("invalid absorber_list in add_absorbing_layer");

  alist->push_back( myabsorber );
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void set_materials_from_geometry(meep::structure *s,
                                 geometric_object_list g,
                                 bool use_anisotropic_averaging,
                                 double tol,
                                 int maxeval,
                                 bool _ensure_periodicity,
                                 bool verbose,
                                 material_type _default_material,
                                 absorber_list alist,
                                 material_type_list extra_materials)
{
  geom_epsilon::verbose=verbose;

  // set global variables in libctlgeom based on data fields in s
  geom_initialize();
  if (_default_material->which_subclass != material_data::MATERIAL_USER &&
      _default_material->which_subclass != material_data::PERFECT_METAL) {
      check_offdiag(&_default_material->medium);
  }
  default_material     = _default_material;
  ensure_periodicity   = _ensure_periodicity;
  meep::grid_volume gv = s->gv;
  double resolution    = gv.a;

  dimensions=3;
  vector3 size = {0.0,0.0,0.0};
  switch (s->user_volume.dim)
   { case meep::D1:   dimensions=1;
                      size.z = s->user_volume.nz()/resolution;
                      break;
     case meep::D2:   dimensions=2;
                      size.x = s->user_volume.nx()/resolution;
                      size.y = s->user_volume.ny()/resolution;
                      break;
     case meep::D3:   dimensions=3;
                      size.x = s->user_volume.nx()/resolution;
                      size.y = s->user_volume.ny()/resolution;
                      size.z = s->user_volume.nz()/resolution;
                      break;
     case meep::Dcyl: dimensions= CYLINDRICAL;
                      size.x = s->user_volume.nr()/resolution;
                      size.z = s->user_volume.nz()/resolution;
                      break;
   };

  set_dimensions(dimensions);

  geometry_lattice.size = size;
  geometry_edge = vector3_to_vec(size) * 0.5;

  master_printf("Working in %s dimensions.\n",
                 meep::dimension_name(s->gv.dim));
  master_printf("Computational cell is %g x %g x %g with resolution %g\n",
                size.x, size.y, size.z, resolution);

  geom_epsilon geps(g, extra_materials, gv.pad().surroundings());

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  if (alist)
   { for(absorber_list_type::iterator layer=alist->begin(); layer!=alist->end(); layer++)
      { LOOP_OVER_DIRECTIONS(gv.dim,d)
         { if (layer->direction!=ALL_DIRECTIONS && layer->direction!=d) continue;
           FOR_SIDES(b)
            { if (layer->side!=ALL_SIDES && layer->side!=b ) continue;
              pml_profile_thunk mythunk;
              mythunk.func      = layer->pml_profile;
              mythunk.func_data = layer->pml_profile_data;
              geps.set_cond_profile(d,b,layer->thickness, gv.inva*0.5,
                                    pml_profile_wrapper, (void *)&mythunk,
                                    layer->R_asymptotic);
            };
         };
      };
   };

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  s->set_materials(geps, use_anisotropic_averaging, tol, maxeval);
  s->remove_susceptibilities();
  geps.add_susceptibilities(s);

  master_printf("-----------\n");
}

/***************************************************************/
/* convenience routines for creating materials of various types*/
/***************************************************************/
material_type make_dielectric(double epsilon)
{
  material_data *md = (material_data *)malloc(sizeof(*md));
  md->which_subclass=material_data::MEDIUM;
  md->user_func=0;
  md->user_data=0;
  md->medium = medium_struct();
  md->medium.epsilon_diag.x=epsilon;
  md->medium.epsilon_diag.y=epsilon;
  md->medium.epsilon_diag.z=epsilon;
  return md;
}

material_type make_user_material(user_material_func user_func,
                                 void *user_data)
{
  material_data *md = (material_data *)malloc(sizeof(*md));
  md->which_subclass=material_data::MATERIAL_USER;
  md->user_func=user_func;
  md->user_data=user_data;
  md->medium = medium_struct();
  return md;
}

// this routine subsumes the content of the old
// 'read_epsilon_file' routine
material_type make_file_material(const char *eps_input_file)
{
  material_data *md = (material_data *)malloc(sizeof(*md));
  md->which_subclass=material_data::MATERIAL_FILE;

  md->epsilon_dims[0] = md->epsilon_dims[1] = md->epsilon_dims[2] = 1;
  if (eps_input_file && eps_input_file[0]) { // file specified
    char *fname = new char[strlen(eps_input_file)+1];
    strcpy(fname, eps_input_file);
    // parse epsilon-input-file as "fname.h5:dataname"
    char *dataname = strrchr(fname, ':');
    if (dataname) *(dataname++) = 0;
    meep::h5file eps_file(fname, meep::h5file::READONLY, false);
    int rank; // ignored since rank < 3 is equivalent to singleton dims
    md->epsilon_data = eps_file.read(dataname, &rank, md->epsilon_dims, 3);
    master_printf("read in %zdx%zdx%zd epsilon-input-file \"%s\"\n",
		  md->epsilon_dims[0],
                  md->epsilon_dims[1],
                  md->epsilon_dims[2],
		  eps_input_file);
  }

  md->medium = medium_struct();

  return md;
}

/******************************************************************************/
/* Helpers from  libctl/utils/geom.c                                          */
/******************************************************************************/
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define MIN(a,b) ((a) < (b) ? (a) : (b))

static void geom_box_intersection(geom_box *bi, const geom_box *b1, const geom_box *b2) {
  bi->low.x = MAX(b1->low.x, b2->low.x);
  bi->low.y = MAX(b1->low.y, b2->low.y);
  bi->low.z = MAX(b1->low.z, b2->low.z);
  bi->high.x = MIN(b1->high.x, b2->high.x);
  bi->high.y = MIN(b1->high.y, b2->high.y);
  bi->high.z = MIN(b1->high.z, b2->high.z);
}

static int geom_boxes_intersect(const geom_box *b1, const geom_box *b2) {
#define BETWEEN(x, low, high) ((x) >= (low) && (x) <= (high))
  /* true if the x, y, and z ranges all intersect. */
  return ((BETWEEN(b1->low.x, b2->low.x, b2->high.x) ||
           BETWEEN(b1->high.x, b2->low.x, b2->high.x) ||
           BETWEEN(b2->low.x, b1->low.x, b1->high.x)) &&
          (BETWEEN(b1->low.y, b2->low.y, b2->high.y) ||
           BETWEEN(b1->high.y, b2->low.y, b2->high.y) ||
           BETWEEN(b2->low.y, b1->low.y, b1->high.y)) &&
          (BETWEEN(b1->low.z, b2->low.z, b2->high.z) ||
           BETWEEN(b1->high.z, b2->low.z, b2->high.z) ||
           BETWEEN(b2->low.z, b1->low.z, b1->high.z)));
}

/******************************************************************************/
/* Fragment Statistics                                                        */
/******************************************************************************/

double fragment_stats::tol = 0;
int fragment_stats::maxeval = 0;
int fragment_stats::resolution = 0;
meep::ndim fragment_stats::dims = meep::D1;

inline static bool is_edge_box(double pt, double half_cell, double box_size, double edge_size) {
  return (pt == -half_cell && edge_size != 0) || pt + box_size > half_cell;
}

static std::vector<geom_box> split_cell_1d(double box_size, vector3 cell_size) {
  double half_box = box_size / 2;
  double half_z = cell_size.z / 2;
  double edge_size_z = fmod(half_z + half_box, box_size);
  std::vector<geom_box> boxes;

  for (double z = -half_z; z < half_z;) {
    double z_increment = is_edge_box(z, half_z, box_size, edge_size_z) ? edge_size_z : box_size;
    vector3 low = {0.0, 0.0, z};
    vector3 high = {0.0, 0.0, z + z_increment};
    geom_box b = {low, high};
    boxes.push_back(b);
    z += z_increment;
  }
  return boxes;
}

static std::vector<geom_box> split_cell_2d(double box_size, vector3 cell_size) {
  double half_box = box_size / 2;
  double half_x = cell_size.x /2;
  double half_y = cell_size.y /2;
  double edge_size_x = fmod(half_x + half_box, box_size);
  double edge_size_y = fmod(half_y + half_box, box_size);
  std::vector<geom_box> boxes;

  for (double x = -half_x; x < half_x;) {
    double x_increment = is_edge_box(x, half_x, box_size, edge_size_x) ? edge_size_x : box_size;
    for (double y = -half_y; y < half_y;) {
      double y_increment = is_edge_box(y, half_y, box_size, edge_size_y) ? edge_size_y : box_size;
      vector3 low = {x, y, 0.0};
      vector3 high = {x + x_increment, y + y_increment, 0.0};
      geom_box b = {low, high};
      boxes.push_back(b);
      y += y_increment;
    }
    x += x_increment;
  }
  return boxes;
}

static std::vector<geom_box> split_cell_3d(double box_size, vector3 cell_size) {
  double half_box = box_size / 2;
  double half_x = cell_size.x /2;
  double half_y = cell_size.y /2;
  double half_z = cell_size.z /2;
  double edge_size_x = fmod(half_x + half_box, box_size);
  double edge_size_y = fmod(half_y + half_box, box_size);
  double edge_size_z = fmod(half_z + half_box, box_size);
  std::vector<geom_box> boxes;

  for (double x = -half_x; x < half_x;) {
    double x_increment = is_edge_box(x, half_x, box_size, edge_size_x) ? edge_size_x : box_size;
    for (double y = -half_y; y < half_y;) {
      double y_increment = is_edge_box(y, half_y, box_size, edge_size_y) ? edge_size_y : box_size;
      for (double z = -half_z; z < half_z;) {
        double z_increment = is_edge_box(z, half_z, box_size, edge_size_z) ? edge_size_z : box_size;
        vector3 low = {x, y, z};
        vector3 high = {x + x_increment, y + y_increment, z + z_increment};
        geom_box b = {low, high};
        boxes.push_back(b);
        z += z_increment;
      }
      y += y_increment;
    }
    x += x_increment;
  }
  return boxes;
}

static std::vector<geom_box> split_cell_cyl(double box_size, vector3 cell_size) {
  double half_box = box_size / 2;
  double half_x = cell_size.x /2;
  double half_z = cell_size.z /2;
  double edge_size_x = fmod(half_x + half_box, box_size);
  double edge_size_z = fmod(half_z + half_box, box_size);
  std::vector<geom_box> boxes;

  for (double x = -half_x; x < half_x;) {
    double x_increment = is_edge_box(x, half_x, box_size, edge_size_x) ? edge_size_x : box_size;
    for (double z = -half_z; z < half_z;) {
      double z_increment = is_edge_box(z, half_z, box_size, edge_size_z) ? edge_size_z : box_size;
      vector3 low = {x, 0.0, z};
      vector3 high = {x + x_increment, 0.0, z + z_increment};
      geom_box b = {low, high};
      boxes.push_back(b);
      z += z_increment;
    }
    x += x_increment;
  }
  return boxes;
}

static std::vector<geom_box> split_cell_into_boxes(meep::grid_volume *gv,
                                                   double box_size,
                                                   vector3 cell_size) {
  switch (gv->dim) {
  case meep::D1:
    return split_cell_1d(box_size, cell_size);
  case meep::D2:
    return split_cell_2d(box_size, cell_size);
  case meep::D3:
    return split_cell_3d(box_size, cell_size);
  case meep::Dcyl:
    return split_cell_cyl(box_size, cell_size);
    break;
  default:
    return std::vector<geom_box>();
  }
}

static size_t get_pixels_in_box(geom_box *b, int empty_pixel=1) {
  int empty_x = b->low.x == b->high.x;
  int empty_y = b->low.y == b->high.y;
  int empty_z = b->low.z == b->high.z;

  double v = ((empty_x ? empty_pixel : b->high.x - b->low.x) *
              (empty_y ? empty_pixel : b->high.y - b->low.y) *
              (empty_z ? empty_pixel : b->high.z - b->low.z));

  return v == 1 ? 0 : (size_t)ceil(v * fragment_stats::resolution);
}

static void center_box(geom_box *b) {
  b->low = vector3_plus(geometry_center, b->low);
  b->high = vector3_plus(geometry_center, b->high);
}

static std::vector<fragment_stats> init_fragments(std::vector<geom_box>& boxes, double tol, int maxeval,
                                                  meep::grid_volume *gv) {
  std::vector<fragment_stats> fragments;
  fragment_stats::tol = tol;
  fragment_stats::maxeval = maxeval;
  fragment_stats::resolution = gv->a;
  fragment_stats::dims = gv->dim;

  for (size_t i = 0; i < boxes.size(); ++i) {
    center_box(&boxes[i]);
    size_t pixels_in_box = get_pixels_in_box(&boxes[i]);
    fragments.push_back(fragment_stats(boxes[i], pixels_in_box));
  }
  return fragments;
}

static void init_libctl(material_type default_mat, bool ensure_per, meep::grid_volume *gv,
                        vector3 cell_size, vector3 cell_center, geometric_object_list *geom) {
  geom_initialize();
  default_material = default_mat;
  ensure_periodicity = ensure_per;
  geometry_center = cell_center;
  dimensions = meep::number_of_directions(gv->dim);
  geometry_lattice.size = cell_size;
  geom_fix_objects0(*geom);
}

std::vector<fragment_stats> compute_fragment_stats(geometric_object_list geom,
                                                   meep::grid_volume *gv,
                                                   vector3 cell_size,
                                                   vector3 cell_center,
                                                   material_type default_mat,
                                                   std::vector<dft_data> dft_data_list,
                                                   double tol,
                                                   int maxeval,
                                                   bool ensure_per,
                                                   double box_size) {

  init_libctl(default_mat, ensure_per, gv, cell_size, cell_center, &geom);
  std::vector<geom_box> boxes = split_cell_into_boxes(gv, box_size, cell_size);
  std::vector<fragment_stats> fragments = init_fragments(boxes, tol, maxeval, gv);

  for (size_t i = 0; i < fragments.size(); ++i) {
    fragments[i].compute_stats(&geom);
    fragments[i].compute_dft_stats(&dft_data_list);
  }
  return fragments;
}

fragment_stats::fragment_stats(geom_box& bx, size_t pixels):
  num_anisotropic_eps_pixels(0),
  num_anisotropic_mu_pixels(0),
  num_nonlinear_pixels(0),
  num_susceptibility_pixels(0),
  num_nonzero_conductivity_pixels(0),
  num_dft_pixels(0),
  num_pixels_in_box(pixels),
  box(bx) {

}

void fragment_stats::update_stats_from_material(material_type mat, size_t pixels) {
  switch (mat->which_subclass) {
  case material_data::MEDIUM: {
    medium_struct *med = &mat->medium;
    count_anisotropic_pixels(med, pixels);
    count_nonlinear_pixels(med, pixels);
    count_susceptibility_pixels(med, pixels);
    count_nonzero_conductivity_pixels(med, pixels);
    break;
  }
  default:
    // Only Medium is supported
    return;
  }
}

void fragment_stats::compute_stats(geometric_object_list *geom) {

  for (int i = 0; i < geom->num_items; ++i) {
    geometric_object *go = &geom->items[i];
    double overlap = box_overlap_with_object(box, *go, tol, maxeval);

    // Count contributions from material of object
    size_t pixels = (size_t)ceil(overlap * num_pixels_in_box);
    if (pixels > 0) {
      material_type mat = (material_type)go->material;
      update_stats_from_material(mat, pixels);
    }

    // Count contributions from default_material
    size_t default_material_pixels = num_pixels_in_box - pixels;
    if (default_material_pixels > 0) {
      update_stats_from_material((material_type)default_material, default_material_pixels);
    }
  }
}

void fragment_stats::count_anisotropic_pixels(medium_struct *med, size_t pixels) {
  size_t eps_offdiag_elements = 0;
  size_t mu_offdiag_elements = 0;

  if (med->epsilon_offdiag.x.re != 0) {
    eps_offdiag_elements++;
  }
  if (med->epsilon_offdiag.y.re != 0) {
    eps_offdiag_elements++;
  }
  if (med->epsilon_offdiag.z.re != 0) {
    eps_offdiag_elements++;
  }
  if (med->mu_offdiag.x.re != 0) {
    mu_offdiag_elements++;
  }
  if (med->mu_offdiag.y.re != 0) {
    mu_offdiag_elements++;
  }
  if (med->mu_offdiag.z.re != 0) {
    mu_offdiag_elements++;
  }

  num_anisotropic_eps_pixels += eps_offdiag_elements * pixels;
  num_anisotropic_mu_pixels += mu_offdiag_elements * pixels;
}

void fragment_stats::count_nonlinear_pixels(medium_struct *med, size_t pixels) {
  size_t nonzero_chi_elements = 0;

  if (med->E_chi2_diag.x != 0) {
    nonzero_chi_elements++;
  }
  if (med->E_chi2_diag.y != 0) {
    nonzero_chi_elements++;
  }
  if (med->E_chi2_diag.z != 0) {
    nonzero_chi_elements++;
  }
  if (med->E_chi3_diag.x != 0) {
    nonzero_chi_elements++;
  }
  if (med->E_chi3_diag.y != 0) {
    nonzero_chi_elements++;
  }
  if (med->E_chi3_diag.z != 0) {
    nonzero_chi_elements++;
  }
  if (med->H_chi2_diag.x != 0) {
    nonzero_chi_elements++;
  }
  if (med->H_chi2_diag.y != 0) {
    nonzero_chi_elements++;
  }
  if (med->H_chi2_diag.z != 0) {
    nonzero_chi_elements++;
  }
  if (med->H_chi3_diag.x != 0) {
    nonzero_chi_elements++;
  }
  if (med->H_chi3_diag.y != 0) {
    nonzero_chi_elements++;
  }
  if (med->H_chi3_diag.z != 0) {
    nonzero_chi_elements++;
  }

  num_nonlinear_pixels += nonzero_chi_elements * pixels;
}

void fragment_stats::count_susceptibility_pixels(medium_struct *med, size_t pixels) {
  num_susceptibility_pixels += med->E_susceptibilities.num_items * pixels;
  num_susceptibility_pixels += med->H_susceptibilities.num_items * pixels;
}

void fragment_stats::count_nonzero_conductivity_pixels(medium_struct *med, size_t pixels) {
  size_t nonzero_conductivity_elements = 0;

  if (med->D_conductivity_diag.x != 0) {
    nonzero_conductivity_elements++;
  }
  if (med->D_conductivity_diag.y != 0) {
    nonzero_conductivity_elements++;
  }
  if (med->D_conductivity_diag.z != 0) {
    nonzero_conductivity_elements++;
  }
  if (med->B_conductivity_diag.x != 0) {
    nonzero_conductivity_elements++;
  }
  if (med->B_conductivity_diag.y != 0) {
    nonzero_conductivity_elements++;
  }
  if (med->B_conductivity_diag.z != 0) {
    nonzero_conductivity_elements++;
  }

  num_nonzero_conductivity_pixels += nonzero_conductivity_elements * pixels;
}

void fragment_stats::compute_dft_stats(std::vector<dft_data> *dft_data_list) {

  for (size_t i = 0; i < dft_data_list->size(); ++i) {
    for (size_t j = 0; j < (*dft_data_list)[i].vols.size(); ++j) {
      geom_box dft_box = gv2box((*dft_data_list)[i].vols[j]);

      if (geom_boxes_intersect(&dft_box, &box)) {
        geom_box overlap_box;
        geom_box_intersection(&overlap_box, &dft_box, &box);

        // Note: Since geom_boxes_intersect returns true if two planes share a line or two volumes
        // share a line or plane, there are cases where some pixels are counted multiple times.
        size_t overlap_pixels = get_pixels_in_box(&overlap_box, 2);
        num_dft_pixels += overlap_pixels * (*dft_data_list)[i].num_freqs * (*dft_data_list)[i].num_components;
      }
    }
  }
}

void fragment_stats::print_stats() {
  master_printf("Fragment stats\n");
  master_printf("  num_anisotropic_eps_pixels: %zd\n", num_anisotropic_eps_pixels);
  master_printf("  num_anisotropic_mu_pixels: %zd\n", num_anisotropic_mu_pixels);
  master_printf("  num_nonlinear_pixels: %zd\n", num_nonlinear_pixels);
  master_printf("  num_susceptibility_pixels: %zd\n", num_susceptibility_pixels);
  master_printf("  num_nonzero_conductivity_pixels: %zd\n", num_nonzero_conductivity_pixels);
  master_printf("  num_dft_pixels: %zd\n", num_dft_pixels);
  master_printf("  num_pixels_in_box: %zd\n", num_pixels_in_box);
  master_printf("  box.low:  {%f, %f, %f}\n", box.low.x, box.low.y, box.low.z);
  master_printf("  box.high: {%f, %f, %f}\n\n", box.high.x, box.high.y, box.high.z);
}

dft_data::dft_data(int freqs, int components, std::vector<meep::volume> volumes):
  num_freqs(freqs),
  num_components(components),
  vols(volumes) {

}
} // namespace meep_geom
