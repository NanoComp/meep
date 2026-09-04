/* Copyright (C) 2005-2026 Massachusetts Institute of Technology
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

#include <algorithm>
#include <set>
#include <vector>
#include "meepgeom.hpp"
#include "meep_internals.hpp"

namespace meep_geom {

double _adjdbg_sum[2] = {0, 0};
long _adjdbg_n[2] = {0, 0};


#define master_printf meep::master_printf

/***************************************************************/
/* global variables for default material                       */
/***************************************************************/
material_data vacuum_material_data;
material_type vacuum = &vacuum_material_data;

void set_default_material(material_type _default_material) {
  if (default_material != NULL) {
    if (default_material == _default_material) return;
    material_free((material_type)default_material);
    default_material = NULL;
  }

  if (_default_material != NULL) {
    material_type new_material = new material_data();
    new_material->copy_from(*_default_material);
    default_material = (void *)new_material;
  }
}

void unset_default_material(void) {
  if (default_material != NULL) {
    material_free((material_type)default_material);
    default_material = NULL;
  }
}

bool susceptibility_equal(const susceptibility &s1, const susceptibility &s2) {
  return (vector3_equal(s1.sigma_diag, s2.sigma_diag) &&
          vector3_equal(s1.sigma_offdiag, s2.sigma_offdiag) && vector3_equal(s1.bias, s2.bias) &&
          s1.frequency == s2.frequency && s1.gamma == s2.gamma && s1.alpha == s2.alpha &&
          s1.noise_amp == s2.noise_amp && s1.drude == s2.drude &&
          s1.saturated_gyrotropy == s2.saturated_gyrotropy && s1.is_file == s2.is_file);
}

bool susceptibility_list_equal(const susceptibility_list &s1, const susceptibility_list &s2) {
  if (s1.size() != s2.size()) return false;
  for (size_t i = 0; i < s1.size(); ++i) {
    if (!susceptibility_equal(s1[i], s2[i])) return false;
  }
  return true;
}

bool medium_struct_equal(const medium_struct *m1, const medium_struct *m2) {
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

bool material_grid_equal(const material_data *m1, const material_data *m2) {
  // a rigorous method for comapring material grids
  int n1, n2;
  n1 = m1->grid_size.x * m1->grid_size.y * m1->grid_size.z;
  n2 = m2->grid_size.x * m2->grid_size.y * m2->grid_size.z;
  if (n1 != n2) return false;
  for (int i = 0; i < n1; i++)
    if (m1->epsilon_data[i] != m2->epsilon_data[i]) return false;

  return (medium_struct_equal(&(m1->medium), &(m2->medium)) &&
          medium_struct_equal(&(m1->medium_1), &(m2->medium_1)) &&
          medium_struct_equal(&(m1->medium_2), &(m2->medium_2)));
}

// garbage collection for material structures: called to deallocate memory
// allocated for susceptibilities in user-defined materials.
// TODO
void material_gc(material_type m) {
  if (!m || m->which_subclass != material_data::MATERIAL_USER) return;
  m->medium.E_susceptibilities.clear();
  m->medium.H_susceptibilities.clear();
  m->medium_1.E_susceptibilities.clear();
  m->medium_1.H_susceptibilities.clear();
  m->medium_2.E_susceptibilities.clear();
  m->medium_2.H_susceptibilities.clear();
}

void material_free(material_type m) {
  if (!m) return;

  m->medium.E_susceptibilities.clear();
  m->medium.H_susceptibilities.clear();
  m->medium_1.E_susceptibilities.clear();
  m->medium_1.H_susceptibilities.clear();
  m->medium_2.E_susceptibilities.clear();
  m->medium_2.H_susceptibilities.clear();

  // NOTE: We do not delete the user_data field here since it is an opaque/void
  // object so will assume that the caller keeps track of its lifetime.
  delete[] m->epsilon_data;
  m->epsilon_data = NULL;

  delete[] m->weights;
  m->weights = NULL;
  delete m;
}

bool material_type_equal(const material_type m1, const material_type m2) {
  if (m1 == m2) return true;
  if (m1->which_subclass != m2->which_subclass) return false;
  switch (m1->which_subclass) {
    case material_data::MATERIAL_FILE:
    case material_data::PERFECT_METAL: return true;
    case material_data::MATERIAL_USER:
      return m1->user_func == m2->user_func && m1->user_data == m2->user_data;
    case material_data::MATERIAL_GRID:
    case material_data::MEDIUM: return medium_struct_equal(&(m1->medium), &(m2->medium));
    default: return false;
  }
}

/***************************************************************/
/***************************************************************/
/***************************************************************/

/* rotate A by a unitary (real) rotation matrix R:
      RAR = transpose(R) * A * R
*/
void sym_matrix_rotate(symm_matrix *RAR, const symm_matrix *A_, const double R[3][3]) {
  int i, j;
  double A[3][3], AR[3][3];
  A[0][0] = A_->m00;
  A[1][1] = A_->m11;
  A[2][2] = A_->m22;
  A[0][1] = A[1][0] = A_->m01;
  A[0][2] = A[2][0] = A_->m02;
  A[1][2] = A[2][1] = A_->m12;
  for (i = 0; i < 3; ++i)
    for (j = 0; j < 3; ++j)
      AR[i][j] = A[i][0] * R[0][j] + A[i][1] * R[1][j] + A[i][2] * R[2][j];
  for (i = 0; i < 3; ++i)
    for (j = i; j < 3; ++j)
      A[i][j] = R[0][i] * AR[0][j] + R[1][i] * AR[1][j] + R[2][i] * AR[2][j];
  RAR->m00 = A[0][0];
  RAR->m11 = A[1][1];
  RAR->m22 = A[2][2];
  RAR->m01 = A[0][1];
  RAR->m02 = A[0][2];
  RAR->m12 = A[1][2];
}

/* Set Vinv = inverse of V, where both V and Vinv are real-symmetric matrices.*/
void sym_matrix_invert(symm_matrix *Vinv, const symm_matrix *V) {
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
    detinv = m00 * m11 * m22 - m02 * m11 * m02 + 2.0 * m01 * m12 * m02 - m01 * m01 * m22 -
             m12 * m12 * m00;

    if (detinv == 0.0) meep::abort("singular 3x3 matrix");

    detinv = 1.0 / detinv;

    Vinv->m00 = detinv * (m11 * m22 - m12 * m12);
    Vinv->m11 = detinv * (m00 * m22 - m02 * m02);
    Vinv->m22 = detinv * (m11 * m00 - m01 * m01);

    Vinv->m02 = detinv * (m01 * m12 - m11 * m02);
    Vinv->m01 = detinv * (m12 * m02 - m01 * m22);
    Vinv->m12 = detinv * (m01 * m02 - m00 * m12);
  }
}

/* Returns whether or not V is positive-definite. */
int sym_matrix_positive_definite(symm_matrix *V) {
  double det2, det3;
  double m00 = V->m00, m11 = V->m11, m22 = V->m22;

#if defined(WITH_HERMITIAN_EPSILON)
  scalar_complex m01 = V->m01, m02 = V->m02, m12 = V->m12;

  det2 = m00 * m11 - CSCALAR_NORMSQR(m01);
  det3 = det2 * m22 - m11 * CSCALAR_NORMSQR(m02) - CSCALAR_NORMSQR(m12) * m00 +
         2.0 * ((m01.re * m12.re - m01.im * m12.im) * m02.re +
                (m01.re * m12.im + m01.im * m12.re) * m02.im);
#else  /* real matrix */
  double m01 = V->m01, m02 = V->m02, m12 = V->m12;

  det2 = m00 * m11 - m01 * m01;
  det3 = det2 * m22 - m02 * m11 * m02 + 2.0 * m01 * m12 * m02 - m12 * m12 * m00;
#endif /* real matrix */

  return (m00 > 0.0 && det2 > 0.0 && det3 > 0.0);
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
static meep::ndim dim = meep::D3;
void set_dimensions(int dims) {
  if (dims == CYLINDRICAL) { dim = meep::Dcyl; }
  else { dim = meep::ndim(dims - 1); }
}

vector3 vec_to_vector3(const meep::vec &pt) {
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

meep::vec vector3_to_vec(const vector3 v3) {
  switch (dim) {
    case meep::D1: return meep::vec(v3.z);
    case meep::D2: return meep::vec(v3.x, v3.y);
    case meep::D3: return meep::vec(v3.x, v3.y, v3.z);
    case meep::Dcyl: return meep::veccyl(v3.x, v3.z);
    default: meep::abort("unknown dimensionality in vector3_to_vec");
  }
}

geom_box gv2box(const meep::volume &v) {
  geom_box box;
  box.low = vec_to_vector3(v.get_min_corner());
  box.high = vec_to_vector3(v.get_max_corner());
  return box;
}

bool is_material_grid(material_type mt) {
  return (mt->which_subclass == material_data::MATERIAL_GRID);
}
bool is_material_grid(void *md) { return is_material_grid((material_type)md); }

bool is_variable(material_type mt) {
  return (mt->which_subclass == material_data::MATERIAL_USER) ||
         (mt->which_subclass == material_data::MATERIAL_GRID);
}
bool is_variable(void *md) { return is_variable((material_type)md); }

bool is_file(material_type md) { return (md->which_subclass == material_data::MATERIAL_FILE); }
bool is_file(void *md) { return is_file((material_type)md); }

bool is_medium(material_type md, medium_struct **m) {
  if (md->which_subclass == material_data::MEDIUM) {
    *m = &(md->medium);
    return true;
  };
  return false;
}

bool is_medium(void *md, medium_struct **m) { return is_medium((material_type)md, m); }

bool is_metal(meep::field_type ft, const material_type *material) {
  material_data *md = *material;
  if (ft == meep::E_stuff) switch (md->which_subclass) {
      case material_data::MEDIUM:
      case material_data::MATERIAL_GRID:
        return (md->medium.epsilon_diag.x < 0 || md->medium.epsilon_diag.y < 0 ||
                md->medium.epsilon_diag.z < 0);
      case material_data::PERFECT_METAL: return true;
      default: meep::abort("unknown material type"); return false;
    }
  else
    switch (md->which_subclass) {
      case material_data::MEDIUM:
      case material_data::MATERIAL_GRID:
        return (md->medium.mu_diag.x < 0 || md->medium.mu_diag.y < 0 || md->medium.mu_diag.z < 0);
      case material_data::PERFECT_METAL:
        return false; // is an electric conductor, but not a magnetic conductor
      default: meep::abort("unknown material type"); return false;
    }
}

bool has_offdiag(const medium_struct *material) {
  if ((material->epsilon_offdiag.x.re != 0) || /* account for offdiagonal components */
      (material->epsilon_offdiag.y.re != 0) || (material->epsilon_offdiag.z.re != 0) ||
      (material->epsilon_offdiag.x.im != 0) || (material->epsilon_offdiag.y.im != 0) ||
      (material->epsilon_offdiag.z.im != 0))
    return true;
  else
    return false;
}

// computes the vector-Jacobian product of the gradient of the matgrid_val function v
// with the Jacobian of the to_geom_box_coords function for geometric_object o
vector3 to_geom_object_coords_VJP(vector3 v, const geometric_object *o) {
  if (!o) { meep::abort("must pass a geometric_object to to_geom_object_coords_VJP.\n"); }

  switch (o->which_subclass) {
    default: {
      vector3 po = {0, 0, 0};
      return po;
    }
    case geometric_object::SPHERE: {
      number radius = o->subclass.sphere_data->radius;
      return vector3_scale(0.5 / radius, v);
    }
    /* case geometric_object::CYLINDER:
       NOT YET IMPLEMENTED */
    case geometric_object::BLOCK: {
      vector3 size = o->subclass.block_data->size;
      if (size.x != 0.0) v.x /= size.x;
      if (size.y != 0.0) v.y /= size.y;
      if (size.z != 0.0) v.z /= size.z;
      return matrix3x3_transpose_vector3_mult(o->subclass.block_data->projection_matrix, v);
    }
      /* case geometric_object::PRISM:
         NOT YET IMPLEMENTED */
  }
}

meep::vec material_grid_grad(vector3 p, material_data *md, const geometric_object *o) {
  if (!is_material_grid(md)) { meep::abort("Invalid material grid detected.\n"); }

  meep::vec gradient(zero_vec(dim));
  double *data = md->weights;
  int nx = md->grid_size.x;
  int ny = md->grid_size.y;
  int nz = md->grid_size.z;
  double rx = p.x;
  double ry = p.y;
  double rz = p.z;
  int stride = 1;
  int x1, y1, z1, x2, y2, z2;
  double dx, dy, dz;
  bool signflip_dx = false, signflip_dy = false, signflip_dz = false;

  meep::map_coordinates(rx, ry, rz, nx, ny, nz, x1, y1, z1, x2, y2, z2, dx, dy, dz,
                        false /* do_fabs */);

  if (dx != fabs(dx)) {
    dx = fabs(dx);
    signflip_dx = true;
  }
  if (dy != fabs(dy)) {
    dy = fabs(dy);
    signflip_dy = true;
  }
  if (dz != fabs(dz)) {
    dz = fabs(dz);
    signflip_dz = true;
  }

  /* define a macro to give us data(x,y,z) on the grid,
     in row-major order: */
#define D(x, y, z) (data[(((x)*ny + (y)) * nz + (z)) * stride])

  double du_dx =
      (signflip_dx ? -1.0 : 1.0) *
      (((-D(x1, y1, z1) + D(x2, y1, z1)) * (1.0 - dy) + (-D(x1, y2, z1) + D(x2, y2, z1)) * dy) *
           (1.0 - dz) +
       ((-D(x1, y1, z2) + D(x2, y1, z2)) * (1.0 - dy) + (-D(x1, y2, z2) + D(x2, y2, z2)) * dy) *
           dz);
  double du_dy = (signflip_dy ? -1.0 : 1.0) * ((-(D(x1, y1, z1) * (1.0 - dx) + D(x2, y1, z1) * dx) +
                                                (D(x1, y2, z1) * (1.0 - dx) + D(x2, y2, z1) * dx)) *
                                                   (1.0 - dz) +
                                               (-(D(x1, y1, z2) * (1.0 - dx) + D(x2, y1, z2) * dx) +
                                                (D(x1, y2, z2) * (1.0 - dx) + D(x2, y2, z2) * dx)) *
                                                   dz);
  double du_dz = (signflip_dz ? -1.0 : 1.0) *
                 (-((D(x1, y1, z1) * (1.0 - dx) + D(x2, y1, z1) * dx) * (1.0 - dy) +
                    (D(x1, y2, z1) * (1.0 - dx) + D(x2, y2, z1) * dx) * dy) +
                  ((D(x1, y1, z2) * (1.0 - dx) + D(x2, y1, z2) * dx) * (1.0 - dy) +
                   (D(x1, y2, z2) * (1.0 - dx) + D(x2, y2, z2) * dx) * dy));

#undef D

  // [du_dx,du_dy,du_dz] is the gradient ∇u with respect to the transformed coordinate
  // r1 of the matgrid_val function but what we want is the gradient of u(g(r2)) with
  // respect to r2 where g(r2) is the to_geom_object_coords function (in libctl/utils/geom.c).
  // computing this quantity involves using the chain rule and thus the vector-Jacobian product
  // ∇u J where J is the Jacobian matrix of g.
  vector3 grad_u;
  grad_u.x = du_dx * nx;
  grad_u.y = du_dy * ny;
  grad_u.z = du_dz * nz;
  if (o != NULL) {
    vector3 grad_u_J = to_geom_object_coords_VJP(grad_u, o);
    gradient.set_direction(meep::X, grad_u_J.x);
    gradient.set_direction(meep::Y, grad_u_J.y);
    gradient.set_direction(meep::Z, grad_u_J.z);
  }
  else {
    gradient.set_direction(meep::X,
                           geometry_lattice.size.x == 0 ? 0 : grad_u.x / geometry_lattice.size.x);
    gradient.set_direction(meep::Y,
                           geometry_lattice.size.y == 0 ? 0 : grad_u.y / geometry_lattice.size.y);
    gradient.set_direction(meep::Z,
                           geometry_lattice.size.z == 0 ? 0 : grad_u.z / geometry_lattice.size.z);
  }

  return gradient;
}

void map_lattice_coordinates(double &px, double &py, double &pz) {
  px = geometry_lattice.size.x == 0 ? 0 : 0.5 + (px - geometry_center.x) / geometry_lattice.size.x;
  py = geometry_lattice.size.y == 0 ? 0 : 0.5 + (py - geometry_center.y) / geometry_lattice.size.y;
  pz = geometry_lattice.size.z == 0 ? 0 : 0.5 + (pz - geometry_center.z) / geometry_lattice.size.z;
}

meep::vec matgrid_grad(vector3 p, geom_box_tree tp, int oi, material_data *md) {
  if (md->material_grid_kinds == material_data::U_MIN ||
      md->material_grid_kinds == material_data::U_PROD)
    meep::abort("%s:%i:matgrid_grad does not support overlapping grids with U_MIN or U_PROD\n",
                __FILE__, __LINE__);

  meep::vec gradient(zero_vec(dim));
  int matgrid_val_count = 0;

  // iterate through object tree at current point
  if (tp) {
    do {
      gradient +=
          material_grid_grad(to_geom_box_coords(p, &tp->objects[oi]),
                             (material_data *)tp->objects[oi].o->material, tp->objects[oi].o);
      if (md->material_grid_kinds == material_data::U_DEFAULT) break;
      ++matgrid_val_count;
      tp = geom_tree_search_next(p, tp, &oi);
    } while (tp && is_material_grid((material_data *)tp->objects[oi].o->material));
  }
  // perhaps there is no object tree and the default material is a material grid
  if (!tp && is_material_grid(default_material)) {
    map_lattice_coordinates(p.x, p.y, p.z);
    gradient =
        material_grid_grad(p, (material_data *)default_material, NULL /* geometric_object *o */);
    ++matgrid_val_count;
  }

  if (md->material_grid_kinds == material_data::U_MEAN)
    gradient = gradient * 1.0 / matgrid_val_count;

  return gradient;
}

double material_grid_val(vector3 p, material_data *md) {
  // given the relative location, p, interpolate the material grid point.

  if (!is_material_grid(md)) { meep::abort("Invalid material grid detected.\n"); }
  return meep::linear_interpolate(p.x, p.y, p.z, md->weights, md->grid_size.x, md->grid_size.y,
                                  md->grid_size.z, 1);
}

static double tanh_projection(double u, double beta, double eta) {
  if (beta == 0) return u;
  if (u == eta) return 0.5; // avoid NaN when beta is Inf
  double tanh_beta_eta = tanh(beta * eta);
  return (tanh_beta_eta + tanh(beta * (u - eta))) / (tanh_beta_eta + tanh(beta * (1 - eta)));
}

double matgrid_val(vector3 p, geom_box_tree tp, int oi, material_data *md) {
  double uprod = 1.0, umin = 1.0, usum = 0.0, udefault = 0.0, u;
  int matgrid_val_count = 0;

  // iterate through object tree at current point
  if (tp) {
    do {
      u = material_grid_val(to_geom_box_coords(p, &tp->objects[oi]),
                            (material_data *)tp->objects[oi].o->material);
      if (md->material_grid_kinds == material_data::U_DEFAULT) {
        udefault = u;
        break;
      }
      if (u < umin) umin = u;
      uprod *= u;
      usum += u;
      ++matgrid_val_count;
      tp = geom_tree_search_next(p, tp, &oi);
    } while (tp && is_material_grid((material_data *)tp->objects[oi].o->material));
  }
  // perhaps there is no object tree and the default material is a material grid
  if (!tp && is_material_grid(default_material)) {
    map_lattice_coordinates(p.x, p.y, p.z);
    u = material_grid_val(p, (material_data *)default_material);
    if (matgrid_val_count == 0) udefault = u;
    if (u < umin) umin = u;
    uprod *= u;
    usum += u;
    ++matgrid_val_count;
  }

  return (md->material_grid_kinds == material_data::U_MIN
              ? umin
              : (md->material_grid_kinds == material_data::U_PROD
                     ? uprod
                     : (md->material_grid_kinds == material_data::U_MEAN ? usum / matgrid_val_count
                                                                         : udefault)));
}
static void cinterp_tensors(vector3 diag_in_1, cvector3 offdiag_in_1, vector3 diag_in_2,
                            cvector3 offdiag_in_2, vector3 *diag_out, cvector3 *offdiag_out,
                            double u) {
  /* convenience routine to interpolate material tensors with real and imaginary components */
  diag_out->x = diag_in_1.x + u * (diag_in_2.x - diag_in_1.x);
  diag_out->y = diag_in_1.y + u * (diag_in_2.y - diag_in_1.y);
  diag_out->z = diag_in_1.z + u * (diag_in_2.z - diag_in_1.z);
  offdiag_out->x.re = offdiag_in_1.x.re + u * (offdiag_in_2.x.re - offdiag_in_1.x.re);
  offdiag_out->x.im = offdiag_in_1.x.im + u * (offdiag_in_2.x.im - offdiag_in_1.x.im);
  offdiag_out->y.re = offdiag_in_1.y.re + u * (offdiag_in_2.y.re - offdiag_in_1.y.re);
  offdiag_out->y.im = offdiag_in_1.y.im + u * (offdiag_in_2.y.im - offdiag_in_1.y.im);
  offdiag_out->z.re = offdiag_in_1.z.re + u * (offdiag_in_2.z.re - offdiag_in_1.z.re);
  offdiag_out->z.im = offdiag_in_1.z.im + u * (offdiag_in_2.z.im - offdiag_in_1.z.im);
}

static void interp_tensors(vector3 diag_in_1, vector3 offdiag_in_1, vector3 diag_in_2,
                           vector3 offdiag_in_2, vector3 *diag_out, vector3 *offdiag_out,
                           double u) {
  /* convenience routine to interpolate material tensors with all real components */
  diag_out->x = diag_in_1.x + u * (diag_in_2.x - diag_in_1.x);
  diag_out->y = diag_in_1.y + u * (diag_in_2.y - diag_in_1.y);
  diag_out->z = diag_in_1.z + u * (diag_in_2.z - diag_in_1.z);
  offdiag_out->x = offdiag_in_1.x + u * (offdiag_in_2.x - offdiag_in_1.x);
  offdiag_out->y = offdiag_in_1.y + u * (offdiag_in_2.y - offdiag_in_1.y);
  offdiag_out->z = offdiag_in_1.z + u * (offdiag_in_2.z - offdiag_in_1.z);
}
// return material of the point p from the material grid
void epsilon_material_grid(material_data *md, double u) {
  // NOTE: assume p lies on normalized grid within (0,1)

  if (!(md->weights)) meep::abort("material params were not initialized!");

  medium_struct *mm = &(md->medium);
  medium_struct *m1 = &(md->medium_1);
  medium_struct *m2 = &(md->medium_2);

  // Linearly interpolate dc epsilon values
  cinterp_tensors(m1->epsilon_diag, m1->epsilon_offdiag, m2->epsilon_diag, m2->epsilon_offdiag,
                  &mm->epsilon_diag, &mm->epsilon_offdiag, u);

  // Interpolate resonant strength from d.p.
  //
  // mm->E_susceptibilities is the concatenation of the two constituents' pole
  // lists, assembled when the grid is created (python/typemap_utils.cpp).  The
  // indexing below relies on that and never resizes, so check it rather than
  // writing past the end for a grid built through some other path.
  if (mm->E_susceptibilities.size() !=
      m1->E_susceptibilities.size() + m2->E_susceptibilities.size())
    meep::abort("material grid susceptibility list was not initialized!");

  vector3 zero_vec;
  zero_vec.x = zero_vec.y = zero_vec.z = 0;
  for (size_t i = 0; i < m1->E_susceptibilities.size(); i++) {
    // iterate through medium1 sus list first
    interp_tensors(zero_vec, zero_vec, m1->E_susceptibilities[i].sigma_diag,
                   m1->E_susceptibilities[i].sigma_offdiag, &mm->E_susceptibilities[i].sigma_diag,
                   &mm->E_susceptibilities[i].sigma_offdiag, (1 - u));
  }
  for (size_t i = 0; i < m2->E_susceptibilities.size(); i++) {
    // iterate through medium2 sus list next
    size_t j = i + m1->E_susceptibilities.size();
    interp_tensors(zero_vec, zero_vec, m2->E_susceptibilities[i].sigma_diag,
                   m2->E_susceptibilities[i].sigma_offdiag, &mm->E_susceptibilities[j].sigma_diag,
                   &mm->E_susceptibilities[j].sigma_offdiag, u);
  }

  // Linearly interpolate electric conductivity
  vector3 zero_offdiag;
  interp_tensors(m1->D_conductivity_diag, zero_vec, m2->D_conductivity_diag, zero_vec,
                 &mm->D_conductivity_diag, &zero_offdiag, u);

  // Add damping factor if we have dispersion.
  // This prevents instabilities when interpolating between sus. profiles.
  if ((m1->E_susceptibilities.size() + m2->E_susceptibilities.size()) > 0.0) {
    // calculate mean harmonic frequency
    double omega_mean = 0;
    for (const susceptibility &m1_sus : m1->E_susceptibilities) {
      omega_mean += m1_sus.frequency;
    }
    for (const susceptibility &m2_sus : m2->E_susceptibilities) {
      omega_mean += m2_sus.frequency;
    }
    omega_mean = omega_mean / (m1->E_susceptibilities.size() + m2->E_susceptibilities.size());

    // assign interpolated, nondimensionalized conductivity term
    // TODO: dampen the lorentzians to improve stability
    // mm->D_conductivity_diag.x = mm->D_conductivity_diag.y = mm->D_conductivity_diag.z = u*(1-u) *
    // omega_mean;
    md->trivial = false;
  }
  double fake_damping = u * (1 - u) * (md->damping);
  mm->D_conductivity_diag.x += fake_damping;
  mm->D_conductivity_diag.y += fake_damping;
  mm->D_conductivity_diag.z += fake_damping;

  // set the trivial flag
  if (md->damping != 0) md->trivial = false;
}

// return material of the point p from the file (assumed already read)
void epsilon_file_material(material_data *md, vector3 p) {
  set_default_material(md);

  if (md->which_subclass != material_data::MATERIAL_FILE)
    meep::abort("epsilon-input-file only works with a type=file default-material");

  if (!(md->epsilon_data)) return;
  medium_struct *mm = &(md->medium);
  double rx =
      geometry_lattice.size.x == 0 ? 0 : 0.5 + (p.x - geometry_center.x) / geometry_lattice.size.x;
  double ry =
      geometry_lattice.size.y == 0 ? 0 : 0.5 + (p.y - geometry_center.y) / geometry_lattice.size.y;
  double rz =
      geometry_lattice.size.z == 0 ? 0 : 0.5 + (p.z - geometry_center.z) / geometry_lattice.size.z;
  mm->epsilon_diag.x = mm->epsilon_diag.y = mm->epsilon_diag.z =
      meep::linear_interpolate(rx, ry, rz, md->epsilon_data, md->epsilon_dims[0],
                               md->epsilon_dims[1], md->epsilon_dims[2], 1);
  mm->epsilon_offdiag.x.re = mm->epsilon_offdiag.y.re = mm->epsilon_offdiag.z.re = 0;
}

/***********************************************************************/

geom_epsilon::geom_epsilon(geometric_object_list g, material_type_list mlist,
                           const meep::volume &v) {
  // Copy the geometry
  int length = g.num_items;
  geometry.num_items = length;
  geometry.items = new geometric_object[length];
  has_user_materials = default_material && ((material_type)default_material)->which_subclass ==
                                               material_data::MATERIAL_USER;
  for (int i = 0; i < length; i++) {
    geometric_object_copy(&g.items[i], &geometry.items[i]);
    geometry.items[i].material = new material_data();
    static_cast<material_data *>(geometry.items[i].material)
        ->copy_from(*(material_data *)(g.items[i].material));
    if (static_cast<material_data *>(g.items[i].material)->which_subclass ==
        material_data::MATERIAL_USER)
      has_user_materials = true;
  }

  extra_materials = mlist;
  current_pol = NULL;

  FOR_DIRECTIONS(d) FOR_SIDES(b) { cond[d][b].prof = NULL; }

  if (meep::am_master()) {
    int num_print = meep::verbosity > 2
                        ? geometry.num_items
                        : std::min(geometry.num_items, meep::verbosity > 0 ? 10 : 0);
    for (int i = 0; i < geometry.num_items; ++i) {

      if (i < num_print) display_geometric_object_info(5, geometry.items[i]);

      medium_struct *mm;
      if (is_medium(geometry.items[i].material, &mm)) {
        mm->check_offdiag_im_zero_or_abort();
        if (i < num_print)
          master_printf("%*sdielectric constant epsilon diagonal "
                        "= (%g,%g,%g)\n",
                        5 + 5, "", mm->epsilon_diag.x, mm->epsilon_diag.y, mm->epsilon_diag.z);
      }
    }
    if (num_print < geometry.num_items && meep::verbosity > 0)
      master_printf("%*s...(+ %d objects not shown)...\n", 5, "", geometry.num_items - num_print);
  }
  geom_fix_object_list(geometry);
  geom_box box = gv2box(v);
  geometry_tree = create_geom_box_tree0(geometry, box);
  if (meep::verbosity > 2 && meep::am_master()) {
    master_printf("Geometric-object bounding-box tree:\n");
    display_geom_box_tree(5, geometry_tree);

    int tree_depth, tree_nobjects;
    geom_box_tree_stats(geometry_tree, &tree_depth, &tree_nobjects);
    master_printf("Geometric object tree has depth %d "
                  "and %d object nodes (vs. %d actual objects)\n",
                  tree_depth, tree_nobjects, geometry.num_items);
  }

  restricted_tree = geometry_tree;
}

// copy constructor
geom_epsilon::geom_epsilon(const geom_epsilon &geps1) {
  // Copy the geometry
  int length = geps1.geometry.num_items;
  geometry.num_items = length;
  geometry.items = new geometric_object[length];
  has_user_materials = geps1.has_user_materials;
  for (int i = 0; i < length; i++) {
    geometric_object_copy(&geps1.geometry.items[i], &geometry.items[i]);
    geometry.items[i].material = new material_data();
    static_cast<material_data *>(geometry.items[i].material)
        ->copy_from(*(material_data *)(geps1.geometry.items[i].material));
  }

  geometry_tree = geps1.geometry_tree;
  restricted_tree = geps1.restricted_tree;
  extra_materials = geps1.extra_materials;
  current_pol = NULL;

  // Parameters stashed for the gradient calculation.  Silently reverting these
  // to their defaults in a copy would not fail, it would just quietly compute
  // against a different operator than the original was built with.
  u_p = geps1.u_p;
  tol = geps1.tol;
  maxeval = geps1.maxeval;
  dt = geps1.dt;

  FOR_DIRECTIONS(d) FOR_SIDES(b) { cond[d][b].prof = geps1.cond[d][b].prof; }
}
geom_epsilon::~geom_epsilon() {
  int length = geometry.num_items;
  for (int i = 0; i < length; i++) {
    material_free((material_type)geometry.items[i].material);
    geometric_object_destroy(geometry.items[i]);
  }
  delete[] geometry.items;
  unset_volume();
  destroy_geom_box_tree(geometry_tree);
  FOR_DIRECTIONS(d) FOR_SIDES(b) {
    if (cond[d][b].prof) delete[] cond[d][b].prof;
  }
}

void geom_epsilon::set_cond_profile(meep::direction dir, meep::boundary_side side, double L,
                                    double dx, double (*P)(int, double *, void *), void *data,
                                    double R) {
  if (cond[dir][side].prof) delete[] cond[dir][side].prof;

  int N = int(L / dx + 0.5);
  cond[dir][side].L = L;
  cond[dir][side].N = N;
  double *prof = cond[dir][side].prof = new double[N + 1];

  double umin = 0, umax = 1, esterr;
  int errflag;
  double prof_int =
      adaptive_integration(P, &umin, &umax, 1, data, 1e-9, 1e-4, 50000, &esterr, &errflag);

  double prefac = (-log(R)) / (4 * L * prof_int);
  for (int i = 0; i <= N; ++i) {
    double u = double(i) / N;
    prof[i] = prefac * P(1, &u, data);
  }
}

void geom_epsilon::unset_volume(void) {
  if (restricted_tree != geometry_tree) {
    destroy_geom_box_tree(restricted_tree);
    restricted_tree = geometry_tree;
  }
}

void geom_epsilon::set_volume(const meep::volume &v) {
  unset_volume();

  geom_box box = gv2box(v);
  if (!restricted_tree) restricted_tree = create_geom_box_tree0(geometry, box);
}

static void material_epsmu(meep::field_type ft, material_type material, symm_matrix *epsmu,
                           symm_matrix *epsmu_inv) {

  material_data *md = material;
  if (ft == meep::E_stuff || ft == meep::D_stuff) switch (md->which_subclass) {

      case material_data::MEDIUM:
      case material_data::MATERIAL_FILE:
      case material_data::MATERIAL_USER:
      case material_data::MATERIAL_GRID:
        epsmu->m00 = md->medium.epsilon_diag.x;
        epsmu->m11 = md->medium.epsilon_diag.y;
        epsmu->m22 = md->medium.epsilon_diag.z;
        epsmu->m01 = md->medium.epsilon_offdiag.x.re;
        epsmu->m02 = md->medium.epsilon_offdiag.y.re;
        epsmu->m12 = md->medium.epsilon_offdiag.z.re;
        sym_matrix_invert(epsmu_inv, epsmu);
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

      default: meep::abort("unknown material type");
    }
  else
    switch (md->which_subclass) {
      case material_data::MEDIUM:
      case material_data::MATERIAL_FILE:
      case material_data::MATERIAL_USER:
      case material_data::MATERIAL_GRID:
        epsmu->m00 = md->medium.mu_diag.x;
        epsmu->m11 = md->medium.mu_diag.y;
        epsmu->m22 = md->medium.mu_diag.z;
        epsmu->m01 = md->medium.mu_offdiag.x.re;
        epsmu->m02 = md->medium.mu_offdiag.y.re;
        epsmu->m12 = md->medium.mu_offdiag.z.re;
        sym_matrix_invert(epsmu_inv, epsmu);
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

      default: meep::abort("unknown material type");
    }
}

// the goal of this routine is to fill in the 'medium' field
// within the material structure as appropriate for the
// material properties at r.
void geom_epsilon::get_material_pt(material_type &material, const meep::vec &r) {
  vector3 p = vec_to_vector3(r);
  boolean inobject;
  material =
      (material_type)material_of_unshifted_point_in_tree_inobject(p, restricted_tree, &inobject);
  material_data *md = material;

  switch (md->which_subclass) {
    // material grid: interpolate onto user specified material grid to get properties at r
    case material_data::MATERIAL_GRID:
      double u;
      int oi;
      geom_box_tree tp;

      tp = geom_tree_search(p, restricted_tree, &oi);
      // interpolate and project onto material grid
      u = tanh_projection(matgrid_val(p, tp, oi, md) + this->u_p, md->beta, md->eta);
      // interpolate material from material grid point
      epsilon_material_grid(md, u);

      return;
    // material read from file: interpolate to get properties at r
    case material_data::MATERIAL_FILE:
      if (md->epsilon_data)
        epsilon_file_material(md, p);
      else
        material = (material_type)default_material;
      return;

    // material specified by user-supplied function: call user
    // function to get properties at r.
    // Note that we initialize the medium to vacuum, so that
    // the user's function only needs to fill in whatever is
    // different from vacuum.
    case material_data::MATERIAL_USER:
      md->medium = medium_struct();
      md->user_func(p, md->user_data, &(md->medium));
      md->medium.check_offdiag_im_zero_or_abort();
      return;

    // position-independent material or metal: there is nothing to do
    case material_data::MEDIUM:
    case material_data::PERFECT_METAL: return;

    default: meep::abort("unknown material type");
  };
}

// returns trace of the tensor diagonal
double geom_epsilon::chi1p1(meep::field_type ft, const meep::vec &r) {
  symm_matrix chi1p1, chi1p1_inv;

#ifdef DEBUG
  vector3 p = vec_to_vector3(r);
  if (p.x < restricted_tree->b.low.x || p.y < restricted_tree->b.low.y ||
      p.z < restricted_tree->b.low.z || p.x > restricted_tree->b.high.x ||
      p.y > restricted_tree->b.high.y || p.z > restricted_tree->b.high.z)
    meep::abort("invalid point (%g,%g,%g)\n", p.x, p.y, p.z);
#endif

  material_type material;
  get_material_pt(material, r);
  material_epsmu(ft, material, &chi1p1, &chi1p1_inv);
  material_gc(material);

  return (chi1p1.m00 + chi1p1.m11 + chi1p1.m22) / 3;
}

/* Find frontmost object in v, along with the constant material behind it.
   Returns false if material behind the object is not constant.

   Requires moderately horrifying logic to figure things out properly,
   stolen from MPB. */
static bool get_front_object(const meep::volume &v, geom_box_tree geometry_tree, vector3 &pcenter,
                             const geometric_object **o_front, vector3 &shiftby_front,
                             material_type &mat_front, material_type &mat_behind) {
  vector3 p;
  const geometric_object *o1 = 0, *o2 = 0;
  vector3 shiftby1 = {0, 0, 0}, shiftby2 = {0, 0, 0};
  geom_box pixel;
  material_type mat1 = vacuum, mat2 = vacuum;
  int id1 = -1, id2 = -1;
  const int num_neighbors[3] = {3, 5, 9};
  const int neighbors[3][9][3] = {{{0, 0, 0},
                                   {0, 0, -1},
                                   {0, 0, 1},
                                   {0, 0, 0},
                                   {0, 0, 0},
                                   {0, 0, 0},
                                   {0, 0, 0},
                                   {0, 0, 0},
                                   {0, 0, 0}},
                                  {{0, 0, 0},
                                   {-1, -1, 0},
                                   {1, 1, 0},
                                   {-1, 1, 0},
                                   {1, -1, 0},
                                   {0, 0, 0},
                                   {0, 0, 0},
                                   {0, 0, 0},
                                   {0, 0, 0}},
                                  {{0, 0, 0},
                                   {1, 1, 1},
                                   {1, 1, -1},
                                   {1, -1, 1},
                                   {1, -1, -1},
                                   {-1, 1, 1},
                                   {-1, 1, -1},
                                   {-1, -1, 1},
                                   {-1, -1, -1}}};
  pixel = gv2box(v);
  pcenter = p = vec_to_vector3(v.center());
  double d1, d2, d3;
  d1 = (pixel.high.x - pixel.low.x) * 0.5;
  d2 = (pixel.high.y - pixel.low.y) * 0.5;
  d3 = (pixel.high.z - pixel.low.z) * 0.5;
  int dimension_index = meep::number_of_directions(dim) - 1;
  for (int i = 0; i < num_neighbors[dimension_index]; ++i) {
    const geometric_object *o;
    material_type mat;
    vector3 q, shiftby;
    int id;
    q.x = p.x + neighbors[dimension_index][i][0] * d1;
    q.y = p.y + neighbors[dimension_index][i][1] * d2;
    q.z = p.z + neighbors[dimension_index][i][2] * d3;
    o = object_of_point_in_tree(q, geometry_tree, &shiftby, &id);
    if ((id == id1 && vector3_equal(shiftby, shiftby1)) ||
        (id == id2 && vector3_equal(shiftby, shiftby2)))
      continue;

    mat = (material_type)default_material;
    if (o) {
      material_data *md = (material_data *)o->material;
      if (md->which_subclass != material_data::MATERIAL_FILE) mat = md;
    }
    if (id1 == -1) {
      o1 = o;
      shiftby1 = shiftby;
      id1 = id;
      mat1 = mat;
    }
    else if (id2 == -1 ||
             ((id >= id1 && id >= id2) && (id1 == id2 || material_type_equal(mat1, mat2)))) {
      o2 = o;
      shiftby2 = shiftby;
      id2 = id;
      mat2 = mat;
    }
    else if (!(id1 < id2 && (id1 == id || material_type_equal(mat1, mat))) &&
             !(id2 < id1 && (id2 == id || material_type_equal(mat2, mat))))
      return false;
  }

  // CHECK(id1 > -1, "bug in object_of_point_in_tree?");
  if (id2 == -1) { /* only one nearby object/material */
    id2 = id1;
    o2 = o1;
    mat2 = mat1;
    shiftby2 = shiftby1;
  }

  if ((o1 && is_variable(o1->material)) || (o2 && is_variable(o2->material)) ||
      ((is_variable(default_material) || is_file(default_material)) &&
       (!o1 || is_file(o1->material) || !o2 || is_file(o2->material))))
    return false;

  if (id1 >= id2) {
    *o_front = o1;
    shiftby_front = shiftby1;
    mat_front = mat1;
    if (id1 == id2)
      mat_behind = mat1;
    else
      mat_behind = mat2;
  }
  if (id2 > id1) {
    *o_front = o2;
    shiftby_front = shiftby2;
    mat_front = mat2;
    mat_behind = mat1;
  }
  return true;
}

void geom_epsilon::eff_chi1inv_row(meep::component c, double chi1inv_row[3], const meep::volume &v,
                                   double tol, int maxeval) {
  symm_matrix meps_inv;
  bool fallback;
  eff_chi1inv_matrix(c, &meps_inv, v, tol, maxeval, fallback);

  if (fallback) { fallback_chi1inv_row(c, chi1inv_row, v, tol, maxeval); }
  else {
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
      case meep::NO_DIRECTION: chi1inv_row[0] = chi1inv_row[1] = chi1inv_row[2] = 0; break;
    }
  }
}

void geom_epsilon::eff_chi1inv_matrix(meep::component c, symm_matrix *chi1inv_matrix,
                                      const meep::volume &v, double tol, int maxeval,
                                      bool &fallback, double fill_override) {
  const geometric_object *o;
  material_type mat, mat_behind;
  symm_matrix meps;
  vector3 p, shiftby, normal;
  fallback = false;

  if (maxeval == 0) {
  noavg:
    get_material_pt(mat, v.center());
  trivial:
    material_epsmu(meep::type(c), mat, &meps, chi1inv_matrix);
    material_gc(mat);
    return;
  }

  if (!get_front_object(v, geometry_tree, p, &o, shiftby, mat, mat_behind)) {
    get_material_pt(mat, v.center());
    if (mat &&
        (mat->which_subclass == material_data::MATERIAL_USER ||
         mat->which_subclass == material_data::MATERIAL_GRID) &&
        mat->do_averaging) {
      fallback = true;
      return;
    }
    else { goto trivial; }
  }

  /* check for trivial case of only one object/material */
  if (material_type_equal(mat, mat_behind)) goto trivial;

  // it doesn't make sense to average metals (electric or magnetic)
  if (is_metal(meep::type(c), &mat) || is_metal(meep::type(c), &mat_behind)) goto noavg;

  normal = unit_vector3(normal_to_fixed_object(vector3_minus(p, shiftby), *o));
  if (normal.x == 0 && normal.y == 0 && normal.z == 0)
    goto noavg; // couldn't get normal vector for this point, punt
  geom_box pixel = gv2box(v);
  pixel.low = vector3_minus(pixel.low, shiftby);
  pixel.high = vector3_minus(pixel.high, shiftby);

  double fill =
      fill_override >= 0.0 ? fill_override : box_overlap_with_object(pixel, *o, tol, maxeval);

  material_epsmu(meep::type(c), mat, &meps, chi1inv_matrix);
  symm_matrix eps2, epsinv2;
  symm_matrix eps1, delta;
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
    double s = Rot[0][2] * Rot[0][2] + Rot[1][2] * Rot[1][2] + Rot[2][2] * Rot[2][2];
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

#define AVG (fill * (EXPR(eps1)) + (1 - fill) * (EXPR(eps2)))
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

  meps.m00 = -1 / delta.m00;
  meps.m11 = delta.m11 - SQR(delta.m01) / delta.m00;
  meps.m22 = delta.m22 - SQR(delta.m02) / delta.m00;
  meps.m01 = -delta.m01 / delta.m00;
  meps.m02 = -delta.m02 / delta.m00;
  meps.m12 = delta.m12 - (delta.m02 * delta.m01) / delta.m00;

#undef SQR

#define SWAP(a, b)                                                                                 \
  {                                                                                                \
    double xxx = a;                                                                                \
    a = b;                                                                                         \
    b = xxx;                                                                                       \
  }
  /* invert rotation matrix = transpose */
  SWAP(Rot[0][1], Rot[1][0]);
  SWAP(Rot[0][2], Rot[2][0]);
  SWAP(Rot[2][1], Rot[1][2]);
  sym_matrix_rotate(&meps, &meps, Rot); /* rotate back */
#undef SWAP

#ifdef DEBUG
  if (!sym_matrix_positive_definite(&meps))
    meep::abort("negative mean epsilon from Kottke algorithm");
#endif

  sym_matrix_invert(chi1inv_matrix, &meps);
}

/* The filling fraction this pixel would smooth with. Returns false when the
   pixel does not straddle an interface between two distinct materials, in
   which case the shape derivative vanishes there and the caller can skip it --
   which is what restricts the gradient loop to the boundary shell. */
bool geom_epsilon::interface_fill(const meep::volume &v, double tol, int maxeval, double &fill,
                                  const geometric_object **which, vector3 *shift,
                                  material_type *front_mat, material_type *behind_mat) {
  const geometric_object *o;
  material_type mat, mat_behind;
  vector3 p, shiftby;
  fill = -1.0;

  if (!get_front_object(v, geometry_tree, p, &o, shiftby, mat, mat_behind)) return false;
  if (material_type_equal(mat, mat_behind)) return false;

  geom_box pixel = gv2box(v);
  pixel.low = vector3_minus(pixel.low, shiftby);
  pixel.high = vector3_minus(pixel.high, shiftby);
  fill = box_overlap_with_object(pixel, *o, tol, maxeval);
  if (which) *which = o;
  if (shift) *shift = shiftby;
  if (front_mat) *front_mat = mat;
  if (behind_mat) *behind_mat = mat_behind;
  return true;
}

static int eps_ever_negative = 0;
static meep::field_type func_ft = meep::E_stuff;

struct matgrid_volavg {
  meep::ndim dim;   // dimensionality of voxel
  double rad;       // (spherical) voxel radius
  double uval;      // bilinearly-interpolated weight at voxel center
  double ugrad_abs; // magnitude of gradient of uval
  double beta;      // thresholding bias
  double eta;       // thresholding erosion/dilation
  double eps1;      // trace of epsilon tensor from medium 1
  double eps2;      // trace of epsilon tensor from medium 2
};

static void get_uproj_w(const matgrid_volavg *mgva, double x0, double &u_proj, double &w) {
  // use a linear approximation for the material grid weights around the Yee grid point
  u_proj = tanh_projection(mgva->uval + mgva->ugrad_abs * x0, mgva->beta, mgva->eta);
  if (mgva->dim == meep::D1)
    w = 1 / (2 * mgva->rad);
  else if (mgva->dim == meep::D2 || mgva->dim == meep::Dcyl)
    w = 2 * sqrt(mgva->rad * mgva->rad - x0 * x0) / (meep::pi * mgva->rad * mgva->rad);
  else if (mgva->dim == meep::D3)
    // 4.0 / 3.0, not 4 / 3: the latter is integer division, which drops the
    // denominator to pi*rad^3 instead of the sphere volume (4/3)*pi*rad^3 and
    // leaves the kernel integrating to 4/3 rather than 1. The 1d and 2d branches
    // above are correctly normalized, so this only ever affected 3d.
    w = meep::pi * (mgva->rad * mgva->rad - x0 * x0) /
        (4.0 / 3.0 * meep::pi * mgva->rad * mgva->rad * mgva->rad);
}

#ifdef CTL_HAS_COMPLEX_INTEGRATION
static cnumber matgrid_ceps_func(int n, number *x, void *mgva_) {
  (void)n; // unused
  double u_proj = 0, w = 0;
  matgrid_volavg *mgva = (matgrid_volavg *)mgva_;
  get_uproj_w(mgva, x[0], u_proj, w);
  cnumber ret;
  ret.re = (1 - u_proj) * mgva->eps1 + u_proj * mgva->eps2;
  ret.im = (1 - u_proj) / mgva->eps1 + u_proj / mgva->eps2;
  return ret * w;
}
#else
static number matgrid_eps_func(int n, number *x, void *mgva_) {
  (void)n; // unused
  double u_proj = 0, w = 0;
  matgrid_volavg *mgva = (matgrid_volavg *)mgva_;
  get_uproj_w(mgva, x[0], u_proj, w);
  return w * ((1 - u_proj) * mgva->eps1 + u_proj * mgva->eps2);
}
static number matgrid_inveps_func(int n, number *x, void *mgva_) {
  (void)n; // unused
  double u_proj = 0, w = 0;
  matgrid_volavg *mgva = (matgrid_volavg *)mgva_;
  get_uproj_w(mgva, x[0], u_proj, w);
  return w * ((1 - u_proj) / mgva->eps1 + u_proj / mgva->eps2);
}
#endif

#ifdef CTL_HAS_COMPLEX_INTEGRATION
static cnumber ceps_func(int n, number *x, void *geomeps_) {
  geom_epsilon *geomeps = (geom_epsilon *)geomeps_;
  vector3 p = {0, 0, 0};
  p.x = x[0];
  p.y = n > 1 ? x[1] : 0;
  p.z = n > 2 ? x[2] : 0;
  double s = 1;
  if (dim == meep::Dcyl) {
    double py = p.y;
    p.y = p.z;
    p.z = py;
    s = p.x;
  }
  cnumber ret;
  double ep = geomeps->chi1p1(func_ft, vector3_to_vec(p));
  if (ep < 0) eps_ever_negative = 1;
  ret.re = ep * s;
  ret.im = s / ep;
  return ret;
}
#else
static number eps_func(int n, number *x, void *geomeps_) {
  geom_epsilon *geomeps = (geom_epsilon *)geomeps_;
  vector3 p = {0, 0, 0};
  double s = 1;
  p.x = x[0];
  p.y = n > 1 ? x[1] : 0;
  p.z = n > 2 ? x[2] : 0;
  if (dim == meep::Dcyl) {
    double py = p.y;
    p.y = p.z;
    p.z = py;
    s = p.x;
  }
  double ep = geomeps->chi1p1(func_ft, vector3_to_vec(p));
  if (ep < 0) eps_ever_negative = 1;
  return ep * s;
}
static number inveps_func(int n, number *x, void *geomeps_) {
  geom_epsilon *geomeps = (geom_epsilon *)geomeps_;
  vector3 p = {0, 0, 0};
  double s = 1;
  p.x = x[0];
  p.y = n > 1 ? x[1] : 0;
  p.z = n > 2 ? x[2] : 0;
  if (dim == meep::Dcyl) {
    double py = p.y;
    p.y = p.z;
    p.z = py;
    s = p.x;
  }
  double ep = geomeps->chi1p1(func_ft, vector3_to_vec(p));
  if (ep < 0) eps_ever_negative = 1;
  return s / ep;
}
#endif

// fallback meaneps using libctl's adaptive cubature routine
void geom_epsilon::fallback_chi1inv_row(meep::component c, double chi1inv_row[3],
                                        const meep::volume &v, double tol, int maxeval) {

  symm_matrix chi1p1, chi1p1_inv;
  material_type material;
  vector3 p = vec_to_vector3(v.center());
  boolean inobject;
  material =
      (material_type)material_of_unshifted_point_in_tree_inobject(p, restricted_tree, &inobject);
  material_data *md = material;
  meep::vec gradient(zero_vec(v.dim));
  double uval = 0;

  if (md->which_subclass == material_data::MATERIAL_GRID) {
    geom_box_tree tp;
    int oi;
    tp = geom_tree_search(p, restricted_tree, &oi);
    gradient = matgrid_grad(p, tp, oi, md);
    uval = matgrid_val(p, tp, oi, md) + this->u_p;
  }
  else { gradient = normal_vector(meep::type(c), v); }

  get_material_pt(material, v.center());
  material_epsmu(meep::type(c), material, &chi1p1, &chi1p1_inv);
  material_gc(material);
  if (chi1p1.m01 != 0 || chi1p1.m02 != 0 || chi1p1.m12 != 0 || chi1p1.m00 != chi1p1.m11 ||
      chi1p1.m11 != chi1p1.m22 || chi1p1.m00 != chi1p1.m22 || meep::abs(gradient) < 1e-8) {
    int rownum = meep::component_direction(c) % 3;
    if (rownum == 0) {
      chi1inv_row[0] = chi1p1_inv.m00;
      chi1inv_row[1] = chi1p1_inv.m01;
      chi1inv_row[2] = chi1p1_inv.m02;
    }
    else if (rownum == 1) {
      chi1inv_row[0] = chi1p1_inv.m01;
      chi1inv_row[1] = chi1p1_inv.m11;
      chi1inv_row[2] = chi1p1_inv.m12;
    }
    else {
      chi1inv_row[0] = chi1p1_inv.m02;
      chi1inv_row[1] = chi1p1_inv.m12;
      chi1inv_row[2] = chi1p1_inv.m22;
    }
    return;
  }
  number esterr;
  integer errflag;
  double meps, minveps;

  if (md->which_subclass == material_data::MATERIAL_GRID) {
    number xmin[1], xmax[1];
    matgrid_volavg mgva;
    mgva.dim = v.dim;
    mgva.ugrad_abs = meep::abs(gradient);
    mgva.uval = uval;
    mgva.rad = v.diameter() / 2;
    mgva.beta = md->beta;
    mgva.eta = md->eta;
    mgva.eps1 =
        (md->medium_1.epsilon_diag.x + md->medium_1.epsilon_diag.y + md->medium_1.epsilon_diag.z) /
        3;
    mgva.eps2 =
        (md->medium_2.epsilon_diag.x + md->medium_2.epsilon_diag.y + md->medium_2.epsilon_diag.z) /
        3;
    xmin[0] = -v.diameter() / 2;
    xmax[0] = v.diameter() / 2;
#ifdef CTL_HAS_COMPLEX_INTEGRATION
    cnumber ret = cadaptive_integration(matgrid_ceps_func, xmin, xmax, 1, (void *)&mgva, 0, tol,
                                        maxeval, &esterr, &errflag);
    meps = ret.re;
    minveps = ret.im;
#else
    meps = adaptive_integration(matgrid_eps_func, xmin, xmax, 1, (void *)&mgva, 0, tol, maxeval,
                                &esterr, &errflag);
    minveps = adaptive_integration(matgrid_inveps_func, xmin, xmax, 1, (void *)&mgva, 0, tol,
                                   maxeval, &esterr, &errflag);
#endif
  }
  else {
    integer n;
    number xmin[3], xmax[3];
    vector3 gvmin, gvmax;
    gvmin = vec_to_vector3(v.get_min_corner());
    gvmax = vec_to_vector3(v.get_max_corner());
    xmin[0] = gvmin.x;
    xmax[0] = gvmax.x;
    if (dim == meep::Dcyl) {
      xmin[1] = gvmin.z;
      xmin[2] = gvmin.y;
      xmax[1] = gvmax.z;
      xmax[2] = gvmax.y;
    }
    else {
      xmin[1] = gvmin.y;
      xmin[2] = gvmin.z;
      xmax[1] = gvmax.y;
      xmax[2] = gvmax.z;
    }
    if (xmin[2] == xmax[2])
      n = xmin[1] == xmax[1] ? 1 : 2;
    else
      n = 3;
    double vol = 1;
    for (int i = 0; i < n; ++i)
      vol *= xmax[i] - xmin[i];
    if (dim == meep::Dcyl) vol *= (xmin[0] + xmax[0]) * 0.5;
    eps_ever_negative = 0;
    func_ft = meep::type(c);
#ifdef CTL_HAS_COMPLEX_INTEGRATION
    cnumber ret = cadaptive_integration(ceps_func, xmin, xmax, n, (void *)this, 0, tol, maxeval,
                                        &esterr, &errflag);
    meps = ret.re / vol;
    minveps = ret.im / vol;
#else
    meps = adaptive_integration(eps_func, xmin, xmax, n, (void *)this, 0, tol, maxeval, &esterr,
                                &errflag) /
           vol;
    minveps = adaptive_integration(inveps_func, xmin, xmax, n, (void *)this, 0, tol, maxeval,
                                   &esterr, &errflag) /
              vol;
#endif
  }
  if (eps_ever_negative) // averaging negative eps causes instability
    minveps = 1.0 / (meps = eps(v.center()));
  {
    double n[3] = {0, 0, 0};
    double nabsinv = 1.0 / meep::abs(gradient);
    LOOP_OVER_DIRECTIONS(gradient.dim, k) { n[k % 3] = gradient.in_direction(k) * nabsinv; }
    int rownum = meep::component_direction(c) % 3;
    for (int i = 0; i < 3; ++i)
      chi1inv_row[i] = n[rownum] * n[i] * (minveps - 1 / meps);
    chi1inv_row[rownum] += 1 / meps;
  }
}

static double get_chi3(meep::component c, const medium_struct *m) {
  switch (c) {
    case meep::Er:
    case meep::Ex: return m->E_chi3_diag.x;
    case meep::Ep:
    case meep::Ey: return m->E_chi3_diag.y;
    case meep::Ez: return m->E_chi3_diag.z;
    case meep::Hr:
    case meep::Hx: return m->H_chi3_diag.x;
    case meep::Hp:
    case meep::Hy: return m->H_chi3_diag.y;
    case meep::Hz: return m->H_chi3_diag.z;
    default: return 0;
  }
}

static double get_chi2(meep::component c, const medium_struct *m) {
  switch (c) {
    case meep::Er:
    case meep::Ex: return m->E_chi2_diag.x;
    case meep::Ep:
    case meep::Ey: return m->E_chi2_diag.y;
    case meep::Ez: return m->E_chi2_diag.z;
    case meep::Hr:
    case meep::Hx: return m->H_chi2_diag.x;
    case meep::Hp:
    case meep::Hy: return m->H_chi2_diag.y;
    case meep::Hz: return m->H_chi2_diag.z;
    default: return 0;
  }
}

static double get_chi(meep::component c, const medium_struct *m, int p) {
  return ((p == 2) ? get_chi2(c, m) : get_chi3(c, m));
}

// the following routines consolidate the former has_chi2, has_chi3 and chi2, chi3
// p=2,3 for chi2,chi3
bool geom_epsilon::has_chi(meep::component c, int p) {
  medium_struct *mm;

  for (int i = 0; i < geometry.num_items; ++i)
    if (is_medium(geometry.items[i].material, &mm))
      if (get_chi(c, mm, p) != 0) return true;

  for (int i = 0; i < extra_materials.num_items; ++i)
    if (is_medium(extra_materials.items[i], &mm))
      if (get_chi(c, mm, p) != 0) return true;

  return (is_medium(default_material, &mm) && get_chi(c, mm, p) != 0);
}

bool geom_epsilon::has_chi3(meep::component c) { return has_chi(c, 3); }

bool geom_epsilon::has_chi2(meep::component c) { return has_chi(c, 2); }

double geom_epsilon::chi(meep::component c, const meep::vec &r, int p) {
  material_type material;
  get_material_pt(material, r);

  double chi_val;

  material_data *md = material;
  switch (md->which_subclass) {
    case material_data::MEDIUM:
    case material_data::MATERIAL_GRID:
    case material_data::MATERIAL_USER: chi_val = get_chi(c, &(md->medium), p); break;

    default: chi_val = 0;
  };

  material_gc(material);
  return chi_val;
}

double geom_epsilon::chi3(meep::component c, const meep::vec &r) { return chi(c, r, 3); }

double geom_epsilon::chi2(meep::component c, const meep::vec &r) { return chi(c, r, 2); }

static bool mu_not_1(material_type m) {
  medium_struct *mm;
  return (is_medium(m, &mm) &&
          (mm->mu_diag.x != 1 || mm->mu_diag.y != 1 || mm->mu_diag.z != 1 ||
           mm->mu_offdiag.x.re != 0 || mm->mu_offdiag.y.re != 0 || mm->mu_offdiag.z.re != 0));
}
static bool mu_not_1(void *m) { return mu_not_1((material_type)m); }

bool geom_epsilon::has_mu() {
  for (int i = 0; i < geometry.num_items; ++i)
    if (mu_not_1(geometry.items[i].material)) return true;

  for (int i = 0; i < extra_materials.num_items; ++i)
    if (mu_not_1(extra_materials.items[i])) return true;

  return (mu_not_1(default_material));
}

/* a global scalar conductivity to add to all materials; this
   is mostly for the convenience of Casimir calculations where
   the global conductivity corresponds to a rotation to
   complex frequencies */
static double global_D_conductivity = 0, global_B_conductivity = 0;

static double get_cnd(meep::component c, const medium_struct *m) {
  switch (c) {
    case meep::Dr:
    case meep::Dx: return m->D_conductivity_diag.x + global_D_conductivity;
    case meep::Dp:
    case meep::Dy: return m->D_conductivity_diag.y + global_D_conductivity;
    case meep::Dz: return m->D_conductivity_diag.z + global_D_conductivity;
    case meep::Br:
    case meep::Bx: return m->B_conductivity_diag.x + global_B_conductivity;
    case meep::Bp:
    case meep::By: return m->B_conductivity_diag.y + global_B_conductivity;
    case meep::Bz: return m->B_conductivity_diag.z + global_B_conductivity;
    default: return 0;
  }
}

static bool has_conductivity(const material_type &md, meep::component c) {
  medium_struct *mm;
  if (is_medium(md, &mm) && get_cnd(c, mm)) return true;
  if (md->which_subclass == material_data::MATERIAL_GRID &&
      (get_cnd(c, &md->medium_1) || get_cnd(c, &md->medium_2) || md->damping != 0))
    return true;
  return false;
}

bool geom_epsilon::has_conductivity(meep::component c) {
  FOR_DIRECTIONS(d) FOR_SIDES(b) {
    if (cond[d][b].prof) return true;
  }
  for (int i = 0; i < geometry.num_items; ++i)
    if (meep_geom::has_conductivity((material_type)geometry.items[i].material, c)) return true;
  for (int i = 0; i < extra_materials.num_items; ++i)
    if (meep_geom::has_conductivity((material_type)extra_materials.items[i], c)) return true;
  return meep_geom::has_conductivity((material_type)default_material, c);
}

static meep::vec geometry_edge; // geometry_lattice.size / 2
// TODO
double geom_epsilon::conductivity(meep::component c, const meep::vec &r) {
  material_type material;
  get_material_pt(material, r);

  double cond_val;
  material_data *md = material;
  switch (md->which_subclass) {
    case material_data::MEDIUM:
    case material_data::MATERIAL_GRID:
    case material_data::MATERIAL_USER: cond_val = get_cnd(c, &(md->medium)); break;
    default: cond_val = 0;
  }
  material_gc(material);

  // if the user specified scalar absorbing layers, add their conductivities
  // to cond_val (isotropically, for both magnetic and electric conductivity).
  LOOP_OVER_DIRECTIONS(r.dim, d) {
    double x = r.in_direction(d);
    double edge = geometry_edge.in_direction(d) - cond[d][meep::High].L;
    if (cond[d][meep::High].prof && x >= edge) {
      int N = cond[d][meep::High].N;
      double ui = N * (x - edge) / cond[d][meep::High].L;
      int i = int(ui);
      if (i >= N)
        cond_val += cond[d][meep::High].prof[N];
      else {
        double di = ui - i;
        cond_val += cond[d][meep::High].prof[i] * (1 - di) + cond[d][meep::High].prof[i + 1] * di;
      }
    }
    edge = cond[d][meep::Low].L - geometry_edge.in_direction(d);
    if (cond[d][meep::Low].prof && x <= edge) {
      int N = cond[d][meep::Low].N;
      double ui = N * (edge - x) / cond[d][meep::Low].L;
      int i = int(ui);
      if (i >= N)
        cond_val += cond[d][meep::Low].prof[N];
      else {
        double di = ui - i;
        cond_val += cond[d][meep::Low].prof[i] * (1 - di) + cond[d][meep::Low].prof[i + 1] * di;
      }
    }
  }

  return cond_val;
}

/* like susceptibility_equal in ctl-io.cpp, but ignores sigma and id
   (must be updated manually, re-copying from ctl-io.cpp), if we
   add new susceptibility subclasses) */
static bool susceptibility_equiv(const susceptibility &o0, const susceptibility &o) {
  if (!vector3_equal(o0.bias, o.bias)) return false;
  if (o0.frequency != o.frequency) return false;
  if (o0.gamma != o.gamma) return false;
  if (o0.alpha != o.alpha) return false;
  if (o0.noise_amp != o.noise_amp) return false;
  if (o0.drude != o.drude) return false;
  if (o0.saturated_gyrotropy != o.saturated_gyrotropy) return false;
  if (o0.is_file != o.is_file) return false;

  if (o0.transitions != o.transitions) return false;
  if (o0.initial_populations != o.initial_populations) return false;

  return true;
}

/* Sum the amplitudes of `mat`'s poles equivalent to the one being registered,
   into the c'th row of its sigma tensor.  Factored out of sigma_row() so that
   eff_sigma_row() can ask the same question of two materials at once. */
void geom_epsilon::material_sigma_row(meep::component c, double sigrow[3], material_type mat) {
  sigrow[0] = sigrow[1] = sigrow[2] = 0.0;
  if (!current_pol) return;

  if (mat->which_subclass != material_data::MATERIAL_USER &&
      mat->which_subclass != material_data::MATERIAL_GRID &&
      mat->which_subclass != material_data::MEDIUM)
    return;

  const susceptibility_list &slist =
      type(c) == meep::E_stuff ? mat->medium.E_susceptibilities : mat->medium.H_susceptibilities;
  /* Accumulate rather than stopping at the first match.  A material grid's
     list is the concatenation of its two constituents' poles, so when both
     media carry the same pole -- same frequency and damping, different
     strength -- two entries are equivalent to the one being registered and
     the mixture's amplitude is their sum.  This is also the right answer for
     a plain medium that happens to list a pole twice, since identical poles
     superpose. */
  for (const susceptibility &susc : slist) {
    if (!susceptibility_equiv(susc, current_pol->user_s)) continue;
    switch (meep::component_index(c)) { // which row of the sigma tensor to return
      case 0:
        sigrow[0] += susc.sigma_diag.x;
        sigrow[1] += susc.sigma_offdiag.x;
        sigrow[2] += susc.sigma_offdiag.y;
        break;
      case 1:
        sigrow[0] += susc.sigma_offdiag.x;
        sigrow[1] += susc.sigma_diag.y;
        sigrow[2] += susc.sigma_offdiag.z;
        break;
      default: // case 2:
        sigrow[0] += susc.sigma_offdiag.y;
        sigrow[1] += susc.sigma_offdiag.z;
        sigrow[2] += susc.sigma_diag.z;
        break;
    }
  }
}

/* Volume-average sigma across a geometric interface.

   The permittivity has always been subpixel-smoothed while sigma was point
   sampled, so a dispersive object's response switched on a whole pixel at a
   time as it moved.  Taper each medium's pole amplitude by the pixel's filling
   fraction instead, using the same front object, shift and fill that
   eff_chi1inv_matrix uses, so the two describe the same interface.

   Linear in fill rather than Kottke, deliberately.  Kottke's normal component
   is a harmonic average of a permittivity; there is no analogous exact result
   for a resonant amplitude, and a harmonic average against a medium without
   the pole would be zero in every boundary pixel.  Scaling each medium's whole
   sigma tensor by a non-negative weight and adding keeps the result positive
   semidefinite -- which is what passivity depends on -- and makes the
   amplitude exactly linear in fill, so the shape derivative is a constant.

   Material grids and user materials are point sampled as before: their
   material varies continuously in space already, and mixing a level set with a
   geometric filling fraction is a separate question. */
void geom_epsilon::eff_sigma_row(meep::component c, double sigrow[3], const meep::volume &v) {
  const geometric_object *o;
  material_type mat, mat_behind;
  vector3 p, shiftby;

  // No smoothing requested: match what the permittivity is doing.
  if (maxeval == 0 || !current_pol) return sigma_row(c, sigrow, v.center());

  /* Out of scope, point sample instead.  A multilevel atom carries per-pixel
     level populations with their own nonlinear dynamics, so there is no reason
     to think scaling an amplitude by a filling fraction gives the response of
     a partially filled pixel the way it does for a linear pole.  A gyrotropic
     sigma is antisymmetric, so the argument that a non-negative combination
     stays positive semidefinite does not apply to it either. */
  const susceptibility &pole = current_pol->user_s;
  if (!pole.transitions.empty() || !pole.initial_populations.empty() || pole.saturated_gyrotropy ||
      pole.bias.x != 0 || pole.bias.y != 0 || pole.bias.z != 0)
    return sigma_row(c, sigrow, v.center());

  if (!get_front_object(v, geometry_tree, p, &o, shiftby, mat, mat_behind))
    return sigma_row(c, sigrow, v.center());

  if (material_type_equal(mat, mat_behind) || is_variable(mat) || is_variable(mat_behind)) {
    material_gc(mat);
    return sigma_row(c, sigrow, v.center());
  }

  geom_box pixel = gv2box(v);
  pixel.low = vector3_minus(pixel.low, shiftby);
  pixel.high = vector3_minus(pixel.high, shiftby);
  double fill = box_overlap_with_object(pixel, *o, tol, maxeval);

  double front[3], behind[3];
  material_sigma_row(c, front, mat);
  material_sigma_row(c, behind, mat_behind);
  for (int i = 0; i < 3; ++i)
    sigrow[i] = fill * front[i] + (1 - fill) * behind[i];

  material_gc(mat);
}

void geom_epsilon::sigma_row(meep::component c, double sigrow[3], const meep::vec &r) {

  vector3 p = vec_to_vector3(r);

  boolean inobject;
  material_type mat =
      (material_type)material_of_unshifted_point_in_tree_inobject(p, restricted_tree, &inobject);

  if (mat->which_subclass == material_data::MATERIAL_USER) {
    mat->medium = medium_struct();
    mat->user_func(p, mat->user_data, &(mat->medium));
    mat->medium.check_offdiag_im_zero_or_abort();
  }

  if (mat->which_subclass == material_data::MATERIAL_GRID) {
    double u;
    int oi;
    geom_box_tree tp;

    tp = geom_tree_search(p, restricted_tree, &oi);
    // interpolate and project onto material grid
    u = tanh_projection(matgrid_val(p, tp, oi, mat) + this->u_p, mat->beta, mat->eta);
    epsilon_material_grid(mat, u); // interpolate material from material grid point
    mat->medium.check_offdiag_im_zero_or_abort();
  }

  material_sigma_row(c, sigrow, mat);
  material_gc(mat);
}

/* make multilevel_susceptibility from python input data */
static meep::susceptibility *make_multilevel_sus(const susceptibility_struct *d) {
  if (!d || d->transitions.size() == 0) return NULL;

  // the user can number the levels however she wants, but we
  // will renumber them to 0...(L-1)
  int minlev = d->transitions[0].to_level;
  int maxlev = minlev;
  for (size_t t = 0; t < d->transitions.size(); ++t) {
    if (minlev > d->transitions[t].from_level) minlev = d->transitions[t].from_level;
    if (minlev > d->transitions[t].to_level) minlev = d->transitions[t].to_level;
    if (maxlev < d->transitions[t].from_level) maxlev = d->transitions[t].from_level;
    if (maxlev < d->transitions[t].to_level) maxlev = d->transitions[t].to_level;
  }
  size_t L = maxlev - minlev + 1; // number of atom levels

  // count number of radiative transitions
  int T = 0;
  for (size_t t = 0; t < d->transitions.size(); ++t)
    if (d->transitions[t].frequency != 0) ++T;
  if (T == 0) return NULL; // don't bother if there is no radiative coupling

  // non-radiative transition-rate matrix Gamma
  meep::realnum *Gamma = new meep::realnum[L * L];
  memset(Gamma, 0, sizeof(meep::realnum) * (L * L));
  for (size_t t = 0; t < d->transitions.size(); ++t) {
    int i = d->transitions[t].from_level - minlev;
    int j = d->transitions[t].to_level - minlev;
    Gamma[i * L + i] += +d->transitions[t].transition_rate + d->transitions[t].pumping_rate;
    Gamma[j * L + i] -= +d->transitions[t].transition_rate + d->transitions[t].pumping_rate;
  }

  // initial populations of each level
  meep::realnum *N0 = new meep::realnum[L];
  memset(N0, 0, sizeof(meep::realnum) * L);
  for (size_t p = 0; p < d->initial_populations.size() && p < L; ++p)
    N0[p] = d->initial_populations[p];

  meep::realnum *alpha = new meep::realnum[L * T];
  memset(alpha, 0, sizeof(meep::realnum) * (L * T));
  meep::realnum *omega = new meep::realnum[T];
  meep::realnum *gamma = new meep::realnum[T];
  meep::realnum *sigmat = new meep::realnum[T * 5];

  const double pi = 3.14159265358979323846264338327950288; // need pi below.

  for (size_t t = 0, tr = 0; t < d->transitions.size(); ++t)
    if (d->transitions[t].frequency != 0) {
      omega[tr] = d->transitions[t].frequency; // no 2*pi here
      gamma[tr] = d->transitions[t].gamma;
      if (dim == meep::Dcyl) {
        sigmat[5 * tr + meep::R] = d->transitions[t].sigma_diag.x;
        sigmat[5 * tr + meep::P] = d->transitions[t].sigma_diag.y;
        sigmat[5 * tr + meep::Z] = d->transitions[t].sigma_diag.z;
      }
      else {
        sigmat[5 * tr + meep::X] = d->transitions[t].sigma_diag.x;
        sigmat[5 * tr + meep::Y] = d->transitions[t].sigma_diag.y;
        sigmat[5 * tr + meep::Z] = d->transitions[t].sigma_diag.z;
      }
      int i = d->transitions[t].from_level - minlev;
      int j = d->transitions[t].to_level - minlev;
      alpha[i * T + tr] = -1.0 / (2 * pi * omega[tr]); // but we *do* need the 2*pi here. -- AWC
      alpha[j * T + tr] = +1.0 / (2 * pi * omega[tr]);
      ++tr;
    }

  meep::multilevel_susceptibility *s =
      new meep::multilevel_susceptibility(L, T, Gamma, N0, alpha, omega, gamma, sigmat);

  delete[] Gamma;
  delete[] N0;
  delete[] alpha;
  delete[] omega;
  delete[] gamma;
  delete[] sigmat;

  return s;
}

// add a polarization to the list if it is not already there
static pol *add_pol(pol *pols, const susceptibility &user_s) {
  struct pol *p = pols;
  while (p && !susceptibility_equiv(user_s, p->user_s))
    p = p->next;
  if (!p) {
    p = new pol;
    p->user_s = user_s;
    p->next = pols;
    pols = p;
  }
  return pols;
}

static pol *add_pols(pol *pols, const susceptibility_list &slist) {
  for (const susceptibility &susc : slist) {
    pols = add_pol(pols, susc);
  }
  return pols;
}

/* Collect the susceptibilities a material can contribute.

   is_medium() is true only for a plain MEDIUM, but a material grid resolves at
   each point to a mixture of its two constituents (see epsilon_material_grid),
   so both of their pole lists have to be registered here.  Without this a
   dispersive material grid is silently non-dispersive: nothing ever calls
   sigma_row() for its poles, so the interpolated amplitudes go nowhere, even
   though the adjoint side reads medium.E_susceptibilities directly and does
   see them. */
static pol *add_material_pols(pol *pols, void *m, meep::field_type ft) {
  medium_struct *mm;
  if (is_medium(m, &mm))
    return add_pols(pols, ft == meep::E_stuff ? mm->E_susceptibilities : mm->H_susceptibilities);

  material_type mat = (material_type)m;
  if (is_material_grid(mat) && ft == meep::E_stuff) {
    /* Only the electric side: epsilon_material_grid does not interpolate the
       magnetic response, and MaterialGrid.__init__ refuses a constituent that
       carries one, so registering H poles here would allocate polarizations
       that sigma_row can only ever report as zero. */
    pols = add_pols(pols, mat->medium_1.E_susceptibilities);
    pols = add_pols(pols, mat->medium_2.E_susceptibilities);
  }
  return pols;
}

void geom_epsilon::add_susceptibilities(meep::structure *s) {
  add_susceptibilities(meep::E_stuff, s);
  add_susceptibilities(meep::H_stuff, s);
}

void geom_epsilon::add_susceptibilities(meep::field_type ft, meep::structure *s) {
  pol *pols = 0;
  // construct a list of the unique susceptibilities in the geometry:
  for (int i = 0; i < geometry.num_items; ++i)
    pols = add_material_pols(pols, geometry.items[i].material, ft);

  for (int i = 0; i < extra_materials.num_items; ++i)
    pols = add_material_pols(pols, extra_materials.items[i], ft);

  pols = add_material_pols(pols, default_material, ft);

  for (struct pol *p = pols; p; p = p->next) {

    susceptibility *ss = &(p->user_s);
    if (ss->is_file) meep::abort("unknown susceptibility");
    bool noisy = (ss->noise_amp != 0.0);
    bool gyrotropic =
        (ss->saturated_gyrotropy || ss->bias.x != 0.0 || ss->bias.y != 0.0 || ss->bias.z != 0.0);
    meep::susceptibility *sus;

    if (ss->transitions.size() != 0 || ss->initial_populations.size() != 0) {
      // multilevel atom
      sus = make_multilevel_sus(ss);
      if (meep::verbosity > 0) master_printf("multilevel atom susceptibility\n");
    }
    else {
      if (noisy) {
        sus = new meep::noisy_lorentzian_susceptibility(ss->noise_amp, ss->frequency, ss->gamma,
                                                        ss->drude);
      }
      else if (gyrotropic) {
        meep::gyrotropy_model model = ss->saturated_gyrotropy ? meep::GYROTROPIC_SATURATED
                                      : ss->drude             ? meep::GYROTROPIC_DRUDE
                                                              : meep::GYROTROPIC_LORENTZIAN;
        sus = new meep::gyrotropic_susceptibility(meep::vec(ss->bias.x, ss->bias.y, ss->bias.z),
                                                  ss->frequency, ss->gamma, ss->alpha, model);
      }
      else { sus = new meep::lorentzian_susceptibility(ss->frequency, ss->gamma, ss->drude); }
      if (meep::verbosity > 0) {
        master_printf("%s%s susceptibility: frequency=%g, gamma=%g",
                      noisy        ? "noisy "
                      : gyrotropic ? "gyrotropic "
                                   : "",
                      ss->saturated_gyrotropy ? "Landau-Lifshitz-Gilbert-type"
                      : ss->drude             ? "drude"
                                              : "lorentzian",
                      ss->frequency, ss->gamma);
        if (noisy) master_printf(", amp=%g ", ss->noise_amp);
        if (gyrotropic) {
          if (ss->saturated_gyrotropy) master_printf(", alpha=%g", ss->alpha);
          master_printf(", bias=(%g,%g,%g)", ss->bias.x, ss->bias.y, ss->bias.z);
        }
        master_printf("\n");
      }
    }

    current_pol = p;
    if (sus) {
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

typedef struct pml_profile_thunk {
  meep::pml_profile_func func;
  void *func_data;
} pml_profile_thunk;

double pml_profile_wrapper(int dim, double *u, void *user_data) {
  (void)dim; // unused
  pml_profile_thunk *mythunk = (pml_profile_thunk *)user_data;
  return mythunk->func(u[0], mythunk->func_data);
}

/***************************************************************/
/* mechanism for allowing users to specify non-PML absorbing   */
/* layers.                                                     */
/***************************************************************/
absorber_list create_absorber_list() {
  absorber_list alist = new absorber_list_type;
  return alist;
}

void destroy_absorber_list(absorber_list alist) { delete alist; }

void add_absorbing_layer(absorber_list alist, double thickness, int direction, int side,
                         double R_asymptotic, double mean_stretch, meep::pml_profile_func func,
                         void *func_data) {
  absorber myabsorber;
  myabsorber.thickness = thickness;
  myabsorber.direction = direction;
  myabsorber.side = side;
  myabsorber.R_asymptotic = R_asymptotic;
  myabsorber.mean_stretch = mean_stretch;
  myabsorber.pml_profile = func;
  myabsorber.pml_profile_data = func_data;

  if (alist == 0) meep::abort("invalid absorber_list in add_absorbing_layer");

  alist->push_back(myabsorber);
}

/***************************************************************/
/***************************************************************/
/***************************************************************/

/* create a geom_epsilon object that can persist
if needed */
geom_epsilon *make_geom_epsilon(meep::structure *s, geometric_object_list *g, vector3 center,
                                bool _ensure_periodicity, material_type _default_material,
                                material_type_list extra_materials) {
  // set global variables in libctlgeom based on data fields in s
  geom_initialize();
  geometry_center = center;

  if (_default_material->which_subclass != material_data::MATERIAL_USER &&
      _default_material->which_subclass != material_data::PERFECT_METAL) {
    _default_material->medium.check_offdiag_im_zero_or_abort();
  }
  set_default_material(_default_material);
  ensure_periodicity = _ensure_periodicity;
  meep::grid_volume gv = s->gv;
  double resolution = gv.a;

  int sim_dims = 3;
  vector3 size = {0.0, 0.0, 0.0};
  switch (s->user_volume.dim) {
    case meep::D1:
      sim_dims = 1;
      size.z = s->user_volume.nz() / resolution;
      break;
    case meep::D2:
      sim_dims = 2;
      size.x = s->user_volume.nx() / resolution;
      size.y = s->user_volume.ny() / resolution;
      break;
    case meep::D3:
      sim_dims = 3;
      size.x = s->user_volume.nx() / resolution;
      size.y = s->user_volume.ny() / resolution;
      size.z = s->user_volume.nz() / resolution;
      break;
    case meep::Dcyl:
      sim_dims = CYLINDRICAL;
      size.x = s->user_volume.nr() / resolution;
      size.z = s->user_volume.nz() / resolution;
      break;
  };

  set_dimensions(sim_dims);

  geometry_lattice.size = size;
  geometry_edge = vector3_to_vec(size) * 0.5;

  if (meep::verbosity > 0) {
    master_printf("Working in %s dimensions.\n", meep::dimension_name(s->gv.dim));
    master_printf("Computational cell is %g x %g x %g with resolution %g\n", size.x, size.y, size.z,
                  resolution);
  }

  geom_epsilon *geps = new geom_epsilon(*g, extra_materials, gv.pad().surroundings());
  return geps;
}

/* set the materials without previously creating
a geom_eps object */
void set_materials_from_geometry(meep::structure *s, geometric_object_list g, vector3 center,
                                 bool use_anisotropic_averaging, double tol, int maxeval,
                                 bool _ensure_periodicity, material_type _default_material,
                                 absorber_list alist, material_type_list extra_materials) {
  meep_geom::geom_epsilon *geps = meep_geom::make_geom_epsilon(s, &g, center, _ensure_periodicity,
                                                               _default_material, extra_materials);
  set_materials_from_geom_epsilon(s, geps, use_anisotropic_averaging, tol, maxeval, alist);
  delete geps;
}

/* from a previously created geom_epsilon object,
set the materials as specified */
void set_materials_from_geom_epsilon(meep::structure *s, geom_epsilon *geps,
                                     bool use_anisotropic_averaging, double tol, int maxeval,
                                     absorber_list alist) {

  // Store for later use in gradient calculations. These must mirror what
  // structure_chunk::set_chi1inv() actually used, and that routine zeroes maxeval
  // when averaging is off (see anisotropic_averaging.cpp). Passing the raw maxeval
  // through would make get_material_gradient() finite-difference a *smoothed*
  // operator that the forward solve never saw, so the adjoint gradient would not
  // be the derivative of the simulation being run.
  geps->dt = s->dt;
  geps->tol = tol;
  geps->maxeval = use_anisotropic_averaging ? maxeval : 0;

  meep::grid_volume gv = s->gv;
  if (alist) {
    for (absorber_list_type::iterator layer = alist->begin(); layer != alist->end(); layer++) {
      LOOP_OVER_DIRECTIONS(gv.dim, d) {
        if (layer->direction != ALL_DIRECTIONS && layer->direction != d) continue;
        FOR_SIDES(b) {
          if (layer->side != ALL_SIDES && layer->side != b) continue;
          pml_profile_thunk mythunk;
          mythunk.func = layer->pml_profile;
          mythunk.func_data = layer->pml_profile_data;
          geps->set_cond_profile(d, b, layer->thickness, gv.inva * 0.5, pml_profile_wrapper,
                                 (void *)&mythunk, layer->R_asymptotic);
        }
      }
    }
  }
  s->set_materials(*geps, use_anisotropic_averaging, tol, maxeval);
  s->remove_susceptibilities();
  geps->add_susceptibilities(s);

  if (meep::verbosity > 0) master_printf("-----------\n");
}

/***************************************************************/
/* convenience routines for creating materials of various types*/
/***************************************************************/
material_type make_dielectric(double epsilon) {
  material_data *md = new material_data();
  md->medium.epsilon_diag.x = epsilon;
  md->medium.epsilon_diag.y = epsilon;
  md->medium.epsilon_diag.z = epsilon;
  return md;
}

material_type make_user_material(user_material_func user_func, void *user_data, bool do_averaging) {
  material_data *md = new material_data();
  md->which_subclass = material_data::MATERIAL_USER;
  md->user_func = user_func;
  md->user_data = user_data;
  md->do_averaging = do_averaging;
  return md;
}

// this routine subsumes the content of the old
// 'read_epsilon_file' routine
material_type make_file_material(const char *eps_input_file) {
  material_data *md = new material_data();
  md->which_subclass = material_data::MATERIAL_FILE;

  md->epsilon_dims[0] = md->epsilon_dims[1] = md->epsilon_dims[2] = 1;
  if (eps_input_file && eps_input_file[0]) { // file specified
    char *fname = new char[strlen(eps_input_file) + 1];
    strcpy(fname, eps_input_file);
    // parse epsilon-input-file as "fname.h5:dataname"
    char *dataname = strrchr(fname, ':');
    if (dataname) *(dataname++) = 0;
    meep::h5file eps_file(fname, meep::h5file::READONLY, false);
    int rank; // ignored since rank < 3 is equivalent to singleton dims
    md->epsilon_data = (double *)eps_file.read(dataname, &rank, md->epsilon_dims, 3, false);
    if (meep::verbosity > 0)
      master_printf("read in %zdx%zdx%zd epsilon-input-file \"%s\"\n", md->epsilon_dims[0],
                    md->epsilon_dims[1], md->epsilon_dims[2], eps_input_file);
    delete[] fname;
  }

  return md;
}

/******************************************************************************/
/* Material grid functions                                                    */
/******************************************************************************/
material_type make_material_grid(bool do_averaging, double beta, double eta, double damping) {
  material_data *md = new material_data();
  md->which_subclass = material_data::MATERIAL_GRID;
  md->do_averaging = do_averaging;
  md->beta = beta;
  md->eta = eta;
  md->damping = damping;
  return md;
}

void update_weights(material_type matgrid, double *weights) {
  size_t N = matgrid->grid_size.x * matgrid->grid_size.y * matgrid->grid_size.z;
  memcpy(matgrid->weights, weights, N * sizeof(double));
}

/******************************************************************************/
/* Helpers from  libctl/utils/geom.c                                          */
/******************************************************************************/
#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define MIN(a, b) ((a) < (b) ? (a) : (b))

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
  return (
      (BETWEEN(b1->low.x, b2->low.x, b2->high.x) || BETWEEN(b1->high.x, b2->low.x, b2->high.x) ||
       BETWEEN(b2->low.x, b1->low.x, b1->high.x)) &&
      (BETWEEN(b1->low.y, b2->low.y, b2->high.y) || BETWEEN(b1->high.y, b2->low.y, b2->high.y) ||
       BETWEEN(b2->low.y, b1->low.y, b1->high.y)) &&
      (BETWEEN(b1->low.z, b2->low.z, b2->high.z) || BETWEEN(b1->high.z, b2->low.z, b2->high.z) ||
       BETWEEN(b2->low.z, b1->low.z, b1->high.z)));
}

/******************************************************************************/
/* Fragment Statistics                                                        */
/******************************************************************************/

double fragment_stats::tol = 0;
int fragment_stats::maxeval = 0;
int fragment_stats::resolution = 0;
meep::ndim fragment_stats::dims = meep::D1;
geometric_object_list fragment_stats::geom = {};
std::vector<dft_data> fragment_stats::dft_data_list;
std::vector<meep::volume> fragment_stats::pml_1d_vols;
std::vector<meep::volume> fragment_stats::pml_2d_vols;
std::vector<meep::volume> fragment_stats::pml_3d_vols;
std::vector<meep::volume> fragment_stats::absorber_vols;
material_type_list fragment_stats::extra_materials = material_type_list();
bool fragment_stats::split_chunks_evenly = false;
bool fragment_stats::eps_averaging = false;

static geom_box make_box_from_cell(vector3 cell_size) {
  double edgex = cell_size.x / 2;
  double edgey = cell_size.y / 2;
  double edgez = cell_size.z / 2;
  vector3 low = {-edgex, -edgey, -edgez};
  vector3 high = {edgex, edgey, edgez};
  geom_box result = {low, high};

  return result;
}

static size_t get_pixels_in_box(geom_box *b, int empty_pixel = 1) {
  int empty_x = b->low.x == b->high.x;
  int empty_y = b->low.y == b->high.y;
  int empty_z = b->low.z == b->high.z;

  double total_pixels =
      ((empty_x ? empty_pixel : (b->high.x - b->low.x) * fragment_stats::resolution) *
       (empty_y ? empty_pixel : (b->high.y - b->low.y) * fragment_stats::resolution) *
       (empty_z ? empty_pixel : (b->high.z - b->low.z) * fragment_stats::resolution));

  return (size_t)ceil(total_pixels);
}

static void center_box(geom_box *b) {
  b->low = vector3_plus(geometry_center, b->low);
  b->high = vector3_plus(geometry_center, b->high);
}

static fragment_stats init_stats(geom_box &box, double tol, int maxeval, meep::grid_volume *gv) {
  fragment_stats::tol = tol;
  fragment_stats::maxeval = maxeval;
  fragment_stats::resolution = gv->a;
  fragment_stats::dims = gv->dim;

  center_box(&box);
  fragment_stats result(box);
  return result;
}

fragment_stats compute_fragment_stats(
    geometric_object_list geom_, meep::grid_volume *gv, vector3 cell_size, vector3 cell_center,
    material_type default_mat, std::vector<dft_data> dft_data_list_,
    std::vector<meep::volume> pml_1d_vols_, std::vector<meep::volume> pml_2d_vols_,
    std::vector<meep::volume> pml_3d_vols_, std::vector<meep::volume> absorber_vols_,
    material_type_list extra_materials_, double tol, int maxeval, bool ensure_per,
    bool eps_averaging) {

  fragment_stats::geom = geom_;
  fragment_stats::dft_data_list = dft_data_list_;
  fragment_stats::pml_1d_vols = pml_1d_vols_;
  fragment_stats::pml_2d_vols = pml_2d_vols_;
  fragment_stats::pml_3d_vols = pml_3d_vols_;
  fragment_stats::absorber_vols = absorber_vols_;
  fragment_stats::extra_materials = extra_materials_;
  fragment_stats::eps_averaging = eps_averaging;

  init_libctl(default_mat, ensure_per, gv, cell_size, cell_center, &geom_);
  geom_box box = make_box_from_cell(cell_size);
  fragment_stats stats = init_stats(box, tol, maxeval, gv);
  stats.compute();
  return stats;
}

fragment_stats::fragment_stats(geom_box &bx)
    : num_anisotropic_eps_pixels(0), num_anisotropic_mu_pixels(0), num_nonlinear_pixels(0),
      num_susceptibility_pixels(0), num_nonzero_conductivity_pixels(0), num_1d_pml_pixels(0),
      num_2d_pml_pixels(0), num_3d_pml_pixels(0), num_dft_pixels(0), num_pixels_in_box(0), box(bx) {

  num_pixels_in_box = get_pixels_in_box(&bx);
}

void init_libctl(material_type default_mat, bool ensure_per, meep::grid_volume *gv,
                 vector3 cell_size, vector3 cell_center, geometric_object_list *geom_) {
  geom_initialize();
  set_default_material(default_mat);
  ensure_periodicity = ensure_per;
  geometry_center = cell_center;
  dimensions = meep::number_of_directions(gv->dim);
  geometry_lattice.size = cell_size;
  geom_fix_object_list(*geom_);
}

bool fragment_stats::has_non_medium_material() {
  for (int i = 0; i < geom.num_items; ++i) {
    material_type mat = (material_type)geom.items[i].material;
    if (mat->which_subclass != material_data::MEDIUM) { return true; }
  }
  if (((material_type)default_material)->which_subclass != material_data::MEDIUM) { return true; }
  return false;
}

void fragment_stats::update_stats_from_material(material_type mat, size_t pixels,
                                                bool anisotropic_pixels_already_added) {
  switch (mat->which_subclass) {
    case material_data::MEDIUM: {
      medium_struct *med = &mat->medium;
      if (!anisotropic_pixels_already_added) { count_anisotropic_pixels(med, pixels); }
      count_nonlinear_pixels(med, pixels);
      count_susceptibility_pixels(med, pixels);
      count_nonzero_conductivity_pixels(med, pixels);
      break;
    }
    case material_data::MATERIAL_USER: {
      bool anisotropic_pixels_extmat_already_added = false;
      bool nonlinear_pixels_extmat_already_added = false;
      bool susceptibility_pixels_extmat_already_added = false;
      bool nonzero_conductivity_pixels_extmat_already_added = false;
      for (int i = 0; i < extra_materials.num_items; ++i) {
        medium_struct *med = &extra_materials.items[i]->medium;
        if (!anisotropic_pixels_already_added && !anisotropic_pixels_extmat_already_added) {
          anisotropic_pixels_extmat_already_added = count_anisotropic_pixels(med, pixels);
        }
        if (!nonlinear_pixels_extmat_already_added) {
          nonlinear_pixels_extmat_already_added = count_nonlinear_pixels(med, pixels);
        }
        if (!susceptibility_pixels_extmat_already_added) {
          susceptibility_pixels_extmat_already_added = count_susceptibility_pixels(med, pixels);
        }
        if (!nonzero_conductivity_pixels_extmat_already_added) {
          nonzero_conductivity_pixels_extmat_already_added =
              count_nonzero_conductivity_pixels(med, pixels);
        }
      }
      break;
    }
    default:
      // Only Medium and Material Function are supported
      return;
  }
}

void fragment_stats::compute_stats() {

  if (geom.num_items == 0) {
    // If there is no geometry, count the default material for the whole fragment
    update_stats_from_material((material_type)default_material, num_pixels_in_box);
  }

  for (int i = 0; i < geom.num_items; ++i) {
    geometric_object *go = &geom.items[i];
    // tolerance and max number of function evaluations of numerical quadrature are increased and
    // decreased from default values of 0.0001 and 100000, respectively, to obtain fast, approximate
    // result
    double overlap = box_overlap_with_object(box, *go, 0.05 /* tol */, 1000 /* maxeval */);

    bool anisotropic_pixels_already_added = false;
    if (eps_averaging) {
      // If the object doesn't overlap the entire box, that implies that
      // an object interface intercepts the box, which means we treat
      // the entire box as anisotropic. This method could give some false
      // positives if there is another object with the same material behind
      // the current object, but in practice it is probably reasonable to
      // assume that there is a material interface somwhere in the box so
      // we won't worry about fancier edge-detection methods for now.
      if (overlap != 1.0) {
        anisotropic_pixels_already_added = true;
        num_anisotropic_eps_pixels += num_pixels_in_box;
        if (mu_not_1(go->material)) { num_anisotropic_mu_pixels += num_pixels_in_box; }
      }
    }

    // Count contributions from material of object
    size_t pixels = (size_t)ceil(overlap * num_pixels_in_box);
    if (pixels > 0) {
      material_type mat = (material_type)go->material;
      update_stats_from_material(mat, pixels, anisotropic_pixels_already_added);
    }

    // Count contributions from default_material
    size_t default_material_pixels = num_pixels_in_box - pixels;
    if (default_material_pixels > 0) {
      update_stats_from_material((material_type)default_material, default_material_pixels,
                                 anisotropic_pixels_already_added);
    }
  }
}

bool fragment_stats::count_anisotropic_pixels(medium_struct *med, size_t pixels) {
  size_t eps_offdiag_elements = 0;
  size_t mu_offdiag_elements = 0;

  if (med->epsilon_offdiag.x.re != 0) { eps_offdiag_elements++; }
  if (med->epsilon_offdiag.y.re != 0) { eps_offdiag_elements++; }
  if (med->epsilon_offdiag.z.re != 0) { eps_offdiag_elements++; }
  if (med->mu_offdiag.x.re != 0) { mu_offdiag_elements++; }
  if (med->mu_offdiag.y.re != 0) { mu_offdiag_elements++; }
  if (med->mu_offdiag.z.re != 0) { mu_offdiag_elements++; }

  num_anisotropic_eps_pixels += eps_offdiag_elements * pixels;
  num_anisotropic_mu_pixels += mu_offdiag_elements * pixels;
  return (eps_offdiag_elements != 0) || (mu_offdiag_elements != 0);
}

bool fragment_stats::count_nonlinear_pixels(medium_struct *med, size_t pixels) {
  size_t nonzero_chi_elements = 0;

  if (med->E_chi2_diag.x != 0) { nonzero_chi_elements++; }
  if (med->E_chi2_diag.y != 0) { nonzero_chi_elements++; }
  if (med->E_chi2_diag.z != 0) { nonzero_chi_elements++; }
  if (med->E_chi3_diag.x != 0) { nonzero_chi_elements++; }
  if (med->E_chi3_diag.y != 0) { nonzero_chi_elements++; }
  if (med->E_chi3_diag.z != 0) { nonzero_chi_elements++; }
  if (med->H_chi2_diag.x != 0) { nonzero_chi_elements++; }
  if (med->H_chi2_diag.y != 0) { nonzero_chi_elements++; }
  if (med->H_chi2_diag.z != 0) { nonzero_chi_elements++; }
  if (med->H_chi3_diag.x != 0) { nonzero_chi_elements++; }
  if (med->H_chi3_diag.y != 0) { nonzero_chi_elements++; }
  if (med->H_chi3_diag.z != 0) { nonzero_chi_elements++; }

  num_nonlinear_pixels += nonzero_chi_elements * pixels;
  return nonzero_chi_elements != 0;
}

bool fragment_stats::count_susceptibility_pixels(medium_struct *med, size_t pixels) {
  num_susceptibility_pixels += med->E_susceptibilities.size() * pixels;
  num_susceptibility_pixels += med->H_susceptibilities.size() * pixels;
  return (med->E_susceptibilities.size() != 0) || (med->H_susceptibilities.size() != 0);
}

bool fragment_stats::count_nonzero_conductivity_pixels(medium_struct *med, size_t pixels) {
  size_t nonzero_conductivity_elements = 0;

  if (med->D_conductivity_diag.x != 0) { nonzero_conductivity_elements++; }
  if (med->D_conductivity_diag.y != 0) { nonzero_conductivity_elements++; }
  if (med->D_conductivity_diag.z != 0) { nonzero_conductivity_elements++; }
  if (med->B_conductivity_diag.x != 0) { nonzero_conductivity_elements++; }
  if (med->B_conductivity_diag.y != 0) { nonzero_conductivity_elements++; }
  if (med->B_conductivity_diag.z != 0) { nonzero_conductivity_elements++; }

  num_nonzero_conductivity_pixels += nonzero_conductivity_elements * pixels;
  return nonzero_conductivity_elements != 0;
}

void fragment_stats::compute_dft_stats() {

  for (size_t i = 0; i < dft_data_list.size(); ++i) {
    for (size_t j = 0; j < dft_data_list[i].vols.size(); ++j) {
      geom_box dft_box = gv2box(dft_data_list[i].vols[j]);

      if (geom_boxes_intersect(&dft_box, &box)) {
        geom_box overlap_box;
        geom_box_intersection(&overlap_box, &dft_box, &box);

        // Note: Since geom_boxes_intersect returns true if two planes share a line or two volumes
        // share a line or plane, there are cases where some pixels are counted multiple times.
        size_t overlap_pixels = get_pixels_in_box(&overlap_box, 2);
        const size_t num_freqs = dft_data_list[i].num_freqs;
        const size_t pade_samples = dft_data_list[i].pade_samples;
        // A Padé-enabled DFT chunk stores a bounded N*pade_samples history. Its
        // corrected value and previous-diagnostic caches can each add N*num_freqs
        // entries. Count all possible allocations so pre-run chunk balancing and
        // memory estimates do not depend on when or how often diagnostics are read.
        const size_t entries_per_component =
            num_freqs + (pade_samples ? pade_samples + 2 * num_freqs : 0);
        num_dft_pixels += overlap_pixels * entries_per_component * dft_data_list[i].num_components;
      }
    }
  }
}

void fragment_stats::compute_pml_stats() {

  const std::vector<meep::volume> *pml_vols[] = {&pml_1d_vols, &pml_2d_vols, &pml_3d_vols};
  size_t *pml_pixels[] = {&num_1d_pml_pixels, &num_2d_pml_pixels, &num_3d_pml_pixels};

  for (int j = 0; j < 3; ++j) {
    for (size_t i = 0; i < pml_vols[j]->size(); ++i) {
      geom_box pml_box = gv2box((*pml_vols[j])[i]);

      if (geom_boxes_intersect(&pml_box, &box)) {
        geom_box overlap_box;
        geom_box_intersection(&overlap_box, &pml_box, &box);
        size_t overlap_pixels = get_pixels_in_box(&overlap_box, 1);
        *pml_pixels[j] += overlap_pixels;
      }
    }
  }
}

void fragment_stats::compute_absorber_stats() {

  for (size_t i = 0; i < absorber_vols.size(); ++i) {
    geom_box absorber_box = gv2box(absorber_vols[i]);

    if (geom_boxes_intersect(&absorber_box, &box)) {
      geom_box overlap_box;
      geom_box_intersection(&overlap_box, &absorber_box, &box);
      size_t overlap_pixels = get_pixels_in_box(&overlap_box, 1);
      num_nonzero_conductivity_pixels += overlap_pixels;
    }
  }
}

void fragment_stats::compute() {
  compute_stats();
  compute_dft_stats();
  compute_pml_stats();
  compute_absorber_stats();
}

// Return the estimated time in seconds this fragment will take to run
// based on a cost function obtained via linear regression on a dataset
// of random simulations.
double fragment_stats::cost() const {
  return (num_anisotropic_eps_pixels * 1.15061674e-04 + num_anisotropic_mu_pixels * 1.26843801e-04 +
          num_nonlinear_pixels * 1.67029547e-04 + num_susceptibility_pixels * 2.24790864e-04 +
          num_nonzero_conductivity_pixels * 4.61260934e-05 + num_dft_pixels * 1.47283950e-04 +
          num_1d_pml_pixels * 9.92955372e-05 + num_2d_pml_pixels * 1.36901107e-03 +
          num_3d_pml_pixels * 6.63939607e-04 + num_pixels_in_box * 3.46518274e-04);
}

void fragment_stats::print_stats() const {
  master_printf("Fragment stats\n");
  master_printf("  anisotropic_eps: %zu\n", num_anisotropic_eps_pixels);
  master_printf("  anisotropic_mu: %zu\n", num_anisotropic_mu_pixels);
  master_printf("  nonlinear: %zu\n", num_nonlinear_pixels);
  master_printf("  susceptibility: %zu\n", num_susceptibility_pixels);
  master_printf("  conductivity: %zu\n", num_nonzero_conductivity_pixels);
  master_printf("  pml_1d: %zu\n", num_1d_pml_pixels);
  master_printf("  pml_2d: %zu\n", num_2d_pml_pixels);
  master_printf("  pml_3d: %zu\n", num_3d_pml_pixels);
  master_printf("  dft: %zu\n", num_dft_pixels);
  master_printf("  pixels_in_box: %zu\n", num_pixels_in_box);
  master_printf("  box.low:  {%f, %f, %f}\n", box.low.x, box.low.y, box.low.z);
  master_printf("  box.high: {%f, %f, %f}\n\n", box.high.x, box.high.y, box.high.z);
}

dft_data::dft_data(int freqs, int components, std::vector<meep::volume> volumes)
    : dft_data(freqs, components, volumes, 0) {}

dft_data::dft_data(int freqs, int components, std::vector<meep::volume> volumes,
                   size_t pade_samples)
    : num_freqs(freqs), num_components(components), pade_samples(pade_samples), vols(volumes) {}

/***************************************************************/
// Gradient calculation routines needed for material grid
/***************************************************************/

geom_box_tree calculate_tree(const meep::volume &v, geometric_object_list g) {
  geom_fix_object_list(g);
  geom_box box = gv2box(v);
  geom_box_tree geometry_tree = create_geom_box_tree0(g, box);
  return geometry_tree;
}

/* convenience routine to get element of material tensors */
static std::complex<double> cvec_to_value(vector3 diag, cvector3 offdiag, int idx) {
  std::complex<double> val = std::complex<double>(0, 0);
  switch (idx) {
    case 0: val = std::complex<double>(diag.x, 0); break;
    case 1: val = std::complex<double>(offdiag.x.re, offdiag.x.im); break;
    case 2: val = std::complex<double>(offdiag.y.re, offdiag.y.im); break;
    case 3: val = std::complex<double>(offdiag.x.re, -offdiag.x.im); break;
    case 4: val = std::complex<double>(diag.y, 0); break;
    case 5: val = std::complex<double>(offdiag.z.re, offdiag.z.im); break;
    case 6: val = std::complex<double>(offdiag.y.re, -offdiag.y.im); break;
    case 7: val = std::complex<double>(offdiag.z.re, -offdiag.z.im); break;
    case 8: val = std::complex<double>(diag.z, 0); break;
    default: meep::abort("Invalid value in switch statement.");
  }
  return val;
}

/* convenience routine to get element of material tensors */
double vec_to_value(vector3 diag, vector3 offdiag, int idx) {
  double val = 0.0;
  switch (idx) {
    case 0: val = diag.x; break;
    case 1: val = offdiag.x; break;
    case 2: val = offdiag.y; break;
    case 3: val = offdiag.x; break;
    case 4: val = diag.y; break;
    case 5: val = offdiag.z; break;
    case 6: val = offdiag.y; break;
    case 7: val = offdiag.z; break;
    case 8: val = diag.z; break;
    default: meep::abort("Invalid value in switch statement.");
  }
  return val;
}

void invert_tensor(std::complex<double> t_inv[9], std::complex<double> t[9]) {

#define m(x, y) t[x * 3 + y]
#define minv(x, y) t_inv[x * 3 + y]
  std::complex<double> det = m(0, 0) * (m(1, 1) * m(2, 2) - m(2, 1) * m(1, 2)) -
                             m(0, 1) * (m(1, 0) * m(2, 2) - m(1, 2) * m(2, 0)) +
                             m(0, 2) * (m(1, 0) * m(2, 1) - m(1, 1) * m(2, 0));
  std::complex<double> invdet = 1.0 / det;
  minv(0, 0) = (m(1, 1) * m(2, 2) - m(2, 1) * m(1, 2)) * invdet;
  minv(0, 1) = (m(0, 2) * m(2, 1) - m(0, 1) * m(2, 2)) * invdet;
  minv(0, 2) = (m(0, 1) * m(1, 2) - m(0, 2) * m(1, 1)) * invdet;
  minv(1, 0) = (m(1, 2) * m(2, 0) - m(1, 0) * m(2, 2)) * invdet;
  minv(1, 1) = (m(0, 0) * m(2, 2) - m(0, 2) * m(2, 0)) * invdet;
  minv(1, 2) = (m(1, 0) * m(0, 2) - m(0, 0) * m(1, 2)) * invdet;
  minv(2, 0) = (m(1, 0) * m(2, 1) - m(2, 0) * m(1, 1)) * invdet;
  minv(2, 1) = (m(2, 0) * m(0, 1) - m(0, 0) * m(2, 1)) * invdet;
  minv(2, 2) = (m(0, 0) * m(1, 1) - m(1, 0) * m(0, 1)) * invdet;
#undef m
#undef minv
}

/* The factor by which D-conductivity scales the permittivity at `freq`.

   Continuum: Ampere's law reads D(-i w + sigma) = curl, so factoring out the
   time derivative gives 1 + i sigma / w.

   The D update actually timestepped is

     D' = ((1 - sigma dt/2) D + dcurl) / (1 + sigma dt/2)

   (step_db.cpp, with condinv from structure.cpp).  Substituting
   D ~ exp(-i w t) and dividing by exp(-i w dt/2) turns the bracket into

     -(2/dt) sin(w dt/2) i + sigma cos(w dt/2)

   i.e. w -> (2/dt) sin(w dt/2) and sigma -> sigma cos(w dt/2), so the factor
   becomes 1 + i sigma (dt/2) cot(pi freq dt).  That differs from the continuum
   form at O((freq*dt)^2); `dt = 0` selects the continuum form. */
static std::complex<double> conductivity_factor(double sigma, double freq, double dt) {
  if (sigma == 0) return std::complex<double>(1.0, 0.0);

  double inv_omega;
  if (dt > 0) {
    double half_phase = meep::pi * freq * dt;
    inv_omega = 0.5 * dt / std::tan(half_phase);
  }
  else { inv_omega = 1.0 / (2 * meep::pi * freq); }
  return std::complex<double>(1.0, sigma * inv_omega);
}

void get_chi1_tensor_disp(std::complex<double> tensor[9], const meep::vec &r, double freq,
                          geom_epsilon *geps) {
  // locate the proper material
  material_type md;
  geps->get_material_pt(md, r);
  const medium_struct *mm = &(md->medium);

  // loop over all the tensor components
  for (int i = 0; i < 9; i++) {
    std::complex<double> a, b;
    // compute first part containing conductivity
    vector3 dummy;
    dummy.x = dummy.y = dummy.z = 0.0;
    double conductivityCur = vec_to_value(mm->D_conductivity_diag, dummy, i);
    a = conductivity_factor(conductivityCur, freq, geps->dt);

    // compute lorentzian component including the instantaneous ε
    b = cvec_to_value(mm->epsilon_diag, mm->epsilon_offdiag, i);
    for (const auto &mm_susc : mm->E_susceptibilities) {
      meep::lorentzian_susceptibility sus =
          meep::lorentzian_susceptibility(mm_susc.frequency, mm_susc.gamma, mm_susc.drude);
      double sigma = vec_to_value(mm_susc.sigma_diag, mm_susc.sigma_offdiag, i);
      b += sus.chi1_discrete(freq, sigma, geps->dt);
    }

    // elementwise multiply
    tensor[i] = a * b;
  }
}

/* Does this material contribute anything the non-dispersive tensor misses? */
bool is_dispersive(material_type m) {
  if (!m) return false;
  const medium_struct *mm = &m->medium;
  return !mm->E_susceptibilities.empty() || mm->D_conductivity_diag.x != 0 ||
         mm->D_conductivity_diag.y != 0 || mm->D_conductivity_diag.z != 0;
}

/* The effective dispersive tensor of a pixel straddling an interface, at a
   filling fraction supplied by the caller.

   This is the dispersive counterpart of the `fill_override` on
   eff_chi1inv_matrix, and it exists for the same reason: a shape derivative
   needs d(chi1inv)/d(fill), and getting it from meep's own assembly is better
   than re-deriving it.  eff_chi1inv_row_disp cannot serve, because it resolves
   the material at a *point* and so knows nothing about a filling fraction.

   The construction matches what the forward solve is actually running.  The
   instantaneous part is Kottke-smoothed, exactly as structure_chunk::set_chi1inv
   does it; each pole's amplitude is taken linearly in fill from whichever
   medium carries it, exactly as eff_sigma_row now writes into the sigma arrays.
   Because chi1 is linear in sigma, a pole contributes L(freq) * weight * sigma
   and the two media's pole lists never have to be matched against each other --
   the same reason epsilon_material_grid can concatenate rather than blend. */
static bool eff_chi1inv_row_disp_fill(meep::component c, std::complex<double> chi1inv_row[3],
                                      const meep::volume &v, double freq, geom_epsilon *geps,
                                      double fill, material_type front, material_type behind) {
  bool fallback = false;
  symm_matrix chi1inv_inf, eps_inf;
  geps->eff_chi1inv_matrix(c, &chi1inv_inf, v, geps->tol, geps->maxeval, fallback, fill);
  if (fallback) return false;
  sym_matrix_invert(&eps_inf, &chi1inv_inf);

  std::complex<double> tensor[9], tensor_inv[9];
  const double inf[9] = {eps_inf.m00, eps_inf.m01, eps_inf.m02, eps_inf.m01, eps_inf.m11,
                         eps_inf.m12, eps_inf.m02, eps_inf.m12, eps_inf.m22};
  for (int i = 0; i < 9; i++)
    tensor[i] = std::complex<double>(inf[i], 0);

  const struct {
    material_type m;
    double w;
  } side[2] = {{front, fill}, {behind, 1.0 - fill}};

  for (int sd = 0; sd < 2; sd++) {
    if (!side[sd].m) continue;
    for (const auto &su : side[sd].m->medium.E_susceptibilities) {
      meep::lorentzian_susceptibility sus(su.frequency, su.gamma, su.drude);
      // chi1 is linear in sigma, so evaluate the lineshape once at unit
      // amplitude and scale.
      std::complex<double> lineshape = sus.chi1_discrete(freq, 1.0, geps->dt);
      for (int i = 0; i < 9; i++)
        tensor[i] += lineshape * side[sd].w * vec_to_value(su.sigma_diag, su.sigma_offdiag, i);
    }
  }

  vector3 zero = {0.0, 0.0, 0.0};
  for (int i = 0; i < 9; i++) {
    double cond = 0;
    for (int sd = 0; sd < 2; sd++)
      if (side[sd].m)
        cond += side[sd].w * vec_to_value(side[sd].m->medium.D_conductivity_diag, zero, i);
    tensor[i] *= conductivity_factor(cond, freq, geps->dt);
  }

  invert_tensor(tensor_inv, tensor);
  int row = 0;
  switch (component_direction(c)) {
    case meep::X:
    case meep::R: row = 0; break;
    case meep::Y:
    case meep::P: row = 1; break;
    case meep::Z: row = 2; break;
    case meep::NO_DIRECTION: return false;
  }
  for (int i = 0; i < 3; i++)
    chi1inv_row[i] = tensor_inv[3 * row + i];
  return true;
}

void eff_chi1inv_row_disp(meep::component c, std::complex<double> chi1inv_row[3],
                          const meep::vec &r, double freq, geom_epsilon *geps) {
  std::complex<double> tensor[9], tensor_inv[9];
  get_chi1_tensor_disp(tensor, r, freq, geps);
  // invert the matrix
  invert_tensor(tensor_inv, tensor);

  // get the row we care about
  switch (component_direction(c)) {
    case meep::X:
    case meep::R:
      chi1inv_row[0] = tensor_inv[0];
      chi1inv_row[1] = tensor_inv[1];
      chi1inv_row[2] = tensor_inv[2];
      break;
    case meep::Y:
    case meep::P:
      chi1inv_row[0] = tensor_inv[3];
      chi1inv_row[1] = tensor_inv[4];
      chi1inv_row[2] = tensor_inv[5];
      break;
    case meep::Z:
      chi1inv_row[0] = tensor_inv[6];
      chi1inv_row[1] = tensor_inv[7];
      chi1inv_row[2] = tensor_inv[8];
      break;
    case meep::NO_DIRECTION: chi1inv_row[0] = chi1inv_row[1] = chi1inv_row[2] = 0; break;
  }
}

std::complex<double> cond_cmp(meep::component c, const meep::vec &r, double freq,
                              geom_epsilon *geps) {
  // locate the proper material
  material_type md;
  geps->get_material_pt(md, r);
  const medium_struct *mm = &(md->medium);

  // get the row we care about
  switch (component_direction(c)) {
    case meep::X:
    case meep::R: return conductivity_factor(mm->D_conductivity_diag.x, freq, geps->dt);
    case meep::Y:
    case meep::P: return conductivity_factor(mm->D_conductivity_diag.y, freq, geps->dt);
    case meep::Z: return conductivity_factor(mm->D_conductivity_diag.z, freq, geps->dt);
    case meep::NO_DIRECTION: meep::abort("Invalid adjoint field component");
  }
}

std::complex<double>
get_material_gradient(const meep::vec &r,              // current point
                      const meep::component adjoint_c, // adjoint field component
                      const meep::component forward_c, // forward field component
                      std::complex<double> fields_f,   // forward field at current point
                      double freq,                     // frequency
                      geom_epsilon *geps,              // material
                      meep::grid_volume &gv,           // simulation grid volume
                      double du,                       // step size
                      const std::vector<double *> &us, // matgrid weights, one per
                                                       // copy of this design grid
                      int idx                          // matgrid index
) {
  /*Compute the Aᵤx product from the -λᵀAᵤx calculation.
  The current adjoint (λ) field component (adjoint_c)
  determines which row of Aᵤ we care about.
  The current forward (x) field component (forward_c)
  determines which column of Aᵤ we care about.

  There are two ways we can compute the required row/
  column of Aᵤ:
    1. If we want to incorporate subpixel smoothing,
    we use the eff_chi1inv_row() routine. This neglects
    conductivities, susceptibilities, etc.
    2. We use eff_chi1inv_row_disp() for all other cases
    (at the expense of not accounting for subpixel smoothing,
    if there were any).

  For now we do a finite difference approach to estimate the
  gradient of the system matrix A since it's fairly accurate,
  cheap, and easy to generalize.*/

  // locate the proper material
  material_type md;
  geps->get_material_pt(md, r);

  // get the tensor column index corresponding to the forward component
  int dir_idx = 0;
  switch (meep::component_direction(forward_c)) {
    case meep::X:
    case meep::R: dir_idx = 0; break;
    case meep::Y:
    case meep::P: dir_idx = 1; break;
    case meep::Z: dir_idx = 2; break;
    case meep::NO_DIRECTION: meep::abort("Invalid forward component!\n");
  }

  // materials are non-dispersive
  if (md->trivial) {
    const double sd = 1.0; // FIXME: make user-changable?
    meep::volume v(r);
    LOOP_OVER_DIRECTIONS(dim, d) {
      v.set_direction_min(d, r.in_direction(d) - 0.5 * gv.inva * sd);
      v.set_direction_max(d, r.in_direction(d) + 0.5 * gv.inva * sd);
    }
    double row_1[3], row_2[3];
    /* One design variable may back several copies of the same grid (the way a
       symmetry is imposed), so it has to be perturbed in every copy at once --
       perturbing one copy answers a different question. */
    std::vector<double> orig(us.size());
    for (size_t k = 0; k < us.size(); ++k) {
      orig[k] = us[k][idx];
      us[k][idx] -= du;
    }
    geps->eff_chi1inv_row(adjoint_c, row_1, v, geps->tol, geps->maxeval);
    for (size_t k = 0; k < us.size(); ++k)
      us[k][idx] = orig[k] + du;
    geps->eff_chi1inv_row(adjoint_c, row_2, v, geps->tol, geps->maxeval);
    for (size_t k = 0; k < us.size(); ++k)
      us[k][idx] = orig[k];
    return fields_f * (row_1[dir_idx] - row_2[dir_idx]) / (2 * du);
  }
  // materials have some dispersion
  else {
    std::vector<double> orig(us.size());
    std::complex<double> row_1[3], row_2[3], dA_du[3];
    for (size_t k = 0; k < us.size(); ++k) {
      orig[k] = us[k][idx];
      us[k][idx] -= du;
    }
    eff_chi1inv_row_disp(adjoint_c, row_1, r, freq, geps);
    for (size_t k = 0; k < us.size(); ++k)
      us[k][idx] = orig[k] + du;
    eff_chi1inv_row_disp(adjoint_c, row_2, r, freq, geps);
    for (size_t k = 0; k < us.size(); ++k)
      us[k][idx] = orig[k];

    return fields_f * (row_1[dir_idx] - row_2[dir_idx]) / (2 * du) *
           cond_cmp(forward_c, r, freq, geps);
  }
}

/* A brute force approach to calculating Aᵤ using finite differences.
Past versions of the code only calculated dA/dε using a finite
difference, and then multiplied by the analytic vJp (dε/du).
With the addition of subpixel smoothing, however, the vJp became
much more complicated and it is easier to calculate the entire gradient
using finite differences (at the cost of slightly less accurate gradients
due to floating-point roundoff errors). */
void add_interpolate_weights(const std::vector<double *> &udatas, const std::vector<vector3> &pbs,
                             double *data, int nx, int ny, int nz, int stride, double scaleby,
                             int ukind, double uval, meep::vec r, geom_epsilon *geps,
                             meep::component adjoint_c, meep::component forward_c,
                             std::complex<double> fwd, std::complex<double> adj, double freq,
                             meep::grid_volume &gv, double du) {
  /* `udatas`/`pbs` are the copies of *one* design grid that cover this point,
     with the box coordinates of each. A design variable is shared across all of
     them, so it is differentiated once, with every copy perturbed together, and
     the union of the copies' interpolation stencils is the set of variables this
     point touches. */

#define IDX(x, y, z) (((x)*ny + (y)) * nz + (z)) * stride

  std::set<int> touched;
  for (size_t k = 0; k < pbs.size(); ++k) {
    int x1, y1, z1, x2, y2, z2;
    double dx, dy, dz;
    meep::map_coordinates(pbs[k].x, pbs[k].y, pbs[k].z, nx, ny, nz, x1, y1, z1, x2, y2, z2, dx, dy,
                          dz);
    int x_list[2] = {x1, x2}, y_list[2] = {y1, y2}, z_list[2] = {z1, z2};
    int lx = (x1 == x2) ? 1 : 2;
    int ly = (y1 == y2) ? 1 : 2;
    int lz = (z1 == z2) ? 1 : 2;

    if (k == 0 && (ukind == material_data::U_MIN || ukind == material_data::U_PROD)) {
      /* U_MIN and U_PROD carry adjustments defined in terms of the interpolated
         value at this point. Overlapping grids of these kinds are rejected by
         the caller, so there is exactly one copy and this is unchanged. */
      const double *U = udatas[k];
      double u = (((U[IDX(x1, y1, z1)] * (1.0 - dx) + U[IDX(x2, y1, z1)] * dx) * (1.0 - dy) +
                   (U[IDX(x1, y2, z1)] * (1.0 - dx) + U[IDX(x2, y2, z1)] * dx) * dy) *
                      (1.0 - dz) +
                  ((U[IDX(x1, y1, z2)] * (1.0 - dx) + U[IDX(x2, y1, z2)] * dx) * (1.0 - dy) +
                   (U[IDX(x1, y2, z2)] * (1.0 - dx) + U[IDX(x2, y2, z2)] * dx) * dy) *
                      dz);
      if (ukind == material_data::U_MIN && u != uval) return; // TODO look into this
      if (ukind == material_data::U_PROD) scaleby *= uval / u;
    }

    for (int xi = 0; xi < lx; xi++)
      for (int yi = 0; yi < ly; yi++)
        for (int zi = 0; zi < lz; zi++)
          touched.insert(IDX(x_list[xi], y_list[yi], z_list[zi]));
  }

  for (std::set<int>::const_iterator it = touched.begin(); it != touched.end(); ++it) {
    std::complex<double> prod =
        adj * get_material_gradient(r, adjoint_c, forward_c, fwd, freq, geps, gv, du, udatas, *it);
    data[*it] += prod.real() * scaleby;
  }

#undef IDX
}

void material_grids_addgradient_point(double *v, vector3 p, double scalegrad, geom_epsilon *geps,
                                      meep::component adjoint_c, meep::component forward_c,
                                      std::complex<double> fwd, std::complex<double> adj,
                                      double freq, meep::grid_volume &gv, double tol,
                                      long owner_grid_id) {
  geom_box_tree tp;
  int oi, ois;
  material_data *mg, *mg_sum;
  double uval;
  int kind;
  tp = geom_tree_search(p, geps->geometry_tree, &oi);

  if (tp &&
      ((material_type)tp->objects[oi].o->material)->which_subclass == material_data::MATERIAL_GRID)
    mg = (material_data *)tp->objects[oi].o->material;
  else if (!tp && ((material_type)default_material)->which_subclass == material_data::MATERIAL_GRID)
    mg = (material_data *)default_material;
  else
    return; /* no material grids at this point */

  // Calculate the number of material grids if we are averaging values
  if ((tp) && ((kind = mg->material_grid_kinds) == material_data::U_MEAN)) {
    int matgrid_val_count = 0;
    geom_box_tree tp_sum;
    tp_sum = geom_tree_search(p, geps->geometry_tree, &ois);
    mg_sum = (material_data *)tp_sum->objects[ois].o->material;
    do {
      tp_sum = geom_tree_search_next(p, tp_sum, &ois);
      ++matgrid_val_count;
      if (tp_sum) mg_sum = (material_data *)tp_sum->objects[ois].o->material;
    } while (tp_sum && is_material_grid(mg_sum));
    /* No scalegrad /= matgrid_val_count here: the finite difference in
       get_material_gradient() runs through matgrid_val(), so the 1/N of the
       U_MEAN average is already contained in the result. */
    (void)matgrid_val_count;
  }
  else if ((tp) && ((mg->material_grid_kinds == material_data::U_MIN) ||
                    (mg->material_grid_kinds == material_data::U_PROD))) {
    meep::abort("%s:%i:material_grids_addgradient_point does not support overlapping "
                "MATERIAL_GRIDs with U_MIN or U_PROD.\n",
                __FILE__, __LINE__);
  }

  // Iterate through grids and add weights as needed
  if (tp) {
    /*NOTE in the future, it may be nice to be able to have *different*
    material grids in a particular design region. This would require checking each
    material grid here and iterating to the next spot in a large array of
    design parameters (see MPB).

    For now, it's actually easier just to assign each "unique" material grid its
    own design region. Unlike MPB, we are only iterating over the grid points inside
    that design region. Note that we can still have multiple copies of
    the same design grid (i.e. for symmetries). Thats why we are looping in this
    fashion. Since we aren't checking if each design grid is unique, however,
    it's up to the user to only have one unique design grid in this volume.*/
    vector3 sz = mg->grid_size;
    uval = tanh_projection(matgrid_val(p, tp, oi, mg), mg->beta, mg->eta);

    /* Collect the copies of this design region's own grid that cover the point.

       Only the owner's copies: v[] is sized and indexed by the owner's
       parameterization, so a contribution from any other grid is mis-attributed,
       and an out-of-bounds write when that grid is larger. That happens whenever
       two grids abut, because the design region's DFT monitor is snapped outward
       to enclosing grid nodes and so covers nodes shared with the neighbour. The
       neighbouring design region visits them on its own pass, so nothing is
       lost. All copies of one grid share a grid_id, so a grid overlapped with
       transformed copies of itself -- the documented way to impose a symmetry --
       still contributes through every copy. */
    std::vector<double *> udatas;
    std::vector<vector3> pbs;
    geom_box_tree tpc = tp;
    int oic = oi;
    do {
      material_data *mg_cur = (material_data *)tpc->objects[oic].o->material;
      if (owner_grid_id < 0 || mg_cur->grid_id == owner_grid_id) {
        if (udatas.empty()) sz = mg_cur->grid_size;
        udatas.push_back(mg_cur->weights);
        pbs.push_back(to_geom_box_coords(p, &tpc->objects[oic]));
      }
      if (kind == material_data::U_DEFAULT) break;
      tpc = geom_tree_search_next(p, tpc, &oic);
    } while (tpc && is_material_grid((material_data *)tpc->objects[oic].o->material));

    if (!udatas.empty())
      add_interpolate_weights(udatas, pbs, v, sz.x, sz.y, sz.z, 1, scalegrad, kind, uval,
                              vector3_to_vec(p), geps, adjoint_c, forward_c, fwd, adj, freq, gv,
                              tol);
  }
  // no object tree -- the whole domain is the material grid
  if (!tp && is_material_grid(default_material)) {
    map_lattice_coordinates(p.x, p.y, p.z);
    vector3 sz = mg->grid_size;
    uval = tanh_projection(material_grid_val(p, mg), mg->beta, mg->eta);
    std::vector<double *> udatas(1, mg->weights);
    std::vector<vector3> pbs(1, p);
    add_interpolate_weights(udatas, pbs, v, sz.x, sz.y, sz.z, 1, scalegrad, kind, uval,
                            vector3_to_vec(p), geps, adjoint_c, forward_c, fwd, adj, freq, gv, tol);
  }
}

/* ------------------------------------------------------------------ */
/* Gradients with respect to a geometric object's centre and size.      */
/* ------------------------------------------------------------------ */

/* The smoothed permittivity depends on the geometry only through two things:
   the filling fraction of each pixel and the interface normal. So

       d(chi1inv)/d(parameter) = d(chi1inv)/d(fill) * d(fill)/d(parameter)

   with the second factor analytic for an axis-aligned block, since the block
   is an intersection of three slabs and the overlap volume factorizes:

       |pixel ∩ block| = prod_i overlap_i(c_i, s_i)
       overlap_i = clamp(min(c_i + s_i/2, hi_i) - max(c_i - s_i/2, lo_i), 0, dx)

   d(overlap_i)/d(c_i) is -1, 0 or +1 depending on which face lies inside the
   pixel, and d(overlap_i)/d(s_i) is 0 or +/-1/2. Everything else is a product
   of the other axes' overlaps.

   The first factor comes from varying `fill` through eff_chi1inv_matrix rather
   than re-deriving Kottke's algebra here, which keeps this consistent with
   whatever meep actually does. That variation is pure local algebra: `delta`
   is linear in fill by construction, and only the final delta -> chi1inv map
   is nonlinear.

   Perturbing the geometry and re-differencing instead -- which is what this
   routine used to do -- means differencing box_overlap_with_object, an
   adaptive quadrature, so its tolerance is amplified by 1/step. That produced
   a gradient with no stable step regime at all.

   Only pixels straddling the boundary contribute, since d(fill)/d(parameter)
   vanishes wherever a pixel is wholly inside or wholly outside. That is the
   discrete form of a shape derivative being a surface integral, and it is why
   the loop below skips interior pixels. */

enum geom_param {
  GEOM_CENTER_X = 0,
  GEOM_CENTER_Y,
  GEOM_CENTER_Z,
  GEOM_SIZE_X,
  GEOM_SIZE_Y,
  GEOM_SIZE_Z
};

/* Overlap of [lo, hi] with the block's extent along one axis, and its
   derivatives with respect to that axis's centre and half-size. */
static void axis_overlap(double lo, double hi, double c, double s, double &overlap,
                         double &d_dcenter, double &d_dsize) {
  const double blo = c - 0.5 * s, bhi = c + 0.5 * s;
  const double left = std::max(lo, blo), right = std::min(hi, bhi);
  overlap = right - left;
  if (overlap <= 0) {
    overlap = 0;
    d_dcenter = d_dsize = 0;
    return;
  }
  /* Moving the centre moves both faces together; growing the size moves them
     apart by half each. A face contributes only in the pixel that contains it.

     The comparisons are half-open on purpose. With strict inequalities on both
     sides, a face lying exactly on a pixel boundary belongs to neither
     neighbour and its contribution vanishes -- silently, and for every pixel
     along that face, which zeroes the whole derivative for that axis. Round
     geometry on a round grid hits this constantly: a block of width 1.0 at
     resolution 20 has faces exactly 10 pixels from its centre.

     Half-open assigns such a face to exactly one of the two pixels, so nothing
     is dropped and nothing is double counted. The derivative is then one-sided
     there, which is the truth: the smoothed permittivity has a kink at a pixel
     boundary and no two-sided derivative exists. */
  const double lower_inside = (blo >= lo && blo < hi) ? 1.0 : 0.0;
  const double upper_inside = (bhi > lo && bhi <= hi) ? 1.0 : 0.0;
  d_dcenter = upper_inside - lower_inside;
  d_dsize = 0.5 * (upper_inside + lower_inside);
}

/* Is the point p inside the box [lo,hi] in *every* direction?

   grid_volume::index() is a flat inner product of the per-direction offsets
   with the strides, so a point that leaves the box in a single direction can
   still produce an index inside [0,N): in 3d, stepping one pixel off the end
   in y lands on an unrelated voxel of the next z slab. Guarding only the flat
   index therefore silently substitutes a neighboring field value. Check each
   direction separately instead. */
static bool ivec_in_box(const meep::ivec &p, const meep::ivec &lo, const meep::ivec &hi,
                        meep::ndim dim) {
  LOOP_OVER_DIRECTIONS(dim, d) {
    if (p.in_direction(d) < lo.in_direction(d) || p.in_direction(d) > hi.in_direction(d))
      return false;
  }
  return true;
}

/* Find the dft_chunk in `chunks` covering the same region as `adj`.

   loop_in_chunks() emits one dft_chunk per (fields_chunk, symmetry image,
   lattice shift) triple, so that triple identifies a chunk uniquely. It is
   *not* safe to pair chunks by position in the next_in_dft list: a fields
   chunk contributes a dft_chunk for a given component only if the monitor
   overlaps that component's owned grid there, and the yee shifts differ
   between components, so the per-component lists can have different lengths
   and different orderings near the edges of the monitor. Returns NULL when
   this fields chunk holds no data for the requested component. */
static meep::dft_chunk *matching_dft_chunk(const std::vector<meep::dft_chunk *> &chunks,
                                           const meep::dft_chunk *adj) {
  for (size_t i = 0; i < chunks.size(); ++i)
    if (chunks[i]->fc == adj->fc && chunks[i]->sn == adj->sn && chunks[i]->shift == adj->shift)
      return chunks[i];
  return NULL;
}

/* Look up the forward DFT at point p, or zero if p is outside this chunk.

   The restriction stencil reaches half a pixel (unit_a*2 cancels against the yee
   shift), so it wants the four forward nodes around the adjoint node. With the
   persist pad clamped to the component's own grid (see dft.cpp) those all lie
   inside this chunk's box, the outermost of them being the chunk's ghost layer,
   which step_boundaries() keeps current, so a lookup should only fall outside
   it at the cell boundary, where zero is the right answer. */
static std::complex<meep::realnum> forward_dft_value(const meep::dft_chunk *ch,
                                                     const std::complex<meep::realnum> *values,
                                                     const meep::grid_volume &gv_fwd,
                                                     meep::grid_volume &gv, const meep::ivec &p,
                                                     size_t nf, size_t f_i) {
  if (ivec_in_box(p, ch->is, ch->ie, gv.dim)) {
    ptrdiff_t i = gv_fwd.index(ch->c, p);
    if (i >= 0 && (size_t)i < ch->N) return values[nf * i + f_i];
  }
  return 0;
}

void geometry_addgradient(double *v, size_t nparams, size_t nf,
                          std::vector<meep::dft_fields *> fields_a,
                          std::vector<meep::dft_fields *> fields_f, double *frequencies,
                          double scalegrad, meep::grid_volume &gv, geom_epsilon *geps,
                          int object_index, int *params, double du) {
  (void)du;
  if (object_index < 0 || object_index >= geps->geometry.num_items)
    meep::abort("geometry_addgradient: object index %d out of range (%d objects)", object_index,
                geps->geometry.num_items);
  geometric_object *obj = &geps->geometry.items[object_index];
  if (obj->which_subclass != geometric_object::BLOCK)
    meep::abort("geometry_addgradient: only blocks are supported");

  const vector3 bc = obj->center;
  const vector3 bs = obj->subclass.block_data->size;
  const double centers[3] = {bc.x, bc.y, bc.z};
  const double sizes[3] = {bs.x, bs.y, bs.z};

  std::vector<meep::dft_chunk *> adjoint_chunks[3], forward_chunks[3];
  for (int i = 0; i < 3; i++) {
    for (meep::dft_chunk *c = fields_a[i]->chunks; c; c = c->next_in_dft)
      adjoint_chunks[i].push_back(c);
    for (meep::dft_chunk *c = fields_f[i]->chunks; c; c = c->next_in_dft)
      forward_chunks[i].push_back(c);
  }

  std::vector<double> local(nf * nparams, 0.0);
  /* `fill` is dimensionless and O(1) and `delta` is linear in it, so this step
     is well conditioned -- unlike a step in a length, which has to be compared
     against the pixel. */
  const double dfill = 1e-6;

  for (size_t f_i = 0; f_i < nf; f_i++) {
    for (int ci_adjoint = 0; ci_adjoint < 3; ci_adjoint++) {
      int num_chunks = adjoint_chunks[ci_adjoint].size();
      if (num_chunks == 0) continue;

      for (int cur = 0; cur < num_chunks; cur++) {
        meep::dft_chunk *adj_chunk = adjoint_chunks[ci_adjoint][cur];
        meep::component adjoint_c = adj_chunk->c;
        meep::grid_volume gv_adj = gv.subvolume(adj_chunk->is, adj_chunk->ie, adjoint_c);

        for (int ci_forward = 0; ci_forward < 3; ci_forward++) {
          /* Pair the forward chunk with this adjoint chunk *spatially*.
             Indexing both lists by `cur` assumes they are the same length and
             in the same order, and they need not be -- which is why
             matching_dft_chunk exists a few functions up, for exactly this
             purpose in material_grids_addgradient.  A foreign chunk here does
             not merely mix the gradient up: it then gets indexed with this
             chunk's offsets and read past the end of its array. */
          if (forward_chunks[ci_forward].empty()) continue;
          meep::dft_chunk *fwd_chunk =
              matching_dft_chunk(forward_chunks[ci_forward], adj_chunk);
          if (!fwd_chunk) continue;
          meep::component forward_c = fwd_chunk->c;
          meep::grid_volume gv_fwd = gv.subvolume(fwd_chunk->is, fwd_chunk->ie, forward_c);

          int dir_idx;
          switch (meep::component_direction(forward_c)) {
            case meep::X:
            case meep::R: dir_idx = 0; break;
            case meep::Y:
            case meep::P: dir_idx = 1; break;
            case meep::Z: dir_idx = 2; break;
            default: continue;
          }

          meep::ivec loop_is = adj_chunk->persist ? adj_chunk->is_old : adj_chunk->is;
          meep::ivec loop_ie = adj_chunk->persist ? adj_chunk->ie_old : adj_chunk->ie;

          LOOP_OVER_IVECS(gv_adj, loop_is, loop_ie, idx_adj) {
            IVEC_LOOP_ILOC(gv_adj, ip);
            IVEC_LOOP_LOC(gv_adj, p);
            std::complex<meep::realnum> adj = adj_chunk->dft[nf * idx_adj + f_i];
            if (adj == 0.0) continue;

            /* The smoothed tensor is diagonal in Cartesian axes only where the
               interface normal lies along an axis, which for a block is true on
               a face and false on an edge or corner. Those pixels are O(N) of
               the O(N^2) on the boundary, so their share falls with resolution
               -- which is the signature of the residual error here. Contracting
               them needs the forward field of another component, restricted to
               the two epsilon nodes between the pair. */
            std::complex<meep::realnum> fwd;
            meep::vec eps_at = p;
            double node_weight = 1.0;
            int num_nodes = 1;
            meep::vec node_pos[2];
            std::complex<meep::realnum> node_fwd[2];
            if (forward_c == adjoint_c) {
              node_pos[0] = p;
              /* `idx_adj` is an offset into gv_adj, so it cannot be used to
                 index the forward chunk.  Look the value up by location, which
                 bounds-checks it as the two lookups below already did. */
              node_fwd[0] =
                  forward_dft_value(fwd_chunk, fwd_chunk->dft, gv_fwd, gv, ip, nf, f_i);
            }
            else {
              num_nodes = 2;
              node_weight = 0.5;
              meep::ivec fwd_p = ip + gv.iyee_shift(forward_c) - gv.iyee_shift(adjoint_c);
              meep::ivec unit_a = unit_ivec(gv.dim, component_direction(adjoint_c));
              meep::ivec unit_f = unit_ivec(gv.dim, component_direction(forward_c));
              meep::ivec pl[2] = {fwd_p, fwd_p + unit_a * 2};
              meep::ivec pr[2] = {fwd_p - unit_f * 2, fwd_p + unit_a * 2 - unit_f * 2};
              for (int nd = 0; nd < 2; nd++) {
                ptrdiff_t i1 = gv_fwd.index(forward_c, pl[nd]);
                ptrdiff_t i2 = gv_fwd.index(forward_c, pr[nd]);
                std::complex<meep::realnum> f1 =
                    ((i1 >= fwd_chunk->N) || (i1 < 0)) ? 0 : fwd_chunk->dft[nf * i1 + f_i];
                std::complex<meep::realnum> f2 =
                    ((i2 >= fwd_chunk->N) || (i2 < 0)) ? 0 : fwd_chunk->dft[nf * i2 + f_i];
                node_fwd[nd] = std::complex<meep::realnum>(0.5, 0) * (f1 + f2);
                node_pos[nd] = gv[(pl[nd] + pr[nd]) / 2];
              }
            }

            for (int node = 0; node < num_nodes; node++) {
              fwd = node_fwd[node];
              if (fwd == 0.0) continue;
              const meep::vec p_node = node_pos[node];

              meep::volume voxel(p_node);
              LOOP_OVER_DIRECTIONS(gv.dim, d) {
                voxel.set_direction_min(d, p_node.in_direction(d) - 0.5 * gv.inva);
                voxel.set_direction_max(d, p_node.in_direction(d) + 0.5 * gv.inva);
              }

              /* Skip pixels that do not straddle the boundary: d(fill)/d(param)
                 is zero there, so they contribute nothing. */
              /* `fill` is the fraction of the pixel inside whichever object
                 get_front_object selected. If that is not the object being
                 differentiated, d(fill)/d(our parameters) is not what this pixel
                 responds to -- and taking it anyway gets the sign wrong wherever
                 the front object is the background instead. */
              double fill;
              const geometric_object *front = NULL;
              vector3 shiftby = {0, 0, 0};
              material_type mat_front = NULL, mat_behind = NULL;
              if (!geps->interface_fill(voxel, geps->tol, geps->maxeval, fill, &front, &shiftby,
                                        &mat_front, &mat_behind))
                continue;
              if (front != obj) continue;
              if (fill <= 0.0 || fill >= 1.0) continue;

              /* d(fill)/d(parameter), analytic and separable. */
              double overlap[3], d_dc[3], d_ds[3], extent[3];
              bool degenerate = false;
              for (int ax = 0; ax < 3; ax++) {
                meep::direction dd = (ax == 0) ? meep::X : (ax == 1) ? meep::Y : meep::Z;
                const bool resolved = (gv.dim == meep::D3) || (gv.dim == meep::D2 && ax < 2) ||
                                      (gv.dim == meep::D1 && ax == 2);
                if (!resolved) {
                  /* A dimension the simulation does not resolve. The block is
                     effectively infinite along it, so it contributes a factor of
                     one to the overlap and nothing to the derivative -- note the
                     block's own size along such an axis is typically zero, so
                     using it here would make every pixel degenerate. */
                  extent[ax] = 1.0;
                  overlap[ax] = 1.0;
                  d_dc[ax] = d_ds[ax] = 0.0;
                  continue;
                }
                const double sh = (ax == 0) ? shiftby.x : (ax == 1) ? shiftby.y : shiftby.z;
                const double lo = voxel.in_direction_min(dd) - sh;
                const double hi = voxel.in_direction_max(dd) - sh;
                extent[ax] = hi - lo;
                axis_overlap(lo, hi, centers[ax], sizes[ax], overlap[ax], d_dc[ax], d_ds[ax]);
                if (overlap[ax] <= 0) degenerate = true;
              }
              if (degenerate) continue;

              double pixel_volume = 1.0;
              for (int ax = 0; ax < 3; ax++)
                pixel_volume *= extent[ax];

              const double fill_lo = std::max(0.0, fill - dfill);
              const double fill_hi = std::min(1.0, fill + dfill);
              const double actual_dfill = fill_hi - fill_lo;
              if (actual_dfill <= 0) continue;

              int row_of = 0;
              switch (meep::component_direction(adjoint_c)) {
                case meep::X:
                case meep::R: row_of = 0; break;
                case meep::Y:
                case meep::P: row_of = 1; break;
                case meep::Z: row_of = 2; break;
                default: continue;
              }

              /* d(chi1inv)/d(fill), from meep's own tensor assembly.

                 Which assembly depends on the materials at *this* interface,
                 not on md->trivial of the front object: a plain dielectric
                 sitting in front of a metal would otherwise take the
                 non-dispersive branch and drop the entire contribution at
                 exactly the pixels that carry the derivative. */
              std::complex<double> dchi_dfill;
              if (is_dispersive(mat_front) || is_dispersive(mat_behind)) {
                std::complex<double> lo_d[3], hi_d[3];
                if (!eff_chi1inv_row_disp_fill(adjoint_c, lo_d, voxel, frequencies[f_i], geps,
                                               fill_lo, mat_front, mat_behind))
                  continue;
                if (!eff_chi1inv_row_disp_fill(adjoint_c, hi_d, voxel, frequencies[f_i], geps,
                                               fill_hi, mat_front, mat_behind))
                  continue;
                dchi_dfill = (hi_d[dir_idx] - lo_d[dir_idx]) / actual_dfill;
              }
              else {
                bool fb_lo = false, fb_hi = false;
                symm_matrix m_lo, m_hi;
                geps->eff_chi1inv_matrix(adjoint_c, &m_lo, voxel, geps->tol, geps->maxeval, fb_lo,
                                         fill_lo);
                geps->eff_chi1inv_matrix(adjoint_c, &m_hi, voxel, geps->tol, geps->maxeval, fb_hi,
                                         fill_hi);
                if (fb_lo || fb_hi) continue;

                const double lo_row[3] = {m_lo.m00, m_lo.m01, m_lo.m02};
                const double hi_row[3] = {m_hi.m00, m_hi.m01, m_hi.m02};
                const double lo_row1[3] = {m_lo.m01, m_lo.m11, m_lo.m12};
                const double hi_row1[3] = {m_hi.m01, m_hi.m11, m_hi.m12};
                const double lo_row2[3] = {m_lo.m02, m_lo.m12, m_lo.m22};
                const double hi_row2[3] = {m_hi.m02, m_hi.m12, m_hi.m22};
                const double *lo_r = (row_of == 0) ? lo_row : (row_of == 1) ? lo_row1 : lo_row2;
                const double *hi_r = (row_of == 0) ? hi_row : (row_of == 1) ? hi_row1 : hi_row2;
                dchi_dfill = (hi_r[dir_idx] - lo_r[dir_idx]) / actual_dfill;
              }

              const std::complex<double> pair =
                  std::complex<double>(double(adj.real()), double(adj.imag())) *
                  std::complex<double>(double(fwd.real()), double(fwd.imag()));
              const double cyl_scale = (gv.dim == meep::Dcyl) ? 2 * p_node.r() : 1;

              for (size_t ip = 0; ip < nparams; ip++) {
                const int which = params[ip];
                const int ax = which % 3;
                double dfill_dp = (which < 3) ? d_dc[ax] : d_ds[ax];
                if (dfill_dp == 0.0) continue;
                for (int other = 0; other < 3; other++)
                  if (other != ax) dfill_dp *= overlap[other];
                dfill_dp /= pixel_volume;

                /* the leading minus matches get_material_gradient's convention,
                   which returns -(d row/d parameter) */
                local[nparams * f_i + ip] -=
                    node_weight * scalegrad * cyl_scale * dfill_dp * std::real(dchi_dfill * pair);
              }
            } // node
          }
        }
      }
    }
  }

  meep::sum_to_all(local.data(), v, int(nf * nparams));
}

void material_grids_addgradient(double *v, size_t ng, size_t nf,
                                std::vector<meep::dft_fields *> fields_a,
                                std::vector<meep::dft_fields *> fields_f, double *frequencies,
                                double scalegrad, meep::grid_volume &gv, geom_epsilon *geps,
                                long owner_grid_id, double du) {
  /* ------------------------------------------------------------ */
  // initialize local gradient array
  /* ------------------------------------------------------------ */
  double *v_local = new double[ng * nf];
  for (int i = 0; i < ng * nf; i++) {
    v_local[i] = 0;
  }

  /* ------------------------------------------------------------ */
  // store chunk info in vectors for simplicity
  /* ------------------------------------------------------------ */
  std::vector<std::vector<meep::dft_chunk *> > adjoint_dft_chunks;
  std::vector<std::vector<meep::dft_chunk *> > forward_dft_chunks;
  std::vector<std::vector<const std::complex<meep::realnum> *> > adjoint_dft_values;
  for (int i = 0; i < 3; i++) {
    std::vector<meep::dft_chunk *> c_adjoint_dft_chunks;
    std::vector<meep::dft_chunk *> c_forward_dft_chunks;
    meep::dft_chunk *current_adjoint_chunk = fields_a[i]->chunks;
    meep::dft_chunk *current_forward_chunk = fields_f[i]->chunks;
    while (current_adjoint_chunk) {
      if (current_adjoint_chunk->omega.size() != nf)
        meep::abort("Supplied frequencies %d don't match dft frequencies %d\n", nf,
                    current_adjoint_chunk->omega.size());
      c_adjoint_dft_chunks.push_back(current_adjoint_chunk);
      current_adjoint_chunk = current_adjoint_chunk->next_in_dft;
    }
    while (current_forward_chunk) {
      if (current_forward_chunk->omega.size() != nf)
        meep::abort("Supplied frequencies %d don't match dft frequencies %d\n", nf,
                    current_forward_chunk->omega.size());
      c_forward_dft_chunks.push_back(current_forward_chunk);
      current_forward_chunk = current_forward_chunk->next_in_dft;
    }
    /* NOTE: the adjoint and forward chunk counts may legitimately differ, both
       from each other and between components -- see matching_dft_chunk(). The
       chunks are paired spatially below, so this is not an error. */
    adjoint_dft_chunks.push_back(c_adjoint_dft_chunks);
    forward_dft_chunks.push_back(c_forward_dft_chunks);

    std::vector<const std::complex<meep::realnum> *> c_adjoint_dft_values;
    for (const meep::dft_chunk *chunk : c_adjoint_dft_chunks)
      c_adjoint_dft_values.push_back(chunk->dft_values());
    adjoint_dft_values.push_back(c_adjoint_dft_values);
  }

  /* ------------------------------------------------------------ */
  // Begin looping
  /* ------------------------------------------------------------ */

  // loop over frequency
  for (size_t f_i = 0; f_i < nf; f_i++) {

    // loop over adjoint components
    for (int ci_adjoint = 0; ci_adjoint < 3; ci_adjoint++) {
      int num_chunks = adjoint_dft_chunks[ci_adjoint].size();
      if (num_chunks == 0) continue;

      // loop over each chunk
      for (int cur_chunk = 0; cur_chunk < num_chunks; cur_chunk++) {
        meep::dft_chunk *adj_chunk = adjoint_dft_chunks[ci_adjoint][cur_chunk];
        const std::complex<meep::realnum> *adj_values = adjoint_dft_values[ci_adjoint][cur_chunk];
        meep::component adjoint_c = adj_chunk->c;
        meep::grid_volume gv_adj = gv.subvolume(adj_chunk->is, adj_chunk->ie, adjoint_c);

        // loop over forward components
        for (int ci_forward = 0; ci_forward < 3; ci_forward++) {
          /* pair the forward chunk with this adjoint chunk *spatially*, not by
             position in the per-component chunk list */
          meep::dft_chunk *fwd_chunk =
              matching_dft_chunk(forward_dft_chunks[ci_forward], adj_chunk);
          if (!fwd_chunk) continue;
          const std::complex<meep::realnum> *fwd_values = fwd_chunk->dft_values();
          meep::component forward_c = fwd_chunk->c;
          meep::grid_volume gv_fwd = gv.subvolume(fwd_chunk->is, fwd_chunk->ie, forward_c);

          // loop over each point of interest
          LOOP_OVER_IVECS(gv_adj, adj_chunk->is_old, adj_chunk->ie_old, idx_adj) {
            double cyl_scale;
            IVEC_LOOP_ILOC(gv_adj, ip);
            IVEC_LOOP_LOC(gv_adj, p);
            std::complex<meep::realnum> adj = adj_values[nf * idx_adj + f_i];
            material_type md;
            geps->get_material_pt(md, p);
            /* if we have conductivities (e.g. for damping)
            then we need to make sure we correctly account
            for that here */
            if (!md->trivial) adj *= cond_cmp(adjoint_c, p, frequencies[f_i], geps);
            if (getenv("MEEP_ADJDBG")) {
              extern double _adjdbg_sum[2]; extern long _adjdbg_n[2];
              int b = adj_chunk->sn ? 1 : 0;
              _adjdbg_sum[b] += std::abs(adj); _adjdbg_n[b]++;
            }

            /**************************************/
            /*            Main Routine            */
            /**************************************/

            /********* compute -λᵀAᵤx *************/

            /* trivial case, no interpolation/restriction needed        */
            if (forward_c == adjoint_c) {
              /* index the forward chunk with its *own* index for this point;
                 idx_adj is an index into gv_adj and is only valid here because
                 the paired chunks share a fields chunk */
              std::complex<meep::realnum> fwd =
                  forward_dft_value(fwd_chunk, fwd_values, gv_fwd, gv, ip, nf, f_i);
              cyl_scale = (gv.dim == meep::Dcyl) ? 2 * p.r()
                                                 : 1; // the pi is already factored in near2far.cpp
              material_grids_addgradient_point(v_local + ng * f_i, vec_to_vector3(p),
                                               scalegrad * cyl_scale, geps, adjoint_c, forward_c,
                                               fwd, adj, frequencies[f_i], gv, du, owner_grid_id);
              /* more complicated case requires interpolation/restriction */
            }
            else if ((md->do_averaging) ||             /* account for subpixel smoothing     */
                     (!is_material_grid(md)) ||        /* account for edge effects of mg     */
                     (has_offdiag(&(md->medium_1))) || /* account for offdiagonal components */
                     (has_offdiag(&(md->medium_2)))) {
              /* we need to restrict the adjoint fields to
              the two nodes of interest (which requires a factor
              of 0.5 to scale), and interpolate the forward fields
              to the same two nodes (which requires another factor of 0.5).
              Then we perform our inner product at these nodes.
              */
              std::complex<meep::realnum> fwd_avg, fwd1, fwd2;

              // identify the first corner of the forward fields
              meep::ivec fwd_p = ip + gv.iyee_shift(forward_c) - gv.iyee_shift(adjoint_c);
              // identify the other three corners
              meep::ivec unit_a = unit_ivec(gv.dim, component_direction(adjoint_c));
              meep::ivec unit_f = unit_ivec(gv.dim, component_direction(forward_c));
              meep::ivec fwd_pa = (fwd_p + unit_a * 2);
              meep::ivec fwd_pf = (fwd_p - unit_f * 2);
              meep::ivec fwd_paf = (fwd_p + unit_a * 2 - unit_f * 2);

              // store in vector for convenience
              std::vector<meep::ivec> fwd_pl = {fwd_p, fwd_pa};
              std::vector<meep::ivec> fwd_pr = {fwd_pf, fwd_paf};

              // identify the two eps points
              std::vector<meep::ivec> ieps = {(fwd_p + fwd_pf) / 2, (fwd_pa + fwd_paf) / 2};

// operate on the each eps node
#pragma unroll
              for (int node = 0; node < 2; node++) { // two nodes
                fwd1 = forward_dft_value(fwd_chunk, fwd_values, gv_fwd, gv, fwd_pl[node], nf, f_i);
                fwd2 = forward_dft_value(fwd_chunk, fwd_values, gv_fwd, gv, fwd_pr[node], nf, f_i);
                fwd_avg = std::complex<meep::realnum>(0.5, 0) * (fwd1 + fwd2);
                meep::vec eps1 = gv[ieps[node]];
                cyl_scale = (gv.dim == meep::Dcyl) ? eps1.r() : 1;
                material_grids_addgradient_point(v_local + ng * f_i, vec_to_vector3(eps1),
                                                 scalegrad * cyl_scale, geps, adjoint_c, forward_c,
                                                 fwd_avg, std::complex<meep::realnum>(0.5, 0) * adj,
                                                 frequencies[f_i], gv, du, owner_grid_id);
              }
            }
            /********* compute λᵀbᵤ ***************/
            /* not yet implemented/needed */
            /**************************************/
          }
        }
      }
    }
  }

  /* ------------------------------------------------------------ */
  // Broadcast results
  /* ------------------------------------------------------------ */
  if (getenv("MEEP_ADJDBG")) {
    extern double _adjdbg_sum[2]; extern long _adjdbg_n[2];
    master_printf("ADJDBG sn0 n=%ld sum|adj|=%.8e | sn1 n=%ld sum|adj|=%.8e\n",
                  _adjdbg_n[0], _adjdbg_sum[0], _adjdbg_n[1], _adjdbg_sum[1]);
    _adjdbg_sum[0] = _adjdbg_sum[1] = 0; _adjdbg_n[0] = _adjdbg_n[1] = 0;
  }
  meep::sum_to_all(v_local, v, ng * nf);

  /* ------------------------------------------------------------ */
  // cleanup
  /* ------------------------------------------------------------ */
  // clear the array used for local sum to all
  delete[] v_local;

  // material_grids_addgradient consumes the adjoint DFT chunks.  Remove them
  // through their owning dft_fields objects so the owners' chunk heads are
  // cleared as well.  Deleting the copied pointers above left fields_a with
  // dangling chunk lists, causing a subsequent dft_fields::remove() (as used
  // by streamed Padé adjoints) to delete the same chunks a second time.
  for (int i = 0; i < 3; i++) {
    fields_a[i]->remove();
  }

} // material_grids_addgradient

static void find_array_min_max(int n, const double *data, double &min_val, double &max_val) {
  min_val = data[0];
  max_val = data[0];
  for (int i = 1; i < n; ++i) {
    if (data[i] < min_val) min_val = data[i];
    if (data[i] > max_val) max_val = data[i];
  }
  return;
}

void get_epsilon_grid(geometric_object_list gobj_list, material_type_list mlist,
                      material_type _default_material, bool _ensure_periodicity,
                      meep::grid_volume gv, vector3 cell_size, vector3 cell_center, int nx,
                      const double *x, int ny, const double *y, int nz, const double *z,
                      std::complex<double> *grid_vals, double frequency) {
  double min_val[3], max_val[3];
  for (int n = 0; n < 3; ++n) {
    int ndir = (n == 0) ? nx : ((n == 1) ? ny : nz);
    if (ndir < 1) meep::abort("get_epsilon_grid: ndir < 1.");
    const double *adir = (n == 0) ? x : ((n == 1) ? y : z);
    find_array_min_max(ndir, adir, min_val[n], max_val[n]);
  }
  const meep::volume vol(meep::vec(min_val[0], min_val[1], min_val[2]),
                         meep::vec(max_val[0], max_val[1], max_val[2]));
  init_libctl(_default_material, _ensure_periodicity, &gv, cell_size, cell_center, &gobj_list);
  dim = gv.dim;
  geom_epsilon geps(gobj_list, mlist, vol);
  for (int i = 0; i < nx; ++i)
    for (int j = 0; j < ny; ++j)
      for (int k = 0; k < nz; ++k) {
        /* obtain the trace of the ε tensor (dispersive or non) for each
           grid point in row-major order (the order used by NumPy) */
        if (frequency == 0)
          grid_vals[k + nz * (j + ny * i)] =
              geps.chi1p1(meep::E_stuff, meep::vec(x[i], y[j], z[k]));
        else {
          std::complex<double> tensor[9];
          get_chi1_tensor_disp(tensor, meep::vec(x[i], y[j], z[k]), frequency, &geps);
          grid_vals[k + nz * (j + ny * i)] = (tensor[0] + tensor[4] + tensor[8]) / 3.0;
        }
      }
}

} // namespace meep_geom
