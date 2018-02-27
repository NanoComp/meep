#include <algorithm>
#include <complex>
#include <cstdlib>
#include <cstddef>
#include <iostream>

#include "config.h"
#include "pympb.hpp"
#include "mpb/scalar.h"
#include "meep/mympi.hpp"

// xyz_loop.h
#ifndef HAVE_MPI
  #define LOOP_XYZ(md) { \
    int n1 = md->nx, n2 = md->ny, n3 = md->nz, i1, i2, i3; \
    for (i1 = 0; i1 < n1; ++i1) \
      for (i2 = 0; i2 < n2; ++i2) \
        for (i3 = 0; i3 < n3; ++i3) { \
          int xyz_index = ((i1 * n2 + i2) * n3 + i3);
#else /* HAVE_MPI */
  /* first two dimensions are transposed in MPI output: */
  #define LOOP_XYZ(md) { \
    int n1 = md->nx, n2 = md->ny, n3 = md->nz, i1, i2_, i3; \
    int local_n2 = md->local_ny, local_y_start = md->local_y_start; \
    for (i2_ = 0; i2_ < local_n2; ++i2_) \
      for (i1 = 0; i1 < n1; ++i1) \
        for (i3 = 0; i3 < n3; ++i3) { \
          int i2 = i2_ + local_y_start; \
          int xyz_index = ((i2_ * n1 + i1) * n3 + i3);
#  endif /* HAVE_MPI */

// TODO: Support MPI
#define mpi_allreduce(sb, rb, n, ctype, t, op, comm) { \
     CHECK((sb) != (rb), "MPI_Allreduce doesn't work for sendbuf == recvbuf");\
     memcpy((rb), (sb), (n) * sizeof(ctype)); \
}

/* "in-place" Allreduce wrapper for reducing a single value */
#define mpi_allreduce_1(b, ctype, t, op, comm) { \
     ctype bbbb = *(b); \
     mpi_allreduce(&bbbb, (b), 1, ctype, t, op, comm); \
}

#ifdef CHECK_DISABLE
  #define CHECK(cond, s) // Do nothing
#else
  #define CHECK(cond, s) if (!(cond)){meep::abort(s "\n");}
#endif

namespace py_mpb {

// TODO: Placeholder
int mpb_comm;

const double inf = 1.0e20;

// This is the function passed to `set_maxwell_dielectric`
static void dielectric_function(symmetric_matrix *eps, symmetric_matrix *eps_inv,
                                const mpb_real r[3], void *epsilon_data) {

  mode_solver *ms = static_cast<mode_solver *>(epsilon_data);
  meep_geom::material_type mat;
  vector3 p;

  // p needs to be in the lattice *unit* vector basis, while r is in the lattice
  // vector basis.  Also, shift origin to the center of the grid.
  p.x = (r[0] - 0.5) * geometry_lattice.size.x;
  p.y = (r[1] - 0.5) * geometry_lattice.size.y;
  p.z = (r[2] - 0.5) * geometry_lattice.size.z;

  // p = shift_to_unit_cell(p);

  ms->get_material_pt(mat, p);
  ms->material_epsmu(mat, eps, eps_inv);
}

static int mean_epsilon_func(symmetric_matrix* meps, symmetric_matrix *meps_inv,
                             mpb_real n[3], mpb_real d1, mpb_real d2, mpb_real d3,
                             mpb_real tol, const mpb_real r[3], void *edata) {

  mode_solver *ms = static_cast<mode_solver *>(edata);
  meep_geom::material_type mat;
  vector3 p;

  // p needs to be in the lattice *unit* vector basis, while r is in the lattice
  // vector basis.  Also, shift origin to the center of the grid.
  p.x = (r[0] - 0.5) * geometry_lattice.size.x;
  p.y = (r[1] - 0.5) * geometry_lattice.size.y;
  p.z = (r[2] - 0.5) * geometry_lattice.size.z;

  mpb_real adjusted_tol = tol > 0.01 ? 0.01 : tol;
  mpb_real d[3] = {d1, d2, d3};
  ms->eff_chi1inv_matrix(meps_inv, d, adjusted_tol, 100 / adjusted_tol, true);

  return 1;
}

/****** utils ******/

/* a couple of utilities to convert libctl data types to the data
   types of the eigensolver & maxwell routines: */

void vector3_to_arr(mpb_real arr[3], vector3 v)
{
     arr[0] = v.x;
     arr[1] = v.y;
     arr[2] = v.z;
}

void matrix3x3_to_arr(mpb_real arr[3][3], matrix3x3 m)
{
     vector3_to_arr(arr[0], m.c0);
     vector3_to_arr(arr[1], m.c1);
     vector3_to_arr(arr[2], m.c2);
}

// Return a string describing the current parity, used for frequency and filename
// prefixes
const char *parity_string(maxwell_data *d) {
  static char s[128];
  strcpy(s, "");
  if (d->parity & EVEN_Z_PARITY) {
    strcat(s, (d->nz == 1) ? "te" : "zeven");
  } else if (d->parity & ODD_Z_PARITY) {
    strcat(s, (d->nz == 1) ? "tm" : "zodd");
  }
  if (d->parity & EVEN_Y_PARITY) {
    strcat(s, "yeven");
  } else if (d->parity & ODD_Y_PARITY) {
    strcat(s, "yodd");
  }
  return s;
}

/* Extract the mean epsilon from the effective inverse dielectric tensor,
   which contains two eigenvalues that correspond to the mean epsilon,
   and one which corresponds to the harmonic mean. */
mpb_real mean_medium_from_matrix(const symmetric_matrix *eps_inv) {
  mpb_real eps_eigs[3];
  maxwell_sym_matrix_eigs(eps_eigs, eps_inv);
  /* the harmonic mean should be the largest eigenvalue (smallest
     epsilon), so we'll ignore it and average the other two: */
  return 2.0 / (eps_eigs[0] + eps_eigs[1]);
}

/* When we are solving for a few bands at a time, we solve for the
   upper bands by "deflation"--by continually orthogonalizing them
   against the already-computed lower bands.  (This constraint
   commutes with the eigen-operator, of course, so all is well.) */

typedef struct {
  evectmatrix Y;  /* the vectors to orthogonalize against; Y must
                     itself be normalized (Yt B Y = 1) */
  evectmatrix BY;  /* B * Y */
  int p;  /* the number of columns of Y to orthogonalize against */
  scalar *S;  /* a matrix for storing the dot products; should have
                 at least p * X.p elements (see below for X) */
  scalar *S2; /* a scratch matrix the same size as S */
} deflation_data;

extern "C" {
void blasglue_gemm(char transa, char transb, int m, int n, int k, mpb_real a, scalar *A,
                   int fdA, scalar *B, int fdB, mpb_real b, scalar *C, int fdC);
}

static void deflation_constraint(evectmatrix X, void *data) {
  deflation_data *d = (deflation_data *) data;

  CHECK(X.n == d->BY.n && d->BY.p >= d->p && d->Y.p >= d->p, "invalid dimensions");

  /* compute (1 - Y (BY)t) X = (1 - Y Yt B) X
      = projection of X so that Yt B X = 0 */

  /* (Sigh...call the BLAS functions directly since we are not
     using all the columns of BY...evectmatrix is not set up for
     this case.) */

  /* compute S = Xt BY (i.e. all the dot products): */
  blasglue_gemm('C', 'N', X.p, d->p, X.n, 1.0, X.data, X.p, d->BY.data, d->BY.p,
                0.0, d->S2, d->p);
// TODO
// #if HAVE_MPI
//   MPI_Allreduce(d->S2, d->S, d->p * X.p * SCALAR_NUMVALS, SCALAR_MPI_TYPE,
//                 MPI_SUM, mpb_comm);
// #else
  memcpy(d->S, d->S2, sizeof(mpb_real) * d->p * X.p * SCALAR_NUMVALS);
// #endif

  /* compute X = X - Y*St = (1 - BY Yt B) X */
  blasglue_gemm('N', 'C', X.n, X.p, d->p, -1.0, d->Y.data, d->Y.p, d->S, d->p,
                1.0, X.data, X.p);
}

/******* mode_solver *******/

mode_solver::mode_solver(int num_bands,
                         int parity,
                         double resolution[3],
                         lattice lat,
                         double tolerance,
                         int mesh_size,
                         meep_geom::material_data *_default_material,
                         geometric_object_list geom,
                         bool reset_fields,
                         bool deterministic,
                         double target_freq,
                         int dims,
                         bool verbose,
                         bool periodicity,
                         double flops,
                         bool negative_epsilon_ok):
  num_bands(num_bands),
  parity(parity),
  target_freq(target_freq),
  tolerance(tolerance),
  mesh_size(mesh_size),
  negative_epsilon_ok(negative_epsilon_ok),
  eigensolver_nwork(3),
  eigensolver_block_size(-11),
  last_parity(-2),
  iterations(0),
  eigensolver_flops(flops),
  vol(0),
  mdata(NULL),
  mtdata(NULL),
  curfield_band(0),
  freqs(num_bands),
  verbose(verbose),
  deterministic(deterministic),
  kpoint_index(0),
  curfield(NULL),
  curfield_type('-') {

  this->lat = lat;

  geometry_lattice = lat;
  dimensions = dims;
  ensure_periodicity = periodicity;

  for (int i = 0; i < 3; ++i) {
    this->resolution[i] = resolution[i];
    for (int j = 0; j < 3; ++j) {
      R[i][j] = 0.0;
      G[i][j] = 0.0;
    }
  }

  default_material = _default_material;
  geometry = geom;

  // `init` is called in the constructor to avoid the need to copy the
  // geometric objects and default material. They can then be safely freed by
  // typemaps once the mode_solver constructor call returns to python.
  init(parity, reset_fields);
}

mode_solver::~mode_solver() {
  destroy_maxwell_data(mdata);
  destroy_maxwell_target_data(mtdata);
  destroy_geom_box_tree(geometry_tree);
  destroy_evectmatrix(H);

  for (int i = 0; i < nwork_alloc; ++i) {
    destroy_evectmatrix(W[i]);
  }

  if (Hblock.data != H.data) {
    destroy_evectmatrix(Hblock);
  }

  if (muinvH.data != H.data) {
    destroy_evectmatrix(muinvH);
  }
}

bool mode_solver::get_front_object(mpb_real v[3], vector3 &pcenter,
                                   const geometric_object **o_front, vector3 &shiftby_front,
                                   meep_geom::material_type &mat_front,
                                   meep_geom::material_type &mat_behind) {
  vector3 p;
  const geometric_object *o1 = 0, *o2 = 0;
  vector3 shiftby1 = {0,0,0}, shiftby2 = {0,0,0};
  geom_box pixel;
  meep_geom::material_type mat1, mat2;
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

  // TODO: this should be r, not v. Same?
  p.x = (v[0] - 0.5) * geometry_lattice.size.x;
  p.y = (v[1] - 0.5) * geometry_lattice.size.y;
  p.z = (v[2] - 0.5) * geometry_lattice.size.z;

  pixel.low.x = p.x - v[0];
  pixel.high.x = p.x + v[0];
  pixel.low.y = p.y - v[1];
  pixel.high.y = p.y + v[1];
  pixel.low.z = p.z - v[2];
  pixel.high.z = p.z + v[2];

  pcenter = p;
  double d1, d2, d3;
  d1 = (pixel.high.x - pixel.low.x) * 0.5;
  d2 = (pixel.high.y - pixel.low.y) * 0.5;
  d3 = (pixel.high.z - pixel.low.z) * 0.5;
  for (int i = 0; i < num_neighbors[dimensions - 1]; ++i) {
    const geometric_object *o;
    meep_geom::material_type mat;
    vector3 q, shiftby;
    int id;
    q.x = p.x + neighbors[dimensions - 1][i][0] * d1;
    q.y = p.y + neighbors[dimensions - 1][i][1] * d2;
    q.z = p.z + neighbors[dimensions - 1][i][2] * d3;
    o = object_of_point_in_tree(q, geometry_tree, &shiftby, &id);
    if ((id == id1 && vector3_equal(shiftby, shiftby1)) ||
        (id == id2 && vector3_equal(shiftby, shiftby2))) {
      continue;
    }

    mat = (meep_geom::material_type) default_material;

    if (o) {
      meep_geom::material_data *md = (meep_geom::material_data *)o->material;
      if (md->which_subclass != meep_geom::material_data::MATERIAL_FILE) {
        mat = md;
      }
    }

    if (id1 == -1) {
      o1 = o;
      shiftby1 = shiftby;
      id1 = id;
      mat1 = mat;
    }
    else if (id2 == -1 || ( (id >= id1 && id >= id2) &&
          (id1 == id2 || material_type_equal(mat1,mat2)))) {
      o2 = o;
      shiftby2 = shiftby;
      id2 = id;
      mat2 = mat;
    }
    else if (!(id1 < id2 && (id1 == id || material_type_equal(mat1,mat))) &&
             !(id2 < id1 && (id2 == id || material_type_equal(mat2,mat)))) {
      return false;
    }
  }

  CHECK(id1 > -1, "bug in object_of_point_in_tree?");
  if (id2 == -1) { /* only one nearby object/material */
    id2 = id1;
    o2 = o1;
    mat2 = mat1;
    shiftby2 = shiftby1;
  }

  if ((o1 && is_variable(o1->material)) ||
      (o2 && is_variable(o2->material)) ||
      ((is_variable(default_material) || is_file(default_material)) &&
        (!o1 || is_file(o1->material) || !o2 || is_file(o2->material)))) {
    return false;
  }

  if (id1 >= id2) {
    *o_front = o1;
    shiftby_front = shiftby1;
    mat_front = mat1;
    if (id1 == id2) {
      mat_behind = mat1;
    }
    else {
      mat_behind = mat2;
    }
  }

  if (id2 > id1) {
    *o_front = o2;
    shiftby_front = shiftby2;
    mat_front = mat2;
    mat_behind = mat1;
  }
  return true;
}

void mode_solver::eff_chi1inv_matrix(symmetric_matrix *chi1inv_matrix, mpb_real d[3],
                                     double tol, int maxeval, bool eps) {
  const geometric_object *o;
  meep_geom::material_type mat, mat_behind;
  symmetric_matrix meps;
  vector3 p, shiftby, normal;
  vector3 center = {d[0] / 2, d[1] / 2, d[2] / 2};

  if (maxeval == 0 || !get_front_object(d, geometry_tree, p, &o, shiftby, mat, mat_behind)) {
  noavg:
    get_material_pt(mat, center);
  trivial:
    material_epsmu(mat, &meps, chi1inv_matrix, eps);
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
  if (material_type_equal(mat, mat_behind)) {
    goto trivial;
  }

  // it doesn't make sense to average metals (electric or magnetic)
  if (is_metal(&mat, eps) || is_metal(&mat_behind, eps)) {
    goto noavg;
  }

  normal = unit_vector3(normal_to_fixed_object(vector3_minus(p, shiftby), *o));
  if (normal.x == 0 && normal.y == 0 && normal.z == 0)
    goto noavg; // couldn't get normal vector for this point, punt
  geom_box pixel = gv2box(v);
  pixel.low = vector3_minus(pixel.low, shiftby);
  pixel.high = vector3_minus(pixel.high, shiftby);

  double fill = box_overlap_with_object(pixel, *o, tol, maxeval);

  material_epsmu(mat, &meps, chi1inv_matrix, eps);
  symmetric_matrix eps2, epsinv2;
  symmetric_matrix eps1, delta;
  double Rot[3][3];
  material_epsmu(mat_behind, &eps2, &epsinv2, eps);
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

void mode_solver::material_epsmu(meep_geom::material_type material, symmetric_matrix *epsmu,
                                 symmetric_matrix *epsmu_inv, bool eps) {

  meep_geom::material_data *md = material;

  if (eps) {
    switch (md->which_subclass) {
      case meep_geom::material_data::MEDIUM:
      case meep_geom::material_data::MATERIAL_FILE:
      case meep_geom::material_data::MATERIAL_USER:
        epsmu->m00 = md->medium.epsilon_diag.x;
        epsmu->m11 = md->medium.epsilon_diag.y;
        epsmu->m22 = md->medium.epsilon_diag.z;
        epsmu->m01 = md->medium.epsilon_offdiag.x;
        epsmu->m02 = md->medium.epsilon_offdiag.y;
        epsmu->m12 = md->medium.epsilon_offdiag.z;
        maxwell_sym_matrix_invert(epsmu_inv, epsmu);
        break;
      case meep_geom::material_data::PERFECT_METAL:
        epsmu->m00 = -inf;
        epsmu->m11 = -inf;
        epsmu->m22 = -inf;
        epsmu->m01 = 0.0;
        epsmu->m02 = 0.0;
        epsmu->m12 = 0.0;
        epsmu_inv->m00 = -0.0;
        epsmu_inv->m11 = -0.0;
        epsmu_inv->m22 = -0.0;
        epsmu_inv->m01 = 0.0;
        epsmu_inv->m02 = 0.0;
        epsmu_inv->m12 = 0.0;
        break;
      default:
        meep::abort("Unknown material type");
    }
  }
  else {
    switch (md->which_subclass) {
      case meep_geom::material_data::MEDIUM:
      case meep_geom::material_data::MATERIAL_FILE:
      case meep_geom::material_data::MATERIAL_USER:
        epsmu->m00 = md->medium.mu_diag.x;
        epsmu->m11 = md->medium.mu_diag.y;
        epsmu->m22 = md->medium.mu_diag.z;
        epsmu->m01 = md->medium.mu_offdiag.x;
        epsmu->m02 = md->medium.mu_offdiag.y;
        epsmu->m12 = md->medium.mu_offdiag.z;
        maxwell_sym_matrix_invert(epsmu_inv, epsmu);
        break;
      case meep_geom::material_data::PERFECT_METAL:
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
}

void mode_solver::get_material_pt(meep_geom::material_type &material, vector3 p) {
  boolean inobject;
  material = (meep_geom::material_type)material_of_unshifted_point_in_tree_inobject(p, geometry_tree, &inobject);
  meep_geom::material_data *md = material;

  switch(md->which_subclass) {
    // material read from file: interpolate to get properties at r
    case meep_geom::material_data::MATERIAL_FILE:
      if (md->epsilon_data) {
        meep_geom::epsilon_file_material(md, p);
      }
      else {
        material = (meep_geom::material_type) default_material;
      }
      return;

    // material specified by user-supplied function: call user
    // function to get properties at r.
    // Note that we initialize the medium to vacuum, so that
    // the user's function only needs to fill in whatever is
    // different from vacuum.
    case meep_geom::material_data::MATERIAL_USER:
      md->medium = meep_geom::medium_struct();
      md->user_func(p, md->user_data, &(md->medium));
      // TODO: update this to allow user's function to set
      //       position-dependent susceptibilities. For now
      //       it's an error if the user's function creates
      //       any.
      if ((md->medium.E_susceptibilities.num_items>0) ||
          (md->medium.H_susceptibilities.num_items>0)) {
        meep::abort("susceptibilities in user-defined-materials not yet supported");
      }
      return;

    // position-independent material or metal: there is nothing to do
    case meep_geom::material_data::MEDIUM:
    case meep_geom::material_data::PERFECT_METAL:
      return;
    default:
      meep::abort("unknown material type");
   }
}

bool mode_solver::using_mu() {
  return mdata && mdata->mu_inv != NULL;
}

void mode_solver::init(int p, bool reset_fields) {
  int have_old_fields = 0;

  n[0] = std::max(resolution[0] * std::ceil(geometry_lattice.size.x), 1.0);
  n[1] = std::max(resolution[1] * std::ceil(geometry_lattice.size.y), 1.0);
  n[2] = std::max(resolution[2] * std::ceil(geometry_lattice.size.z), 1.0);

  if (target_freq != 0.0) {
    meep::master_printf("Target frequency is %g\n", target_freq);
  }

  int true_rank = n[2] > 1 ? 3 : (n[1] > 1 ? 2 : 1);
  if (true_rank < dimensions) {
    dimensions = true_rank;
  } else if (true_rank > dimensions) {
    meep::master_printf("WARNING: rank of grid is > dimensions.\n"
                        "         setting extra grid dims. to 1.\n");
    // force extra dims to be 1
    if (dimensions <= 2) {
      n[2] = 1;
    }
    if (dimensions <= 1) {
      n[1] = 1;
    }
  }

  meep::master_printf("Working in %d dimensions.\n", dimensions);
  meep::master_printf("Grid size is %d x %d x %d.\n", n[0], n[1], n[2]);

  int block_size;

  if (eigensolver_block_size != 0 && eigensolver_block_size < num_bands) {
    block_size = eigensolver_block_size;
    if (block_size < 0) {
      // Guess a block_size near -block_size, chosen so that all blocks are nearly equal in size
      block_size = (num_bands - block_size - 1) / (-block_size);
      block_size = (num_bands + block_size - 1) / block_size;
    }
    meep::master_printf("Solving for %d bands at a time.\n", block_size);
  } else {
    block_size = num_bands;
  }

  if (deterministic) {
    // seed should be the same for each run, although
    // it should be different for each process.
    // TODO: MPI
    // int rank = meep::my_rank();
    srand(314159); // * (rank + 1));
  }

  meep::master_printf("Creating Maxwell data...\n");
  mdata = create_maxwell_data(n[0], n[1], n[2], &local_N, &N_start, &alloc_N, block_size, NUM_FFT_BANDS);

  if (target_freq != 0.0) {
    mtdata = create_maxwell_target_data(mdata, target_freq);
  }

  init_epsilon();

  if (check_maxwell_dielectric(mdata, 0)) {
    meep::abort("invalid dielectric function for MPB");
  }

  if (!have_old_fields) {
    meep::master_printf("Allocating fields...\n");

    int N = n[0] * n[1] * n[2];
    int c = 2;

    H = create_evectmatrix(N, c, num_bands, local_N, N_start, alloc_N);
    nwork_alloc = eigensolver_nwork + (mdata->mu_inv != NULL);


    for (int i = 0; i < nwork_alloc; ++i) {
      W[i] = create_evectmatrix(N, c, block_size, local_N, N_start, alloc_N);
    }

    if (block_size < num_bands) {
      Hblock = create_evectmatrix(N, c, block_size, local_N, N_start, alloc_N);
    } else {
      Hblock = H;
    }

    if (using_mu() && block_size < num_bands) {
      muinvH = create_evectmatrix(N, c, num_bands, local_N, N_start, alloc_N);
    } else {
      muinvH = H;
    }
  }

  set_parity(p);

  if (!have_old_fields || reset_fields) {
    randomize_fields();
  }

  evectmatrix_flops = eigensolver_flops;
}

void mode_solver::init_epsilon() {
  int no_size_x = geometry_lattice.size.x == 0 ? 1 : geometry_lattice.size.x;
  int no_size_y = geometry_lattice.size.y == 0 ? 1 : geometry_lattice.size.y;
  int no_size_z = geometry_lattice.size.z == 0 ? 1 : geometry_lattice.size.z;

  meep::master_printf("Mesh size is %d.\n", mesh_size);

  Rm.c0 = vector3_scale(no_size_x, geometry_lattice.basis.c0);
  Rm.c1 = vector3_scale(no_size_y, geometry_lattice.basis.c1);
  Rm.c2 = vector3_scale(no_size_z, geometry_lattice.basis.c2);

  meep::master_printf("Lattice vectors:\n");
  meep::master_printf("     (%g, %g, %g)\n", Rm.c0.x, Rm.c0.y, Rm.c0.z);
  meep::master_printf("     (%g, %g, %g)\n", Rm.c1.x, Rm.c1.y, Rm.c1.z);
  meep::master_printf("     (%g, %g, %g)\n", Rm.c2.x, Rm.c2.y, Rm.c2.z);

  vol = fabs(matrix3x3_determinant(Rm));
  meep::master_printf("Cell volume = %g\n", vol);

  Gm = matrix3x3_inverse(matrix3x3_transpose(Rm));
  meep::master_printf("Reciprocal lattice vectors (/ 2 pi):\n");
  meep::master_printf("     (%g, %g, %g)\n", Gm.c0.x, Gm.c0.y, Gm.c0.z);
  meep::master_printf("     (%g, %g, %g)\n", Gm.c1.x, Gm.c1.y, Gm.c1.z);
  meep::master_printf("     (%g, %g, %g)\n", Gm.c2.x, Gm.c2.y, Gm.c2.z);

  matrix3x3_to_arr(R, Rm);
  matrix3x3_to_arr(G, Gm);

  geom_fix_objects0(geometry);

  meep::master_printf("Geometric objects:\n");
  if (meep::am_master()) {
    for (int i = 0; i < geometry.num_items; ++i) {
      display_geometric_object_info(5, geometry.items[i]);

      // meep_geom::medium_struct *mm;
      // if (meep_geom::is_medium(geometry.items[i].material, &mm)) {
      //   printf("%*sdielectric constant epsilon diagonal = (%g,%g,%g)\n", 5 + 5, "",
      //          mm->epsilon_diag.x, mm->epsilon_diag.y, mm->epsilon_diag.z);
      // }
    }
  }

  {
    // Replace 0 with 1e-20 for no size
    vector3 tmp_size;
    tmp_size.x = geometry_lattice.size.x == 0 ? 1e-20 : geometry_lattice.size.x;
    tmp_size.y = geometry_lattice.size.y == 0 ? 1e-20 : geometry_lattice.size.y;
    tmp_size.z = geometry_lattice.size.z == 0 ? 1e-20 : geometry_lattice.size.z;

    geom_box b0;
    b0.low = vector3_plus(geometry_center, vector3_scale(-0.5, tmp_size));
    b0.high = vector3_plus(geometry_center, vector3_scale(0.5, tmp_size));
    /* pad tree boundaries to allow for sub-pixel averaging */
    b0.low.x -= tmp_size.x / mdata->nx;
    b0.low.y -= tmp_size.y / mdata->ny;
    b0.low.z -= tmp_size.z / mdata->nz;
    b0.high.x += tmp_size.x / mdata->nx;
    b0.high.y += tmp_size.y / mdata->ny;
    b0.high.z += tmp_size.z / mdata->nz;
    geometry_tree = create_geom_box_tree0(geometry, b0);
  }

  if (verbose && meep::am_master()) {
    printf("Geometry object bounding box tree:\n");
    display_geom_box_tree(5, geometry_tree);
  }

  int tree_depth;
  int tree_nobjects;
  geom_box_tree_stats(geometry_tree, &tree_depth, &tree_nobjects);
  meep::master_printf("Geometric object tree has depth %d and %d object nodes"
                      " (vs. %d actual objects)\n", tree_depth, tree_nobjects, geometry.num_items);

  // restricted_tree = geometry_tree;

  reset_epsilon();
}

void mode_solver::reset_epsilon() {
  int mesh[3] = {
    mesh_size,
    (dimensions > 1) ? mesh_size : 1,
    (dimensions > 2) ? mesh_size : 1,
  };

  // TODO: Support epsilon_input_file
  // get_epsilon_file_func(epsilon_input_file, &d.epsilon_file_func, &d.epsilon_file_func_data);
  // get_epsilon_file_func(mu_input_file, &d.mu_file_func, &d.mu_file_func_data);
  meep::master_printf("Initializing epsilon function...\n");
  set_maxwell_dielectric(mdata, mesh, R, G, dielectric_function, mean_epsilon_func, static_cast<void *>(this));

  // TODO
  // if (has_mu(&d)) {
  //   mpi_one_printf("Initializing mu function...\n");
  //   set_maxwell_mu(mdata, mesh, R, G, mu_func, mean_mu_func, &d);
  // }

  // destroy_epsilon_file_func_data(d.epsilon_file_func_data);
  // destroy_epsilon_file_func_data(d.mu_file_func_data);
}

void mode_solver::set_parity(integer p) {

  if (!mdata) {
    meep::master_fprintf(stderr, "init must be called before set-parity!\n");
    return;
  }

  if (p == -1) {
    p = last_parity < 0 ? NO_PARITY : last_parity;
  }

  set_maxwell_data_parity(mdata, p);
  if (mdata->parity != p) {
    meep::master_fprintf(stderr, "k vector incompatible with parity\n");
    exit(EXIT_FAILURE);
  }
  meep::master_printf("Solving for band polarization: %s.\n", parity_string(mdata));

  last_parity = p;
  set_kpoint_index(0);  /* reset index */
}

int mode_solver::get_kpoint_index() {
  return kpoint_index;
}

void mode_solver::set_kpoint_index(int i) {
  kpoint_index = i;
}

void mode_solver::randomize_fields() {

  if (!mdata) {
    return;
  }
  meep::master_printf("Initializing fields to random numbers...\n");

  for (int i = 0; i < H.n * H.p; ++i) {
    ASSIGN_SCALAR(H.data[i], rand() * 1.0 / RAND_MAX, rand() * 1.0 / RAND_MAX);
  }
}

void mode_solver::solve_kpoint(vector3 kvector) {

  // if we get too close to singular k==0 point, just set k=0 exploit our
  // special handling of this k
  if (vector3_norm(kvector) < 1e-10) {
    kvector.x = kvector.y = kvector.z = 0;
  }

  meep::master_printf("solve_kpoint (%g,%g,%g):\n", kvector.x, kvector.y, kvector.z);

  curfield_reset();

  if (num_bands == 0) {
    meep::master_printf("  num-bands is zero, not solving for any bands\n");
    return;
  }

  if (!mdata) {
    meep::master_fprintf(stderr, "init must be called before solve_kpoint!\n");
    return;
  }

  // If this is the first k point, print out a header line for the frequency
  // grep data.
  if (!kpoint_index && meep::am_master()) {
    printf("%sfreqs:, k index, k1, k2, k3, kmag/2pi", parity_string(mdata));

    for (int i = 0; i < num_bands; ++i) {
      printf(", %s%sband %d", parity_string(mdata), mdata->parity == NO_PARITY ? "" : " ", i + 1);
    }
    printf("\n");
  }

  cur_kvector = kvector;
  mpb_real k[3];
  vector3_to_arr(k, kvector);

  update_maxwell_data_k(mdata, k, G[0], G[1], G[2]);

  std::vector<mpb_real> eigvals(num_bands);

  // TODO: Get flags from python
  int flags = EIGS_DEFAULT_FLAGS;
  if (verbose) {
   flags |= EIGS_VERBOSE;
  }

  // Constant (zero frequency) bands at k=0 are handled specially, so remove
  // them from the solutions for the eigensolver.
  int ib0;
  if (mdata->zero_k && !mtdata) {
    ib0 = maxwell_zero_k_num_const_bands(H, mdata);
    for (int in = 0; in < H.n; ++in) {
      for (int ip = 0; ip < H.p - ib0; ++ip) {
        H.data[in * H.p + ip] = H.data[in * H.p + ip + ib0];
      }
    }
    evectmatrix_resize(&H, H.p - ib0, 1);
  }
  else {
    ib0 = 0; /* solve for all bands */
  }

  // Set up deflation data.
  deflation_data deflation;
  if (muinvH.data != Hblock.data) {
    deflation.Y = H;
    deflation.BY = muinvH.data != H.data ? muinvH : H;
    deflation.p = 0;
    deflation.S = (scalar *)malloc(sizeof(scalar) * H.p * Hblock.p);
    deflation.S2 = (scalar *)malloc(sizeof(scalar) * H.p * Hblock.p);
  }

  int total_iters = 0;

  for (int ib = ib0; ib < num_bands; ib += Hblock.alloc_p) {
    evectconstraint_chain *constraints;
    int num_iters;

    // Don't solve for too many bands if the block size doesn't divide the number
    // of bands.
    if (ib + mdata->num_bands > num_bands) {
      maxwell_set_num_bands(mdata, num_bands - ib);

      for (int i = 0; i < nwork_alloc; ++i) {
        evectmatrix_resize(&W[i], num_bands - ib, 0);
      }
      evectmatrix_resize(&Hblock, num_bands - ib, 0);
    }

    meep::master_printf("Solving for bands %d to %d...\n", ib + 1, ib + Hblock.p);

    constraints = NULL;
    constraints = evect_add_constraint(constraints, maxwell_parity_constraint, (void *) mdata);

    if (mdata->zero_k) {
      constraints = evect_add_constraint(constraints, maxwell_zero_k_constraint, (void *) mdata);
    }

    if (Hblock.data != H.data) {  /* initialize fields of block from H */
      for (int in = 0; in < Hblock.n; ++in) {
        for (int ip = 0; ip < Hblock.p; ++ip) {
          Hblock.data[in * Hblock.p + ip] = H.data[in * H.p + ip + (ib-ib0)];
        }
      }

      deflation.p = ib-ib0;
      if (deflation.p > 0) {
        if (deflation.BY.data != H.data) {
          evectmatrix_resize(&deflation.BY, deflation.p, 0);
          maxwell_muinv_operator(H, deflation.BY, (void *) mdata, 1, deflation.BY);
        }
        constraints = evect_add_constraint(constraints, deflation_constraint, &deflation);
      }
    }

    if (mtdata) {  /* solving for bands near a target frequency */
      // TODO
      // if (eigensolver_davidsonp) {
      // }
      CHECK(mdata->mu_inv==NULL, "targeted solver doesn't handle mu");
      // TODO: simple_preconditionerp ? maxwell_target_preconditioner : maxwell_target_preconditioner2
      eigensolver(Hblock, eigvals.data() + ib, maxwell_target_operator, (void *)mtdata,
                  NULL, NULL, maxwell_target_preconditioner2, (void *)mtdata,
                  evectconstraint_chain_func, (void *)constraints, W, nwork_alloc,
                  tolerance, &num_iters, flags);

      // now, diagonalize the real Maxwell operator in the solution subspace to
      // get the true eigenvalues and eigenvectors
      CHECK(nwork_alloc >= 2, "not enough workspace");
      eigensolver_get_eigenvals(Hblock, eigvals.data() + ib, maxwell_operator, mdata, W[0], W[1]);
    }
    else {
      // TODO
      // if (eigensolver_davidsonp) {
      // }

      eigensolver(Hblock, eigvals.data() + ib, maxwell_operator, (void *) mdata,
                  mdata->mu_inv ? maxwell_muinv_operator : NULL, (void *) mdata,
                  maxwell_preconditioner2, (void *) mdata, evectconstraint_chain_func,
                  (void *) constraints, W, nwork_alloc, tolerance, &num_iters,
                  flags);
    }

    if (Hblock.data != H.data) {  /* save solutions of current block */
      for (int in = 0; in < Hblock.n; ++in) {
        for (int ip = 0; ip < Hblock.p; ++ip) {
          H.data[in * H.p + ip + (ib-ib0)] = Hblock.data[in * Hblock.p + ip];
        }
      }
    }

    evect_destroy_constraints(constraints);

    meep::master_printf("Finished solving for bands %d to %d after %d iterations.\n",
                        ib + 1, ib + Hblock.p, num_iters);

    total_iters += num_iters * Hblock.p;
  }

  if (num_bands - ib0 > Hblock.alloc_p) {
    meep::master_printf("Finished k-point with %g mean iterations/band.\n",
                        total_iters * 1.0 / num_bands);
  }

  // Manually put in constant (zero-frequency) solutions for k=0.
  if (mdata->zero_k && !mtdata) {
    evectmatrix_resize(&H, H.alloc_p, 1);
    for (int in = 0; in < H.n; ++in) {
      for (int ip = H.p - ib0 - 1; ip >= 0; --ip) {
        H.data[in * H.p + ip + ib0] = H.data[in * H.p + ip];
      }
    }
    maxwell_zero_k_set_const_bands(H, mdata);
    for (int ib = 0; ib < ib0; ++ib) {
      eigvals[ib] = 0;
    }
  }

  /* Reset scratch matrix sizes: */
  evectmatrix_resize(&Hblock, Hblock.alloc_p, 0);
  for (int i = 0; i < nwork_alloc; ++i) {
    evectmatrix_resize(&W[i], W[i].alloc_p, 0);
  }

  maxwell_set_num_bands(mdata, Hblock.alloc_p);

  /* Destroy deflation data: */
  if (H.data != Hblock.data) {
    free(deflation.S2);
    free(deflation.S);
  }

  // TODO
  // if (num_write_output_vars > 0) {
  //   //  clean up from prev. call
  //   destroy_output_vars();
  // }

  iterations = total_iters; /* iterations output variable */

  set_kpoint_index(kpoint_index + 1);

  meep::master_printf("%sfreqs:, %d, %g, %g, %g, %g", parity_string(mdata), kpoint_index, (double)k[0],
                      (double)k[1], (double)k[2], vector3_norm(matrix3x3_vector3_mult(Gm, kvector)));

  for (int i = 0; i < num_bands; ++i) {
    freqs[i] = negative_epsilon_ok ? eigvals[i] : sqrt(eigvals[i]);
    meep::master_printf(", %g", freqs[i]);
  }
  meep::master_printf("\n");

  eigensolver_flops = evectmatrix_flops;
}

/* get the epsilon function, and compute some statistics */
void mode_solver::get_epsilon() {
  mpb_real eps_mean = 0;
  mpb_real eps_inv_mean = 0;
  mpb_real eps_high = -1e20;
  mpb_real eps_low = 1e20;

  int fill_count = 0;

  if (!mdata) {
    meep::master_fprintf(stderr, "init-params must be called before get-epsilon!\n");
    return;
  }

  curfield = (scalar_complex *) mdata->fft_data;
  mpb_real *epsilon = (mpb_real *) curfield;
  curfield_band = 0;
  curfield_type = epsilon_CURFIELD_TYPE;

  /* get epsilon.  Recall that we actually have an inverse
     dielectric tensor at each point; define an average index by
     the inverse of the average eigenvalue of the 1/eps tensor.
     i.e. 3/(trace 1/eps). */

  int N = mdata->fft_output_size;

  for (int i = 0; i < N; ++i) {
    if (mdata->eps_inv == NULL) {
      epsilon[i] = 1.0;
    }
    else {
      epsilon[i] = mean_medium_from_matrix(mdata->eps_inv + i);
    }
    if (epsilon[i] < eps_low) {
      eps_low = epsilon[i];
    }
    if (epsilon[i] > eps_high) {
      eps_high = epsilon[i];
    }
    eps_mean += epsilon[i];
    eps_inv_mean += 1/epsilon[i];
    if (epsilon[i] > 1.0001) {
      ++fill_count;
    }

  }

  mpi_allreduce_1(&eps_mean, mpb_real, SCALAR_MPI_TYPE, MPI_SUM, mpb_comm);
  mpi_allreduce_1(&eps_inv_mean, mpb_real, SCALAR_MPI_TYPE, MPI_SUM, mpb_comm);
  mpi_allreduce_1(&eps_low, mpb_real, SCALAR_MPI_TYPE, MPI_MIN, mpb_comm);
  mpi_allreduce_1(&eps_high, mpb_real, SCALAR_MPI_TYPE, MPI_MAX, mpb_comm);
  mpi_allreduce_1(&fill_count, int, MPI_INT, MPI_SUM, mpb_comm);

  N = mdata->nx * mdata->ny * mdata->nz;
  eps_mean /= N;
  eps_inv_mean = N/eps_inv_mean;

  meep::master_printf("epsilon: %g-%g, mean %g, harm. mean %g, %g%% > 1, %g%% \"fill\"\n",
                      eps_low, eps_high, eps_mean, eps_inv_mean, (100.0 * fill_count) / N,
                      eps_high == eps_low ? 100.0 : 100.0 * (eps_mean-eps_low) / (eps_high-eps_low));
}

void mode_solver::curfield_reset() {
  curfield = NULL;
  curfield_type = '-';
}

/* get the specified component of the dielectric tensor,
   or the inverse tensor if inv != 0 */
void mode_solver::get_epsilon_tensor(int c1, int c2, int imag, int inv) {
  int conj = 0, offset = 0;

  curfield_type = '-'; /* only used internally, for now */
  mpb_real *epsilon = (mpb_real *) mdata->fft_data;
  int N = mdata->fft_output_size;

  switch (c1 * 3 + c2) {
    case 0:
      offset = offsetof(symmetric_matrix, m00);
      break;
    case 1:
      offset = offsetof(symmetric_matrix, m01);
      break;
    case 2:
      offset = offsetof(symmetric_matrix, m02);
      break;
    case 3:
      offset = offsetof(symmetric_matrix, m01); /* = conj(m10) */
      conj = imag;
      break;
    case 4:
      offset = offsetof(symmetric_matrix, m11);
      break;
    case 5:
      offset = offsetof(symmetric_matrix, m12);
      break;
    case 6:
      offset = offsetof(symmetric_matrix, m02); /* = conj(m20) */
      conj = imag;
      break;
    case 7:
      offset = offsetof(symmetric_matrix, m12); /* = conj(m21) */
      conj = imag;
      break;
    case 8:
      offset = offsetof(symmetric_matrix, m22);
      break;
  }

#ifdef WITH_HERMITIAN_EPSILON
  if (c1 != c2 && imag)
    offset += offsetof(scalar_complex, im);
#endif

  for (int i = 0; i < N; ++i) {
    if (inv) {
      epsilon[i] = *((mpb_real *) (((char *) &mdata->eps_inv[i]) + offset));
    }
    else {
      symmetric_matrix eps;
      maxwell_sym_matrix_invert(&eps, &mdata->eps_inv[i]);
      epsilon[i] = *((mpb_real *) (((char *) &eps) + offset));
    }
    if (conj)
      epsilon[i] = -epsilon[i];
  }
}

void mode_solver::load_eigenvectors(char *filename) {
  meep::master_printf("Loading eigenvectors from \"%s\"...\n", filename);
  // TODO: Write in python
  // evectmatrixio_readall_raw(filename, H);
  curfield_reset();
}

std::vector<mpb_real> mode_solver::get_freqs() {
  return freqs;
}

size_t mode_solver::get_field_size() {
  return mdata ? mdata->fft_output_size * 3 : 0;
}

void mode_solver::get_efield(std::complex<mpb_real> *cdata, int size, int band) {

  get_dfield(cdata, size, band);
  get_efield_from_dfield();

  for (int i = 0; i < size; ++i) {
    cdata[i] = std::complex<mpb_real>(curfield[i].re, curfield[i].im);
  }
}

void mode_solver::get_efield_from_dfield() {

  if (!curfield || curfield_type != 'd') {
    meep::master_fprintf(stderr, "get_dfield must be called before get-efield-from-dfield!\n");
    return;
  }

  maxwell_compute_e_from_d(mdata, curfield, 1);
  curfield_type = 'e';
}

void mode_solver::get_dfield(std::complex<mpb_real> *cdata, int size, int band) {

  if (!kpoint_index) {
    meep::master_fprintf(stderr, "solve_kpoint must be called before get_dfield\n");
    return;
  }

  if (band < 1 || band > H.p) {
    meep::master_fprintf(stderr, "Must have 1 <= band index <= num_bands (%d)\n", H.p);
    return;
  }

  curfield = (scalar_complex *)mdata->fft_data;
  curfield_band = band;
  curfield_type = 'd';

  if (mdata->mu_inv == NULL) {
    maxwell_compute_d_from_H(mdata, H, curfield, band - 1, 1);
  }
  else {
    evectmatrix_resize(&W[0], 1, 0);
    maxwell_compute_H_from_B(mdata, H, W[0], curfield, band - 1, 0, 1);
    maxwell_compute_d_from_H(mdata, W[0], curfield, 0, 1);
    evectmatrix_resize(&W[0], W[0].alloc_p, 0);
  }

  // Here, we correct for the fact that compute_d_from_H actually computes just
  // (k+G) x H, whereas the actual D field is i/omega i(k+G) x H...so, there is
  // an added factor of -1/omega.

  // We also divide by the cell volume so that the integral of H*B or of D*E is
  // unity.  (From the eigensolver + FFT, they are initially normalized to sum to
  // nx*ny*nz.)

  double scale;

  if (freqs[band - 1] != 0.0) {
    scale = -1.0 / freqs[band - 1];
  }
  else
    scale = -1.0; /* arbitrary */

  scale /= sqrt(vol);

  for (int i = 0; i < size; ++i) {
    curfield[i].re *= scale;
    curfield[i].im *= scale;
    // Copy curfield into our output array for numpy
    cdata[i] = std::complex<mpb_real>(curfield[i].re, curfield[i].im);
  }
}

void mode_solver::get_hfield(std::complex<mpb_real> *cdata, int size, int band) {
  if (!kpoint_index) {
    meep::master_fprintf(stderr, "solve_kpoint must be called before get_dfield\n");
    return;
  }

  if (band < 1 || band > H.p) {
    meep::master_fprintf(stderr, "Must have 1 <= band index <= num_bands (%d)\n", H.p);
    return;
  }

  curfield = (scalar_complex *)mdata->fft_data;
  curfield_band = band;
  curfield_type = 'h';

  if (mdata->mu_inv == NULL)
    maxwell_compute_h_from_H(mdata, H, curfield, band - 1, 1);
  else {
    evectmatrix_resize(&W[0], 1, 0);
    maxwell_compute_H_from_B(mdata, H, W[0], curfield, band - 1, 0, 1);
    maxwell_compute_h_from_H(mdata, W[0], curfield, 0, 1);
    evectmatrix_resize(&W[0], W[0].alloc_p, 0);
  }

  // Divide by the cell volume so that the integral of H*B or of D*E is unity.
  // (From the eigensolver + FFT, they are initially normalized to sum to
  // nx*ny*nz.)

  double scale;
  scale = 1.0 / sqrt(vol);

  for (int i = 0; i < size; ++i) {
    curfield[i].re *= scale;
    curfield[i].im *= scale;
    // Copy curfield to our output array for numpy
    cdata[i] = std::complex<mpb_real>(curfield[i].re, curfield[i].im);
  }
}

void mode_solver::get_bfield(std::complex<mpb_real> *cdata, int size, int band) {
  if (!kpoint_index) {
    meep::master_fprintf(stderr, "solve_kpoint must be called before get_dfield\n");
    return;
  }

  if (band < 1 || band > H.p) {
    meep::master_fprintf(stderr, "Must have 1 <= band index <= num_bands (%d)\n", H.p);
    return;
  }

  curfield = (scalar_complex *)mdata->fft_data;
  curfield_band = band;
  curfield_type = 'b';
  maxwell_compute_h_from_H(mdata, H, curfield, band - 1, 1);

  // Divide by the cell volume so that the integral of H*B or of D*E is unity.
  // (From the eigensolver + FFT, they are initially normalized to sum to nx*ny*nz.) */

  double scale;
  scale = 1.0 / sqrt(vol);

  for (int i = 0; i < size; ++i) {
    curfield[i].re *= scale;
    curfield[i].im *= scale;
    // Copy curfield to our output array for numpy
    cdata[i] = std::complex<mpb_real>(curfield[i].re, curfield[i].im);
  }
}

char mode_solver::get_curfield_type() {
  return curfield_type;
}

void mode_solver::set_curfield_type(char t) {
  curfield_type = t;
}

std::string mode_solver::get_parity_string() {
  std::string s(parity_string(mdata));
  return s;
}

std::vector<int> mode_solver::get_dims() {
  std::vector<int> dims;

  if (mdata->nx > 1) {
    dims.push_back(mdata->nx);
  }
  if (mdata->ny > 1) {
    dims.push_back(mdata->ny);
  }
  if (mdata->nz > 1) {
    dims.push_back(mdata->nz);
  }

  return dims;
}

void mode_solver::get_curfield(double *data, int size) {
  mpb_real *p = (mpb_real *)curfield;

  for (int i = 0; i < size; ++i) {
    data[i] = p[i];
  }
}

void mode_solver::get_curfield_cmplx(std::complex<mpb_real> *cdata, int size) {
  scalar_complex *p = (scalar_complex *)curfield;

  for (int i = 0; i < size; ++i) {
    cdata[i] = std::complex<mpb_real>(p[i].re, p[i].im);
  }
}

void mode_solver::set_curfield(double *data, int size) {
  mpb_real *p = (mpb_real *)curfield;

  for (int i = 0; i < size; ++i) {
    p[i] = data[i];
  }
}

void mode_solver::set_curfield_cmplx(std::complex<mpb_real> *cdata, int size) {
  scalar_complex *p = (scalar_complex *)curfield;

  for (int i = 0; i < size; ++i) {
    scalar_complex s = {cdata[i].real(), cdata[i].imag()};
    p[i] = s;
  }
}

// internal function for compute_field_energy, below
double mode_solver::compute_field_energy_internal(mpb_real comp_sum[6]) {
  mpb_real comp_sum2[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  mpb_real energy_sum = 0.0;
  mpb_real *energy_density = (mpb_real *) curfield;

  int N = mdata->fft_output_size;

  for (int i = 0; i < N; ++i) {
    scalar_complex field[3];
    mpb_real comp_sqr0, comp_sqr1, comp_sqr2, comp_sqr3, comp_sqr4, comp_sqr5;

    /* energy is either |curfield|^2 / mu or |curfield|^2 / epsilon,
       depending upon whether it is B or D. */
    if (curfield_type == 'd') {
      assign_symmatrix_vector(field, mdata->eps_inv[i], curfield+3*i);
    }
    else if (curfield_type == 'b' && mdata->mu_inv != NULL) {
      assign_symmatrix_vector(field, mdata->mu_inv[i], curfield+3*i);
    }
    else {
      field[0] = curfield[3*i];
      field[1] = curfield[3*i+1];
      field[2] = curfield[3*i+2];
    }

    comp_sum2[0] += comp_sqr0 = field[0].re *   curfield[3*i].re;
    comp_sum2[1] += comp_sqr1 = field[0].im *   curfield[3*i].im;
    comp_sum2[2] += comp_sqr2 = field[1].re * curfield[3*i+1].re;
    comp_sum2[3] += comp_sqr3 = field[1].im * curfield[3*i+1].im;
    comp_sum2[4] += comp_sqr4 = field[2].re * curfield[3*i+2].re;
    comp_sum2[5] += comp_sqr5 = field[2].im * curfield[3*i+2].im;

    /* Note: here, we write to energy_density[i]; this is
       safe, even though energy_density is aliased to curfield,
       since energy_density[i] is guaranteed to come at or before
       curfield[i] (which we are now done with). */
    energy_sum += energy_density[i] = comp_sqr0+comp_sqr1+comp_sqr2+comp_sqr3+comp_sqr4+comp_sqr5;
  }

  mpi_allreduce_1(&energy_sum, mpb_real, SCALAR_MPI_TYPE, MPI_SUM, mpb_comm);
  mpi_allreduce(comp_sum2, comp_sum, 6, mpb_real, SCALAR_MPI_TYPE, MPI_SUM, mpb_comm);

  // remember that we now have energy density; denoted by capital D/H
  curfield_type = toupper(curfield_type);

  return energy_sum;
}

/* Replace curfield (either d or h) with the scalar energy density function,
   normalized to one.  While we're at it, compute some statistics about
   the relative strength of different field components.  Also return
   the integral of the energy density, which should be unity. */
std::vector<mpb_real> mode_solver::compute_field_energy() {
  std::vector<mpb_real> retval;

  if (!curfield || !strchr("dhb", curfield_type)) {
    meep::master_fprintf(stderr, "The D or H field must be loaded first.\n");
    return retval;
  }
  else if (curfield_type == 'h' && mdata->mu_inv != NULL) {
    meep::master_fprintf(stderr, "B, not H, must be loaded if we have mu.\n");
    return retval;
  }

  mpb_real comp_sum[6];
  mpb_real energy_sum = compute_field_energy_internal(comp_sum);

  meep::master_printf("%c-energy-components:, %d, %d", curfield_type, kpoint_index, curfield_band);
  for (int i = 0; i < 6; ++i) {
    comp_sum[i] /= (energy_sum == 0 ? 1 : energy_sum);
    if (i % 2 == 1) {
         meep::master_printf(", %g", comp_sum[i] + comp_sum[i-1]);
    }
  }
  meep::master_printf("\n");

  /* The return value is a list of 7 items: the total energy,
     followed by the 6 elements of the comp_sum array (the fraction
     of the energy in the real/imag. parts of each field component). */

  retval.push_back(energy_sum * vol / H.N);

  for (int i = 0; i < 6; ++i) {
    retval.push_back(comp_sum[i]);
  }

  return retval;
}


std::vector<mpb_real> mode_solver::get_output_k() {
  std::vector<mpb_real> output_k;

  output_k.push_back(R[0][0]*mdata->current_k[0]
                     + R[0][1]*mdata->current_k[1]
                     + R[0][2]*mdata->current_k[2]);

  output_k.push_back(R[1][0]*mdata->current_k[0]
                     + R[1][1]*mdata->current_k[1]
                     + R[1][2]*mdata->current_k[2]);

  output_k.push_back(R[2][0]*mdata->current_k[0]
                     + R[2][1]*mdata->current_k[1]
                     + R[2][2]*mdata->current_k[2]);
  return output_k;
}

void mode_solver::multiply_bloch_phase() {

  std::vector<mpb_real> kvector = get_output_k();

  int dims[] = {mdata->nx, mdata->ny, mdata->nz};
  int local_dims[] = {mdata->local_nx, mdata->ny, mdata->nz};
  int start[] = {mdata->local_x_start, 0, 0};

  mpb_real s[3]; /* the step size between grid points dotted with k */
  std::vector<scalar_complex> phasex(local_dims[0]);
  std::vector<scalar_complex> phasey(local_dims[1]);
  std::vector<scalar_complex> phasez(local_dims[2]);

  for (int i = 0; i < 3; ++i) {
    s[i] = TWOPI * kvector[i] / dims[i];
  }

  /* cache exp(ikx) along each of the directions, for speed */
  for (int i = 0; i < local_dims[0]; ++i) {
    mpb_real phase = s[0] * (i + start[0]);
    phasex[i].re = cos(phase);
    phasex[i].im = sin(phase);
  }
  for (int j = 0; j < local_dims[1]; ++j) {
    mpb_real phase = s[1] * (j + start[1]);
    phasey[j].re = cos(phase);
    phasey[j].im = sin(phase);
  }
  for (int k = 0; k < local_dims[2]; ++k) {
    mpb_real phase = s[2] * (k + start[2]);
    phasez[k].re = cos(phase);
    phasez[k].im = sin(phase);
  }

  /* Now, multiply field by exp(i k*r): */
  for (int i = 0; i < local_dims[0]; ++i) {
    scalar_complex px = phasex[i];

    for (int j = 0; j < local_dims[1]; ++j) {
      scalar_complex py;
      mpb_real re = phasey[j].re;
      mpb_real im = phasey[j].im;
      py.re = px.re * re - px.im * im;
      py.im = px.re * im + px.im * re;

      for (int k = 0; k < local_dims[2]; ++k) {
        int ijk = ((i*local_dims[1] + j)*local_dims[2] + k)*3;
        mpb_real p_re, p_im;
        mpb_real re = phasez[k].re, im = phasez[k].im;

        p_re = py.re * re - py.im * im;
        p_im = py.re * im + py.im * re;

        for (int component = 0; component < 3; ++component) {
          int ijkc = ijk + component;
          re = curfield[ijkc].re; im = curfield[ijkc].im;
          curfield[ijkc].re = re * p_re - im * p_im;
          curfield[ijkc].im = im * p_re + re * p_im;
        }
      }
    }
  }
}

// Replace the current field with its scalar divergence; only works for Bloch fields
void mode_solver::compute_field_divergence()
{
  scalar *field = (scalar *) curfield;
  scalar *field2 = mdata->fft_data == mdata->fft_data2 ? field :
    (field == mdata->fft_data ? mdata->fft_data2 : mdata->fft_data);
  mpb_real scale;

  if (!curfield || !strchr("dhbec", curfield_type)) {
    meep::master_fprintf(stderr, "A Bloch-periodic field must be loaded.\n");
    return;
  }

  /* convert back to Fourier space */
  maxwell_compute_fft(-1, mdata, field, field2, 3, 3, 1);

  /* compute (k+G) dot field */
  for (int i = 0; i < mdata->other_dims; ++i) {
    for (int j = 0; j < mdata->last_dim; ++j) {
      int ij = i * mdata->last_dim_size + j;
      k_data cur_k = mdata->k_plus_G[ij];
      /* k+G = |k+G| (m x n) */
      mpb_real kx = cur_k.kmag * (cur_k.my*cur_k.nz-cur_k.mz*cur_k.ny);
      mpb_real ky = cur_k.kmag * (cur_k.mz*cur_k.nx-cur_k.mx*cur_k.nz);
      mpb_real kz = cur_k.kmag * (cur_k.mx*cur_k.ny-cur_k.my*cur_k.nz);
      ASSIGN_SCALAR(field2[ij], SCALAR_RE(field2[3*ij+0]) * kx +
                    SCALAR_RE(field2[3*ij+1]) * ky + SCALAR_RE(field2[3*ij+2]) * kz,
                    SCALAR_IM(field2[3*ij+0]) * kx + SCALAR_IM(field2[3*ij+1]) * ky +
                    SCALAR_IM(field2[3*ij+2]) * kz);
    }
  }

  /* convert scalar field back to position space */
  maxwell_compute_fft(+1, mdata, field2, field, 1, 1, 1);

  // multiply by i (from divergence) and normalization (from FFT) and 2*pi (from k+G)
  scale = TWOPI / H.N;
  int N = mdata->fft_output_size;
  for (int i = 0; i < N; ++i) {
    CASSIGN_SCALAR(curfield[i], -CSCALAR_IM(curfield[i]) * scale, CSCALAR_RE(curfield[i]) * scale);
  }

  curfield_type = 'C'; // complex (Bloch) scalar field
}

/* Fix the phase of the current field (e/h/b/d) to a canonical value.
   Also changes the phase of the corresponding eigenvector by the
   same amount, so that future calculations will have a consistent
   phase.

   The following procedure is used, derived from a suggestion by Doug
   Allan of Corning: First, choose the phase to maximize the sum of
   the squares of the real parts of the components.  This doesn't fix
   the overall sign, though.  That is done (after incorporating the
   above phase) by: (1) find the largest absolute value of the real
   part, (2) find the point with the greatest spatial array index that
   has |real part| at least half of the largest value, and (3) make
   that point positive.

   In the case of inversion symmetry, on the other hand, the overall phase
   is already fixed, to within a sign, by the choice to make the Fourier
   transform purely real.  So, in that case we simply pick a sign, in
   a manner similar to (2) and (3) above. */
void mode_solver::fix_field_phase()
{
  mpb_real sq_sum2[2] = {0,0};
  mpb_real sq_sum[2];
  mpb_real maxabs = 0.0;

  int maxabs_index = 0;
  int maxabs_sign = 1;
  int i;

  double theta;
  scalar phase;

  if (!curfield || !strchr("dhbecv", curfield_type)) {
    meep::master_fprintf(stderr, "The D/H/E field must be loaded first.\n");
    return;
  }
  int N = mdata->fft_output_size * 3;

  /* Compute the phase that maximizes the sum of the squares of
     the real parts of the components.  Equivalently, maximize
     the real part of the sum of the squares. */
  for (i = 0; i < N; ++i) {
    mpb_real a = curfield[i].re;
    mpb_real b = curfield[i].im;
    sq_sum2[0] += a*a - b*b;
    sq_sum2[1] += 2*a*b;
  }
  mpi_allreduce(sq_sum2, sq_sum, 2, mpb_real, SCALAR_MPI_TYPE, MPI_SUM, mpb_comm);

  /* compute the phase = exp(i*theta) maximizing the real part of
     the sum of the squares.  i.e., maximize:
       cos(2*theta)*sq_sum[0] - sin(2*theta)*sq_sum[1] */
  theta = 0.5 * atan2(-sq_sum[1], sq_sum[0]);
  phase.re = cos(theta);
  phase.im = sin(theta);

  /* Next, fix the overall sign.  We do this by first computing the
     maximum |real part| of the jmax component (after multiplying
     by phase), and then finding the last spatial index at which
     |real part| is at least half of this value.  The sign is then
     chosen to make the real part positive at that point.

        (Note that we can't just make the point of maximum |real part|
         positive, as that would be ambiguous in the common case of an
         oscillating field within the unit cell.)

        In the case of inversion symmetry (!SCALAR_COMPLEX), we work with
        (real part - imag part) instead of (real part), to insure that we
        have something that is nonzero somewhere. */

  for (i = 0; i < N; ++i) {
    mpb_real r = fabs(curfield[i].re * phase.re - curfield[i].im * phase.im);

    if (r > maxabs) {
      maxabs = r;
    }
  }

  mpi_allreduce_1(&maxabs, mpb_real, SCALAR_MPI_TYPE, MPI_MAX, mpb_comm);

  for (i = N - 1; i >= 0; --i) {
    mpb_real r = curfield[i].re * phase.re - curfield[i].im * phase.im;

    if (fabs(r) >= 0.5 * maxabs) {
         maxabs_index = i;
         maxabs_sign = r < 0 ? -1 : 1;
         break;
    }
  }

  if (i >= 0) { /* convert index to global index in distributed array: */
    maxabs_index += mdata->local_y_start * mdata->nx * mdata->nz;
  }

  {
    /* compute maximum index and corresponding sign over all the
       processors, using the MPI_MAXLOC reduction operation: */
    struct twoint_struct {int i; int s;} x;
    x.i = maxabs_index; x.s = maxabs_sign;
    mpi_allreduce_1(&x, struct twoint_struct, MPI_2INT, MPI_MAXLOC, mpb_comm);
    maxabs_index = x.i; maxabs_sign = x.s;
  }

  ASSIGN_SCALAR(phase, SCALAR_RE(phase)*maxabs_sign, SCALAR_IM(phase)*maxabs_sign);

  meep::master_printf("Fixing %c-field (band %d) phase by %g + %gi; "
                      "max ampl. = %g\n", curfield_type, curfield_band,
                      SCALAR_RE(phase), SCALAR_IM(phase), maxabs);

  /* Now, multiply everything by this phase, *including* the
     stored "raw" eigenvector in H, so that any future fields
     that we compute will have a consistent phase: */
  for (i = 0; i < N; ++i) {
    mpb_real a = curfield[i].re;
    mpb_real b = curfield[i].im;
    curfield[i].re = a*SCALAR_RE(phase) - b*SCALAR_IM(phase);
    curfield[i].im = a*SCALAR_IM(phase) + b*SCALAR_RE(phase);
  }
  for (int i = 0; i < H.n; ++i) {
    mpb_real bbbb_re = H.data[i*H.p + curfield_band - 1].re;
    mpb_real bbbb_im = H.data[i*H.p + curfield_band - 1].im;
    mpb_real cccc_re = phase.re;
    mpb_real cccc_im = phase.im;
    H.data[i*H.p + curfield_band - 1].re = bbbb_re * cccc_re - bbbb_im * cccc_im;
    H.data[i*H.p + curfield_band - 1].im = bbbb_re * cccc_im + bbbb_im * cccc_re;
  }
}

void mode_solver::get_lattice(double data[3][3]) {
  matrix3x3_to_arr(data, Rm);
}

double mode_solver::get_eigensolver_flops() {
  return eigensolver_flops;
}

int mode_solver::get_iterations() {
  return iterations;
}

std::vector<mpb_real> mode_solver::compute_zparities() {
  std::vector<mpb_real> z_parity(num_bands);
  double *d = maxwell_zparity(H, mdata);

  for (int i = 0; i < num_bands; ++i) {
    z_parity[i] = d[i];
  }

  free(d);
  return z_parity;
}

std::vector<mpb_real> mode_solver::compute_yparities() {
  std::vector<mpb_real> y_parity(num_bands);
  double *d = maxwell_yparity(H, mdata);

  for (int i = 0; i < num_bands; ++i) {
    y_parity[i] = d[i];
  }

  free(d);
  return y_parity;
}

/* Compute the group velocity dw/dk in the given direction d (where
   the length of d is ignored).  d is in the reciprocal lattice basis.
   Should only be called after solve_kpoint.  Returns a list of the
   group velocities, one for each band, in units of c. */
std::vector<mpb_real> mode_solver::compute_group_velocity_component(vector3 d) {
  curfield_reset(); // has the side effect of overwriting curfield scratch

  if (!mdata) {
    meep::master_fprintf(stderr, "mode_solver.init must be called first!\n");
    return std::vector<mpb_real>(0);
  }
  if (!kpoint_index) {
    meep::master_fprintf(stderr, "mode_solver.solve_kpoint must be called first!\n");
    return std::vector<mpb_real>(0);
  }

  /* convert d to unit vector in Cartesian coords: */
  d = unit_vector3(matrix3x3_vector3_mult(Gm, d));
  mpb_real u[] = {d.x, d.y, d.z};

  std::vector<mpb_real> group_v(num_bands);
  std::vector<mpb_real> gv_scratch(num_bands * 2);

  /* now, compute group_v.items = diag Re <H| curl 1/eps i u x |H>: */

  /* ...we have to do this in blocks of eigensolver_block_size since
     the work matrix W[0] may not have enough space to do it all at once. */

  for (int ib = 0; ib < num_bands; ib += Hblock.alloc_p) {
    if (ib + mdata->num_bands > num_bands) {
      maxwell_set_num_bands(mdata, num_bands - ib);
      evectmatrix_resize(&W[0], num_bands - ib, 0);
      evectmatrix_resize(&Hblock, num_bands - ib, 0);
    }
    maxwell_compute_H_from_B(mdata, H, Hblock, (scalar_complex *) mdata->fft_data, ib, 0, Hblock.p);
    maxwell_ucross_op(Hblock, W[0], mdata, u);
    evectmatrix_XtY_diag_real(Hblock, W[0], gv_scratch.data(), gv_scratch.data() + group_v.size());
    {
      for (int ip = 0; ip < Hblock.p; ++ip)
        group_v[ib + ip] = gv_scratch[ip];
    }
  }

  /* Reset scratch matrix sizes: */
  evectmatrix_resize(&Hblock, Hblock.alloc_p, 0);
  evectmatrix_resize(&W[0], W[0].alloc_p, 0);
  maxwell_set_num_bands(mdata, Hblock.alloc_p);

  /* The group velocity is given by:

    grad_k(omega)*d = grad_k(omega^2)*d / 2*omega
       = grad_k(<H|maxwell_op|H>)*d / 2*omega
       = Re <H| curl 1/eps i u x |H> / omega

    Note that our k is in units of 2*Pi/a, and omega is in
    units of 2*Pi*c/a, so the result will be in units of c. */

  for (int i = 0; i < num_bands; ++i) {
    if (freqs[i] == 0) { /* v is undefined in this case */
         group_v[i] = 0.0;  /* just set to zero */
    }
    else {
      group_v[i] /= negative_epsilon_ok ? sqrt(fabs(freqs[i])) : freqs[i];
    }
  }

  return group_v;
}

bool mode_solver::with_hermitian_epsilon() {
#ifdef WITH_HERMITIAN_EPSILON
  return true;
#else
  return false;
#endif
}

/* For curfield and energy density, compute the fraction of the energy
   that resides inside the given list of geometric objects.   Later
   objects in the list have precedence, just like the ordinary
   geometry list. */
double mode_solver::compute_energy_in_objects(geometric_object_list objects) {

  mpb_real *energy = (mpb_real *)curfield;
  mpb_real energy_sum = 0;

  if (!curfield || !strchr("DHBR", curfield_type)) {
    meep::master_fprintf(stderr, "The D or H energy density must be loaded first.\n");
    return 0.0;
  }

  for (int i = 0; i < objects.num_items; ++i) {
    geom_fix_object(objects.items[i]);
  }

  int n1 = mdata->nx;
  int n2 = mdata->ny;
  int n3 = mdata->nz;

  mpb_real s1 = geometry_lattice.size.x / n1;
  mpb_real s2 = geometry_lattice.size.y / n2;
  mpb_real s3 = geometry_lattice.size.z / n3;
  mpb_real c1 = n1 <= 1 ? 0 : geometry_lattice.size.x * 0.5;
  mpb_real c2 = n2 <= 1 ? 0 : geometry_lattice.size.y * 0.5;
  mpb_real c3 = n3 <= 1 ? 0 : geometry_lattice.size.z * 0.5;

  LOOP_XYZ(mdata) {
    vector3 p;
    int n;
    p.x = i1 * s1 - c1; p.y = i2 * s2 - c2; p.z = i3 * s3 - c3;
    for (n = objects.num_items - 1; n >= 0; --n) {
      if (point_in_periodic_fixed_objectp(p, objects.items[n])) {
        // TODO:
        // if (((meep_geom::material_data *)objects.items[n].material)->which_subclass == MATERIAL_TYPE_SELF) {
        //   break; /* treat as a "nothing" object */
        // }
        energy_sum += energy[xyz_index];
        break;
      }
    }
  }}}

  mpi_allreduce_1(&energy_sum, mpb_real, SCALAR_MPI_TYPE, MPI_SUM, mpb_comm);
  energy_sum *= vol / H.N;
  return energy_sum;
}
} // namespace meep_mpb
