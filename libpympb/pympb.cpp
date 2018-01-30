#include <algorithm>
#include <complex>
#include <cstdlib>
#include <cstddef>
#include <iostream>

#include "config.h"
#include "pympb.hpp"
#include "mpb/scalar.h"
#include "meep/mympi.hpp"

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

namespace py_mpb {

// TODO: Placeholder
int mpb_comm;

const double inf = 1.0e20;

// TODO: Replace this functionality with h5py
#include "matrixio.cpp"

// This is the function passed to `set_maxwell_dielectric`
void dielectric_function(symmetric_matrix *eps, symmetric_matrix *eps_inv,
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

// TODO: Store as class member?
/* return a string describing the current parity, used for frequency
   and filename prefixes */
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

mode_solver::mode_solver(int num_bands, int parity, double resolution, lattice lat,
                         double tolerance, meep_geom::material_data *_default_material,
                         geometric_object_list geom, bool reset_fields, bool deterministic):
  num_bands(num_bands),
  parity(parity),
  resolution(resolution),
  tolerance(tolerance),
  eigensolver_nwork(3),
  eigensolver_block_size(-11),
  mesh_size(3),
  last_parity(-2),
  negative_epsilon_ok(false),
  iterations(0),
  vol(0),
  mdata(NULL),
  mtdata(NULL),
  curfield_band(0),
  freqs(num_bands),
  verbose(true),
  deterministic(deterministic),
  kpoint_index(0),
  curfield(NULL),
  curfield_type('-') {

  this->lat = lat;

  geometry_lattice.size.x = lat.size.x;
  geometry_lattice.size.y = lat.size.y;
  geometry_lattice.size.z = lat.size.z;

  for (int i = 0; i < 3; ++i) {
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
  destroy_geom_box_tree(geometry_tree);
}

void mode_solver::material_epsmu(meep_geom::material_type material, symmetric_matrix *epsmu,
                                 symmetric_matrix *epsmu_inv) {

  meep_geom::material_data *md = material;

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

  // TODO: if (field_type != meep::E_stuff)

  // switch (md->which_subclass) {
  //   case meep_geom::material_data::MEDIUM:
  //   case meep_geom::material_data::MATERIAL_FILE:
  //   case meep_geom::material_data::MATERIAL_USER:
  //     epsmu->m00 = md->medium.mu_diag.x;
  //     epsmu->m11 = md->medium.mu_diag.y;
  //     epsmu->m22 = md->medium.mu_diag.z;
  //     epsmu->m01 = md->medium.mu_offdiag.x;
  //     epsmu->m02 = md->medium.mu_offdiag.y;
  //     epsmu->m12 = md->medium.mu_offdiag.z;
  //     maxwell_sym_matrix_invert(epsmu_inv, epsmu);
  //     break;
  //   case meep_geom::material_data::PERFECT_METAL:
  //     epsmu->m00 = 1.0;
  //     epsmu->m11 = 1.0;
  //     epsmu->m22 = 1.0;
  //     epsmu_inv->m00 = 1.0;
  //     epsmu_inv->m11 = 1.0;
  //     epsmu_inv->m22 = 1.0;
  //     epsmu->m01 = epsmu->m02 = epsmu->m12 = 0.0;
  //     epsmu_inv->m01 = epsmu_inv->m02 = epsmu_inv->m12 = 0.0;
  //     break;
  //   default:
  //     meep::abort("unknown material type");
  // }
}

void mode_solver::get_material_pt(meep_geom::material_type &material, vector3 p) {
  boolean inobject;
  material = (meep_geom::material_type)material_of_unshifted_point_in_tree_inobject(p, restricted_tree, &inobject);
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

  n[0] = std::max(resolution * std::ceil(geometry_lattice.size.x), 1.0);
  n[1] = std::max(resolution * std::ceil(geometry_lattice.size.y), 1.0);
  n[2] = std::max(resolution * std::ceil(geometry_lattice.size.z), 1.0);

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

  if (mdata) {
    // TODO: Clean up if mdata is not NULL
  }

  meep::master_printf("Creating Maxwell data...\n");
  mdata = create_maxwell_data(n[0], n[1], n[2], &local_N, &N_start, &alloc_N, num_bands, num_bands);

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
      if (using_mu() && block_size < num_bands) {
        muinvH = create_evectmatrix(N, c, num_bands, local_N, N_start, alloc_N);
      } else {
        muinvH = H;
      }
    }
  }

  set_parity(p);

  if (!have_old_fields || reset_fields) {
    randomize_fields();
  }

  // TODO
  // evectmatrix_flops = eigensolver_flops;
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
    geom_box b0;
    b0.low = vector3_plus(geometry_center, vector3_scale(-0.5, geometry_lattice.size));
    b0.high = vector3_plus(geometry_center, vector3_scale(0.5, geometry_lattice.size));
    /* pad tree boundaries to allow for sub-pixel averaging */
    b0.low.x -= geometry_lattice.size.x / mdata->nx;
    b0.low.y -= geometry_lattice.size.y / mdata->ny;
    b0.low.z -= geometry_lattice.size.z / mdata->nz;
    b0.high.x += geometry_lattice.size.x / mdata->nx;
    b0.high.y += geometry_lattice.size.y / mdata->ny;
    b0.high.z += geometry_lattice.size.z / mdata->nz;
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

  restricted_tree = geometry_tree;

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
  set_maxwell_dielectric(mdata, mesh, R, G, dielectric_function, NULL, static_cast<void *>(this));

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

void mode_solver::set_kpoint_index(int i) {
  kpoint_index = i;
}

void mode_solver::randomize_fields() {
  int i;

  if (!mdata) {
    return;
  }
  meep::master_printf("Initializing fields to random numbers...\n");

  if (deterministic) {
    // seed should be the same for each run, although
    // it should be different for each process.
    // TODO: MPI
    // int rank;
    // MPI_Comm_rank(mpb_comm, &rank);
    srand(314159); // * (rank + 1));
  }

  for (i = 0; i < H.n * H.p; ++i) {
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

  mpb_real *eigvals = new mpb_real[num_bands];

  // TODO
  // flags = eigensolver_flags;
  // if (verbose) {
  //  flags |= EIGS_VERBOSE;
  // }



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

    // TODO
    // if (mtdata) {  /* solving for bands near a target frequency */
    // }

    // TODO
    // if (eigensolver_davidsonp) {
    // }

    eigensolver(Hblock, eigvals + ib, maxwell_operator, (void *) mdata,
                mdata->mu_inv ? maxwell_muinv_operator : NULL, (void *) mdata,
                maxwell_preconditioner2, (void *) mdata, evectconstraint_chain_func,
                (void *) constraints, W, nwork_alloc, tolerance, &num_iters,
                EIGS_DEFAULT_FLAGS);

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

  // TODO
  // eigensolver_flops = evectmatrix_flops

  delete eigvals;
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
  int last_dim = mdata->last_dim;
  int last_dim_stored = mdata->last_dim_size / (sizeof(scalar_complex) / sizeof(scalar));
  int nx = mdata->nx;
  int nz = mdata->nz;
  int local_y_start = mdata->local_y_start;

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

  (void)last_dim;
  (void)last_dim_stored;
  (void)nx;
  (void)nz;
  (void)local_y_start;

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

/* given the field in curfield, store it to HDF (or whatever) using
   the matrixio (fieldio) routines.  Allow the component to be specified
   (which_component 0/1/2 = x/y/z, -1 = all) for vector fields.
   Also allow the user to specify a prefix string for the filename. */
void mode_solver::output_field_to_file(int which_component, char *filename_prefix) {
  char fname[100];
  char *fname2;
  char description[100];

  int dims[3];
  int local_dims[3];
  int start[3] = {0,0,0};

  matrixio_id file_id = {-1,1};

  int attr_dims[2] = {3, 3};
  mpb_real output_k[3]; /* kvector in reciprocal lattice basis */
  mpb_real output_R[3][3];

  if (!curfield) {
    meep::master_fprintf(stderr, "fields, energy dens., or epsilon must be loaded first.\n");
    return;
  }

// TODO: Support MPI
// #ifdef HAVE_MPI
//   /* The first two dimensions (x and y) of the position-space fields
//      are transposed when we use MPI, so we need to transpose everything. */
//   dims[0] = mdata->ny;
//   local_dims[1] = dims[1] = mdata->nx;
//   local_dims[2] = dims[2] = mdata->nz;
//   local_dims[0] = mdata->local_ny;
//   start[0] = mdata->local_y_start;

//   output_k[0] = R[1][0]*mdata->current_k[0] + R[1][1]*mdata->current_k[1]
//                 + R[1][2]*mdata->current_k[2];
//   output_k[1] = R[0][0]*mdata->current_k[0] + R[0][1]*mdata->current_k[1]
//                 + R[0][2]*mdata->current_k[2];
//   output_k[2] = R[2][0]*mdata->current_k[0] + R[2][1]*mdata->current_k[1]
//                 + R[2][2]*mdata->current_k[2];
//   output_R[0][0]=R[1][0]; output_R[0][1]=R[1][1]; output_R[0][2]=R[1][2];
//   output_R[1][0]=R[0][0]; output_R[1][1]=R[0][1]; output_R[1][2]=R[0][2];
//   output_R[2][0]=R[2][0]; output_R[2][1]=R[2][1]; output_R[2][2]=R[2][2];
// #else /* ! HAVE_MPI */
  dims[0] = mdata->nx;
  local_dims[1] = dims[1] = mdata->ny;
  local_dims[2] = dims[2] = mdata->nz;
  local_dims[0] = mdata->local_nx;
  start[0] = mdata->local_x_start;
  output_k[0] = R[0][0]*mdata->current_k[0] + R[0][1]*mdata->current_k[1]
                + R[0][2]*mdata->current_k[2];
  output_k[1] = R[1][0]*mdata->current_k[0] + R[1][1]*mdata->current_k[1]
                + R[1][2]*mdata->current_k[2];
  output_k[2] = R[2][0]*mdata->current_k[0] + R[2][1]*mdata->current_k[1]
                + R[2][2]*mdata->current_k[2];
  output_R[0][0]=R[0][0]; output_R[0][1]=R[0][1]; output_R[0][2]=R[0][2];
  output_R[1][0]=R[1][0]; output_R[1][1]=R[1][1]; output_R[1][2]=R[1][2];
  output_R[2][0]=R[2][0]; output_R[2][1]=R[2][1]; output_R[2][2]=R[2][2];
// #endif /* ! HAVE_MPI */

  if (strchr("Rv", curfield_type)) /* generic scalar/vector field */
    output_k[0] = output_k[1] = output_k[2] = 0.0; /* don't know k */

  if (strchr("dhbecv", curfield_type)) { /* outputting vector field */
    matrixio_id data_id[6] = {{-1,1},{-1,1},{-1,1},{-1,1},{-1,1},{-1,1}};
    int i;

    sprintf(fname, "%c.k%02d.b%02d", curfield_type, kpoint_index, curfield_band);
    if (which_component >= 0) {
      char comp_str[] = ".x";
      comp_str[1] = 'x' + which_component;
      strcat(fname, comp_str);
    }
    sprintf(description, "%c field, kpoint %d, band %d, freq=%g", curfield_type,
            kpoint_index, curfield_band, freqs[curfield_band - 1]);
    fname2 = fix_fname(fname, filename_prefix, mdata, 1);
    meep::master_printf("Outputting fields to %s...\n", fname2);
    file_id = matrixio_create(fname2);
    free(fname2);
    fieldio_write_complex_field(curfield, 3, dims, local_dims, start, which_component, 3,
                                output_k, file_id, 0, data_id);

    for (i = 0; i < 6; ++i) {
      if (data_id[i].id >= 0) {
        matrixio_close_dataset(data_id[i]);
      }
    }
    matrixio_write_data_attr(file_id, "Bloch wavevector", output_k, 1, attr_dims);
  }
  else if (strchr("C", curfield_type)) { /* outputting cmplx scalar field */
    matrixio_id data_id[2] = {{-1,1},{-1,1}};
    int i;

    sprintf(fname, "%c.k%02d.b%02d", curfield_type, kpoint_index, curfield_band);
    sprintf(description, "%c field, kpoint %d, band %d, freq=%g", curfield_type, kpoint_index,
            curfield_band, freqs[curfield_band - 1]);
    fname2 = fix_fname(fname, filename_prefix, mdata, 1);
    meep::master_printf("Outputting complex scalar field to %s...\n", fname2);
    file_id = matrixio_create(fname2);
    free(fname2);
    fieldio_write_complex_field(curfield, 3, dims, local_dims, start, which_component, 1,
                                output_k, file_id, 0, data_id);

    for (i = 0; i < 2; ++i) {
         if (data_id[i].id >= 0) {
        matrixio_close_dataset(data_id[i]);
      }
    }
    matrixio_write_data_attr(file_id, "Bloch wavevector", output_k, 1, attr_dims);
  }
  else if (strchr("DHBnmR", curfield_type)) { /* scalar field */
    if (curfield_type == 'n') {
      sprintf(fname, "epsilon");
      sprintf(description, "dielectric function, epsilon");
    }
    else if (curfield_type == 'm') {
      sprintf(fname, "mu");
      sprintf(description, "permeability mu");
    }
    else {
      sprintf(fname, "%cpwr.k%02d.b%02d", tolower(curfield_type), kpoint_index, curfield_band);
      sprintf(description, "%c field energy density, kpoint %d, band %d, freq=%g",
              curfield_type, kpoint_index, curfield_band, freqs[curfield_band - 1]);
    }
    fname2 = fix_fname(fname, filename_prefix, mdata,
             /* no parity suffix for epsilon: */
             curfield_type != 'n' && curfield_type != 'm');
    meep::master_printf("Outputting %s...\n", fname2);
    file_id = matrixio_create(fname2);
    free(fname2);

    output_scalarfield((mpb_real *) curfield, dims, local_dims, start, file_id, "data");

    if (curfield_type == 'n') {
      int c1, c2, inv;
      char dataname[100];

      for (inv = 0; inv < 2; ++inv)
        for (c1 = 0; c1 < 3; ++c1)
          for (c2 = c1; c2 < 3; ++c2) {
            get_epsilon_tensor(c1,c2, 0, inv);
            sprintf(dataname, "%s.%c%c", inv ? "epsilon_inverse" : "epsilon",
                    c1 + 'x', c2 + 'x');
            output_scalarfield((mpb_real *) curfield, dims, local_dims, start, file_id, dataname);

#if defined(WITH_HERMITIAN_EPSILON)
            if (c1 != c2) {
              get_epsilon_tensor(c1,c2, 1, inv);
              strcat(dataname, ".i");
              output_scalarfield((mpt_real *) curfield, dims, local_dims, start, file_id, dataname);
            }
#endif
          }
    }
  }
  else {
    meep::master_fprintf(stderr, "unknown field type!\n");
  }

  if (file_id.id >= 0) {
    matrixio_write_data_attr(file_id, "lattice vectors", &output_R[0][0], 2, attr_dims);
    matrixio_write_string_attr(file_id, "description", description);
    matrixio_close(file_id);
  }

  /* We have destroyed curfield (by multiplying it by phases,
     and/or reorganizing in the case of real-amplitude fields). */
  curfield_reset();
}

void mode_solver::output_scalarfield(mpb_real *vals,
                                     const int dims[3],
                                     const int local_dims[3],
                                     const int start[3],
                                     matrixio_id file_id,
                                     const char *dataname) {

  matrixio_id data_id = {-1, 1};

  fieldio_write_real_vals(vals, 3, dims, local_dims, start, file_id, 0, dataname, &data_id);

  if (data_id.id >= 0)
    matrixio_close_dataset(data_id);
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

/* Prepend the prefix to the fname, and (if parity_suffix is true)
   append a parity specifier (if any) (e.g. ".te"), returning a new
   string, which should be deallocated with free().  fname or prefix
   may be NULL, in which case they are treated as the empty string. */
char *mode_solver::fix_fname(const char *fname, const char *prefix, maxwell_data *d, int parity_suffix)
{
  int fname_len = fname ? strlen(fname) : 0;
  int prefix_len = prefix ? strlen(prefix) : 0;

  char *s = (char *)malloc(sizeof(char) * (fname_len + prefix_len + 20));
  strcpy(s, prefix ? prefix : "");
  strcat(s, fname ? fname : "");
  if (parity_suffix && d->parity != NO_PARITY) {
    /* assumes parity suffix is less than 20 characters;
       currently it is less than 12 */
    strcat(s, ".");
    strcat(s, parity_string(d));
  }
  return s;
}

void mode_solver::load_eigenvectors(char *filename) {
  meep::master_printf("Loading eigenvectors from \"%s\"...\n", filename);
  evectmatrixio_readall_raw(filename, H);
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

// void mode_solver::get_curfield(std::complex<mpb_real> *cdata, int size) {
//   for (int i = 0; i < size; ++i) {
//     cdata[i] = curfield[i];
//   }
// }


} // namespace meep_mpb
