#include <algorithm>
#include <complex>
#include <cstdlib>
#include <iostream>

#include "pympb.hpp"
#include "meep/mympi.hpp"

namespace py_mpb {

const double inf = 1.0e20;

// TODO: Temporary matrixio stuff
#if defined(HAVE_HDF5)
/* don't use new HDF5 1.8 API (which isn't even fully documented yet, grrr) */
#  define H5_USE_16_API 1
#  include <hdf5.h>
typedef hid_t matrixio_id_;
/* HDF5 changed this datatype in their interfaces starting in version 1.6.4 */
#  if H5_VERS_MAJOR > 1 \
     || (H5_VERS_MAJOR == 1 && H5_VERS_MINOR > 6) \
     || (H5_VERS_MAJOR == 1 && H5_VERS_MINOR == 6 && H5_VERS_RELEASE > 3)
typedef hsize_t start_t;
#  else
typedef hssize_t start_t;
#  endif
#else /* no HDF */
typedef int matrixio_id_; /* dummy */
#endif

typedef struct {
     matrixio_id_ id;
     int parallel;
} matrixio_id;

// This is the function passed to `set_maxwell_dielectric`
void dielectric_function(symmetric_matrix *eps, symmetric_matrix *eps_inv,
                         const mpb_real r[3], void *epsilon_data) {

  mode_solver *ms = static_cast<mode_solver *>(epsilon_data);
  meep_geom::material_type mat;
  vector3 p = {r[0], r[1], r[2]};
  ms->get_material_pt(mat, p);
  ms->material_epsmu(mat, eps, eps_inv);
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

/******* mode_solver *******/

mode_solver::mode_solver(int num_bands, int parity, double resolution, lattice lat,
                         double tolerance, meep_geom::material_data *_default_material,
                         geometric_object_list geom, bool reset_fields):
  num_bands(num_bands),
  parity(parity),
  resolution(resolution),
  tolerance(tolerance),
  eigensolver_nwork(3),
  eigensolver_block_size(-11),
  last_parity(-2),
  negative_epsilon_ok(false),
  mdata(NULL),
  mtdata(NULL),
  curfield(NULL),
  curfield_band(0),
  curfield_type('-'),
  kpoint_index(0) {

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
  // typemaps once this call returns to python.
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

  // TODO: Support mu
  // switch (md->which_subclass) {
  // case material_data::MEDIUM:
  // case material_data::MATERIAL_FILE:
  // case material_data::MATERIAL_USER:
  //   epsmu->m00 = md->medium.mu_diag.x;
  //   epsmu->m11 = md->medium.mu_diag.y;
  //   epsmu->m22 = md->medium.mu_diag.z;
  //   epsmu->m01 = md->medium.mu_offdiag.x;
  //   epsmu->m02 = md->medium.mu_offdiag.y;
  //   epsmu->m12 = md->medium.mu_offdiag.z;
  //   sym_matrix_invert(epsmu_inv,epsmu);
  //   break;

  // case material_data::PERFECT_METAL:
  //   epsmu->m00 = 1.0;
  //   epsmu->m11 = 1.0;
  //   epsmu->m22 = 1.0;
  //   epsmu_inv->m00 = 1.0;
  //   epsmu_inv->m11 = 1.0;
  //   epsmu_inv->m22 = 1.0;
  //   epsmu->m01 = epsmu->m02 = epsmu->m12 = 0.0;
  //   epsmu_inv->m01 = epsmu_inv->m02 = epsmu_inv->m12 = 0.0;
  //   break;
  // default:
  //   meep::abort("unknown material type");
  }
}

// return material of the point p from the file (assumed already read)
void mode_solver::epsilon_file_material(meep_geom::material_data *md, vector3 p)
{
  default_material = (void*) md;

  if (md->which_subclass != meep_geom::material_data::MATERIAL_FILE) {
    meep::abort("epsilon-input-file only works with a type=file default-material");
  }

  if (!(md->epsilon_data)) {
    return;
  }

  meep_geom::medium_struct *mm = &(md->medium);

  double rx = geometry_lattice.size.x == 0
    ? 0 : 0.5 + (p.x-geometry_center.x) / geometry_lattice.size.x;
  double ry = geometry_lattice.size.y == 0
    ? 0 : 0.5 + (p.y-geometry_center.y) / geometry_lattice.size.y;
  double rz = geometry_lattice.size.z == 0
    ? 0 : 0.5 + (p.z-geometry_center.z) / geometry_lattice.size.z;

  double interpolate_result = meep_geom::linear_interpolate(rx, ry, rz, md->epsilon_data,
                                                            md->epsilon_dims[0],
                                                            md->epsilon_dims[1],
                                                            md->epsilon_dims[2], 1);

  mm->epsilon_diag.x = interpolate_result;
  mm->epsilon_diag.y = interpolate_result;
  mm->epsilon_diag.z = interpolate_result;

  mm->epsilon_offdiag.x = 0;
  mm->epsilon_offdiag.y = 0;
  mm->epsilon_offdiag.z = 0;
}

void mode_solver::get_material_pt(meep_geom::material_type &material, vector3 p) {
  boolean inobject;
  material = (meep_geom::material_type)material_of_unshifted_point_in_tree_inobject(p, restricted_tree, &inobject);
  meep_geom::material_data *md = material;

  switch(md->which_subclass) {
    // material read from file: interpolate to get properties at r
    case meep_geom::material_data::MATERIAL_FILE:
      if (md->epsilon_data) {
        epsilon_file_material(md, p);
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

  meep::master_printf("Creating Maxwell data\n");
  mdata = create_maxwell_data(n[0], n[1], n[2], &local_N, &N_start, &alloc_N, num_bands, num_bands);

  double s[] = {
    lat.size.x == 0 ? 1 : lat.size.x,
    lat.size.y == 0 ? 1 : lat.size.y,
    lat.size.z == 0 ? 1 : lat.size.z
  };

  // TODO: Currently getting R and G twice.
  for (int i = 0; i < 3; ++i) {
    R[i][i] = s[i];
    G[i][i] = 1 / R[i][i]; // recip. latt. vectors / 2 pi
  }

  int mesh_size[] = {3, 3, 3};

  // init_epsilon
  meep::master_printf("Mesh size is %d.\n", mesh_size[0]);

  // matrix3x3 version of R.
  matrix3x3 Rm;

  Rm.c0 = vector3_scale(s[0], geometry_lattice.basis.c0);
  Rm.c1 = vector3_scale(s[1], geometry_lattice.basis.c1);
  Rm.c2 = vector3_scale(s[2], geometry_lattice.basis.c2);

  meep::master_printf("Lattice vectors:\n");
  meep::master_printf("     (%g, %g, %g)\n", Rm.c0.x, Rm.c0.y, Rm.c0.z);
  meep::master_printf("     (%g, %g, %g)\n", Rm.c1.x, Rm.c1.y, Rm.c1.z);
  meep::master_printf("     (%g, %g, %g)\n", Rm.c2.x, Rm.c2.y, Rm.c1.z);

  mpb_real vol = fabs(matrix3x3_determinant(Rm));
  meep::master_printf("Cell volume = %g\n", vol);

  Gm = matrix3x3_inverse(matrix3x3_transpose(Rm));
  meep::master_printf("Reciprocal lattice vectors (/ 2 pi):\n");
  meep::master_printf("     (%g, %g, %g)\n", Gm.c0.x, Gm.c0.y, Gm.c0.z);
  meep::master_printf("     (%g, %g, %g)\n", Gm.c1.x, Gm.c1.y, Gm.c1.z);
  meep::master_printf("     (%g, %g, %g)\n", Gm.c2.x, Gm.c2.y, Gm.c2.z);

  geom_fix_objects0(geometry);

  meep::master_printf("Geometric objects:\n");
  if (meep::am_master()) {
    for (int i = 0; i < geometry.num_items; ++i) {
      display_geometric_object_info(5, geometry.items[i]);

      // meep_geom::material_type m = (meep_geom::material_type)geometry.items[i].material;
      // if (m->which_subclass == meep_geom::material_data::MEDIUM) {
      //   printf("%*sepsilon = %g, mu = %g\n", 5 + 5, "", m->medium.epsilon_diag.x,
      //          m->medium.mu_diag.x);
      // }
    }
  }

  // TODO: Need to destroy tree from previous runs?
  // destroy_geom_box_tree(geometry_tree);

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

  int tree_depth;
  int tree_nobjects;
  geom_box_tree_stats(geometry_tree, &tree_depth, &tree_nobjects);
  meep::master_printf("Geometric object tree has depth %d and %d object nodes"
                      " (vs. %d actual objects)\n", tree_depth, tree_nobjects, geometry.num_items);

  // TODO
  // reset_epsilon();
  meep::master_printf("Initializing epsilon function...\n");
  set_maxwell_dielectric(mdata, mesh_size, R, G, dielectric_function, NULL, static_cast<void *>(this));

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

  for (i = 0; i < H.n * H.p; ++i) {
    ASSIGN_SCALAR(H.data[i], rand() * 1.0 / RAND_MAX, rand() * 1.0 / RAND_MAX);
  }
}

void mode_solver::solve_kpoint(vector3 kpoint) {
  mpb_real k[] = {
    kpoint.x,
    kpoint.y,
    kpoint.z
  };

  meep::master_printf("solve_kpoint (%g,%g,%g):\n", k[0], k[1], k[2]);

  if (!kpoint_index && meep::am_master()) {
    printf("%sfreqs:, k index, k1, k2, k3, kmag/2pi", parity_string(mdata));

    for (int i = 0; i < num_bands; ++i) {
      printf(", %s%sband %d", parity_string(mdata), mdata->parity == NO_PARITY ? "" : " ", i + 1);
    }
    printf("\n");
  }

  update_maxwell_data_k(mdata, k, G[0], G[1], G[2]);

  // TODO
  // if (mtdata) {  /* solving for bands near a target frequency */
  // }

  // TODO
  // if (eigensolver_davidsonp) {
  // }

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
    ib0 = 0;
  }

  mpb_real *eigvals = new mpb_real[num_bands];
  int total_iters = 0;

  for (int ib = ib0; ib < num_bands; ib += Hblock.alloc_p) {
    evectconstraint_chain *constraints;
    int num_iters;

    /* don't solve for too many bands if the block size doesn't divide
       the number of bands: */
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

      // deflation.p = ib-ib0;
      // if (deflation.p > 0) {
      //   if (deflation.BY.data != H.data) {
      //     evectmatrix_resize(&deflation.BY, deflation.p, 0);
      //     maxwell_muinv_operator(H, deflation.BY, (void *) mdata, 1, deflation.BY);
      //   }
      //   constraints = evect_add_constraint(constraints, deflation_constraint, &deflation);
      // }
    }

    eigensolver(Hblock, eigvals + ib, maxwell_operator, (void *) mdata, NULL, NULL, maxwell_preconditioner2,
                (void *) mdata, evectconstraint_chain_func, (void *) constraints, W, 3, tolerance, &num_iters,
                0); // EIGS_DEFAULT_FLAGS | (meep::am_master() && !quiet ? EIGS_VERBOSE : 0));

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

  // TODO: Get this to python
  mpb_real *freqs = new mpb_real[num_bands];

  kpoint_index += 1;

  meep::master_printf("%sfreqs:, %d, %g, %g, %g, %g", parity_string(mdata), kpoint_index, (double)k[0],
                      (double)k[1], (double)k[2], vector3_norm(matrix3x3_vector3_mult(Gm, kpoint)));

  for (int i = 0; i < num_bands; ++i) {
    freqs[i] = negative_epsilon_ok ? eigvals[i] : sqrt(eigvals[i]);
    meep::master_printf(", %g", freqs[i]);
  }
  meep::master_printf("\n");

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

#ifndef SCALAR_COMPLEX
  /* most points need to be counted twice, by rfftw output symmetry: */
    {
      int last_index;
      int last_dim_stored = mdata->last_dim_size / (sizeof(scalar_complex)/sizeof(scalar));
#ifdef HAVE_MPI
      if (mdata->nz == 1) { /* 2d calculation: 1st dim. is truncated one */
        last_index = i / mdata->nx + mdata->local_y_start;
      }
      else {
        last_index = i % last_dim_stored;
      }
#else
      last_index = i % last_dim_stored;
#endif
      if (last_index != 0 && 2*last_index != mdata->last_dim) {
        eps_mean += epsilon[i];
        eps_inv_mean += 1/epsilon[i];
        if (epsilon[i] > 1.0001) {
           ++fill_count;
        }
      }
    }
#endif
  }

  // TODO
  // mpi_allreduce_1(&eps_mean, mpb_real, SCALAR_MPI_TYPE, MPI_SUM, mpb_comm);
  // mpi_allreduce_1(&eps_inv_mean, mpb_real, SCALAR_MPI_TYPE, MPI_SUM, mpb_comm);
  // mpi_allreduce_1(&eps_low, mpb_real, SCALAR_MPI_TYPE, MPI_MIN, mpb_comm);
  // mpi_allreduce_1(&eps_high, mpb_real, SCALAR_MPI_TYPE, MPI_MAX, mpb_comm);
  // mpi_allreduce_1(&fill_count, int, MPI_INT, MPI_SUM, mpb_comm);
  N = mdata->nx * mdata->ny * mdata->nz;
  eps_mean /= N;
  eps_inv_mean = N/eps_inv_mean;

  meep::master_printf("epsilon: %g-%g, mean %g, harm. mean %g, %g%% > 1, %g%% \"fill\"\n",
                      eps_low, eps_high, eps_mean, eps_inv_mean, (100.0 * fill_count) / N,
                      eps_high == eps_low ? 100.0 : 100.0 * (eps_mean-eps_low) / (eps_high-eps_low));
}

/* given the field in curfield, store it to HDF (or whatever) using
   the matrixio (fieldio) routines.  Allow the component to be specified
   (which_component 0/1/2 = x/y/z, -1 = all) for vector fields.
   Also allow the user to specify a prefix string for the filename. */
void mode_solver::output_field_to_file(int which_component, string filename_prefix) {
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

  /* where to put "otherhalf" block of output, only used for real scalars */
  int last_dim_index = 0;
  int last_dim_start = 0;
  int last_dim_size = 0;
  int first_dim_start = 0;
  int first_dim_size = 0;
  int write_start0_special = 0;

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
// #ifndef SCALAR_COMPLEX
//   /* Ugh, hairy.  See also maxwell_vectorfield_otherhalf. */
//   if (dims[2] == 1) {
//     last_dim_index = 0;
//     first_dim_size = local_dims[0];
//     first_dim_start = dims[0] - (start[0] + local_dims[0] - 1);

//     if (start[0] == 0)
//       --first_dim_size; /* DC frequency is not in other half */
//     if (start[0] + local_dims[0] == mdata->last_dim_size / 2 && dims[0] % 2 == 0) {
//       --first_dim_size; /* Nyquist frequency is not in other half */
//       ++first_dim_start;
//     }

//     last_dim_start = first_dim_start;
//     last_dim_size = first_dim_size;
//   }
//   else {
//     last_dim_index = 2;
//     local_dims[last_dim_index] = mdata->last_dim_size / 2;
//     if (start[0] == 0) {
//       first_dim_size = local_dims[0] - 1;
//       first_dim_start = dims[0] - first_dim_size;
//       write_start0_special = 1;
//     }
//     else {
//       first_dim_start = dims[0] - (start[0] + local_dims[0] - 1);
//       first_dim_size = local_dims[0];
//     }
//     last_dim_start = local_dims[last_dim_index];
//     last_dim_size = dims[last_dim_index] - local_dims[last_dim_index];
//   }
// #endif /* ! SCALAR_COMPLEX */
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
#ifndef SCALAR_COMPLEX
  last_dim_index = dims[2] == 1 ? (dims[1] == 1 ? 0 : 1) : 2;
  local_dims[last_dim_index] = mdata->last_dim_size / 2;
  last_dim_start = local_dims[last_dim_index];
  last_dim_size = dims[last_dim_index] - local_dims[last_dim_index];
  first_dim_start = last_dim_index ? 0 : last_dim_start;
  first_dim_size = last_dim_index ? local_dims[0] : last_dim_size;
#endif
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
            kpoint_index, curfield_band, freqs.items[curfield_band - 1]);
    fname2 = fix_fname(fname, filename_prefix, mdata, 1);
    mpi_one_printf("Outputting fields to %s...\n", fname2);
    file_id = matrixio_create(fname2);
    free(fname2);
    fieldio_write_complex_field(curfield, 3, dims, local_dims, start, which_component, 3,
                                output_k, file_id, 0, data_id);

#ifndef SCALAR_COMPLEX
    /* Here's where it gets hairy. */
    maxwell_vectorfield_otherhalf(mdata, curfield, output_k[0], output_k[1], output_k[2]);
    start[last_dim_index] = last_dim_start;
    local_dims[last_dim_index] = last_dim_size;
    start[0] = first_dim_start;
    local_dims[0] = first_dim_size;
    if (write_start0_special) {
      /* The conjugated array half may be discontiguous.
         First, write the part not containing start[0], and
         then write the start[0] slab. */
      fieldio_write_complex_field(curfield + 3 * local_dims[1] * local_dims[2],
                                  3, dims, local_dims, start, which_component, 3, NULL,
                                  file_id, 1, data_id);
      local_dims[0] = 1;
      start[0] = 0;
      fieldio_write_complex_field(curfield, 3, dims,local_dims,start, which_component, 3,
                                  NULL, file_id, 1, data_id);
    }
    else {
      fieldio_write_complex_field(curfield, 3, dims,local_dims,start, which_component, 3,
                                  NULL, file_id, 1, data_id);
    }
#endif

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
            curfield_band, freqs.items[curfield_band - 1]);
    fname2 = fix_fname(fname, filename_prefix, mdata, 1);
    mpi_one_printf("Outputting complex scalar field to %s...\n", fname2);
    file_id = matrixio_create(fname2);
    free(fname2);
    fieldio_write_complex_field(curfield, 3, dims, local_dims, start, which_component, 1,
                                output_k, file_id, 0, data_id);

#ifndef SCALAR_COMPLEX
    /* Here's where it gets hairy. */
    maxwell_cscalarfield_otherhalf(mdata, curfield, output_k[0], output_k[1], output_k[2]);
    start[last_dim_index] = last_dim_start;
    local_dims[last_dim_index] = last_dim_size;
    start[0] = first_dim_start;
    local_dims[0] = first_dim_size;
    if (write_start0_special) {
      /* The conjugated array half may be discontiguous.
         First, write the part not containing start[0], and
         then write the start[0] slab. */
      fieldio_write_complex_field(curfield + local_dims[1] * local_dims[2], 3, dims,
                                  local_dims, start, which_component, 1, NULL, file_id,
                                  1, data_id);
      local_dims[0] = 1;
      start[0] = 0;
      fieldio_write_complex_field(curfield, 3, dims,local_dims,start, which_component, 1, NULL,
                                  file_id, 1, data_id);
    }
    else {
      fieldio_write_complex_field(curfield, 3, dims,local_dims,start, which_component, 1, NULL,
                                  file_id, 1, data_id);
    }
#endif

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
              curfield_type, kpoint_index, curfield_band, freqs.items[curfield_band - 1]);
    }
    fname2 = fix_fname(fname, filename_prefix, mdata,
             /* no parity suffix for epsilon: */
             curfield_type != 'n' && curfield_type != 'm');
    mpi_one_printf("Outputting %s...\n", fname2);
    file_id = matrixio_create(fname2);
    free(fname2);

    output_scalarfield((real *) curfield, dims, local_dims, start, file_id, "data",
                       last_dim_index, last_dim_start, last_dim_size, first_dim_start,
                       first_dim_size, write_start0_special);

    if (curfield_type == 'n') {
      int c1, c2, inv;
      char dataname[100];

      for (inv = 0; inv < 2; ++inv)
        for (c1 = 0; c1 < 3; ++c1)
          for (c2 = c1; c2 < 3; ++c2) {
            get_epsilon_tensor(c1,c2, 0, inv);
            sprintf(dataname, "%s.%c%c", inv ? "epsilon_inverse" : "epsilon",
                    c1 + 'x', c2 + 'x');
            output_scalarfield((real *) curfield, dims, local_dims, start, file_id, dataname,
                               last_dim_index, last_dim_start, last_dim_size, first_dim_start,
                               first_dim_size, write_start0_special);

#if defined(WITH_HERMITIAN_EPSILON)
            if (c1 != c2) {
              get_epsilon_tensor(c1,c2, 1, inv);
              strcat(dataname, ".i");
#ifndef SCALAR_COMPLEX /* scalarfield_otherhalf isn't right */
              strcat(dataname, ".screwy");
#endif
              output_scalarfield((real *) curfield, dims, local_dims, start, file_id,
                                 dataname, last_dim_index, last_dim_start, last_dim_size,
                                 first_dim_start, first_dim_size, write_start0_special);
            }
#endif
          }
    }
  }
  else {
    mpi_one_fprintf(stderr, "unknown field type!\n");
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
} // namespace meep_mpb
