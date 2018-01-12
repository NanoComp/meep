#include <algorithm>
#include <complex>
#include <cstdlib>
#include <iostream>

#include "pympb.hpp"
#include "meep/mympi.hpp"

namespace py_mpb {

const double inf = 1.0e20;

// TODO: Remove globals
geom_box_tree geometry_tree;
geom_box_tree restricted_tree;

// TODO: medium_struct should have these values by default. Create constructor.
meep_geom::medium_struct vacuum_medium = {
  {1.0, 1.0, 1.0}, /* epsilon_diag    */
  {0.0, 0.0, 0.0}, /* epsilon_offdiag */
  {1.0, 1.0, 1.0}, /* mu_diag         */
  {0.0, 0.0, 0.0}, /* mu_offdiag      */
  {0, 0},          /* E_susceptibilities */
  {0, 0},          /* H_susceptibilities */
  {0.0, 0.0, 0.0}, /* E_chi2_diag     */
  {0.0, 0.0, 0.0}, /* E_chi3_diag     */
  {0.0, 0.0, 0.0}, /* H_chi2_diag     */
  {0.0, 0.0, 0.0}, /* H_chi3_diag     */
  {0.0, 0.0, 0.0}, /* D_conductivity_diag  */
  {0.0, 0.0, 0.0}  /* B_conductivity_diag  */
};

// TODO: What do we need in here?
struct eps_data {
  void *placeholder;
};

// /* When we are solving for a few bands at a time, we solve for the
//    upper bands by "deflation"--by continually orthogonalizing them
//    against the already-computed lower bands.  (This constraint
//    commutes with the eigen-operator, of course, so all is well.) */

// typedef struct {
//   evectmatrix Y;   the vectors to orthogonalize against; Y must itself be normalized (Yt B Y = 1) 
//   evectmatrix BY;  /* B * Y */
//   int p;  /* the number of columns of Y to orthogonalize against */
//   scalar *S;  /* a matrix for storing the dot products; should have at least p * X.p elements (see below for X) */
//   scalar *S2; /* a scratch matrix the same size as S */
// } deflation_data;

/* Linearly interpolate a given point in a 3d grid of data.  The point
   coordinates should be in the range [0,1], or at the very least [-1,2]
   ... anything outside [0,1] is *mirror* reflected into [0,1] */
static mpb_real linear_interpolate(mpb_real rx, mpb_real ry, mpb_real rz, mpb_real *data,
                                   int nx, int ny, int nz, int stride) {
     int x, y, z, x2, y2, z2;
     mpb_real dx, dy, dz;

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

static void material_epsmu(meep_geom::material_type material, symmetric_matrix *epsmu,
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
static void epsilon_file_material(meep_geom::material_data *md, vector3 p)
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

  double interpolate_result = linear_interpolate(rx, ry, rz, md->epsilon_data,
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

void get_material_pt(meep_geom::material_type &material, vector3 p) {
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
      md->medium = vacuum_medium;
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

// Pass this to `set_maxwell_dielectric`
void dielectric_function(symmetric_matrix *eps, symmetric_matrix *eps_inv,
                         const mpb_real r[3], void *epsilon_data) {

  // TODO: What should epsilon_data contain?
  (void)epsilon_data;
  meep_geom::material_type mat;
  vector3 p = {r[0], r[1], r[2]};
  get_material_pt(mat, p);
  material_epsmu(mat, eps, eps_inv);
}

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

mode_solver::mode_solver(int num_bands, int parity, double resolution,
                         lattice lat, double tolerance, meep_geom::material_data *_default_material,
                         geometric_object_list geom):
  num_bands(num_bands),
  parity(parity),
  last_parity(-2),
  kpoint_index(0),
  negative_epsilon_ok(false),
  resolution(resolution),
  tolerance(tolerance),
  eigensolver_nwork(3),
  eigensolver_block_size(-11),
  mdata(NULL) {

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
  geometry.num_items = geom.num_items;
  geometry.items = new geometric_object[geometry.num_items];

  // TODO: Avoid the need to copy the geometric objects, or write classes with copy constructors
  for (int i = 0; i < geometry.num_items; ++i) {
    geometric_object_copy(&geom.items[i], &geometry.items[i]);

    geometry.items[i].material = new meep_geom::material_data();
    memcpy(geometry.items[i].material, geom.items[i].material, sizeof(meep_geom::material_data));

    meep_geom::material_type m1 = (meep_geom::material_type)geometry.items[i].material;
    meep_geom::material_type m2 = (meep_geom::material_type)geom.items[i].material;

    int num_E_suceptibilites = m1->medium.E_susceptibilities.num_items;
    if (num_E_suceptibilites > 0) {
      m1->medium.E_susceptibilities.items = new meep_geom::susceptibility_struct[num_E_suceptibilites];
      memcpy(m1->medium.E_susceptibilities.items, m2->medium.E_susceptibilities.items,
             num_E_suceptibilites * sizeof(meep_geom::susceptibility_struct));
    }

    int num_H_suceptibilites = m1->medium.H_susceptibilities.num_items;
    if (num_H_suceptibilites > 0) {
      m1->medium.H_susceptibilities.items = new meep_geom::susceptibility_struct[num_H_suceptibilites];
      memcpy(m1->medium.H_susceptibilities.items, m2->medium.H_susceptibilities.items,
             num_H_suceptibilites * sizeof(meep_geom::susceptibility_struct));
    }
  }
}

mode_solver::~mode_solver() {

  destroy_maxwell_data(mdata);
  destroy_geom_box_tree(geometry_tree);

  for (int i = 0; i < geometry.num_items; ++i) {
    geometric_object_destroy(geometry.items[i]);
    delete[] ((meep_geom::material_type)geometry.items[i].material)->medium.E_susceptibilities.items;
    delete[] ((meep_geom::material_type)geometry.items[i].material)->medium.H_susceptibilities.items;
    delete geometry.items[i].material;
  }
  delete[] geometry.items;
}

bool mode_solver::using_mup() {
  return mdata && mdata->mu_inv != NULL;
}

void mode_solver::init(int p, bool reset_fields) {
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

  int s[] = {
    lat.size.x == 0 ? 1 : lat.size.x,
    lat.size.y == 0 ? 1 : lat.size.y,
    lat.size.z == 0 ? 1 : lat.size.z
  };

  // TODO: Currently getting R and G twice.
  for (int i = 0; i < 3; ++i) {
    R[i][i] = s[i];
    G[i][i] = 1 / R[i][i]; // recip. latt. vectors / 2 pi
  }

  eps_data ed;
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

  destroy_geom_box_tree(geometry_tree);

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

  set_maxwell_dielectric(mdata, mesh_size, R, G, dielectric_function, NULL, &ed);

  if (check_maxwell_dielectric(mdata, 0)) {
    meep::abort("invalid dielectric function for MPB");
  }

  // if (!have_old_fields) {
  meep::master_printf("Allocating fields...\n");
  H = create_evectmatrix(n[0] * n[1] * n[2], 2, num_bands, local_N, N_start, alloc_N);
  nwork_alloc = eigensolver_nwork + (mdata->mu_inv != NULL);

  for (int i = 0; i < nwork_alloc; ++i) {
    W[i] = create_evectmatrix(n[0] * n[1] * n[2], 2, block_size, local_N, N_start, alloc_N);
  }

  if (block_size < num_bands) {
    Hblock = create_evectmatrix(n[0] * n[1] * n[2], 2, block_size, local_N, N_start, alloc_N);
  } else {
    Hblock = H;
    if (using_mup() && block_size < num_bands) {
      muinvH = create_evectmatrix(n[0] * n[1] * n[2], 2, num_bands, local_N, N_start, alloc_N);
    } else {
      muinvH = H;
    }
  }
  // }

  set_parity(p);

  // if (!have_old_fields || reset_fields) {
  randomize_fields();
  // }

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
  kpoint_index = 0;  /* reset index */
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

  // TODO
  int ib0 = 0;

  mpb_real *eigvals = new mpb_real[num_bands];

  // deflation_data deflation;

  // /* Set up deflation data: */
  // if (muinvH.data != Hblock.data) {
  //   deflation.Y = H;
  //   deflation.BY = muinvH.data != H.data ? muinvH : H;
  //   deflation.p = 0;
  //   deflation.S = new scalar[H.p * Hblock.p];
  //   deflation.S2 = new scalar[H.p * Hblock.p];
  // }

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
      int in, ip;
      for (in = 0; in < Hblock.n; ++in) {
        for (ip = 0; ip < Hblock.p; ++ip) {
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
  }

  // TODO: Get this to python
  mpb_real *freqs = new mpb_real[num_bands];

  kpoint_index += 1;

  meep::master_printf("%sfreqs:, %d, %g, %g, %g, %g", parity, kpoint_index, (double)k[0],
                      (double)k[1], (double)k[2], vector3_norm(matrix3x3_vector3_mult(Gm, kpoint)));

  for (int i = 0; i < num_bands; ++i) {
    freqs[i] = negative_epsilon_ok ? eigvals[i] : sqrt(eigvals[i]);
    meep::master_printf(", %g", freqs[i]);
  }
  meep::master_printf("\n");

  delete eigvals;
}


// void add_eigenmode_source(int band_num, const vector3 &kpoint, bool match_frequency,
//                      int parity, double resolution, double tolerance) {
//                      // std::complex<double> amp, std::complex<double> A(const vec &)) {

//   // if (resolution <= 0) {
//   //   resolution = 2 * gv.a; // default to twice resolution
//   // }

//   int n[3];
//   int local_N;
//   int N_start;
//   int alloc_N;
//   int mesh_size[3] = {1,1,1};

//   mpb_real k[3] = {0,0,0};
//   mpb_real kcart[3] = {0,0,0};

//   double s[3] = {0,0,0};
//   double o[3] = {0,0,0};

//   mpb_real R[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
//   mpb_real G[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
//   mpb_real kdir[3] = {0,0,0};

//   // double omega_src = real(src.frequency()), kscale = 1.0;
//   double match_tol = tolerance * 10;

  // if (d == NO_DIRECTION || coordinate_mismatch(gv.dim, d))
  //   abort("invalid direction in add_eigenmode_source");
  // if (where.dim != gv.dim || eig_vol.dim != gv.dim)
  //   abort("invalid volume dimensionality in add_eigenmode_source");

  // if (!eig_vol.contains(where))
  //   abort("invalid grid_volume in add_eigenmode_source (WHERE must be in EIG_VOL)");

  // switch (gv.dim) {
  // case D3:
  //   o[0] = eig_vol.in_direction_min(X);
  //   o[1] = eig_vol.in_direction_min(Y);
  //   o[2] = eig_vol.in_direction_min(Z);
  //   s[0] = eig_vol.in_direction(X);
  //   s[1] = eig_vol.in_direction(Y);
  //   s[2] = eig_vol.in_direction(Z);
  //   k[0] = kpoint.in_direction(X);
  //   k[1] = kpoint.in_direction(Y);
  //   k[2] = kpoint.in_direction(Z);
  //   break;
  // case D2:
  //   o[0] = eig_vol.in_direction_min(X);
  //   o[1] = eig_vol.in_direction_min(Y);
  //   s[0] = eig_vol.in_direction(X);
  //   s[1] = eig_vol.in_direction(Y);
  //   k[0] = kpoint.in_direction(X);
  //   k[1] = kpoint.in_direction(Y);
  //   break;
  // case D1:
  //   o[2] = eig_vol.in_direction_min(Z);
  //   s[2] = eig_vol.in_direction(Z);
  //   k[2] = kpoint.in_direction(Z);
  //   break;
  // default:
  //   abort("unsupported dimensionality in add_eigenmode_source");
  // }

  // master_printf("KPOINT: %g, %g, %g\n", k[0], k[1], k[2]);

  // if match_frequency is true, all we need is a direction for k
  // and a crude guess for its value; we must supply this if k==0.
  // if (match_frequency && k[0] == 0 && k[1] == 0 && k[2] == 0) {
  //   k[d-X] = omega_src * sqrt(get_eps(eig_vol.center()));
  //   master_printf("NEW KPOINT: %g, %g, %g\n", k[0], k[1], k[2]);
  //   if (s[d-X] > 0) {
  //     k[d-X] *= s[d-X]; // put k in G basis (inverted when we compute kcart)
  //     if (fabs(k[d-X]) > 0.4)  // ensure k is well inside the Brillouin zone
  //   k[d-X] = k[d-X] > 0 ? 0.4 : -0.4;
  //     master_printf("NEWER KPOINT: %g, %g, %g\n", k[0], k[1], k[2]);
  //   }
  // }

  // for (int i = 0; i < 3; ++i) {
  //   n[i] = int(resolution * s[i] + 0.5);

  //   if (n[i] == 0) {
  //     n[i] = 1;
  //   }

  //   R[i][i] = s[i] = s[i] == 0 ? 1 : s[i];
  //   G[i][i] = 1 / R[i][i]; // recip. latt. vectors / 2 pi
  // }

  // for (int i = 0; i < 3; ++i)
  //   for (int j = 0; j < 3; ++j)
  //     kcart[i] += G[j][i] * k[j];

  // double klen0 = sqrt(k[0]*k[0]+k[1]*k[1]+k[2]*k[2]);
  // double klen = sqrt(kcart[0]*kcart[0]+kcart[1]*kcart[1]+kcart[2]*kcart[2]);

  // if (klen == 0.0) {
  //   if (match_frequency) {
  //     // abort("need nonzero kpoint guess to match frequency");
  //     abort();
  //   }
  //   klen = 1;
  // }

  // kdir[0] = kcart[0] / klen;
  // kdir[1] = kcart[1] / klen;
  // kdir[2] = kcart[2] / klen;

  // maxwell_data *mdata = create_maxwell_data(n[0], n[1], n[2], &local_N, &N_start,
  //                                           &alloc_N, band_num, band_num);

  // if (local_N != n[0] * n[1] * n[2])
  //   abort("MPI version of MPB library not supported");

  // set_maxwell_data_parity(mdata, parity);
  // update_maxwell_data_k(mdata, k, G[0], G[1], G[2]);

  // if (k[0] == 0 && k[1] == 0 && k[2] == 0) {
  //   evectmatrix H; H.p = band_num; H.c = 2;
  //   band_num -= maxwell_zero_k_num_const_bands(H, mdata);
  //   if (band_num == 0)
  //     abort("zero-frequency bands at k=0 are ill-defined");
  // }

  // eps_data ed;
  // set_maxwell_dielectric(mdata, mesh_size, R, G, dielectric_function, NULL, &ed);

//   if (check_maxwell_dielectric(mdata, 0))
//     abort("invalid dielectric function for MPB");

//   evectmatrix H = create_evectmatrix(n[0] * n[1] * n[2], 2, band_num, local_N,
//                                      N_start, alloc_N);

//   for (int i = 0; i < H.n * H.p; ++i) {
//     ASSIGN_SCALAR(H.data[i], rand() * 1.0/RAND_MAX, rand() * 1.0/RAND_MAX);
//   }

//   mpb_real *eigvals = new mpb_real[band_num];
//   int num_iters;
//   evectmatrix W[3];
//   for (int i = 0; i < 3; ++i)
//     W[i] = create_evectmatrix(n[0] * n[1] * n[2], 2, band_num, local_N, N_start, alloc_N);

//   evectconstraint_chain *constraints = NULL;
//   constraints = evect_add_constraint(constraints, maxwell_parity_constraint, (void *) mdata);

//   // if (k[0] == 0 && k[1] == 0 && k[2] == 0)
//   //   constraints = evect_add_constraint(constraints, maxwell_zero_k_constraint, (void *) mdata);

//   mpb_real knew[3];
//   for (int i = 0; i < 3; ++i)
//     knew[i] = k[i];

//   do {
//     eigensolver(H, eigvals, maxwell_operator, (void *) mdata,
// #if MPB_VERSION_MAJOR > 1 || (MPB_VERSION_MAJOR == 1 && MPB_VERSION_MINOR >= 6)
//                 NULL, NULL, /* eventually, we can support mu here */
// #endif
//         maxwell_preconditioner2, (void *) mdata,
//         evectconstraint_chain_func,
//         (void *) constraints,
//         W, 3,
//         tolerance, &num_iters,
//         EIGS_DEFAULT_FLAGS |
//         (am_master() && !quiet ? EIGS_VERBOSE : 0));

//     if (!quiet) {
//       master_printf("MPB solved for omega_%d(%g,%g,%g) = %g after %d iters\n",
//                     band_num, knew[0],knew[1],knew[2],
//                     sqrt(eigvals[band_num-1]), num_iters);
//     }

//     if (match_frequency) {
//       // copy desired single eigenvector into scratch arrays
//       evectmatrix_resize(&W[0], 1, 0);
//       evectmatrix_resize(&W[1], 1, 0);

//       for (int i = 0; i < H.n; ++i) {
//         W[0].data[i] = H.data[H.p-1 + i * H.p];
//       }

//       // compute the group velocity in the k direction
//       maxwell_ucross_op(W[0], W[1], mdata, kdir); // W[1] = (dTheta/dk) W[0]
//       mpb_real v, vscratch; // v = Re( W[0]* (dTheta/dk) W[0] ) = g. velocity
//       evectmatrix_XtY_diag_real(W[0], W[1], &v, &vscratch);
//       v /= sqrt(eigvals[band_num - 1]);

//       // return to original size
//       evectmatrix_resize(&W[0], band_num, 0);
//       evectmatrix_resize(&W[1], band_num, 0);

//       // update k via Newton step
//       kscale = kscale - (sqrt(eigvals[band_num - 1]) - omega_src) / (v*klen0);

//       if (!quiet) {
//         master_printf("Newton step: group velocity v=%g, kscale=%g\n", v, kscale);
//       }

//       if (kscale < 0 || kscale > 100) {
//         abort("Newton solver not converging -- need a better starting kpoint");
//       }

//       for (int i = 0; i < 3; ++i)
//         knew[i] = k[i] * kscale;

//       update_maxwell_data_k(mdata, knew, G[0], G[1], G[2]);
//     }
//   } while (match_frequency &&
//            fabs(sqrt(eigvals[band_num - 1]) - omega_src) > omega_src * match_tol);

//   evect_destroy_constraints(constraints);

//   for (int i = 0; i < 3; ++i)
//     destroy_evectmatrix(W[i]);

//   // src_time *src_mpb = src.clone();

//   if (!match_frequency) {
//     src_mpb->set_frequency(omega_src = sqrt(eigvals[band_num - 1]));
//   }

//   complex<mpb_real> *cdata = (complex<mpb_real> *) mdata->fft_data;
//   meep_mpb_A_s = s;
//   meep_mpb_A_n = n;
//   meep_mpb_A_data = cdata;
//   // meep_mpb_A_center = eig_vol.center() - where.center();
//   meep_mpb_A_A = A ? A : one;

//   maxwell_compute_h_from_H(mdata, H, (scalar_complex*)cdata, band_num - 1, 1);

//   /* choose deterministic phase, maximizing power in real part;
//      see fix_field_phase routine in MPB.*/
//   {
//     int i;
//     int N = mdata->fft_output_size * 3;

//     double sq_sum0 = 0;
//     double sq_sum1 = 0;
//     double maxabs = 0.0;
//     double theta;

//     for (i = 0; i < N; ++i) {
//       double a = real(cdata[i]);
//       double b = imag(cdata[i]);
//       sq_sum0 += a*a - b*b;
//       sq_sum1 += 2*a*b;
//     }

//     theta = 0.5 * atan2(-sq_sum1, sq_sum0);
//     complex<mpb_real> phase(cos(theta), sin(theta));

//     for (i = 0; i < N; ++i) {
//       double r = fabs(real(cdata[i] * phase));
//       if (r > maxabs) {
//         maxabs = r;
//       }
//     }

//     for (i = N-1; i >= 0 && fabs(real(cdata[i] * phase)) < 0.5 * maxabs; --i)
//       ;

//     if (real(cdata[i] * phase) < 0) {
//       phase = -phase;
//     }

//     for (i = 0; i < N; ++i)
//       cdata[i] *= phase;

//     complex<mpb_real> *hdata = (complex<mpb_real> *) H.data;

//     for (i = 0; i < H.n; ++i)
//       hdata[i*H.p + (band_num-1)] *= phase;
//   }

//   if (is_D(c0)) {
//     c0 = direction_component(Ex, component_direction(c0));
//   }
//   if (is_B(c0)) {
//     c0 = direction_component(Hx, component_direction(c0));
//   }

//   // use principle of equivalence to obtain equivalent currents
//   FOR_ELECTRIC_COMPONENTS(c)
//     if (gv.has_field(c) && (c0 == Centered || c0 == c)
//         && component_direction(c) != d
//         && (gv.dim != D2 || !(parity & (EVEN_Z_PARITY | ODD_Z_PARITY))
//         || ((parity & EVEN_Z_PARITY) && !is_tm(c))
//         || ((parity & ODD_Z_PARITY) && is_tm(c)))) {
//       // E current source = d x (eigenmode H)
//       if ((d + 1) % 3 == component_direction(c) % 3) {
//         meep_mpb_A_component = (d + 2) % 3;
//         add_volume_source(c, *src_mpb, where, meep_mpb_A, -amp);
//       }
//       else {
//         meep_mpb_A_component = (d + 1) % 3;
//         add_volume_source(c, *src_mpb, where, meep_mpb_A, amp);
//       }
//     }

//   maxwell_compute_d_from_H(mdata, H, (scalar_complex*)cdata, band_num - 1, 1);

//   { // d_from_H actually computes -omega*D (see mpb/src/maxwell/maxwell_op.c)
//     double scale = -1.0 / omega_src;
//     int N = mdata->fft_output_size * 3;

//     for (int i = 0; i < N; ++i)
//       cdata[i] *= scale;
//   }

//   maxwell_compute_e_from_d(mdata, (scalar_complex*)cdata, 1);

//   // use principle of equivalence to obtain equivalent currents
//   FOR_MAGNETIC_COMPONENTS(c)
//     if (gv.has_field(c) && (c0 == Centered || c0 == c)
//         && component_direction(c) != d
//         && (gv.dim != D2 || !(parity & (EVEN_Z_PARITY | ODD_Z_PARITY))
//         || ((parity & EVEN_Z_PARITY) && !is_tm(c))
//         || ((parity & ODD_Z_PARITY) && is_tm(c)))) {
//       // H current source = - d x (eigenmode E)
//       if ((d + 1) % 3 == component_direction(c) % 3) {
//         meep_mpb_A_component = (d + 2) % 3;
//         add_volume_source(c, *src_mpb, where, meep_mpb_A, amp);
//       }
//       else {
//         meep_mpb_A_component = (d + 1) % 3;
//         add_volume_source(c, *src_mpb, where, meep_mpb_A, -amp);
//       }
//     }

//   delete src_mpb;
//   destroy_evectmatrix(H);
//   delete[] eigvals;
//   destroy_maxwell_data(mdata);

// }
} // namespace meep_mpb