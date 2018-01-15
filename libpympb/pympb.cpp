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
                         geometric_object_list geom):
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

  // TODO: Write class with copy constructor or eliminate the need to copy this
  default_md = new meep_geom::material_data_struct();
  default_md->which_subclass = _default_material->which_subclass;
  // TODO: Assuming the suscebtibility lists are empty
  memcpy( &(default_md->medium), &_default_material->medium, sizeof(meep_geom::medium_struct));
  // TODO: user_func, user_data
  // TODO: epsilon_data, epsilon_dims

  default_material = (void*)default_md;

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
    delete (meep_geom::material_type)geometry.items[i].material;
  }
  delete[] geometry.items;
  // TODO: susceptibilites in default_md->medium
  delete default_md;
}

bool mode_solver::using_mup() {
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
      if (using_mup() && block_size < num_bands) {
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
} // namespace meep_mpb
