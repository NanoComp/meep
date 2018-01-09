#include <complex>
#include <cstdlib>

#include "mpb.h"
#include "meepgeom.hpp"

namespace py_mpb {

const double inf = 1.0e20;

static void material_epsmu(meep_geom::material_type material, symmetric_matrix *epsmu,
                           symmetric_matrix *epsmu_inv) {

  material_data *md = material;

  switch (md->which_subclass) {
    case material_data::MEDIUM:
    case material_data::MATERIAL_FILE:
    case material_data::MATERIAL_USER:
      epsmu->m00 = md->medium.epsilon_diag.x;
      epsmu->m11 = md->medium.epsilon_diag.y;
      epsmu->m22 = md->medium.epsilon_diag.z;
      epsmu->m01 = md->medium.epsilon_offdiag.x;
      epsmu->m02 = md->medium.epsilon_offdiag.y;
      epsmu->m12 = md->medium.epsilon_offdiag.z;
      sym_matrix_invert(epsmu_inv, epsmu);
      break;
    case material_data::PERFECT_METAL:
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
      std::cerr << "Unknown material type" << std::endl;
      abort();

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

void get_material_pt(meep_geom::material_type &material, vector3 r) {
  boolean inobject;
  material = (material_type)material_of_unshifted_point_in_tree_inobject(p, restricted_tree,&inobject);
  material_data *md = material;

  switch(md->which_subclass) {
    // material read from file: interpolate to get properties at r
    case material_data::MATERIAL_FILE:
      if (md->epsilon_data) {
        epsilon_file_material(md, p);
      }
      else {
        material = (material_type) default_material;
      }
      return;

    // material specified by user-supplied function: call user
    // function to get properties at r.
    // Note that we initialize the medium to vacuum, so that
    // the user's function only needs to fill in whatever is
    // different from vacuum.
    case material_data::MATERIAL_USER:
      md->medium = vacuum_medium;
      md->user_func(p, md->user_data, &(md->medium));
      // TODO: update this to allow user's function to set
      //       position-dependent susceptibilities. For now
      //       it's an error if the user's function creates
      //       any.
      if ((md->medium.E_susceptibilities.num_items>0) ||
          (md->medium.H_susceptibilities.num_items>0)) {
        abort();
      }
      return;

    // position-independent material or metal: there is nothing to do
    case material_data::MEDIUM:
    case material_data::PERFECT_METAL:
      return;
    default:
      abort();
   }
}

// Argument to set_maxwell_dielectric
void dielectric_function(symmetric_matrix *eps, symmetric_matrix *eps_inv, const real r[3], void *epsilon_data) {
  material_type mat;
  vector3 p = {r[0], r[1], r[2]};
  get_material_pt(mat, p);
  material_epsmu(mat, eps, eps_inv);
}

// TODO: Name?
void add_eigenmode_source(int band_num, const vector3 &kpoint, bool match_frequency,
                     int parity, double resolution, double eigensolver_tol) {
                     // std::complex<double> amp, std::complex<double> A(const vec &)) {

  // if (resolution <= 0) {
  //   resolution = 2 * gv.a; // default to twice resolution
  // }

  int n[3];
  int local_N;
  int N_start;
  int alloc_N;
  int mesh_size[3] = {1,1,1};

  mpb_real k[3] = {0,0,0};
  mpb_real kcart[3] = {0,0,0};

  double s[3] = {0,0,0};
  double o[3] = {0,0,0};

  mpb_real R[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
  mpb_real G[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
  mpb_real kdir[3] = {0,0,0};

  // double omega_src = real(src.frequency()), kscale = 1.0;
  double match_tol = eigensolver_tol * 10;

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

  for (int i = 0; i < 3; ++i) {
    n[i] = int(resolution * s[i] + 0.5);

    if (n[i] == 0) {
      n[i] = 1;
    }

    R[i][i] = s[i] = s[i] == 0 ? 1 : s[i];
    G[i][i] = 1 / R[i][i]; // recip. latt. vectors / 2 pi
  }

  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      kcart[i] += G[j][i] * k[j];

  double klen0 = sqrt(k[0]*k[0]+k[1]*k[1]+k[2]*k[2]);
  double klen = sqrt(kcart[0]*kcart[0]+kcart[1]*kcart[1]+kcart[2]*kcart[2]);

  if (klen == 0.0) {
    if (match_frequency) {
      // abort("need nonzero kpoint guess to match frequency");
      abort();
    }
    klen = 1;
  }

  kdir[0] = kcart[0] / klen;
  kdir[1] = kcart[1] / klen;
  kdir[2] = kcart[2] / klen;

  maxwell_data *mdata = create_maxwell_data(n[0], n[1], n[2], &local_N, &N_start,
                                            &alloc_N, band_num, band_num);

  // if (local_N != n[0] * n[1] * n[2])
  //   abort("MPI version of MPB library not supported");

  set_maxwell_data_parity(mdata, parity);
  update_maxwell_data_k(mdata, k, G[0], G[1], G[2]);

  // if (k[0] == 0 && k[1] == 0 && k[2] == 0) {
  //   evectmatrix H; H.p = band_num; H.c = 2;
  //   band_num -= maxwell_zero_k_num_const_bands(H, mdata);
  //   if (band_num == 0)
  //     abort("zero-frequency bands at k=0 are ill-defined");
  // }

  eps_data ed;
  // eps_data.s = s; eps_data.o = o; eps_data.dim = gv.dim; eps_data.f = this;
  set_maxwell_dielectric(mdata, mesh_size, R, G, get_material_pt, NULL, &eps_data);

  if (check_maxwell_dielectric(mdata, 0))
    abort("invalid dielectric function for MPB");

  evectmatrix H = create_evectmatrix(n[0] * n[1] * n[2], 2, band_num, local_N,
                                     N_start, alloc_N);

  for (int i = 0; i < H.n * H.p; ++i) {
    ASSIGN_SCALAR(H.data[i], rand() * 1.0/RAND_MAX, rand() * 1.0/RAND_MAX);
  }

  mpb_real *eigvals = new mpb_real[band_num];
  int num_iters;
  evectmatrix W[3];
  for (int i = 0; i < 3; ++i)
    W[i] = create_evectmatrix(n[0] * n[1] * n[2], 2, band_num, local_N, N_start, alloc_N);

  evectconstraint_chain *constraints = NULL;
  constraints = evect_add_constraint(constraints, maxwell_parity_constraint, (void *) mdata);

  // if (k[0] == 0 && k[1] == 0 && k[2] == 0)
  //   constraints = evect_add_constraint(constraints, maxwell_zero_k_constraint, (void *) mdata);

  mpb_real knew[3];
  for (int i = 0; i < 3; ++i)
    knew[i] = k[i];

  do {
    eigensolver(H, eigvals, maxwell_operator, (void *) mdata,
#if MPB_VERSION_MAJOR > 1 || (MPB_VERSION_MAJOR == 1 && MPB_VERSION_MINOR >= 6)
                NULL, NULL, /* eventually, we can support mu here */
#endif
        maxwell_preconditioner2, (void *) mdata,
        evectconstraint_chain_func,
        (void *) constraints,
        W, 3,
        eigensolver_tol, &num_iters,
        EIGS_DEFAULT_FLAGS |
        (am_master() && !quiet ? EIGS_VERBOSE : 0));

    if (!quiet) {
      master_printf("MPB solved for omega_%d(%g,%g,%g) = %g after %d iters\n",
                    band_num, knew[0],knew[1],knew[2],
                    sqrt(eigvals[band_num-1]), num_iters);
    }

    if (match_frequency) {
      // copy desired single eigenvector into scratch arrays
      evectmatrix_resize(&W[0], 1, 0);
      evectmatrix_resize(&W[1], 1, 0);

      for (int i = 0; i < H.n; ++i) {
        W[0].data[i] = H.data[H.p-1 + i * H.p];
      }

      // compute the group velocity in the k direction
      maxwell_ucross_op(W[0], W[1], mdata, kdir); // W[1] = (dTheta/dk) W[0]
      mpb_real v, vscratch; // v = Re( W[0]* (dTheta/dk) W[0] ) = g. velocity
      evectmatrix_XtY_diag_real(W[0], W[1], &v, &vscratch);
      v /= sqrt(eigvals[band_num - 1]);

      // return to original size
      evectmatrix_resize(&W[0], band_num, 0);
      evectmatrix_resize(&W[1], band_num, 0);

      // update k via Newton step
      kscale = kscale - (sqrt(eigvals[band_num - 1]) - omega_src) / (v*klen0);

      if (!quiet) {
        master_printf("Newton step: group velocity v=%g, kscale=%g\n", v, kscale);
      }

      if (kscale < 0 || kscale > 100) {
        abort("Newton solver not converging -- need a better starting kpoint");
      }

      for (int i = 0; i < 3; ++i)
        knew[i] = k[i] * kscale;

      update_maxwell_data_k(mdata, knew, G[0], G[1], G[2]);
    }
  } while (match_frequency &&
           fabs(sqrt(eigvals[band_num - 1]) - omega_src) > omega_src * match_tol);

  evect_destroy_constraints(constraints);

  for (int i = 0; i < 3; ++i)
    destroy_evectmatrix(W[i]);

  // src_time *src_mpb = src.clone();

  if (!match_frequency) {
    src_mpb->set_frequency(omega_src = sqrt(eigvals[band_num - 1]));
  }

  complex<mpb_real> *cdata = (complex<mpb_real> *) mdata->fft_data;
  meep_mpb_A_s = s;
  meep_mpb_A_n = n;
  meep_mpb_A_data = cdata;
  // meep_mpb_A_center = eig_vol.center() - where.center();
  meep_mpb_A_A = A ? A : one;

  maxwell_compute_h_from_H(mdata, H, (scalar_complex*)cdata, band_num - 1, 1);

  /* choose deterministic phase, maximizing power in real part;
     see fix_field_phase routine in MPB.*/
  {
    int i;
    int N = mdata->fft_output_size * 3;

    double sq_sum0 = 0;
    double sq_sum1 = 0;
    double maxabs = 0.0;
    double theta;

    for (i = 0; i < N; ++i) {
      double a = real(cdata[i]);
      double b = imag(cdata[i]);
      sq_sum0 += a*a - b*b;
      sq_sum1 += 2*a*b;
    }

    theta = 0.5 * atan2(-sq_sum1, sq_sum0);
    complex<mpb_real> phase(cos(theta), sin(theta));

    for (i = 0; i < N; ++i) {
      double r = fabs(real(cdata[i] * phase));
      if (r > maxabs) {
        maxabs = r;
      }
    }

    for (i = N-1; i >= 0 && fabs(real(cdata[i] * phase)) < 0.5 * maxabs; --i)
      ;

    if (real(cdata[i] * phase) < 0) {
      phase = -phase;
    }

    for (i = 0; i < N; ++i)
      cdata[i] *= phase;

    complex<mpb_real> *hdata = (complex<mpb_real> *) H.data;

    for (i = 0; i < H.n; ++i)
      hdata[i*H.p + (band_num-1)] *= phase;
  }

  if (is_D(c0)) {
    c0 = direction_component(Ex, component_direction(c0));
  }
  if (is_B(c0)) {
    c0 = direction_component(Hx, component_direction(c0));
  }

  // use principle of equivalence to obtain equivalent currents
  FOR_ELECTRIC_COMPONENTS(c)
    if (gv.has_field(c) && (c0 == Centered || c0 == c)
        && component_direction(c) != d
        && (gv.dim != D2 || !(parity & (EVEN_Z_PARITY | ODD_Z_PARITY))
        || ((parity & EVEN_Z_PARITY) && !is_tm(c))
        || ((parity & ODD_Z_PARITY) && is_tm(c)))) {
      // E current source = d x (eigenmode H)
      if ((d + 1) % 3 == component_direction(c) % 3) {
        meep_mpb_A_component = (d + 2) % 3;
        add_volume_source(c, *src_mpb, where, meep_mpb_A, -amp);
      }
      else {
        meep_mpb_A_component = (d + 1) % 3;
        add_volume_source(c, *src_mpb, where, meep_mpb_A, amp);
      }
    }

  maxwell_compute_d_from_H(mdata, H, (scalar_complex*)cdata, band_num - 1, 1);

  { // d_from_H actually computes -omega*D (see mpb/src/maxwell/maxwell_op.c)
    double scale = -1.0 / omega_src;
    int N = mdata->fft_output_size * 3;

    for (int i = 0; i < N; ++i)
      cdata[i] *= scale;
  }

  maxwell_compute_e_from_d(mdata, (scalar_complex*)cdata, 1);

  // use principle of equivalence to obtain equivalent currents
  FOR_MAGNETIC_COMPONENTS(c)
    if (gv.has_field(c) && (c0 == Centered || c0 == c)
        && component_direction(c) != d
        && (gv.dim != D2 || !(parity & (EVEN_Z_PARITY | ODD_Z_PARITY))
        || ((parity & EVEN_Z_PARITY) && !is_tm(c))
        || ((parity & ODD_Z_PARITY) && is_tm(c)))) {
      // H current source = - d x (eigenmode E)
      if ((d + 1) % 3 == component_direction(c) % 3) {
        meep_mpb_A_component = (d + 2) % 3;
        add_volume_source(c, *src_mpb, where, meep_mpb_A, amp);
      }
      else {
        meep_mpb_A_component = (d + 1) % 3;
        add_volume_source(c, *src_mpb, where, meep_mpb_A, -amp);
      }
    }

  delete src_mpb;
  destroy_evectmatrix(H);
  delete[] eigvals;
  destroy_maxwell_data(mdata);

}
} // namespace meep_mpb