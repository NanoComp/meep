#include <memory>
#include <math.h>
#include <unistd.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>

#include <meep.hpp>
#include "meep_internals.hpp"
#include "config.h"
using namespace meep;
using std::complex;
using std::max;
using std::min;

const double xsize = 2.0;
const double ysize = 2.0;
const double zsize = 0.6;

const double r = 0.5;
const double eps_k = 2 * pi / 1.0;

double funky_eps_2d(const vec &p_) {
  vec p = p_ - vec(xsize / 2, ysize / 2);
  if (fabs(p & p) < r * r) return 1.0;
  return 2.0 + cos(p.x() * eps_k) * cos(p.y() * eps_k);
}

double funky_eps_3d(const vec &p_) {
  vec p = p_ - vec(xsize / 2, ysize / 2, zsize / 2);
  if (fabs(p & p) < r * r) return 1.0;
  return 2.0 + cos(p.x() * eps_k) * cos(p.y() * eps_k) * cos(p.z() * eps_k);
}

symmetry make_identity(const grid_volume &gv) {
  (void)gv; // unused
  return identity();
}

symmetry make_mirrorx(const grid_volume &gv) { return mirror(X, gv); }

symmetry make_mirrory(const grid_volume &gv) { return mirror(Y, gv); }

symmetry make_mirrorxy(const grid_volume &gv) { return mirror(X, gv) + mirror(Y, gv); }

symmetry make_rotate4z(const grid_volume &gv) { return rotate4(Z, gv); }

typedef symmetry (*symfunc)(const grid_volume &);

const double tol = sizeof(realnum) == sizeof(float) ? 2.5e-4 : 1e-8;
double compare(double a, double b, const char *nam, size_t i0, size_t i1, size_t i2) {
  if (fabs(a - b) > tol * tol + fabs(b) * tol || b != b) {
    master_printf("%g vs. %g differs by\t%g\n", a, b, fabs(a - b));
    master_printf("This gives a fractional error of %g\n", fabs(a - b) / fabs(b));
    meep::abort("Error in %s at (%zd,%zd,%zd)\n", nam, i0, i1, i2);
  }
  return fabs(a - b);
}

double get_reim(complex<double> x, int reim) { return reim ? imag(x) : real(x); }

bool check_2d(double eps(const vec &), double a, int splitting, symfunc Sf, double kx, double ky,
              component src_c, int file_c, volume file_gv, bool real_fields, int expected_rank,
              const char *name, const char *mydirname) {
  const grid_volume gv = vol2d(xsize, ysize, a);
  structure s(gv, eps, no_pml(), Sf(gv), splitting);
  s.set_output_directory(mydirname);
  fields f(&s);

  f.use_bloch(X, real_fields ? 0.0 : kx);
  f.use_bloch(Y, real_fields ? 0.0 : ky);

  if (real_fields) f.use_real_fields();
  f.add_point_source(src_c, 0.3, 2.0, 0.0, 1.0, gv.center(), 1.0, 1);

  if (file_c >= int(Dielectric)) real_fields = true;

  while (f.time() <= 3.0)
    f.step();

  h5file *file = f.open_h5file(name);
  if (is_derived(file_c))
    f.output_hdf5(derived_component(file_c), file_gv, file);
  else
    f.output_hdf5(component(file_c), file_gv, file);

  file->write("stringtest", "Hello, world!\n");

  delete file;
  all_wait();
  sync();
  file = f.open_h5file(name, h5file::READONLY);

  {
    char *str = file->read("stringtest");
    if (strcmp(str, "Hello, world!\n"))
      meep::abort("Failed to read back string test from %s...", name);
    delete[] str;
  }

  // compute corner coordinate of file data
  vec loc0(file_gv.get_min_corner());
  ivec iloc0(gv.dim);
  LOOP_OVER_DIRECTIONS(gv.dim, d) {
    iloc0.set_direction(d, 1 + 2 * int(floor(loc0.in_direction(d) * a - .5)));
    if (file_gv.in_direction(d) == 0.0 &&
        1. - file_gv.in_direction_min(d) * a + 0.5 * iloc0.in_direction(d) <=
            1. + file_gv.in_direction_max(d) * a - 0.5 * (iloc0.in_direction(d) + 2))
      iloc0.set_direction(d, iloc0.in_direction(d) + 2); // snap to grid
  }
  loc0 = gv[iloc0];

  double data_min = meep::infinity, data_max = -meep::infinity;
  double err_max = 0;
  for (int reim = 0; reim < (real_fields ? 1 : 2); ++reim) {
    int rank;
    size_t dims[2] = {1, 1};

    char dataname[256];
    snprintf(dataname, 256, "%s%s", component_name(file_c),
             reim ? ".i" : (real_fields ? "" : ".r"));

    realnum *h5data =
        (realnum *)file->read(dataname, &rank, dims, 2, sizeof(realnum) == sizeof(float));
    file->prevent_deadlock(); // hackery
    if (!h5data) meep::abort("failed to read dataset %s:%s\n", name, dataname);
    if (rank != expected_rank)
      meep::abort("incorrect rank (%d instead of %d) in %s:%s\n", rank, expected_rank, name,
                  dataname);
    if (expected_rank == 1 && file_gv.in_direction_min(X) == file_gv.in_direction_max(X)) {
      dims[1] = dims[0];
      dims[0] = 1;
    }
    vec loc(loc0.dim);
    for (size_t i0 = 0; i0 < dims[0]; ++i0) {
      for (size_t i1 = 0; i1 < dims[1]; ++i1) {
        loc.set_direction(X, loc0.in_direction(X) + i0 * gv.inva);
        loc.set_direction(Y, loc0.in_direction(Y) + i1 * gv.inva);
        ptrdiff_t idx = i0 * dims[1] + i1;

        /* Ugh, for rotational symmetries (which mix up components etc.),
           we can't guarantee that a component is *exactly* the
           same as its rotated version, and we don't know which one
           was written to the file. */
        int cs = file_c;
        complex<double> ph = 1.0;
        double diff = fabs(get_reim(f.get_field(file_c, loc), reim) - h5data[idx]);
        for (int sn = 1; sn < f.S.multiplicity(); ++sn) {
          vec loc2(f.S.transform(loc, sn));
          int cs2 = f.S.transform(file_c, sn);
          complex<double> ph2 = f.S.phase_shift(cs2, -sn);
          double diff2 = fabs(get_reim(f.get_field(cs2, loc2) * ph2, reim) - h5data[idx]);
          if (diff2 < diff) {
            loc = loc2;
            cs = cs2;
            ph = ph2;
            diff = diff2;
          }
        }

        double err =
            compare(h5data[idx], get_reim(f.get_field(cs, loc) * ph, reim), name, i0, i1, 0);
        err_max = max<double>(err, err_max);
        data_min = min<double>(data_min, h5data[idx]);
        data_max = max<double>(data_max, h5data[idx]);
      }
    }
    delete[] h5data;
  }

  // file->remove();
  delete file;

  master_printf("Passed %s (%g..%g), err=%g\n", name, data_min, data_max,
                err_max / max(fabs(data_min), fabs(data_max)));

  return true;
}

bool check_3d(double eps(const vec &), double a, int splitting, symfunc Sf, component src_c,
              int file_c, volume file_gv, bool real_fields, int expected_rank, const char *name,
              const char *mydirname) {
  const grid_volume gv = vol3d(xsize, ysize, zsize, a);
  structure s(gv, eps, no_pml(), Sf(gv), splitting);
  s.set_output_directory(mydirname);
  fields f(&s);

  if (real_fields) f.use_real_fields();
  f.add_point_source(src_c, 0.3, 2.0, 0.0, 1.0, gv.center(), 1.0, 1);

  if (file_c >= Dielectric) real_fields = true;

  while (f.time() <= 3.0)
    f.step();

  h5file *file = f.open_h5file(name);
  if (is_derived(file_c))
    f.output_hdf5(derived_component(file_c), file_gv, file);
  else
    f.output_hdf5(component(file_c), file_gv, file);

  file->write("stringtest", "Hello, world!\n");

  delete file;
  all_wait();
  sync();
  file = f.open_h5file(name, h5file::READONLY);

  {
    char *str = file->read("stringtest");
    if (strcmp(str, "Hello, world!\n"))
      meep::abort("Failed to read back string test from %s...", name);
    delete[] str;
  }

  // compute corner coordinate of file data
  vec loc0(file_gv.get_min_corner());
  ivec iloc0(gv.dim);
  LOOP_OVER_DIRECTIONS(gv.dim, d) {
    iloc0.set_direction(d, 1 + 2 * int(floor(loc0.in_direction(d) * a - .5)));
    if (file_gv.in_direction(d) == 0.0 &&
        1. - file_gv.in_direction_min(d) * a + 0.5 * iloc0.in_direction(d) <=
            1. + file_gv.in_direction_max(d) * a - 0.5 * (iloc0.in_direction(d) + 2))
      iloc0.set_direction(d, iloc0.in_direction(d) + 2); // snap to grid
  }
  loc0 = gv[iloc0];

  double data_min = meep::infinity, data_max = -meep::infinity;
  double err_max = 0;
  for (int reim = 0; reim < (real_fields ? 1 : 2); ++reim) {
    int rank;
    size_t dims[3] = {1, 1, 1};

    char dataname[256];
    snprintf(dataname, 256, "%s%s", component_name(file_c),
             reim ? ".i" : (real_fields ? "" : ".r"));

    realnum *h5data =
        (realnum *)file->read(dataname, &rank, dims, 3, sizeof(realnum) == sizeof(float));
    file->prevent_deadlock(); // hackery
    if (!h5data) meep::abort("failed to read dataset %s:%s\n", name, dataname);
    if (rank != expected_rank)
      meep::abort("incorrect rank (%d instead of %d) in %s:%s\n", rank, expected_rank, name,
                  dataname);
    vec loc(loc0.dim);
    for (size_t i0 = 0; i0 < dims[0]; ++i0) {
      for (size_t i1 = 0; i1 < dims[1]; ++i1) {
        for (size_t i2 = 0; i2 < dims[2]; ++i2) {
          loc.set_direction(X, loc0.in_direction(X) + i0 * gv.inva);
          loc.set_direction(Y, loc0.in_direction(Y) + i1 * gv.inva);
          loc.set_direction(Z, loc0.in_direction(Z) + i2 * gv.inva);
          ptrdiff_t idx = (i0 * dims[1] + i1) * dims[2] + i2;

          /* Ugh, for rotational symmetries (which mix up components etc.),
             we can't guarantee that a component is *exactly* the
             same as its rotated version, and we don't know which one
             was written to the file. */
          int cs = file_c;
          complex<double> ph = 1.0;
          double diff = fabs(get_reim(f.get_field(file_c, loc), reim) - h5data[idx]);
          for (int sn = 1; sn < f.S.multiplicity(); ++sn) {
            vec loc2(f.S.transform(loc, sn));
            int cs2 = f.S.transform(file_c, sn);
            complex<double> ph2 = f.S.phase_shift(cs2, -sn);
            double diff2 = fabs(get_reim(f.get_field(cs2, loc2) * ph2, reim) - h5data[idx]);
            if (diff2 < diff) {
              loc = loc2;
              cs = cs2;
              ph = ph2;
              diff = diff2;
            }
          }

          double err =
              compare(h5data[idx], get_reim(f.get_field(cs, loc) * ph, reim), name, i0, i1, i2);
          err_max = max(err, err_max);
          data_min = min<double>(data_min, h5data[idx]);
          data_max = max<double>(data_max, h5data[idx]);
        }
      }
    }
    delete[] h5data;
  }

  // file->remove();
  delete file;

  master_printf("Passed %s (%g..%g), err=%g\n", name, data_min, data_max,
                err_max / (max(fabs(data_min), fabs(data_max)) + 1e-16));

  return 1;
}

bool check_2d_monitor(double eps(const vec &), double a, int splitting, symfunc Sf, component src_c,
                      int file_c, const vec &pt, bool real_fields, const char *name,
                      const char *mydirname) {
  const grid_volume gv = vol2d(xsize, ysize, a);
  structure s(gv, eps, no_pml(), Sf(gv), splitting);
  s.set_output_directory(mydirname);
  fields f(&s);

  if (real_fields) f.use_real_fields();
  f.add_point_source(src_c, 0.3, 2.0, 0.0, 1.0, gv.center(), 1.0, 1);

  if (file_c >= Dielectric) real_fields = true;

  h5file *file = f.open_h5file(name);

  // compute pt snapped onto dielectric grid
  ivec iloc0(gv.dim);
  LOOP_OVER_DIRECTIONS(gv.dim, d) {
    iloc0.set_direction(d, 1 + 2 * int(floor(pt.in_direction(d) * a - .5)));
    if (1. - pt.in_direction(d) * a + 0.5 * iloc0.in_direction(d) <=
        1. + pt.in_direction(d) * a - 0.5 * (iloc0.in_direction(d) + 2))
      iloc0.set_direction(d, iloc0.in_direction(d) + 2); // snap to grid
  }
  vec pt0(gv[iloc0]);

  const double T = 3.0;
  int NT = int(T / f.dt) + 2;
  complex<double> *mon = new complex<double>[NT];
  while (f.time() <= T) {
    if (is_derived(file_c))
      f.output_hdf5(derived_component(file_c), volume(pt, pt), file, true);
    else
      f.output_hdf5(component(file_c), volume(pt, pt), file, true);
    mon[f.t] = f.get_field(file_c, pt0);
    f.step();
  }

  delete file;
  all_wait();
  sync();
  file = f.open_h5file(name, h5file::READONLY);

  double data_min = meep::infinity, data_max = -meep::infinity;
  double err_max = 0;

  for (int reim = 0; reim < (real_fields ? 1 : 2); ++reim) {
    int rank;
    size_t dims[1] = {1};

    char dataname[256];
    snprintf(dataname, 256, "%s%s", component_name(file_c),
             reim ? ".i" : (real_fields ? "" : ".r"));

    realnum *h5data =
        (realnum *)file->read(dataname, &rank, dims, 2, sizeof(realnum) == sizeof(float));
    file->prevent_deadlock(); // hackery
    if (!h5data) meep::abort("failed to read dataset %s:%s\n", file->file_name(), dataname);
    if (rank != 1) meep::abort("monitor-point data is not one-dimensional");
    if (dims[0] != size_t(f.t)) meep::abort("incorrect size of monitor-point data");

    for (int i = 0; i < f.t; ++i) {
      double err = compare(h5data[i], get_reim(mon[i], reim), name, i, 0, 0);
      err_max = max(err, err_max);
      data_min = min<double>(data_min, h5data[i]);
      data_max = max<double>(data_max, h5data[i]);
    }
    delete[] h5data;
  }

  delete[] mon;

  // file->remove();
  delete file;

  master_printf("Passed %s (%g..%g), err=%g\n", name, data_min, data_max,
                err_max / max(fabs(data_min), fabs(data_max)));

  return 1;
}

int main(int argc, char **argv) {
  const double a = 10.0;
  initialize mpi(argc, argv);
  int chances;
  verbosity = 0;
  std::unique_ptr<const char[]> temp_dir(make_output_directory());
#ifdef HAVE_HDF5
  const double pad1 = 0.314159, pad2 = 0.27183, pad3 = 0.14142;

  volume gv_2d[4] = {
      volume(vec(pad1, pad2), vec(xsize - pad2, ysize - pad1)),
      volume(vec(-pad1, -pad2), vec(2 * xsize - pad2, 2 * ysize - pad1)),
      volume(vec(pad1, pad2), vec(xsize - pad2, pad2)),
      volume(vec(pad1, pad2), vec(pad1, pad2)),
  };
  char gv_2d_name[4][20] = {"plane", "plane-supercell", "line", "point"};
  int gv_2d_rank[4] = {2, 2, 1, 0};
  int tm_c[5] = {Dielectric, Ez, Hy, Sx, D_EnergyDensity};
  symfunc Sf2[5] = {make_identity, make_mirrorx, make_mirrory, make_mirrorxy, make_rotate4z};
  char Sf2_name[5][32] = {"identity", "mirrorx", "mirrory", "mirrorxy", "rotate4z"};
  double Sf2_kx[5] = {0.3, 0, 0.3, 0, 0};
  double Sf2_ky[5] = {0.2, 0.2, 0, 0, 0};

#if 0
  master_printf("Running initial check...\n");
  if (!check_2d(funky_eps_2d, a, 1,
		Sf2[3], Sf2_kx[3], Sf2_ky[3],
		Ez, tm_c[3], gv_2d[1],
		1, gv_2d_rank[1], "initial check", temp_dir))
    return 1;
#endif

  /* this test takes too long, so only do 1/chances of the cases,
     "randomly" selected */
  srand(314159); /* deterministic "rand" */
  chances = argc > 1 ? atoi(argv[1]) : 5;

  for (int iS = 0; iS < 5; ++iS)
    for (int splitting = 0; splitting < 5; ++splitting)
      for (int igv = 0; igv < 4; ++igv)
        for (int ic = 0; ic < 5; ++ic)
          for (int use_real = 1; use_real >= 0; --use_real)
            if (broadcast(0, rand()) % chances == 0) {
              char name[1024];
              snprintf(name, 1024, "check_2d_tm_%s_%d_%s_%s%s", Sf2_name[iS], splitting,
                       gv_2d_name[igv], component_name(tm_c[ic]), use_real ? "_r" : "");
              master_printf("Checking %s...\n", name);
              if (!check_2d(funky_eps_2d, a, splitting, Sf2[iS], Sf2_kx[iS], Sf2_ky[iS], Ez,
                            tm_c[ic], gv_2d[igv], use_real, gv_2d_rank[igv], name, temp_dir.get()))
                return 1;
            }

  for (int iS = 0; iS < 5; ++iS)
    for (int splitting = 0; splitting < 5; ++splitting)
      for (int ic = 0; ic < 4; ++ic)
        for (int use_real = 1; use_real >= 0; --use_real)
          if (broadcast(0, rand()) % chances == 0) {
            char name[1024];
            snprintf(name, 1024, "check_2d_monitor_tm_%s_%d_%s%s", Sf2_name[iS], splitting,
                     component_name(tm_c[ic]), use_real ? "_r" : "");
            master_printf("Checking %s...\n", name);
            if (!check_2d_monitor(funky_eps_2d, a, splitting, Sf2[iS], Ez, tm_c[ic],
                                  vec(pad1, pad2), use_real, name, temp_dir.get()))
              return 1;
          }

  volume gv_3d[4] = {
      volume(vec(pad1, pad2, pad3), vec(xsize - pad2, ysize - pad1, zsize - pad3)),
      volume(vec(pad1, pad2, pad3), vec(xsize - pad2, ysize - pad1, pad3)),
      volume(vec(pad1, pad2, pad3), vec(xsize - pad2, pad2, pad3)),
      volume(vec(pad1, pad2, pad3), vec(pad1, pad2, pad3)),
  };
  char gv_3d_name[4][10] = {"volume", "plane", "line", "point"};
  int gv_3d_rank[4] = {3, 2, 1, 0};
  int c3d[7] = {Ex, Dielectric, Dy, Ez, Sz, H_EnergyDensity, EnergyDensity};
  symfunc Sf3[3] = {make_identity, make_mirrorxy, make_rotate4z};
  char Sf3_name[3][32] = {"identity", "mirrorxy", "rotate4z"};

  for (int iS = 0; iS < 3; ++iS)
    for (int splitting = 0; splitting < 5; splitting += 3)
      for (int igv = 0; igv < 4; ++igv) {
        for (int ic = 0; ic < 1; ++ic)
          if (broadcast(0, rand()) % chances == 0) {
            bool use_real = true;
            char name[1024];
            snprintf(name, 1024, "check_3d_ezsrc_%s_%d_%s_%s%s", Sf3_name[iS], splitting,
                     gv_3d_name[igv], component_name(c3d[ic]), use_real ? "_r" : "");
            master_printf("Checking %s...\n", name);
            if (!check_3d(funky_eps_3d, a, splitting, Sf3[iS], Ez, c3d[ic], gv_3d[igv], use_real,
                          gv_3d_rank[igv], name, temp_dir.get()))
              return 1;
          }
      }
#endif /* HAVE_HDF5 */

  delete_directory(temp_dir.get());

  return 0;
}
