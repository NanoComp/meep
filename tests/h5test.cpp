#include <math.h>
#include <stdio.h>

#include "meep.h"
#include "meep_internals.h"
#include "h5io.h"
using namespace meep;

const double xsize = 2.0;
const double ysize = 2.0;
const double zsize = 0.6;

const double r = 0.5;
const double eps_k = 2*pi / 1.0;

double funky_eps_2d(const vec &p_) {
  vec p = p_ - vec2d(xsize / 2, ysize / 2);
  if (fabs(p & p) < r * r)
    return 1.0;
  return 2.0 + cos(p.x() * eps_k) * cos(p.y() * eps_k);
}

double funky_eps_3d(const vec &p_) {
  vec p = p_ - vec(xsize / 2, ysize / 2, zsize / 2);
  if (fabs(p & p) < r * r)
    return 1.0;
  return 2.0 + cos(p.x() * eps_k) * cos(p.y() * eps_k) * cos(p.z() * eps_k);
}

symmetry make_identity(const volume &v)
{
  (void) v; // unused
  return identity();
}

symmetry make_mirrorx(const volume &v)
{
  return mirror(X, v);
}

symmetry make_mirrory(const volume &v)
{
  return mirror(Y, v);
}

symmetry make_mirrorxy(const volume &v)
{
  return mirror(X, v) + mirror(Y, v);
}

symmetry make_rotate4z(const volume &v)
{
  return rotate4(Z, v);
}

typedef symmetry (*symfunc)(const volume &);

const double tol = 1e-8;
double compare(double a, double b, const char *nam) {
  if (fabs(a-b) > fabs(b)*tol || b != b) {
    master_printf("%g vs. %g differs by\t%g\n", a, b, fabs(a-b));
    master_printf("This gives a fractional error of %g\n", fabs(a-b)/fabs(b));
    abort("Error in %s\n", nam);
  }
  return fabs(a-b);
}

double get_reim(complex<double> x, int reim)
{
  return reim ? imag(x) : real(x);
}

bool check_2d(double eps(const vec &), double a, int splitting, symfunc Sf,
	      component src_c, component file_c,
	      geometric_volume file_gv,
	      bool real_fields, int expected_rank,
	      const char *name) {
  const volume v = vol2d(xsize, ysize, a);
  structure s(v, eps, splitting, Sf(v));
  fields f(&s);

  f.add_point_source(src_c, 0.3, 2.0, 0.0, 1.0, v.center(), 1.0, 1);
  if (real_fields) f.use_real_fields();

  while (f.time() <= 3.0 && !interrupt)
    f.step();

  char fname[1024];
  snprintf(fname, 1024, "%s.h5", name);
  double res = 1.54321 * a;
  f.output_hdf5(fname, file_gv, res, file_c);

  // compute corner coordinate of file data
  double resinv = 1.0 / res;
  vec loc0(file_gv.get_min_corner());
  LOOP_OVER_DIRECTIONS(loc0.dim, d) {
     int minpt = int(ceil(file_gv.in_direction_min(d) * res));
     int maxpt = int(floor(file_gv.in_direction_max(d) * res));
     if (minpt < maxpt)
	  loc0.set_direction(d, minpt * resinv);
  }

  double data_min = infinity, data_max = -infinity;
  double err_max = 0;
  for (int reim = 0; reim < (real_fields ? 1 : 2); ++reim) {
    int rank, dims[2] = {1, 1};

    char dataname[256];
    snprintf(dataname, 256, "%s.%s", component_name(file_c), reim ? "i" : "r");

    double *h5data = h5io::read(fname, dataname,
				&rank, dims, 2);
    if (!h5data)
	 abort("failed to read dataset %s:%s\n", fname, dataname);
    if (rank != expected_rank)
	 abort("incorrect rank (%d instead of %d) in %s:%s\n",
	       rank, expected_rank, fname, dataname);
    if (expected_rank == 1 && 
	file_gv.in_direction_min(X) == file_gv.in_direction_max(X)) {
      dims[1] = dims[0];
      dims[0] = 1;
    }
    vec loc(loc0.dim);
    for (int i0 = 0; i0 < dims[0]; ++i0) {
      for (int i1 = 0; i1 < dims[1]; ++i1) {
	loc.set_direction(X, loc0.in_direction(X) + i0 * resinv);
	loc.set_direction(Y, loc0.in_direction(Y) + i1 * resinv);
	int idx = i0 * dims[1] + i1;

	/* Ugh, for rotational symmetries (which mix up components etc.),
	   we can't guarantee that a component is *exactly* the
	   same as its rotated version, and we don't know which one
	   was written to the file. */
	component cs = file_c;
	complex<double> ph = 1.0;
	double diff = fabs(get_reim(f.get_field(file_c, loc), reim) -
			   h5data[idx]);
	for (int sn = 1; sn < f.S.multiplicity(); ++sn) {
	  vec loc2(f.S.transform(loc, sn));
	  component cs2 = f.S.transform(file_c, sn);
	  complex<double> ph2 = f.S.phase_shift(cs2, -sn);
	  double diff2 = fabs(get_reim(f.get_field(cs2, loc2)*ph2, reim) -
			      h5data[idx]);
	  if (diff2 < diff) {
	    loc = loc2;
	    cs = cs2;
	    ph = ph2;
	    diff = diff2;
	  }
	}

	double err = compare(h5data[idx],
			     get_reim(f.get_field(cs, loc) * ph, reim),
			     name);
	err_max = max(err, err_max);
	data_min = min(data_min, h5data[idx]);
	data_max = max(data_max, h5data[idx]);
      }
    }
    delete[] h5data;
  }

  remove(fname);

  master_printf("Passed %s (%g..%g), err=%g\n", name,
		data_min, data_max,
		err_max / max(fabs(data_min), fabs(data_max)));

  return 1;
}

bool check_3d(double eps(const vec &), double a, int splitting, symfunc Sf,
	      component src_c, component file_c,
	      geometric_volume file_gv,
	      bool real_fields, int expected_rank,
	      const char *name) {
  const volume v = vol3d(xsize, ysize, zsize, a);
  structure s(v, eps, splitting, Sf(v));
  fields f(&s);

  f.add_point_source(src_c, 0.3, 2.0, 0.0, 1.0, v.center(), 1.0, 1);
  if (real_fields) f.use_real_fields();

  while (f.time() <= 3.0 && !interrupt)
    f.step();

  char fname[1024];
  snprintf(fname, 1024, "%s.h5", name);
  double res = 0.54321 * a;
  f.output_hdf5(fname, file_gv, res, file_c);

  // compute corner coordinate of file data
  double resinv = 1.0 / res;
  vec loc0(file_gv.get_min_corner());
  LOOP_OVER_DIRECTIONS(loc0.dim, d) {
     int minpt = int(ceil(file_gv.in_direction_min(d) * res));
     int maxpt = int(floor(file_gv.in_direction_max(d) * res));
     if (minpt < maxpt)
	  loc0.set_direction(d, minpt * resinv);
  }

  double data_min = infinity, data_max = -infinity;
  double err_max = 0;
  for (int reim = 0; reim < (real_fields ? 1 : 2); ++reim) {
    int rank, dims[3] = {1, 1, 1};

    char dataname[256];
    snprintf(dataname, 256, "%s.%s", component_name(file_c), reim ? "i" : "r");

    double *h5data = h5io::read(fname, dataname,
				&rank, dims, 3);
    if (!h5data)
	 abort("failed to read dataset %s:%s\n", fname, dataname);
    if (rank != expected_rank)
	 abort("incorrect rank (%d instead of %d) in %s:%s\n",
	       rank, expected_rank, fname, dataname);
    vec loc(loc0.dim);
    for (int i0 = 0; i0 < dims[0]; ++i0) {
      for (int i1 = 0; i1 < dims[1]; ++i1) {
	for (int i2 = 0; i2 < dims[2]; ++i2) {
	  loc.set_direction(X, loc0.in_direction(X) + i0 * resinv);
	  loc.set_direction(Y, loc0.in_direction(Y) + i1 * resinv);
	  loc.set_direction(Z, loc0.in_direction(Z) + i2 * resinv);
	  int idx = (i0 * dims[1] + i1) * dims[2] + i2;
	  
	  /* Ugh, for rotational symmetries (which mix up components etc.),
	     we can't guarantee that a component is *exactly* the
	     same as its rotated version, and we don't know which one
	     was written to the file. */
	  component cs = file_c;
	  complex<double> ph = 1.0;
	  double diff = fabs(get_reim(f.get_field(file_c, loc), reim) -
			     h5data[idx]);
	  for (int sn = 1; sn < f.S.multiplicity(); ++sn) {
	    vec loc2(f.S.transform(loc, sn));
	    component cs2 = f.S.transform(file_c, sn);
	    complex<double> ph2 = f.S.phase_shift(cs2, -sn);
	    double diff2 = fabs(get_reim(f.get_field(cs2, loc2)*ph2, reim) -
				h5data[idx]);
	    if (diff2 < diff) {
	      loc = loc2;
	      cs = cs2;
	      ph = ph2;
	      diff = diff2;
	    }
	  }
	  
	  double err = compare(h5data[idx],
			       get_reim(f.get_field(cs, loc)*ph,reim),
			       name);
	  err_max = max(err, err_max);
	  data_min = min(data_min, h5data[idx]);
	  data_max = max(data_max, h5data[idx]);
	}
      }
    }
    delete[] h5data;
  }

  remove(fname);

  master_printf("Passed %s (%g..%g), err=%g\n", name,
		data_min, data_max,
		err_max / (max(fabs(data_min), fabs(data_max)) + 1e-16));

  return 1;
}

bool check_2d_monitor(double eps(const vec &),
		      double a, int splitting, symfunc Sf,
		      component src_c, component file_c,
		      const vec &pt,
		      bool real_fields,
		      const char *name) {
  const volume v = vol2d(xsize, ysize, a);
  structure s(v, eps, splitting, Sf(v));
  fields f(&s);

  f.add_point_source(src_c, 0.3, 2.0, 0.0, 1.0, v.center(), 1.0, 1);
  if (real_fields) f.use_real_fields();

  char fname[1024];
  snprintf(fname, 1024, "%s.h5", name);

  const double T = 3.0;
  int NT = int(T / (f.inva * c)) + 2;
  complex<double> *mon = new complex<double>[NT];
  while (f.time() <= T && !interrupt) {
    f.output_hdf5(fname, geometric_volume(pt, pt), a, file_c,
		  true, f.t);
    mon[f.t] = f.get_field(file_c, pt);
    f.step();
  }

  double data_min = infinity, data_max = -infinity;
  double err_max = 0;

  for (int reim = 0; reim < (real_fields ? 1 : 2); ++reim) {
    int rank, dims[1] = {1};

    char dataname[256];
    snprintf(dataname, 256, "%s.%s", component_name(file_c), reim ? "i" : "r");

    double *h5data = h5io::read(fname, dataname, &rank, dims, 2);
    if (!h5data)
	 abort("failed to read dataset %s:%s\n", fname, dataname);
    if (rank != 1)
      abort("monitor-point data is not one-dimensional");
    if (dims[0] != f.t)
      abort("incorrect size of monitor-point data");

    for (int i = 0; i < f.t; ++i) {
      double err = compare(h5data[i], get_reim(mon[i], reim), name);
      err_max = max(err, err_max);
      data_min = min(data_min, h5data[i]);
      data_max = max(data_max, h5data[i]);
    }
    delete[] h5data;
  }

  delete[] mon;

  remove(fname);

  master_printf("Passed %s (%g..%g), err=%g\n", name,
		data_min, data_max,
		err_max / max(fabs(data_min), fabs(data_max)));

  return 1;
}

int main(int argc, char **argv)
{
  const double a = 10.0;
  initialize mpi(argc, argv);
  const double pad1 = 0.3, pad2 = 0.2, pad3 = 0.1;

  geometric_volume gv_2d[3] = {
       geometric_volume(vec2d(pad1,pad2), vec2d(xsize-pad2,ysize-pad1)),
       geometric_volume(vec2d(pad1,pad2), vec2d(xsize-pad2,pad2)),
       geometric_volume(vec2d(pad1,pad2), vec2d(pad1,pad2)),
  };
  char gv_2d_name[3][10] = {"plane", "line", "point"};
  int gv_2d_rank[3] = {2,1,0};
  component tm_c[4] = {Ez, Dz, Hx, Hy};
  symfunc Sf2[5] = {make_identity, make_mirrorx, make_mirrory, make_mirrorxy,
		   make_rotate4z};
  char Sf2_name[5][32] = {"identity", "mirrorx", "mirrory", "mirrorxy",
			 "rotate4z"};

  for (int iS = 0; iS < 5; ++iS)
    for (int splitting = 0; splitting < 5; ++splitting)
      for (int igv = 0; igv < 3; ++igv)
	for (int ic = 0; ic < 4; ++ic)
	  for (int use_real = 0; use_real <= 1; ++use_real) {
	    char name[1024];
	    snprintf(name, 1024, "check_2d_tm_%s_%d_%s_%s%s",
		     Sf2_name[iS], splitting, gv_2d_name[igv],
		     component_name(tm_c[ic]), use_real ? "_r" : "");
	    master_printf("Checking %s...\n", name);
	    if (!check_2d(funky_eps_2d, a, splitting, Sf2[iS], Ez, tm_c[ic],
			  gv_2d[igv], use_real, gv_2d_rank[igv], name))
	      return 1;
	  }
  
  for (int iS = 0; iS < 5; ++iS)
    for (int splitting = 0; splitting < 5; ++splitting)
      for (int ic = 0; ic < 4; ++ic)
	for (int use_real = 0; use_real <= 1; ++use_real) {
	  char name[1024];
	  snprintf(name, 1024, "check_2d_monitor_tm_%s_%d_%s%s",
		   Sf2_name[iS], splitting,
		   component_name(tm_c[ic]), use_real ? "_r" : "");
	  master_printf("Checking %s...\n", name);
	  if (!check_2d_monitor(funky_eps_2d, a, splitting, Sf2[iS], Ez, 
				tm_c[ic], vec2d(pad1,pad2), use_real, name))
	    return 1;
	}
  
  geometric_volume gv_3d[4] = {
       geometric_volume(vec(pad1,pad2,pad3), vec(xsize-pad2,ysize-pad1,zsize-pad3)),
       geometric_volume(vec(pad1,pad2,pad3), vec(xsize-pad2,ysize-pad1,pad3)),
       geometric_volume(vec(pad1,pad2,pad3), vec(xsize-pad2,pad2,pad3)),
       geometric_volume(vec(pad1,pad2,pad3), vec(pad1,pad2,pad3)),
  };
  char gv_3d_name[4][10] = {"volume", "plane", "line", "point"};
  int gv_3d_rank[4] = {3,2,1,0};
  component c3d[6] = {Ex,Ey,Ez, Hx,Hy,Hz};
  symfunc Sf3[3] = {make_identity, make_mirrorxy, make_rotate4z};
  char Sf3_name[3][32] = {"identity", "mirrorxy", "rotate4z"};

  for (int iS = 0; iS < 3; ++iS)
    for (int splitting = 0; splitting < 5; splitting += 3)
      for (int igv = 0; igv < 4; ++igv) {
	int ic = 0;
	bool use_real = true;
	char name[1024];
	snprintf(name, 1024, "check_3d_ezsrc_%s_%d_%s_%s%s", Sf3_name[iS],
		 splitting, gv_3d_name[igv], component_name(c3d[ic]),
		 use_real ? "_r" : "");
	master_printf("Checking %s...\n", name);
	if (!check_3d(funky_eps_3d, a, splitting, Sf3[iS], Ez, c3d[ic],
		      gv_3d[igv], use_real, gv_3d_rank[igv], name))
	  return 1;
      }

  return 0;
}

