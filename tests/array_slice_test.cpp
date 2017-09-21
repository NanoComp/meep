#include <math.h>
#include <unistd.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include <meep.hpp>
#include "meep_internals.hpp"
#include "config.h"

using namespace meep;
using namespace std;

typedef complex<double> cdouble;

const double xsize = 2.0;
const double ysize = 2.0;
const double zsize = 0.6;

const double r = 0.5;
const double eps_k = 2*pi / 1.0;

double funky_eps_2d(const vec &p_) {
  vec p = p_ - vec(xsize / 2, ysize / 2);
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

symmetry make_identity(const grid_volume &gv)
{
  (void) gv; // unused
  return identity();
}

symmetry make_mirrorx(const grid_volume &gv)
{
  return mirror(X, gv);
}

symmetry make_mirrory(const grid_volume &gv)
{
  return mirror(Y, gv);
}

symmetry make_mirrorxy(const grid_volume &gv)
{
  return mirror(X, gv) + mirror(Y, gv);
}

symmetry make_rotate4z(const grid_volume &gv)
{
  return rotate4(Z, gv);
}

typedef symmetry (*symfunc)(const grid_volume &);

const double tol = sizeof(realnum) == sizeof(float) ? 1e-4 : 1e-8;
double compare(double a, double b, const char *nam, int i0,int i1,int i2) {
  if (fabs(a-b) > tol*tol + fabs(b) * tol || b != b) {
    master_printf("%g vs. %g differs by\t%g\n", a, b, fabs(a-b));
    master_printf("This gives a fractional error of %g\n", fabs(a-b)/fabs(b));
    abort("Error in %s at (%d,%d,%d)\n", nam, i0,i1,i2);
  }
  return fabs(a-b);
}

#if 0
bool check_2d(double eps(const vec &), double a, int splitting, symfunc Sf,
	      double kx, double ky,
	      component src_c, int file_c,
	      volume file_gv,
	      bool real_fields, int expected_rank,
	      const char *name) {
  const grid_volume gv = vol2d(xsize, ysize, a);
  structure s(gv, eps, no_pml(), Sf(gv), splitting);
  fields f(&s);

  f.use_bloch(X, real_fields ? 0.0 : kx);
  f.use_bloch(Y, real_fields ? 0.0 : ky);

  if (real_fields) f.use_real_fields();
  f.add_point_source(src_c, 0.3, 2.0, 0.0, 1.0, gv.center(), 1.0, 1);

  if (file_c >= int(Dielectric)) real_fields = true;

  while (f.time() <= 3.0 && !interrupt)
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

  char *str = file->read("stringtest");
  if (strcmp(str, "Hello, world!\n"))
       abort("Failed to read back string test from %s...", name);

  // compute corner coordinate of file data
  vec loc0(file_gv.get_min_corner());
  ivec iloc0(gv.dim);
  LOOP_OVER_DIRECTIONS(gv.dim, d) {
    iloc0.set_direction(d, 1+2*int(floor(loc0.in_direction(d)*a-.5)));
    if (file_gv.in_direction(d) == 0.0 &&
	1. - file_gv.in_direction_min(d)*a + 0.5*iloc0.in_direction(d)
	<= 1. + file_gv.in_direction_max(d)*a - 0.5*(iloc0.in_direction(d)+2))
      iloc0.set_direction(d, iloc0.in_direction(d) + 2); // snap to grid
  }
  loc0 = gv[iloc0];

  double data_min = meep::infinity, data_max = -meep::infinity;
  double err_max = 0;
  for (int reim = 0; reim < (real_fields ? 1 : 2); ++reim) {
    int rank, dims[2] = {1, 1};

    char dataname[256];
    snprintf(dataname, 256, "%s%s", component_name(file_c),
	     reim ? ".i" : (real_fields ? "" : ".r"));

    realnum *h5data = file->read(dataname, &rank, dims, 2);
    file->prevent_deadlock(); // hackery
    if (!h5data)
	 abort("failed to read dataset %s:%s\n", name, dataname);
    if (rank != expected_rank)
	 abort("incorrect rank (%d instead of %d) in %s:%s\n",
	       rank, expected_rank, name, dataname);
    if (expected_rank == 1 && 
	file_gv.in_direction_min(X) == file_gv.in_direction_max(X)) {
      dims[1] = dims[0];
      dims[0] = 1;
    }
    vec loc(loc0.dim);
    for (int i0 = 0; i0 < dims[0]; ++i0) {
      for (int i1 = 0; i1 < dims[1]; ++i1) {
	loc.set_direction(X, loc0.in_direction(X) + i0 * gv.inva);
	loc.set_direction(Y, loc0.in_direction(Y) + i1 * gv.inva);
	int idx = i0 * dims[1] + i1;

	/* Ugh, for rotational symmetries (which mix up components etc.),
	   we can't guarantee that a component is *exactly* the
	   same as its rotated version, and we don't know which one
	   was written to the file. */
	int cs = file_c;
	complex<double> ph = 1.0;
	double diff = fabs(get_reim(f.get_field(file_c, loc), reim) -
			   h5data[idx]);
	for (int sn = 1; sn < f.S.multiplicity(); ++sn) {
	  vec loc2(f.S.transform(loc, sn));
	  int cs2 = f.S.transform(file_c, sn);
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
			     name, i0,i1,0);
	err_max = max(err, err_max);
	data_min = min(data_min, h5data[idx]);
	data_max = max(data_max, h5data[idx]);
      }
    }
    delete[] h5data;
  }

  //file->remove();
  delete file;

  master_printf("Passed %s (%g..%g), err=%g\n", name,
		data_min, data_max,
		err_max / max(fabs(data_min), fabs(data_max)));

  return true;
}

bool check_3d(double eps(const vec &), double a, int splitting, symfunc Sf,
	      component src_c, int file_c,
	      volume file_gv,
	      bool real_fields, int expected_rank,
	      const char *name) {
  const grid_volume gv = vol3d(xsize, ysize, zsize, a);
  structure s(gv, eps, no_pml(), Sf(gv), splitting);
  fields f(&s);

  if (real_fields) f.use_real_fields();
  f.add_point_source(src_c, 0.3, 2.0, 0.0, 1.0, gv.center(), 1.0, 1);

  if (file_c >= Dielectric) real_fields = true;

  while (f.time() <= 3.0 && !interrupt)
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

  char *str = file->read("stringtest");
  if (strcmp(str, "Hello, world!\n"))
       abort("Failed to read back string test from %s...", name);

  // compute corner coordinate of file data
  vec loc0(file_gv.get_min_corner());
  ivec iloc0(gv.dim);
  LOOP_OVER_DIRECTIONS(gv.dim, d) {
    iloc0.set_direction(d, 1+2*int(floor(loc0.in_direction(d)*a-.5)));
    if (file_gv.in_direction(d) == 0.0 &&
	1. - file_gv.in_direction_min(d)*a + 0.5*iloc0.in_direction(d)
	<= 1. + file_gv.in_direction_max(d)*a - 0.5*(iloc0.in_direction(d)+2))
      iloc0.set_direction(d, iloc0.in_direction(d) + 2); // snap to grid
  }
  loc0 = gv[iloc0];

  double data_min = meep::infinity, data_max = -meep::infinity;
  double err_max = 0;
  for (int reim = 0; reim < (real_fields ? 1 : 2); ++reim) {
    int rank, dims[3] = {1, 1, 1};

    char dataname[256];
    snprintf(dataname, 256, "%s%s", component_name(file_c),
	     reim ? ".i" : (real_fields ? "" : ".r"));

    realnum *h5data = file->read(dataname, &rank, dims, 3);
    file->prevent_deadlock(); // hackery
    if (!h5data)
	 abort("failed to read dataset %s:%s\n", name, dataname);
    if (rank != expected_rank)
	 abort("incorrect rank (%d instead of %d) in %s:%s\n",
	       rank, expected_rank, name, dataname);
    vec loc(loc0.dim);
    for (int i0 = 0; i0 < dims[0]; ++i0) {
      for (int i1 = 0; i1 < dims[1]; ++i1) {
	for (int i2 = 0; i2 < dims[2]; ++i2) {
	  loc.set_direction(X, loc0.in_direction(X) + i0 * gv.inva);
	  loc.set_direction(Y, loc0.in_direction(Y) + i1 * gv.inva);
	  loc.set_direction(Z, loc0.in_direction(Z) + i2 * gv.inva);
	  int idx = (i0 * dims[1] + i1) * dims[2] + i2;
	  
	  /* Ugh, for rotational symmetries (which mix up components etc.),
	     we can't guarantee that a component is *exactly* the
	     same as its rotated version, and we don't know which one
	     was written to the file. */
	  int cs = file_c;
	  complex<double> ph = 1.0;
	  double diff = fabs(get_reim(f.get_field(file_c, loc), reim) -
			     h5data[idx]);
	  for (int sn = 1; sn < f.S.multiplicity(); ++sn) {
	    vec loc2(f.S.transform(loc, sn));
	    int cs2 = f.S.transform(file_c, sn);
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
			       name, i0,i1,i2);
	  err_max = max(err, err_max);
	  data_min = min(data_min, h5data[idx]);
	  data_max = max(data_max, h5data[idx]);
	}
      }
    }
    delete[] h5data;
  }

  //file->remove();
  delete file;

  master_printf("Passed %s (%g..%g), err=%g\n", name,
		data_min, data_max,
		err_max / (max(fabs(data_min), fabs(data_max)) + 1e-16));

  return 1;
}
#endif

/***************************************************************/
/***************************************************************/
/***************************************************************/
int main(int argc, char **argv)
{
  initialize mpi(argc, argv);

  const double a = 10.0; // resolution

  const double pad1 = 0.314159, pad2 = 0.27183, pad3 = 0.14142;

  volume gv_1d( vec(pad1),
                vec(xsize-pad1)
              );

  volume gv_2d( vec(pad1,       pad2),
                vec(xsize-pad1, ysize-pad2)
              );

  volume gv_3d( vec(pad1,      pad2,      pad3),
                vec(xsize-pad1,ysize-pad2,zsize-pad3)
              );

  volume gvs[3] = {gv_1d, gv_2d, gv_3d};

  /***************************************************************/
  /* 1D checks ***************************************************/
  /***************************************************************/
  // TODO

  /***************************************************************/
  /* 2D checks ***************************************************/
  /***************************************************************/
  int tm_c[5] = {Dielectric, Ez, Hy, Sx, D_EnergyDensity};
  symfunc Sf2[5] = {make_identity, make_mirrorx, make_mirrory, make_mirrorxy, make_rotate4z};
  char Sf2_name[5][32] = {"identity", "mirrorx", "mirrory", "mirrorxy", "rotate4z"};
  double Sf2_kx[5] = {0.3, 0, 0.3, 0, 0};
  double Sf2_ky[5] = {0.2, 0.2, 0, 0, 0};

  const grid_volume gv = vol2d(xsize, ysize, a);
  for(int iS=0; iS<1; /*3;*/ iS++)
   for(int splitting=0; splitting<1; /*5;*/ splitting++)
    { 
      symfunc Sf=Sf2[iS];
      structure s(gv, funky_eps_2d, no_pml(), Sf(gv), splitting);
      fields f(&s);
      f.use_bloch(X, 0.3);
      f.use_bloch(Y, 0.2);
      f.add_point_source(Ez, 0.3, 2.0, 0.0, 1.0, gv.center(), 1.0, 1);
      while (f.time() <= 3.0 && !interrupt)
       f.step();

      int dims[3], dirs[3];
      int rank=f.get_array_slice_dimensions(gv_2d, dims, dirs);
      if (rank!=2)
       abort("incorrect rank in get_array_slice_dimensions");
      printf("(iS,s)=(%i,%i): dims={%i,%i}\n",iS,splitting,dims[0],dims[1]);
      
      std::vector<component> components(Ez);
      components.push_back(Hz);
      cdouble *slice=f.get_array_slice(gv_2d, components);
      //delete[] slice;
      printf("Howdage foryaf! \n");
    };

  return 0;

}
