/***************************************************************/
/***************************************************************/
/***************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <complex>
#include <vector>

#include "meep.hpp"

#include "ctl-math.h"
#include "ctlgeom.h"

#include "meepgeom.hpp"

using namespace meep;
using namespace std;

vector3 v3(double x=0.0, double y=0.0, double z=0.0)
{
  vector3 v;
  v.x=x; v.y=y; v.z=z;
  return v;
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void test_array_metadata(meep::fields &f, const volume &where,
                         bool collapse_empty_dimensions)
{
  /***************************************************************/
  /* step 1: get coordinate grids and weights as reported by     */
  /* get_array_metadata                                          */
  /***************************************************************/
  size_t dims[3];
  direction dirs[3];
  int rank=f.get_array_slice_dimensions(where, dims, dirs,
                                        collapse_empty_dimensions);
  printf("Rank=%i, dims=");
  for(int r=0; r<rank; r++)
   printf("%c%zu ",r==0 ? '{' : ',',dims[r]);
  printf("}\n");
  size_t nxyz[3]={0,0,0}, nw=1;
  for(int r=0; r<rank; r++)
   { nxyx[dirs[r]-X] = dims[r];
     nw*=dims[r];
   }

  vector<double> xgrid(nxyx[0],0.0);
  vector<double> ygrid(nxyx[1],0.0);
  vector<double> zgrid(nxyx[2],0.0);
  vector<double> weights(nw,0.0);
  f.get_array_metadata(where, xgrid, ygrid, zgrid, weights, collapse_empty_dimensions);
  size_t stride[3];
  stride[2] = 1;
  stride[1] = ( zgrid.size() > 0 ? zgrid.size() : 1);
  stride[0] = ( ygrid.size() > 0 ? ygrid.size() : 1) * stride[1];

  /***************************************************************/
  /* step 2: initialize loop over grid points in the volume via  */
  /*         standard libmeep looping primitives                 */
  /***************************************************************/
  component cgrid=Centered;
  vec yee_c(gv.yee_shift(Centered) - gv.yee_shift(cgrid));
  ivec iyee_c(gv.iyee_shift(Centered) - gv.iyee_shift(cgrid));
  volume wherec(where + yee_c);
  ivec is(vec2diel_floor(wherec.get_min_corner(), gv.a, zero_ivec(gv.dim)));
  ivec ie(vec2diel_ceil(wherec.get_max_corner(), gv.a, zero_ivec(gv.dim)));
  
  ivec imin=gv.little_corner()+one_ivec(gv.dim), imax=gv.big_corner()-one_ivec(gv.dim);
  LOOP_OVER_DIRECTIONS(gv.dim, d)
   { if (is.in_direction(d) < imin.in_direction(d))
      is.set_direction(d,imin.in_direction(d));
     if (ie.in_direction(d) > imax.in_direction(d))
      ie.set_direction(d,imax.in_direction(d));
   }

  bool snap_empty_dims=true;
  vec s0(gv.dim), e0(gv.dim), s1(gv.dim), e1(gv.dim);
  // this initialization step seems to be necessary here to avoid winding
  // up with zero or undefined integration weights; I don't know why it
  // seems to be unnecessary for loop_in_chunks above.
  FOR_DIRECTIONS(d)
   if (!has_direction(gv.dim,d))
    { s0.set_direction(d,1.0);
      e0.set_direction(d,1.0);
      s1.set_direction(d,1.0);
      e1.set_direction(d,1.0);
    }
  compute_boundary_weights(gv, wherec, is, ie, snap_empty_dims,
                           s0, e0, s1, e1);

  // Determine integration "volumes" dV0 and dV1
  double dV0 = 1.0, dV1 = 0.0;
  LOOP_OVER_DIRECTIONS(gv.dim, d)
   if (wherec.in_direction(d) > 0.0)
    dV0 *= gv.inva;

  /***************************************************************/
  /* step 3: execute the loop and check that coordinates and     */
  /*         weights of each point as determined from the return */
  /*         values of get_array_metadata agree with those       */
  /*         determined by the libmeep loop primitives           */
  /***************************************************************/
  int np=0;
  FILE *f=fopen("test-array-metadata.out","w");
  LOOP_OVER_IVECS(gv, is, ie, idx)
   {
     // correct coordinates and weight for current grid point
     double xyzw_loop[4]={0.0,0.0,0.0,0.0};
     IVEC_LOOP_LOC(gv, loc);
     xyzw_loop[0] = has_direction(gv.dim,X) ? loc.x() : 0.0;
     xyzw_loop[1] = has_direction(gv.dim,Y) ? loc.y() : 0.0;
     xyzw_loop[2] = has_direction(gv.dim,Z) ? loc.z() : 0.0;
     xyzw_loop[3] = IVEC_LOOP_WEIGHT(s0, s1, e0, e1, dV0 + dV1 * loop_i2);

     // coordinates and weight for current grid point according to metadata
     double xyzw_metadata[4]={0.0,0.0,0.0,0.0};
     IVEC_LOOP_ILOC(gv, iloc);
     ivec twon = iloc - is;
     int nx=0, ny=0, nz=0, index=0;
     if (has_direction(gv.dim,X))
      { nx=twon.in_direction(X)/2;
        xyzw_metadata[0]=xgrid[nx];
        index += nx*stride[0];
      }
     if (has_direction(gv.dim,Y))
      { ny=twon.in_direction(Y)/2;
        xyzw_metadata[1]=ygrid[ny];
        index += ny*stride[1];
      }
     if (has_direction(gv.dim,Z))
      { nz=twon.in_direction(Z)/2;
        xyzw_metadata[2]=zgrid[nz];
        index += nz*stride[0];
      }
     xyzw_metadata[3]=weights[index];

     fprintf(f,"%i ",np);
     fprintf(f,"%e %e %e %e ", xyzw_loop[0], xyzw_loop[1], xyzw_loop[2], xyzw_loop[3]);
     fprintf(f,"%e %e %e %e ", xyzw_meta[0], xyzw_meta[1], xyzw_meta[2], xyzw_meta[3]);
     fprintf(f,"\n");

   } // LOOP_OVER_IVECS(gv, is, ie, idx)
  fclose(f);
}


/***************************************************************/
/* dummy material function needed to pass to structure( )      */
/* constructor as a placeholder before we can call             */
/* set_materials_from_geometry                                 */
/***************************************************************/
double dummy_eps(const vec &) { return 1.0; }

/***************************************************************/
/***************************************************************/
/***************************************************************/
int main(int argc, char *argv[])
{
  initialize mpi(argc, argv);

  /*--------------------------------------------------------------*/
  /*- set default geometric parameters ---------------------------*/
  /*--------------------------------------------------------------*/
  // size of computational cell
  double sx = 10.0;
  double sy =  5.0;
  double sz =  0.0;

  // corners of array volume
  double vxmin = -2.5,   vxmax = -2.5;
  double vymin = -1.0,   vymax = +3.0;
  double vzmin =  0.0,   vzmax =  0.0;

  double res=10.0;

  bool collapse_empty_dimensions=false;

  // double-valued command-line parameters
  vector <const char *> parm_name;
  vector <double *> parm_adrs;
  parm_name.push_back("--vxmin"); parm_adrs.push_back(&vxmin);
  parm_name.push_back("--vymin"); parm_adrs.push_back(&vymin);
  parm_name.push_back("--vzmin"); parm_adrs.push_back(&vzmin);
  parm_name.push_back("--vxmax"); parm_adrs.push_back(&vxmax);
  parm_name.push_back("--vymax"); parm_adrs.push_back(&vymax);
  parm_name.push_back("--vzmax"); parm_adrs.push_back(&vzmax);
  parm_name.push_back("--res");   parm_adrs.push_back(&res);
  parm_name.push_back("--sz");    parm_adrs.push_back(&sz);

  /*--------------------------------------------------------------*/
  /*- parse arguments --------------------------------------------*/
  /*--------------------------------------------------------------*/
  for(int narg=1; narg<argc; narg++)
   {
     // process boolean-valued parameters
     if (!strcasecmp(argv[narg],"--collapse"))
      { collapse_empty_dimensions=true;
        printf("Collapsing empty array dimensions.\n");
        continue;
      }

     // process double-valued parameters
     size_t np;
     for(np=0; np<parm_name.size(); np++)
      if (!strcasecmp(argv[narg], parm_name[np]))
       break;
     if (np==parm_name.size())
      meep::abort("unknown command-line option %s",argv[narg]);
     if (narg+1 == argc)
      meep::abort("no option specified for %s",argv[narg]);
     if (1!=sscanf(argv[narg+1],"%le",parm_adrs[np]))
      meep::abort("invalid value %s specified for %s",argv[narg+1],argv[narg]);
     printf("Setting %s=%e.\n",argv[narg],*(parm_adrs[np]));
     narg++;
   }

  /*--------------------------------------------------------------*/
  /*- initialize geometry ----------------------------------------*/
  /*--------------------------------------------------------------*/
  geometry_lattice.size.x=sx;
  geometry_lattice.size.y=sy;
  geometry_lattice.size.z=sz;
  grid_volume gv;
  if (sx==0.0 && sy==0.0) 
   gv=vol1d(sz,res);
  else if (sz==0.0)
   gv=vol2d(sx,sy,res);
  else
   gv=vol3d(sx,sy,sz,res);
  gv.center_origin();
  structure the_structure(gv, dummy_eps, pml(1.0));

  meep_geom::material_type silicon = meep_geom::make_dielectric(12.0);
  geometric_object objects[1];
  vector3 origin   = v3(0.0,    0.0,    0.0);
  vector3 wvg_size = v3(0.5*sx, 0.5*sy, 0.5*sz);
  vector3 xhat   = {1.0, 0.0, 0.0};
  vector3 yhat   = {0.0, 1.0, 0.0};
  vector3 zhat   = {0.0, 0.0, 1.0};
  objects[0] = make_block(silicon, center, xhat, yhat, zhat, wvg_size);
  geometric_object_list g={ 1, objects };
  meep_geom::set_materials_from_geometry(&the_structure, g);

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  fields f(&the_structure);
  double fcen = 0.2;
  double df   = 0.1;
  gaussian_src_time src(fcen, df);
  vec src_point = origin;
  vec src_size  = v3(0.0, 0.0, 0.0);
  f.add_point_source(src_cmpt, src, src_point);
  f.init_fields();

  volume slice( v3(vxmin, vymin, vzmin), v3(vxmax, vymax, vzmax);
  test_array_metadata(f, slice, collapse_empty_dimensions);

}
