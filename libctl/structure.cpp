#include "meep-ctl.hpp"
#include <ctlgeom.h>

using namespace ctlio;

#define master_printf meep::master_printf
#define MTS material_type_struct

typedef struct {
  double m00, m01, m02,
              m11, m12,
                   m22;
} symmetric_matrix;

/* rotate A by a unitary (real) rotation matrix R:
      RAR = transpose(R) * A * R
*/
void sym_matrix_rotate(symmetric_matrix *RAR,
			       const symmetric_matrix *A_,
			       const double R[3][3])
{
     int i,j;
     double A[3][3], AR[3][3];
     A[0][0] = A_->m00;
     A[1][1] = A_->m11;
     A[2][2] = A_->m22;
     A[0][1] = A[1][0] = A_->m01;
     A[0][2] = A[2][0] = A_->m02;
     A[1][2] = A[2][1] = A_->m12;
     for (i = 0; i < 3; ++i) for (j = 0; j < 3; ++j) 
	  AR[i][j] = A[i][0]*R[0][j] + A[i][1]*R[1][j] + A[i][2]*R[2][j];
     for (i = 0; i < 3; ++i) for (j = i; j < 3; ++j) 
	  A[i][j] = R[0][i]*AR[0][j] + R[1][i]*AR[1][j] + R[2][i]*AR[2][j];
     RAR->m00 = A[0][0];
     RAR->m11 = A[1][1];
     RAR->m22 = A[2][2];
     RAR->m01 = A[0][1];
     RAR->m02 = A[0][2];
     RAR->m12 = A[1][2];
}

/* Set Vinv = inverse of V, where both V and Vinv are real-symmetric matrices.*/
void sym_matrix_invert(symmetric_matrix *Vinv, 
			       const symmetric_matrix *V)
{
     double m00 = V->m00, m11 = V->m11, m22 = V->m22;
     double m01 = V->m01, m02 = V->m02, m12 = V->m12;

     if (m01 == 0.0 && m02 == 0.0 && m12 == 0.0) {
	  /* optimize common case of a diagonal matrix: */
	  Vinv->m00 = 1.0 / m00;
	  Vinv->m11 = 1.0 / m11;
	  Vinv->m22 = 1.0 / m22;
	  Vinv->m01 = Vinv->m02 = Vinv->m12 = 0.0;
     }
     else {
	  double detinv;
	  
	  /* compute the determinant: */
	  detinv = m00*m11*m22 - m02*m11*m02 + 2.0 * m01*m12*m02 -
	       m01*m01*m22 - m12*m12*m00;
	  
	  if (detinv == 0.0) meep::abort( "singular 3x3 matrix");
	  
	  detinv = 1.0/detinv;
	  
	  Vinv->m00 = detinv * (m11*m22 - m12*m12);
	  Vinv->m11 = detinv * (m00*m22 - m02*m02);
	  Vinv->m22 = detinv * (m11*m00 - m01*m01);
	  
	  Vinv->m02 = detinv * (m01*m12 - m11*m02);
	  Vinv->m01 = detinv * (m12*m02 - m01*m22);
	  Vinv->m12 = detinv * (m01*m02 - m00*m12);
     }
}

/* Returns whether or not V is positive-definite. */
int sym_matrix_positive_definite(symmetric_matrix *V)
{
     double det2, det3;
     double m00 = V->m00, m11 = V->m11, m22 = V->m22;

#if defined(WITH_HERMITIAN_EPSILON)
     scalar_complex m01 = V->m01, m02 = V->m02, m12 = V->m12;

     det2 = m00*m11 - CSCALAR_NORMSQR(m01);
     det3 = det2*m22 - m11*CSCALAR_NORMSQR(m02) - CSCALAR_NORMSQR(m12)*m00 +
	  2.0 * ((m01.re * m12.re - m01.im * m12.im) * m02.re +
		 (m01.re * m12.im + m01.im * m12.re) * m02.im);
#else /* real matrix */
     double m01 = V->m01, m02 = V->m02, m12 = V->m12;

     det2 = m00*m11 - m01*m01;
     det3 = det2*m22 - m02*m11*m02 + 2.0 * m01*m12*m02 - m12*m12*m00;
#endif /* real matrix */
     
     return (m00 > 0.0 && det2 > 0.0 && det3 > 0.0);
}

static meep::ndim dim = meep::D3;

/***********************************************************************/

void set_dimensions(int dims)
{
  if (dims == CYLINDRICAL) {
    dimensions = 2;
    dim = meep::Dcyl;
  }
  else {
    dimensions = dims;
    dim = meep::ndim(dims - 1);
  }
}

vector3 vec_to_vector3(const meep::vec &v)
{
  vector3 v3;
  
  switch (v.dim) {
  case meep::D1:
    v3.x = 0;
    v3.y = 0;
    v3.z = v.z();
    break;
  case meep::D2:
    v3.x = v.x();
    v3.y = v.y();
    v3.z = 0;
    break;
  case meep::D3:
    v3.x = v.x();
    v3.y = v.y();
    v3.z = v.z();
    break;
  case meep::Dcyl:
    v3.x = v.r();
    v3.y = 0;
    v3.z = v.z();
    break;
  }
  return v3;
}

meep::vec vector3_to_vec(const vector3 v3)
{
  switch (dim) {
  case meep::D1:
    return meep::vec(v3.z);
  case meep::D2:
    return meep::vec(v3.x, v3.y);
  case meep::D3:
    return meep::vec(v3.x, v3.y, v3.z);
  case meep::Dcyl:
    return meep::veccyl(v3.x, v3.z);
  default:
    meep::abort("unknown dimensionality in vector3_to_vec");
  }
}

static geom_box gv2box(const meep::geometric_volume &gv)
{
  geom_box box;
  box.low = vec_to_vector3(gv.get_min_corner());
  box.high = vec_to_vector3(gv.get_max_corner());
  return box;
}

/***********************************************************************/

class geom_epsilon : public meep::material_function {
  geometric_object_list geometry;
  geom_box_tree geometry_tree;
  geom_box_tree restricted_tree;
  
public:
  geom_epsilon(geometric_object_list g, material_type_list mlist,
	       const meep::geometric_volume &gv);
  virtual ~geom_epsilon();
  
  virtual void set_volume(const meep::geometric_volume &gv);
  virtual void unset_volume(void);

  virtual bool has_chi3(meep::component c);
  virtual double chi3(meep::component c, const meep::vec &r);
  virtual bool has_chi2(meep::component c);
  virtual double chi2(meep::component c, const meep::vec &r);

  virtual bool has_mu();

  virtual bool has_conductivity(meep::component c);
  virtual double conductivity(meep::component c, const meep::vec &r);

  virtual double chi1p1(meep::field_type ft, const meep::vec &r);
  virtual void eff_chi1inv_row(meep::component c, double chi1inv_row[3],
			       const meep::geometric_volume &gv, 
			       double tol, int maxeval);

  void fallback_chi1inv_row(meep::component c, double chi1inv_row[3],
			      const meep::geometric_volume &gv,
			      double tol, int maxeval);

  virtual void sigma_row(meep::component c, double sigrow[3],
			 const meep::vec &r);
  void add_polarizabilities(meep::structure *s);
  void add_polarizabilities(meep::field_type ft, meep::structure *s);

private:
  bool get_material_pt(material_type &material, const meep::vec &r);

  material_type_list extra_materials;
};

geom_epsilon::geom_epsilon(geometric_object_list g, material_type_list mlist,
			   const meep::geometric_volume &gv)
{
  geometry = g; // don't bother making a copy, only used in one place
  extra_materials = mlist;

  if (meep::am_master()) {
    for (int i = 0; i < geometry.num_items; ++i) {
      display_geometric_object_info(5, geometry.items[i]);
      
      if (geometry.items[i].material.which_subclass 
	  == MTS::MEDIUM)
	printf("%*sdielectric constant epsilon diagonal = (%g,%g,%g)\n",
	       5 + 5, "",
	       geometry.items[i].material.
	       subclass.medium_data->epsilon_diag.x,
	       geometry.items[i].material.
	       subclass.medium_data->epsilon_diag.y,
	       geometry.items[i].material.
	       subclass.medium_data->epsilon_diag.z);
    }
  }
  
  geom_fix_objects0(geometry);
  geom_box box = gv2box(gv);
  geometry_tree = create_geom_box_tree0(geometry, box);
  if (verbose && meep::am_master()) {
    printf("Geometric-object bounding-box tree:\n");
    display_geom_box_tree(5, geometry_tree);
    
    int tree_depth, tree_nobjects;
    geom_box_tree_stats(geometry_tree, &tree_depth, &tree_nobjects);
    master_printf("Geometric object tree has depth %d "
		  "and %d object nodes (vs. %d actual objects)\n",
		  tree_depth, tree_nobjects, geometry.num_items);
  }
  
  restricted_tree = geometry_tree;
}

geom_epsilon::~geom_epsilon()
{
  unset_volume();
  destroy_geom_box_tree(geometry_tree);
}

void geom_epsilon::unset_volume(void)
{
  if (restricted_tree != geometry_tree) {
    destroy_geom_box_tree(restricted_tree);
    restricted_tree = geometry_tree;
  }
}

void geom_epsilon::set_volume(const meep::geometric_volume &gv)
{
  unset_volume();
  
  geom_box box = gv2box(gv);
  restricted_tree = create_geom_box_tree0(geometry, box);
}

static material_type eval_material_func(function material_func, vector3 p)
{
  SCM pscm = ctl_convert_vector3_to_scm(p);
  material_type material;
  SCM mo;
  
  mo = gh_call1(material_func, pscm);
  material_type_input(mo, &material);
  
  while (material.which_subclass == MTS::MATERIAL_FUNCTION) {
    material_type m;
    
    mo = gh_call1(material.subclass.
		  material_function_data->material_func,
		  pscm);
    material_type_input(mo, &m);
    material_type_destroy(material);
    material = m;
  }
  
  if (material.which_subclass == MTS::MATERIAL_TYPE_SELF) {
    material_type_copy(&default_material, &material);
  }
  CK(material.which_subclass != MTS::MATERIAL_FUNCTION,
     "infinite loop in material functions");
  
  return material;
}

static int variable_material(int which_subclass)
{
     return (which_subclass == MTS::MATERIAL_FUNCTION);
}

static bool is_metal(meep::field_type ft, const material_type *material) {
  if (ft == meep::E_stuff)
    switch (material->which_subclass) {
    case MTS::MEDIUM:
      return (material->subclass.medium_data->epsilon_diag.x < 0 ||
	      material->subclass.medium_data->epsilon_diag.y < 0 ||
	      material->subclass.medium_data->epsilon_diag.z < 0);
    case MTS::PERFECT_METAL:
      return true;
    default:
      meep::abort("unknown material type");
  }
  else
    switch (material->which_subclass) {
    case MTS::MEDIUM:
      return (material->subclass.medium_data->mu_diag.x < 0 ||
	      material->subclass.medium_data->mu_diag.y < 0 ||
	      material->subclass.medium_data->mu_diag.z < 0);
    case MTS::PERFECT_METAL:
      return false; // is an electric conductor, but not a magnetic conductor
    default:
      meep::abort("unknown material type");
  }
}

static void material_epsmu(meep::field_type ft, material_type material, 
		    symmetric_matrix *epsmu, symmetric_matrix *epsmu_inv) {
  if (ft == meep::E_stuff)
    switch (material.which_subclass) {
    case MTS::MEDIUM:
      {
      epsmu->m00 = material.subclass.medium_data->epsilon_diag.x;
      epsmu->m11 = material.subclass.medium_data->epsilon_diag.y;
      epsmu->m22 = material.subclass.medium_data->epsilon_diag.z;
      epsmu->m01 = material.subclass.medium_data->epsilon_offdiag.x;
      epsmu->m02 = material.subclass.medium_data->epsilon_offdiag.y;
      epsmu->m12 = material.subclass.medium_data->epsilon_offdiag.z;
      sym_matrix_invert(epsmu_inv,epsmu);
      break;
      }
    case MTS::PERFECT_METAL:
      {
      epsmu->m00 = -meep::infinity;
      epsmu->m11 = -meep::infinity;
      epsmu->m22 = -meep::infinity;
      epsmu_inv->m00 = -0.0;
      epsmu_inv->m11 = -0.0;
      epsmu_inv->m22 = -0.0;
      break;
      }
    default:
      meep::abort("unknown material type");
  }
  else
    switch (material.which_subclass) {
    case MTS::MEDIUM:
      {
      epsmu->m00 = material.subclass.medium_data->mu_diag.x;
      epsmu->m11 = material.subclass.medium_data->mu_diag.y;
      epsmu->m22 = material.subclass.medium_data->mu_diag.z;
      epsmu->m01 = material.subclass.medium_data->mu_offdiag.x;
      epsmu->m02 = material.subclass.medium_data->mu_offdiag.y;
      epsmu->m12 = material.subclass.medium_data->mu_offdiag.z;
      sym_matrix_invert(epsmu_inv,epsmu);
      break;
      }
    case MTS::PERFECT_METAL:
      {
      epsmu->m00 = 1.0;
      epsmu->m11 = 1.0;
      epsmu->m22 = 1.0;
      epsmu_inv->m00 = 1.0;
      epsmu_inv->m11 = 1.0;
      epsmu_inv->m22 = 1.0;
      break;
      }
    default:
      meep::abort("unknown material type");
  }
}

bool geom_epsilon::get_material_pt(material_type &material, const meep::vec &r)
{
  vector3 p = vec_to_vector3(r);
  boolean inobject;
  material =
    material_of_unshifted_point_in_tree_inobject(p, restricted_tree,&inobject);

  bool destroy_material = false;
  if (material.which_subclass == MTS::MATERIAL_TYPE_SELF) {
    material = default_material;
  }
  if (material.which_subclass == MTS::MATERIAL_FUNCTION) {
    material = eval_material_func(material.subclass.
                                  material_function_data->material_func,
                                  p);
    destroy_material = true;
  }
  return destroy_material;
}

// returns trace of the tensor diagonal
double geom_epsilon::chi1p1(meep::field_type ft, const meep::vec &r)
{
  symmetric_matrix chi1p1, chi1p1_inv;

#ifdef DEBUG
  vector3 p = vec_to_vector3(r);
  if (p.x < restricted_tree->b.low.x ||
      p.y < restricted_tree->b.low.y ||
      p.z < restricted_tree->b.low.z ||
      p.x > restricted_tree->b.high.x ||
      p.y > restricted_tree->b.high.y ||
      p.z > restricted_tree->b.high.z)
    meep::abort("invalid point (%g,%g,%g)\n", p.x,p.y,p.z);
#endif

  material_type material;
  bool destroy_material = get_material_pt(material, r);

  material_epsmu(ft, material, &chi1p1, &chi1p1_inv);  
  
  if (destroy_material)
    material_type_destroy(material);
  
  return (chi1p1.m00 + chi1p1.m11 + chi1p1.m22)/3;
}

/* Find frontmost object in gv, along with the constant material behind it.
   Returns false if material behind the object is not constant.
   
   Requires moderately horrifying logic to figure things out properly,
   stolen from MPB. */
static bool get_front_object(const meep::geometric_volume &gv,
			     geom_box_tree geometry_tree,
			     vector3 &pcenter,
			     const geometric_object **o_front,
			     vector3 &shiftby_front,
			     material_type &mat_front,
			     material_type &mat_behind) {
  vector3 p;
  const geometric_object *o1 = 0, *o2 = 0;
  vector3 shiftby1 = {0,0,0}, shiftby2 = {0,0,0};
  geom_box pixel;
  material_type mat1, mat2;
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
  pixel = gv2box(gv);
  pcenter = p = vec_to_vector3(gv.center());
  double d1, d2, d3;
  d1 = (pixel.high.x - pixel.low.x) * 0.5;
  d2 = (pixel.high.y - pixel.low.y) * 0.5;
  d3 = (pixel.high.z - pixel.low.z) * 0.5;
  for (int i = 0; i < num_neighbors[dimensions - 1]; ++i) {
    const geometric_object *o;
    material_type mat;
    vector3 q, shiftby;
    int id;
    q.x = p.x + neighbors[dimensions - 1][i][0] * d1;
    q.y = p.y + neighbors[dimensions - 1][i][1] * d2;
    q.z = p.z + neighbors[dimensions - 1][i][2] * d3;
    o = object_of_point_in_tree(q, geometry_tree, &shiftby, &id);
    if ((id == id1 && vector3_equal(shiftby, shiftby1)) ||
	(id == id2 && vector3_equal(shiftby, shiftby2)))
      continue;
    mat = (o && o->material.which_subclass != MTS::MATERIAL_TYPE_SELF)
      ? o->material : default_material;
    if (id1 == -1) {
      o1 = o;
      shiftby1 = shiftby;
      id1 = id;
      mat1 = mat;
    }
    else if (id2 == -1 || ((id >= id1 && id >= id2) &&
			   (id1 == id2 
			    || material_type_equal(&mat1,&mat2)))) {
      o2 = o;
      shiftby2 = shiftby;
      id2 = id;
      mat2 = mat;
    }
    else if (!(id1 < id2 && 
	       (id1 == id || material_type_equal(&mat1,&mat))) &&
	     !(id2 < id1 &&
	       (id2 == id || material_type_equal(&mat2,&mat))))
      return false;
  }

  // CHECK(id1 > -1, "bug in object_of_point_in_tree?");
  if (id2 == -1) { /* only one nearby object/material */
    id2 = id1;
    o2 = o1;
    mat2 = mat1;
    shiftby2 = shiftby1;
  }

  if ((o1 && variable_material(o1->material.which_subclass)) ||
      (o2 && variable_material(o2->material.which_subclass)) ||
      (variable_material(default_material.which_subclass)
       && (!o1 || !o2 ||
	   o1->material.which_subclass == MTS::MATERIAL_TYPE_SELF ||
	   o2->material.which_subclass == MTS::MATERIAL_TYPE_SELF)))
    return false;

  if (id1 >= id2) {
    *o_front = o1;
    shiftby_front = shiftby1;
    mat_front = mat1;
    if (id1 == id2) mat_behind = mat1; else mat_behind = mat2;
  }
  if (id2 > id1) {
    *o_front = o2;
    shiftby_front = shiftby2;
    mat_front = mat2;
    mat_behind = mat1;
  }
  return true;
}

void geom_epsilon::eff_chi1inv_row(meep::component c, double chi1inv_row[3],
				   const meep::geometric_volume &gv,
				   double tol, int maxeval) {
  const geometric_object *o;
  material_type mat, mat_behind;
  symmetric_matrix meps, meps_inv;
  vector3 p, shiftby, normal;
  bool destroy_material = false;

  if (maxeval == 0 || !get_front_object(gv, geometry_tree,
					p, &o, shiftby, mat, mat_behind)) {
  noavg:
    destroy_material = get_material_pt(mat, gv.center());
  trivial:    
    material_epsmu(meep::type(c), mat, &meps, &meps_inv);
    switch (component_direction(c)) {
    case meep::X: case meep::R:
      chi1inv_row[0] = meps_inv.m00;
      chi1inv_row[1] = meps_inv.m01;
      chi1inv_row[2] = meps_inv.m02;
      break;
    case meep::Y: case meep::P:
      chi1inv_row[0] = meps_inv.m01;
      chi1inv_row[1] = meps_inv.m11;
      chi1inv_row[2] = meps_inv.m12;
      break;
    case meep::Z:
      chi1inv_row[0] = meps_inv.m02;
      chi1inv_row[1] = meps_inv.m12;
      chi1inv_row[2] = meps_inv.m22;
      break;
    case meep::NO_DIRECTION: chi1inv_row[0] = chi1inv_row[1] = chi1inv_row[2] = 0;
    }
    if (destroy_material) material_type_destroy(mat);
    return;
  }

  // FIXME: reimplement support for fallback integration, without
  //        messing up anisotropic support
  //  if (!get_front_object(gv, geometry_tree,
  //			p, &o, shiftby, mat, mat_behind)) {
  //     fallback_chi1inv_row(c, chi1inv_row, gv, tol, maxeval);
  //     return;
  //  }

  /* check for trivial case of only one object/material */
  if (material_type_equal(&mat, &mat_behind)) goto trivial;

  // it doesn't make sense to average metals (electric or magnetic)
  if (is_metal(meep::type(c), &mat) || is_metal(meep::type(c), &mat_behind))
    goto noavg;

  normal = unit_vector3(normal_to_fixed_object(vector3_minus(p, shiftby), *o));
  geom_box pixel = gv2box(gv);
  pixel.low = vector3_minus(pixel.low, shiftby);
  pixel.high = vector3_minus(pixel.high, shiftby);

  double fill = box_overlap_with_object(pixel, *o, tol, maxeval); 

  material_epsmu(meep::type(c), mat, &meps, &meps_inv);
  symmetric_matrix eps2, epsinv2;
  symmetric_matrix eps1, delta;
  double Rot[3][3];
  material_epsmu(meep::type(c), mat_behind, &eps2, &epsinv2);
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

#  ifdef DEBUG
	  if(!sym_matrix_positive_definite(&meps))
	    meep::abort("negative mean epsilon from Kottke algorithm");
#  endif

  sym_matrix_invert(&meps_inv, &meps);
  switch (component_direction(c)) {
  case meep::X: case meep::R:
    chi1inv_row[0] = meps_inv.m00;
    chi1inv_row[1] = meps_inv.m01;
    chi1inv_row[2] = meps_inv.m02;
    break;
  case meep::Y: case meep::P:
    chi1inv_row[0] = meps_inv.m01;
    chi1inv_row[1] = meps_inv.m11;
    chi1inv_row[2] = meps_inv.m12;
    break;
  case meep::Z:
    chi1inv_row[0] = meps_inv.m02;
    chi1inv_row[1] = meps_inv.m12;
    chi1inv_row[2] = meps_inv.m22;
    break;
  case meep::NO_DIRECTION: chi1inv_row[0] = chi1inv_row[1] = chi1inv_row[2] = 0;
  }
}

static int eps_ever_negative = 0;
static meep::field_type func_ft = meep::E_stuff;

#ifdef CTL_HAS_COMPLEX_INTEGRATION
static cnumber ceps_func(int n, number *x, void *geomeps_)
{
  geom_epsilon *geomeps = (geom_epsilon *) geomeps_;
  vector3 p = {0,0,0};
  p.x = x[0]; p.y = n > 1 ? x[1] : 0; p.z = n > 2 ? x[2] : 0;
  double s = 1;
  if (dim == meep::Dcyl) { double py = p.y; p.y = p.z; p.z = py; s = p.x; }
  cnumber ret;
  double ep = geomeps->chi1p1(func_ft, vector3_to_vec(p));
  if (ep < 0) eps_ever_negative = 1;
  ret.re = ep * s;
  ret.im = s / ep;
  return ret;
}
#else
static number eps_func(int n, number *x, void *geomeps_)
{
  geom_epsilon *geomeps = (geom_epsilon *) geomeps_;
  vector3 p = {0,0,0};
  double s = 1;
  p.x = x[0]; p.y = n > 1 ? x[1] : 0; p.z = n > 2 ? x[2] : 0;
  if (dim == meep::Dcyl) { double py = p.y; p.y = p.z; p.z = py; s = p.x; }
  double ep = geomeps->chi1p1(func_ft, vector3_to_vec(p));
  if (ep < 0) eps_ever_negative = 1;
  return ep * s;
}
static number inveps_func(int n, number *x, void *geomeps_)
{
  geom_epsilon *geomeps = (geom_epsilon *) geomeps_;
  vector3 p = {0,0,0};
  double s = 1;
  p.x = x[0]; p.y = n > 1 ? x[1] : 0; p.z = n > 2 ? x[2] : 0;
  if (dim == meep::Dcyl) { double py = p.y; p.y = p.z; p.z = py; s = p.x; }
  double ep = geomeps->chi1p1(func_ft, vector3_to_vec(p));
  if (ep < 0) eps_ever_negative = 1;
  return s / ep;
}
#endif

// fallback meaneps using libctl's adaptive cubature routine
void geom_epsilon::fallback_chi1inv_row(meep::component c,
					double chi1inv_row[3],
					const meep::geometric_volume &gv,
					double tol, int maxeval)
{

  symmetric_matrix chi1p1, chi1p1_inv;
  material_type material;
  bool destroy_material = get_material_pt(material, gv.center());
  material_epsmu(meep::type(c), material, &chi1p1, &chi1p1_inv);
  if (destroy_material)
    material_type_destroy(material);
  if (chi1p1.m01 != 0 || chi1p1.m02 != 0 || chi1p1.m12 != 0
      || chi1p1.m00 != chi1p1.m11 || chi1p1.m11 != chi1p1.m22 || 
      chi1p1.m00 != chi1p1.m22) {
    int rownum = meep::component_direction(c) % 3;
    if (rownum == 0) {
      chi1inv_row[0] = chi1p1.m00;
      chi1inv_row[1] = chi1p1.m01;
      chi1inv_row[2] = chi1p1.m02;
    } 
    else if (rownum == 1) {
      chi1inv_row[0] = chi1p1.m01;
      chi1inv_row[1] = chi1p1.m11;
      chi1inv_row[2] = chi1p1.m12;
    }
    else {
      chi1inv_row[0] = chi1p1.m02;
      chi1inv_row[1] = chi1p1.m12;
      chi1inv_row[2] = chi1p1.m22;
    }
    return;
      }

  number esterr;
  integer errflag, n;
  number xmin[3], xmax[3];
  vector3 gvmin, gvmax;
  gvmin = vec_to_vector3(gv.get_min_corner());
  gvmax = vec_to_vector3(gv.get_max_corner());
  xmin[0] = gvmin.x; xmax[0] = gvmax.x; 
  if (dim == meep::Dcyl) {
    xmin[1] = gvmin.z; xmin[2] = gvmin.y; xmax[1] = gvmax.z; xmax[2] = gvmax.y;
  }
  else{
    xmin[1] = gvmin.y; xmin[2] = gvmin.z; xmax[1] = gvmax.y; xmax[2] = gvmax.z;
  }
  if (xmin[2] == xmax[2])
    n = xmin[1] == xmax[1] ? 1 : 2;
  else
    n = 3;
  double vol = 1;
  for (int i = 0; i < n; ++i) vol *= xmax[i] - xmin[i];
  if (dim == meep::Dcyl) vol *= (xmin[0] + xmax[0]) * 0.5;
  eps_ever_negative = 0;
  func_ft = meep::type(c);
  double meps, minveps;
#ifdef CTL_HAS_COMPLEX_INTEGRATION
  cnumber ret = cadaptive_integration(ceps_func, xmin, xmax, n, (void*) this,
				      0, tol, maxeval, &esterr, &errflag);
  meps = ret.re / vol;
  minveps = ret.im / vol;
#else
  meps = adaptive_integration(eps_func, xmin, xmax, n, (void*) this,
			      0, tol, maxeval, &esterr, &errflag) / vol;
  minveps = adaptive_integration(inveps_func, xmin, xmax, n, (void*) this,
				 0, tol, maxeval, &esterr, &errflag) / vol;
#endif
  if (eps_ever_negative) // averaging negative eps causes instability
    minveps = 1.0 / (meps = eps(gv.center()));

  {
    meep::vec gradient(normal_vector(meep::type(c), gv));
    double n[3] = {0,0,0};
    double nabsinv = 1.0/meep::abs(gradient);
    LOOP_OVER_DIRECTIONS(gradient.dim, k)
      n[k%3] = gradient.in_direction(k) * nabsinv;
    int rownum = meep::component_direction(c) % 3;
    for (int i=0; i<3; ++i)
      chi1inv_row[i] = n[rownum] * n[i] * (minveps - 1/meps);
    chi1inv_row[rownum] += 1/meps;
  }
}

static double get_chi3(meep::component c, const medium *m) {
  switch (c) {
  case meep::Er: case meep::Ex: return m->E_chi3_diag.x;
  case meep::Ep: case meep::Ey: return m->E_chi3_diag.y;
  case meep::Ez: return m->E_chi3_diag.z;
  case meep::Hr: case meep::Hx: return m->H_chi3_diag.x;
  case meep::Hp: case meep::Hy: return m->H_chi3_diag.y;
  case meep::Hz: return m->H_chi3_diag.z;
  default: return 0;
  }
}

bool geom_epsilon::has_chi3(meep::component c)
{
  for (int i = 0; i < geometry.num_items; ++i) {
    if (geometry.items[i].material.which_subclass == MTS::MEDIUM) {
      if (get_chi3(c, geometry.items[i].material.subclass.medium_data) != 0)
	return true; 
    }
  }
  for (int i = 0; i < extra_materials.num_items; ++i)
    if (extra_materials.items[i].which_subclass == MTS::MEDIUM)
      if (get_chi3(c, extra_materials.items[i].subclass.medium_data) != 0)
        return true;
  return (default_material.which_subclass == MTS::MEDIUM &&
	  get_chi3(c, default_material.subclass.medium_data) != 0);
}

double geom_epsilon::chi3(meep::component c, const meep::vec &r) {
  material_type material;
  bool destroy_material = get_material_pt(material, r);
  
  double chi3_val;
  switch (material.which_subclass) {
  case MTS::MEDIUM:
    chi3_val = get_chi3(c, material.subclass.medium_data);
    break;
  default:
    chi3_val = 0;
  }
  
  if (destroy_material)
    material_type_destroy(material);
  
  return chi3_val;
}

static double get_chi2(meep::component c, const medium *m) {
  switch (c) {
  case meep::Er: case meep::Ex: return m->E_chi2_diag.x;
  case meep::Ep: case meep::Ey: return m->E_chi2_diag.y;
  case meep::Ez: return m->E_chi2_diag.z;
  case meep::Hr: case meep::Hx: return m->H_chi2_diag.x;
  case meep::Hp: case meep::Hy: return m->H_chi2_diag.y;
  case meep::Hz: return m->H_chi2_diag.z;
  default: return 0;
  }
}

bool geom_epsilon::has_chi2(meep::component c)
{
  for (int i = 0; i < geometry.num_items; ++i) {
    if (geometry.items[i].material.which_subclass == MTS::MEDIUM) {
      if (get_chi2(c, geometry.items[i].material.subclass.medium_data) != 0)
	return true; 
    }
  }
  for (int i = 0; i < extra_materials.num_items; ++i)
    if (extra_materials.items[i].which_subclass == MTS::MEDIUM)
      if (get_chi2(c, extra_materials.items[i].subclass.medium_data) != 0)
        return true;
  return (default_material.which_subclass == MTS::MEDIUM &&
	  get_chi2(c, default_material.subclass.medium_data) != 0);
}

double geom_epsilon::chi2(meep::component c, const meep::vec &r) {
  material_type material;
  bool destroy_material = get_material_pt(material, r);
  
  double chi2_val;
  switch (material.which_subclass) {
  case MTS::MEDIUM:
    chi2_val = get_chi2(c, material.subclass.medium_data);
    break;
  default:
    chi2_val = 0;
  }
  
  if (destroy_material)
    material_type_destroy(material);
  
  return chi2_val;
}

static bool mu_not_1(material_type &m)
{
  return (m.which_subclass == MTS::MEDIUM &&
	  (m.subclass.medium_data->mu_diag.x != 1 ||
	   m.subclass.medium_data->mu_diag.y != 1 ||
	   m.subclass.medium_data->mu_diag.z != 1 ||
	   m.subclass.medium_data->mu_offdiag.x != 0 ||
	   m.subclass.medium_data->mu_offdiag.y != 0 ||
	   m.subclass.medium_data->mu_offdiag.z != 0));
}

bool geom_epsilon::has_mu()
{
  for (int i = 0; i < geometry.num_items; ++i) {
    if (mu_not_1(geometry.items[i].material))
      return true; 
  }
  for (int i = 0; i < extra_materials.num_items; ++i)
    if (mu_not_1(extra_materials.items[i]))
      return true;
  return (mu_not_1(default_material));
}

/* a global scalar conductivity to add to all materials; this
   is mostly for the convenience of Casimir calculations where
   the global conductivity corresponds to a rotation to
   complex frequencies */
static double global_D_conductivity = 0, global_B_conductivity = 0;

static double get_cnd(meep::component c, const medium *m) {
  switch (c) {
  case meep::Dr: case meep::Dx: return m->D_conductivity_diag.x + global_D_conductivity;
  case meep::Dp: case meep::Dy: return m->D_conductivity_diag.y + global_D_conductivity;
  case meep::Dz: return m->D_conductivity_diag.z + global_D_conductivity;
  case meep::Br: case meep::Bx: return m->B_conductivity_diag.x + global_B_conductivity;
  case meep::Bp: case meep::By: return m->B_conductivity_diag.y + global_B_conductivity;
  case meep::Bz: return m->B_conductivity_diag.z + global_B_conductivity;
  default: return 0;
  }
}

bool geom_epsilon::has_conductivity(meep::component c)
{
  for (int i = 0; i < geometry.num_items; ++i) {
    if (geometry.items[i].material.which_subclass == MTS::MEDIUM) {
      if (get_cnd(c, geometry.items[i].material.subclass.medium_data) != 0)
	return true; 
    }
  }
  for (int i = 0; i < extra_materials.num_items; ++i)
    if (extra_materials.items[i].which_subclass == MTS::MEDIUM)
      if (get_cnd(c, extra_materials.items[i].subclass.medium_data) != 0)
        return true;
  return (default_material.which_subclass == MTS::MEDIUM &&
	  get_cnd(c, default_material.subclass.medium_data) != 0);
}

double geom_epsilon::conductivity(meep::component c, const meep::vec &r) {
  material_type material;
  bool destroy_material = get_material_pt(material, r);

  double cond_val;
  switch (material.which_subclass) {
  case MTS::MEDIUM:
    cond_val = get_cnd(c, material.subclass.medium_data);
    break;
  default:
    cond_val = 0;
  }
  
  if (destroy_material)
    material_type_destroy(material);
  
  return cond_val;
}

void geom_epsilon::sigma_row(meep::component c, double sigrow[3], 
			     const meep::vec &r) {
  vector3 p = vec_to_vector3(r);

  boolean inobject;
  material_type material =
    material_of_unshifted_point_in_tree_inobject(p, restricted_tree, &inobject);
  
  int destroy_material = 0;
  if (material.which_subclass == MTS::MATERIAL_TYPE_SELF) {
    material = default_material;
  }
  if (material.which_subclass == MTS::MATERIAL_FUNCTION) {
    material = eval_material_func(material.subclass.
				  material_function_data->material_func,
				  p);
    destroy_material = 1;
  }
  
  sigrow[0] = sigrow[1] = sigrow[2] = 0.0;
  if (material.which_subclass == MTS::MEDIUM) {
    polarizability_list plist = 
      pol_ft == meep::E_stuff ? material.subclass.medium_data->E_polarizations
      : material.subclass.medium_data->H_polarizations;
    for (int j = 0; j < plist.num_items; ++j)
      if (plist.items[j].omega == omega &&
	  plist.items[j].gamma == gamma) {
	int ic = meep::component_index(c);
	sigrow[ic] = (ic == 0 ? plist.items[j].sigma_diag.x
		      : (ic == 1 ? plist.items[j].sigma_diag.y
			 : plist.items[j].sigma_diag.z));
	break;
      }
  }
  
  if (destroy_material)
    material_type_destroy(material);
}

struct pol {
  double omega, gamma;
  struct pol *next;
};

// add a polarization to the list if it is not already there
static pol *add_pol(pol *pols,
		    double omega, double gamma)
{
  struct pol *p = pols;
  while (p && !(p->omega == omega && p->gamma == gamma))
    p = p->next;
  if (!p) {
    p = new pol;
    p->omega = omega;
    p->gamma = gamma;
    p->next = pols;
    pols = p;
  }
  return pols;
}

static pol *add_pols(pol *pols, const polarizability_list plist) {
  for (int j = 0; j < plist.num_items; ++j) {
    pols = add_pol(pols,
		   plist.items[j].omega, plist.items[j].gamma);
  }
  return pols;
}

void geom_epsilon::add_polarizabilities(meep::structure *s) {
  add_polarizabilities(meep::E_stuff, s);
  add_polarizabilities(meep::H_stuff, s);
}

void geom_epsilon::add_polarizabilities(meep::field_type ft, 
					meep::structure *s) {
  pol *pols = 0;

  // construct a list of the unique polarizabilities in the geometry:
  for (int i = 0; i < geometry.num_items; ++i) {
    if (geometry.items[i].material.which_subclass == MTS::MEDIUM)
      pols = add_pols(pols, ft == meep::E_stuff
		      ? geometry.items[i].material
		      .subclass.medium_data->E_polarizations
		      : geometry.items[i].material
		      .subclass.medium_data->H_polarizations);
  }
  for (int i = 0; i < extra_materials.num_items; ++i)
    if (extra_materials.items[i].which_subclass == MTS::MEDIUM)
      pols = add_pols(pols, ft == meep::E_stuff
		      ? extra_materials.items[i]
		      .subclass.medium_data->E_polarizations
		      : extra_materials.items[i]
		      .subclass.medium_data->H_polarizations);
  if (default_material.which_subclass == MTS::MEDIUM)
    pols = add_pols(pols, ft == meep::E_stuff
		    ? default_material.subclass.medium_data->E_polarizations
		    : default_material.subclass.medium_data->H_polarizations);
    
  for (struct pol *p = pols; p; p = p->next) {
    master_printf("polarizability: omega=%g, gamma=%g\n",
		  p->omega, p->gamma);
    s->add_polarizability(*this, ft, p->omega, p->gamma);
  }
  
  while (pols) {
    struct pol *p = pols;
    pols = pols->next;
    delete p;
  }
}

/***********************************************************************/

// wrapper around Scheme function for PML profile
static double scm_pml_profile(double u, void *f_)
{
  SCM f = (SCM) f_;
  return ctl_convert_number_to_c(gh_call1(f, ctl_convert_number_to_scm(u)));
}

// for passing to multidimensional integration routine
static double scm_pml_profile2(int dim, double *u, void *f_)
{
  SCM f = (SCM) f_; (void) dim;
  return ctl_convert_number_to_c(gh_call1(f, ctl_convert_number_to_scm(*u)));
}

meep::structure *make_structure(int dims, vector3 size, vector3 center,
				double resolution, bool enable_averaging,
				double subpixel_tol, int subpixel_maxeval,
				bool ensure_periodicity_p,
				geometric_object_list geometry,
				material_type_list extra_materials,
				material_type default_mat,
				pml_list pml_layers,
				symmetry_list symmetries,
				int num_chunks, double Courant,
				double global_D_conductivity_,
				double global_B_conductivity_)
{
  master_printf("-----------\nInitializing structure...\n");
  
  // only cartesian lattices are currently allowed
  geom_initialize();
  geometry_center = center;

  global_D_conductivity = global_D_conductivity_;
  global_B_conductivity = global_B_conductivity_;
  
  number no_size = 2.0 / ctl_get_number("infinity");
  if (size.x <= no_size)
    size.x = 0.0;
  if (size.y <= no_size)
    size.y = 0.0;
  if (size.z <= no_size)
    size.z = 0.0;
  
  set_dimensions(dims);
  
  geometry_lattice.size = size;

  master_printf("Working in %s dimensions.\n", meep::dimension_name(dim));
  master_printf("Computational cell is %g x %g x %g with resolution %g\n",
                size.x, size.y, size.z, resolution);  

  meep::volume v;
  switch (dims) {
  case 0: case 1:
    v = meep::vol1d(size.z, resolution);
    break;
  case 2:
    v = meep::vol2d(size.x, size.y, resolution);
    break;
  case 3:
    v = meep::vol3d(size.x, size.y, size.z, resolution);
    break;
  case CYLINDRICAL:
    v = meep::volcyl(size.x, size.z, resolution);
    break;
  default:
    CK(0, "unsupported dimensionality");
  }
  v.center_origin();
  v.shift_origin(vector3_to_vec(center));
  
  meep::symmetry S;
  for (int i = 0; i < symmetries.num_items; ++i) 
    switch (symmetries.items[i].which_subclass) {
    case symmetry::SYMMETRY_SELF: break; // identity
    case symmetry::MIRROR_SYM:
      S = S + meep::mirror(meep::direction(symmetries.items[i].direction), v)
	* complex<double>(symmetries.items[i].phase.re,
			  symmetries.items[i].phase.im);
      break;
    case symmetry::ROTATE2_SYM:
      S = S + meep::rotate2(meep::direction(symmetries.items[i].direction), v)
	* complex<double>(symmetries.items[i].phase.re,
			  symmetries.items[i].phase.im);
      break;
    case symmetry::ROTATE4_SYM:
      S = S + meep::rotate4(meep::direction(symmetries.items[i].direction), v)
	* complex<double>(symmetries.items[i].phase.re,
			  symmetries.items[i].phase.im);
      break;
    }

  meep::boundary_region br;
  for (int i = 0; i < pml_layers.num_items; ++i) {
    double umin = 0, umax = 1, esterr;
    int errflag;
    using namespace meep;
    if (pml_layers.items[i].direction == -1) {
      LOOP_OVER_DIRECTIONS(v.dim, d) {
	if (pml_layers.items[i].side == -1) {
	  FOR_SIDES(b)
	    br = br + meep::boundary_region
	    (meep::boundary_region::PML,
	     pml_layers.items[i].thickness,
	     pow(pml_layers.items[i].R_asymptotic,
		 pml_layers.items[i].strength),
	     scm_pml_profile, pml_layers.items[i].pml_profile,
	     adaptive_integration(scm_pml_profile2, &umin, &umax, 1,
				  (void*) pml_layers.items[i].pml_profile,
				  1e-9, 1e-4, 50000, &esterr, &errflag),
	     d, b);
	}
	else
	  br = br + meep::boundary_region
	    (meep::boundary_region::PML,
	     pml_layers.items[i].thickness,
	     pow(pml_layers.items[i].R_asymptotic,
		 pml_layers.items[i].strength),
	     scm_pml_profile, pml_layers.items[i].pml_profile,
	     adaptive_integration(scm_pml_profile2, &umin, &umax, 1,
				  (void*) pml_layers.items[i].pml_profile,
				  1e-9, 1e-4, 50000, &esterr, &errflag),
	     d,
	     (meep::boundary_side) pml_layers.items[i].side);
      }
    }
    else {
	if (pml_layers.items[i].side == -1) {
	  FOR_SIDES(b)
	    br = br + meep::boundary_region
	    (meep::boundary_region::PML,
	     pml_layers.items[i].thickness,
	     pow(pml_layers.items[i].R_asymptotic,
		 pml_layers.items[i].strength),
	     scm_pml_profile, pml_layers.items[i].pml_profile,
	     adaptive_integration(scm_pml_profile2, &umin, &umax, 1,
				  (void*) pml_layers.items[i].pml_profile,
				  1e-9, 1e-4, 50000, &esterr, &errflag),
	     (meep::direction) pml_layers.items[i].direction,
	     b);
	}
	else
	  br = br + meep::boundary_region
	    (meep::boundary_region::PML,
	     pml_layers.items[i].thickness,
	     pow(pml_layers.items[i].R_asymptotic,
		 pml_layers.items[i].strength),
	     scm_pml_profile, pml_layers.items[i].pml_profile,
	     adaptive_integration(scm_pml_profile2, &umin, &umax, 1,
				  (void*) pml_layers.items[i].pml_profile,
				  1e-9, 1e-4, 50000, &esterr, &errflag),
	     (meep::direction) pml_layers.items[i].direction,
	     (meep::boundary_side) pml_layers.items[i].side);
    }
  }
  
  ensure_periodicity = ensure_periodicity_p;
  default_material = default_mat;
  geom_epsilon geps(geometry, extra_materials, v.pad().surroundings());

  if (subpixel_maxeval < 0) subpixel_maxeval = 0; // no limit

  meep::structure *s = new meep::structure(v, geps, br, S, 
					   num_chunks, Courant,
					   enable_averaging,
					   subpixel_tol,
					   subpixel_maxeval);

  geps.add_polarizabilities(s);

  master_printf("-----------\n");
  
  return s;
}

/*************************************************************************/
