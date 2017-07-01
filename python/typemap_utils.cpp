/* Copyright (C) 2005-2017 Massachusetts Institute of Technology  
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2, or (at your option)
 *  any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software Foundation,
 *  Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */

// Utility functions for pymeep typemaps

static PyObject* vec2py(const meep::vec &v) {

    double x = 0, y = 0, z = 0;

    switch (v.dim) {
      case meep::D1:
        z = v.z();
        break;
      case meep::D2:
        x = v.x();
        y = v.y();
        break;
      case meep::D3:
        x = v.x();
        y = v.y();
        z = v.z();
        break;
      case meep::Dcyl:
        return Py_BuildValue("(ddd)", v.r(), v.z(), 0);
    }

    return Py_BuildValue("(ddd)", x, y, z);
}

static double py_callback_wrap(const meep::vec &v) {
    PyObject *pyv = vec2py(v);
    PyObject *pyret = PyObject_CallFunctionObjArgs(py_callback, pyv, NULL);
    Py_XDECREF(pyv);
    double ret = PyFloat_AsDouble(pyret);
    Py_XDECREF(pyret);
    return ret;
}


static int pyv3_to_v3(PyObject *po, vector3 *v) {

    PyObject *py_x = PyObject_GetAttrString(po, "x");
    PyObject *py_y = PyObject_GetAttrString(po, "y");
    PyObject *py_z = PyObject_GetAttrString(po, "z");

    if(!py_x || !py_y || !py_z) {
        PyErr_SetString(PyExc_ValueError, "Vector3 is not initialized");
        return 0;
    }

    double x = PyFloat_AsDouble(py_x);
    double y = PyFloat_AsDouble(py_y);
    double z = PyFloat_AsDouble(py_z);
    Py_DECREF(py_x);
    Py_DECREF(py_y);
    Py_DECREF(py_z);

    v->x = x;
    v->y = y;
    v->z = z;

    return 1;
}

static int get_attr_v3(PyObject *py_obj, vector3 *v, const char *name) {
    PyObject *py_attr = PyObject_GetAttrString(py_obj, name);

    if(!py_attr) {
        PyErr_Format(PyExc_ValueError, "Class attribute '%s' is None\n", name);
        return 0;
    }

    if(!pyv3_to_v3(py_attr, v)) {
        return 0;
    }

    Py_XDECREF(py_attr);
    return 1;
}

static int get_attr_dbl(PyObject *py_obj, double *result, const char *name) {
    PyObject *py_attr = PyObject_GetAttrString(py_obj, name);

    if(!py_attr) {
        PyErr_Format(PyExc_ValueError, "Class attribute '%s' is None\n", name);
        return 0;
    }

    *result = PyFloat_AsDouble(py_attr);
    Py_XDECREF(py_attr);
    return 1;
}

static int get_attr_material(PyObject *po, material_type *m) {
    PyObject *py_material = PyObject_GetAttrString(po, "material");

    if(!py_material) {
        PyErr_SetString(PyExc_ValueError, "Object's material is not set\n");
        return 0;
    }

    if(!pymaterial_to_material(py_material, m)) {
        return 0;
    }

    Py_XDECREF(py_material);

    return 1;
}

static int py_susceptibility_to_susceptibility(PyObject *po, susceptibility_struct *s) {
    // TODO(chogan): Accont for subclasses.
    if(!get_attr_v3(po, &s->sigma_diag, "sigma_diag") ||
       !get_attr_v3(po, &s->sigma_offdiag, "sigma_offdiag")) {

        return 0;
    }

    return 1;
}

static int py_list_to_susceptibility_list(PyObject *po, susceptibility_list *sl) {
    if(!PyList_Check(po)) {
        PyErr_SetString(PyExc_TypeError, "Expected a list\n");
        return 0;
    }

    int length = PyList_Size(po);
    sl->num_items = length;
    sl->items = new susceptibility_struct[length];

    for(int i = 0; i < length; i++) {
        // TODO(chogan): Account for subclasses
        susceptibility_struct s;
        if(!py_susceptibility_to_susceptibility(PyList_GetItem(po, i), &s)) {
            return 0;
        }
        sl->items[i].sigma_diag = s.sigma_diag;
        sl->items[i].sigma_offdiag = s.sigma_offdiag;
    }

    return 1;
}

static int pymaterial_to_material(PyObject *po, material_type *mt) {

    material_data *md = new material_data();
    md->which_subclass = material_data::MEDIUM;
    md->user_func = 0;
    md->user_data = 0;
    md->medium = new medium_struct();

    if(!get_attr_v3(po, &md->medium->epsilon_diag, "epsilon_diag") ||
       !get_attr_v3(po, &md->medium->epsilon_offdiag, "epsilon_offdiag") ||
       !get_attr_v3(po, &md->medium->mu_diag, "mu_diag") ||
       !get_attr_v3(po, &md->medium->mu_offdiag, "mu_offdiag")) {

        return 0;
    }

    PyObject *py_e_susceptibilities = PyObject_GetAttrString(po, "E_susceptibilities");
    PyObject *py_h_susceptibilities = PyObject_GetAttrString(po, "H_susceptibilities");

    if(!py_e_susceptibilities || !py_h_susceptibilities) {
        return 0;
    }

    if(!py_list_to_susceptibility_list(py_e_susceptibilities, &md->medium->E_susceptibilities) ||
       !py_list_to_susceptibility_list(py_h_susceptibilities, &md->medium->H_susceptibilities)) {

        return 0;
    }

    Py_XDECREF(py_e_susceptibilities);
    Py_XDECREF(py_h_susceptibilities);

    if(!get_attr_v3(po, &md->medium->E_chi2_diag, "E_chi2_diag") ||
       !get_attr_v3(po, &md->medium->E_chi3_diag, "E_chi3_diag") ||
       !get_attr_v3(po, &md->medium->H_chi2_diag, "H_chi2_diag") ||
       !get_attr_v3(po, &md->medium->H_chi3_diag, "H_chi3_diag") ||
       !get_attr_v3(po, &md->medium->D_conductivity_diag, "D_conductivity_diag") ||
       !get_attr_v3(po, &md->medium->B_conductivity_diag, "B_conductivity_diag")) {

        return 0;
    }

    mt->data = (void *)md;

    return 1;
}

static int pysphere_to_sphere(PyObject *py_sphere, geometric_object *go) {

    material_type material;
    vector3 center;
    double radius;
    
    if(!get_attr_v3(py_sphere, &center, "center") ||
       !get_attr_dbl(py_sphere, &radius, "radius") ||
       !get_attr_material(py_sphere, &material)) {

        return 0;
    }

    *go = make_sphere(material, center, radius);

    return 1;
}

static int pycylinder_to_cylinder(PyObject *py_cyl, geometric_object *cyl) {
    material_type material;
    vector3 center, axis;
    double radius, height;

    if(!get_attr_v3(py_cyl, &center, "center") ||
       !get_attr_v3(py_cyl, &axis, "axis") ||
       !get_attr_dbl(py_cyl, &radius, "radius") ||
       !get_attr_dbl(py_cyl, &height, "height") ||
       !get_attr_material(py_cyl, &material)) {

        return 0;
    }

    *cyl = make_cylinder(material, center, radius, height, axis);

    return 1;
}

static int pywedge_to_wedge(PyObject *py_wedge, geometric_object *wedge) {
    geometric_object cyl;
    if(!pycylinder_to_cylinder(py_wedge, &cyl)) {
        return 0;
    }

    double wedge_angle; 
    vector3 wedge_start;

    if(!get_attr_dbl(py_wedge, &wedge_angle, "wedge_angle") ||
       !get_attr_v3(py_wedge, &wedge_start, "wedge_start")) {

        return 0;
    }

    double radius = cyl.subclass.cylinder_data->radius;
    double height = cyl.subclass.cylinder_data->height;
    vector3 axis = cyl.subclass.cylinder_data->axis;

    *wedge = make_wedge(cyl.material, cyl.center, radius, height, axis, wedge_angle, wedge_start);

    geometric_object_destroy(cyl);

    return 1;
}

static int pycone_to_cone(PyObject *py_cone, geometric_object *cone) {
    geometric_object cyl;
    if(!pycylinder_to_cylinder(py_cone, &cyl)) {
        return 0;
    }

    double radius2;
    if(!get_attr_dbl(py_cone, &radius2, "radius2")) {
        return 0;
    }

    double radius = cyl.subclass.cylinder_data->radius;
    double height = cyl.subclass.cylinder_data->height;
    vector3 axis = cyl.subclass.cylinder_data->axis;

    *cone = make_cone(cyl.material, cyl.center, radius, height, axis, radius2);

    geometric_object_destroy(cyl);

    return 1;
}

static int pyblock_to_block(PyObject *py_blk, geometric_object *blk) {
    material_type material;
    vector3 center, e1, e2, e3, size;

    if(!get_attr_material(py_blk, &material) ||
       !get_attr_v3(py_blk, &center, "center") ||
       !get_attr_v3(py_blk, &e1, "e1") ||
       !get_attr_v3(py_blk, &e2, "e2") ||
       !get_attr_v3(py_blk, &e3, "e3") ||
       !get_attr_v3(py_blk, &size, "size")) {

        return 0;
    }

    *blk = make_block(material, center, e1, e2, e3, size);
    return 1;
}

static int pyellipsoid_to_ellipsoid(PyObject *py_ell, geometric_object *e) {
    geometric_object blk;
    if(!pyblock_to_block(py_ell, &blk)) {
        return 0;
    }

    material_type material = blk.material;
    vector3 center = blk.center;
    vector3 e1 = blk.subclass.block_data->e1;
    vector3 e2 = blk.subclass.block_data->e2;
    vector3 e3 = blk.subclass.block_data->e3;
    vector3 size = blk.subclass.block_data->size;

    *e = make_ellipsoid(material, center, e1, e2, e3, size);

    geometric_object_destroy(blk);

    return 1;
}

static int py_list_to_gobj_list(PyObject *po, geometric_object_list *l) {
    if(!PyList_Check(po)) {
        PyErr_SetString(PyExc_TypeError, "Expected a list");
        return 0;
    }

    int length = PyList_Size(po);

    l->num_items = length;
    l->items = new geometric_object[length];

    for(int i = 0; i < length; i++) {
        PyObject *py_gobj = PyList_GetItem(po, i);
        geometric_object go;
        if(!py_gobj_to_gobj(py_gobj, &go)) {
            return 0;
        }
        geometric_object_copy(&go, &l->items[i]);
    }

    return 1;
}

static int pycgo_to_cgo(PyObject *py_cgo, geometric_object *o) {
    vector3 center;
    material_type material;

    if(!get_attr_v3(py_cgo, &center, "center") || !get_attr_material(py_cgo, &material)) {
        return 0;
    }

    *o = make_geometric_object(material, center);

    o->which_subclass = geometric_object::COMPOUND_GEOMETRIC_OBJECT;
    o->subclass.compound_geometric_object_data = new compound_geometric_object();

    PyObject *py_component_objects = PyObject_GetAttrString(py_cgo, "component_objects");

    if(!py_component_objects) {
        return 0;
    }

    if(!py_list_to_gobj_list(py_component_objects, &o->subclass.compound_geometric_object_data->component_objects)) {
        return 0;
    }

    Py_XDECREF(py_component_objects);

    return 1;
}

static int py_gobj_to_gobj(PyObject *po, geometric_object *o) {
    int success = 0;
    PyObject *py_type = PyObject_Type(po);
    PyObject *name = PyObject_GetAttrString(py_type, "__name__");

    if(!py_type || !name) {
        return 0;
    }
    
    char *bytes;

    if(PyUnicode_Check(name)) {
        // Python 3
        Py_ssize_t size;
        bytes = PyUnicode_AsUTF8AndSize(name, &size);
        if(!size) {
            return 0;
        }
    } else {
        // Python 2
        bytes = PyBytes_AsString(name);
    }

    if(!bytes) {
        return 0;
    }

    std::string go_type(bytes);

    if(go_type == "Sphere") {
        success = pysphere_to_sphere(po, o);
    }
    else if(go_type == "Cylinder") {
        success = pycylinder_to_cylinder(po, o);
    }
    else if(go_type == "Wedge") {
        success = pywedge_to_wedge(po, o);
    }
    else if(go_type == "Cone") {
        success = pycone_to_cone(po, o);
    }
    else if(go_type == "Block") {
        success = pyblock_to_block(po, o);
    }
    else if(go_type == "Ellipsoid") {
        success = pyellipsoid_to_ellipsoid(po, o);
    }
    else if(go_type == "CompoundGeometricObject") {
        success = pycgo_to_cgo(po, o);
    }
    else {
        PyErr_Format(PyExc_TypeError, "Error: %s is not a valid GeometricObject type\n", go_type.c_str());
        return 0;
    }
    Py_XDECREF(py_type);
    Py_XDECREF(name);

    return success;
}
