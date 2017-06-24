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


static vector3 pyv3_to_v3(PyObject *po) {
    PyObject *py_x = PyObject_GetAttrString(po, "x");
    PyObject *py_y = PyObject_GetAttrString(po, "y");
    PyObject *py_z = PyObject_GetAttrString(po, "z");
    double x = PyFloat_AsDouble(py_x);
    double y = PyFloat_AsDouble(py_y);
    double z = PyFloat_AsDouble(py_z);
    Py_DECREF(py_x);
    Py_DECREF(py_y);
    Py_DECREF(py_z);

    vector3 v = {x, y, z};
    return v;
}

static vector3 get_attr_v3(PyObject *py_obj, const char *name) {
    PyObject *py_attr = PyObject_GetAttrString(py_obj, name);

    if(!py_attr) {
        PyErr_Format(PyExc_ValueError, "Class attribute '%s' is None\n", name);
        // TODO(chogan): Pass error to wrapper, then to python
    }

    vector3 result = pyv3_to_v3(py_attr);
    Py_XDECREF(py_attr);
    return result;
}

static double get_attr_dbl(PyObject *py_obj, const char *name) {
    PyObject *py_attr = PyObject_GetAttrString(py_obj, name);

    if(!py_attr) {
        PyErr_Format(PyExc_ValueError, "Class attribute '%s' is None\n", name);
        // TODO(chogan): Pass error to wrapper, then to python
    }

    double result = PyFloat_AsDouble(py_attr);
    Py_XDECREF(py_attr);
    return result;
}

static material_type get_attr_material(PyObject *po) {
    PyObject *py_material = PyObject_GetAttrString(po, "material");

    if(!py_material) {
        PyErr_SetString(PyExc_ValueError, "Object's material is not set\n");
        // TODO(chogan): Pass error to wrapper
    }

    material_type m = pymaterial_to_material(py_material);
    Py_XDECREF(py_material);

    return m;
}

static susceptibility_struct py_susceptibility_to_susceptibility(PyObject *po) {
    // TODO(chogan): Accont for subclasses. Are all subclasses needed? Several are duplicates
    susceptibility_struct s;
    s.sigma_diag = get_attr_v3(po, "sigma_diag");
    s.sigma_offdiag = get_attr_v3(po, "sigma_offdiag");

    return s;
}

static susceptibility_list py_list_to_susceptibility_list(PyObject *po) {
    if(!PyList_Check(po)) {
        PyErr_SetString(PyExc_TypeError, "Expected a list\n");
        // TODO(chogan): Pass error to wrapper, then to python
    }

    susceptibility_list sl;
    sl.num_items = PyList_Size(po);
    sl.items = new susceptibility_struct[sl.num_items];

    for(int i = 0; i < sl.num_items; i++) {
        // TODO(chogan): Account for subclasses
        susceptibility_struct s = py_susceptibility_to_susceptibility(PyList_GetItem(po, i));
        sl.items[i].sigma_diag = s.sigma_diag;
        sl.items[i].sigma_offdiag = s.sigma_offdiag;
    }

    return sl;
}

static material_type pymaterial_to_material(PyObject *po) {

    material_data *md = new material_data();
    md->which_subclass = material_data::MEDIUM;
    md->user_func = 0;
    md->user_data = 0;
    md->medium = new medium_struct();

    md->medium->epsilon_diag = get_attr_v3(po, "epsilon_diag");
    md->medium->epsilon_offdiag = get_attr_v3(po, "epsilon_offdiag");
    md->medium->mu_diag = get_attr_v3(po, "mu_diag");
    md->medium->mu_offdiag = get_attr_v3(po, "mu_offdiag");

    PyObject *py_e_susceptibilities = PyObject_GetAttrString(po, "E_susceptibilities");
    PyObject *py_h_susceptibilities = PyObject_GetAttrString(po, "H_susceptibilities");

    md->medium->E_susceptibilities = py_list_to_susceptibility_list(py_e_susceptibilities);
    md->medium->H_susceptibilities = py_list_to_susceptibility_list(py_h_susceptibilities);

    Py_XDECREF(py_e_susceptibilities);
    Py_XDECREF(py_h_susceptibilities);

    md->medium->E_chi2_diag = get_attr_v3(po, "E_chi2_diag");
    md->medium->E_chi3_diag = get_attr_v3(po, "E_chi3_diag");
    md->medium->H_chi2_diag = get_attr_v3(po, "H_chi2_diag");
    md->medium->H_chi3_diag = get_attr_v3(po, "H_chi3_diag");
    md->medium->D_conductivity_diag = get_attr_v3(po, "D_conductivity_diag");
    md->medium->B_conductivity_diag = get_attr_v3(po, "B_conductivity_diag");

    material_type mt = { (void *)md };

    return mt;
}

static geometric_object pysphere_to_sphere(PyObject *py_sphere) {
    material_type material = get_attr_material(py_sphere);
    vector3 center = get_attr_v3(py_sphere, "center");
    double radius = get_attr_dbl(py_sphere, "radius");

    return make_sphere(material, center, radius);
}

static geometric_object pycylinder_to_cylinder(PyObject *py_cyl) {
    material_type material = get_attr_material(py_cyl);
    vector3 center = get_attr_v3(py_cyl, "center");
    vector3 axis = get_attr_v3(py_cyl, "axis");
    double radius = get_attr_dbl(py_cyl, "radius");
    double height = get_attr_dbl(py_cyl, "height");

    return make_cylinder(material, center, radius, height, axis);
}

static geometric_object pywedge_to_wedge(PyObject *py_wedge) {
    geometric_object cyl = pycylinder_to_cylinder(py_wedge);
    double wedge_angle = get_attr_dbl(py_wedge, "wedge_angle");
    vector3 wedge_start = get_attr_v3(py_wedge, "wedge_start");

    double radius = cyl.subclass.cylinder_data->radius;
    double height = cyl.subclass.cylinder_data->height;
    vector3 axis = cyl.subclass.cylinder_data->axis;

    geometric_object w = make_wedge(cyl.material, cyl.center, radius, height, axis, wedge_angle, wedge_start);

    geometric_object_destroy(cyl);

    return w;
}

static geometric_object pycone_to_cone(PyObject *py_cone) {
    geometric_object cyl = pycylinder_to_cylinder(py_cone);
    double radius2 = get_attr_dbl(py_cone, "radius2");
    double radius = cyl.subclass.cylinder_data->radius;
    double height = cyl.subclass.cylinder_data->height;
    vector3 axis = cyl.subclass.cylinder_data->axis;

    geometric_object c = make_cone(cyl.material, cyl.center, radius, height, axis, radius2);

    geometric_object_destroy(cyl);

    return c;
}

static geometric_object pyblock_to_block(PyObject *py_blk) {
    material_type material = get_attr_material(py_blk);
    vector3 center = get_attr_v3(py_blk, "center");
    vector3 e1 = get_attr_v3(py_blk, "e1");
    vector3 e2 = get_attr_v3(py_blk, "e2");
    vector3 e3 = get_attr_v3(py_blk, "e3");
    vector3 size = get_attr_v3(py_blk, "size");

    return make_block(material, center, e1, e2, e3, size);
}

static geometric_object pyellipsoid_to_ellipsoid(PyObject *py_ell) {
    geometric_object blk = pyblock_to_block(py_ell);

    material_type material = blk.material;
    vector3 center = blk.center;
    vector3 e1 = blk.subclass.block_data->e1;
    vector3 e2 = blk.subclass.block_data->e2;
    vector3 e3 = blk.subclass.block_data->e3;
    vector3 size = blk.subclass.block_data->size;

    geometric_object ellipsoid = make_ellipsoid(material, center, e1, e2, e3, size);

    geometric_object_destroy(blk);

    return ellipsoid;
}

static geometric_object py_gobj_to_gobj(PyObject *po) {
    geometric_object o;
    PyObject *py_type = PyObject_Type(po);
    PyObject *name = PyObject_GetAttrString(py_type, "__name__");

    std::string go_type(PyString_AsString(name));

    if(go_type == "Sphere") {
        o = pysphere_to_sphere(po);
    }
    else if(go_type == "Cylinder") {
        o = pycylinder_to_cylinder(po);
    }
    else if(go_type == "Wedge") {
        o = pywedge_to_wedge(po);
    }
    else if(go_type == "Cone") {
        o = pycone_to_cone(po);
    }
    else if(go_type == "Block") {
        o = pyblock_to_block(po);
    }
    else if(go_type == "Ellipsoid") {
        o = pyellipsoid_to_ellipsoid(po);
    }
    else if(go_type == "CompoundGeometricObject") {
        o = pycgo_to_cgo(po);
    }
    else {
        // TODO(chogan): Exception
        PyErr_Format(PyExc_TypeError, "Error: %s is not a valid GeometricObject type\n", go_type.c_str());
    }
    Py_XDECREF(py_type);
    Py_XDECREF(name);

    return o;
}

static geometric_object_list py_list_to_gobj_list(PyObject *po) {
    if(!PyList_Check(po)) {
        PyErr_SetString(PyExc_TypeError, "Expected a list");
        // TODO(chogan): Pass error to wrapper
    }

    int length = PyList_Size(po);

    geometric_object_list l;
    l.num_items = length;
    l.items = new geometric_object[length];

    for(int i = 0; i < length; i++) {
        PyObject *py_gobj = PyList_GetItem(po, i);
        geometric_object go = py_gobj_to_gobj(py_gobj);
        geometric_object_copy(&go, l.items + i);
    }

    return l;
}

static geometric_object pycgo_to_cgo(PyObject *py_cgo) {
    vector3 center = get_attr_v3(py_cgo, "center");
    material_type material = get_attr_material(py_cgo);

    geometric_object o = make_geometric_object(material, center);

    o.which_subclass = geometric_object::COMPOUND_GEOMETRIC_OBJECT;
    o.subclass.compound_geometric_object_data = new compound_geometric_object();

    PyObject *py_component_objects = PyObject_GetAttrString(py_cgo, "component_objects");
    o.subclass.compound_geometric_object_data->component_objects = py_list_to_gobj_list(py_component_objects);

    Py_XDECREF(py_component_objects);

    return o;
}
