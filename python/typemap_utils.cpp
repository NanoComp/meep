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

#if PY_MAJOR_VERSION >= 3
    #define PyObject_ToCharPtr(n) PyUnicode_AsUTF8(n)
    #define IsPyString(n) PyUnicode_Check(n)
    #define PyInteger_Check(n) PyLong_Check(n)
    #define PyInteger_AsLong(n) PyLong_AsLong(n)
    #define PyInteger_FromLong(n) PyLong_FromLong(n)
#else
    #define PyObject_ToCharPtr(n) py2_string_as_utf8(n)
    #define IsPyString(n) PyString_Check(n) || PyUnicode_Check(n)
    #define PyInteger_Check(n) PyInt_Check(n)
    #define PyInteger_AsLong(n) PyInt_AsLong(n)
    #define PyInteger_FromLong(n) PyInt_FromLong(n)
#endif

PyObject *py_callback = NULL;
PyObject *py_callback_v3 = NULL;
PyObject *py_amp_func = NULL;

static int pymedium_to_medium(PyObject *po, medium_struct *m);
static int pymaterial_to_material(PyObject *po, material_type *mt);

#if PY_MAJOR_VERSION == 2
static char *py2_string_as_utf8(PyObject *po) {
    if (PyString_Check(po)) {
        return PyString_AsString(po);
    }
    else if (PyUnicode_Check(po)) {
        PyObject *s = PyUnicode_AsUTF8String(po);
        char *result = PyString_AsString(s);
        Py_DECREF(s);
        return result;
    }
    else {
        return NULL;
    }
}
#endif

static PyObject *get_geom_mod() {
    static PyObject *geom_mod = NULL;
    if (geom_mod == NULL) {
        geom_mod = PyImport_ImportModule("meep.geom");
    }
    return geom_mod;
}

static PyObject *py_material_object() {
    static PyObject *material_object = NULL;
    if (material_object == NULL) {
        PyObject *geom_mod = get_geom_mod();
        material_object = PyObject_GetAttrString(geom_mod, "Medium");
    }
    return material_object;
}

static PyObject *py_vector3_object() {
    static PyObject *vector3_object = NULL;
    if (vector3_object == NULL) {
        PyObject *geom_mod = get_geom_mod();
        vector3_object = PyObject_GetAttrString(geom_mod, "Vector3");
    }
    return vector3_object;
}

static PyObject* vec2py(const meep::vec &v, bool newobj=false) {

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
        x = v.r();
        z = v.z();
        break;
    }

    if (newobj) {
        PyObject *v3_class = py_vector3_object();
        PyObject *args = Py_BuildValue("(d,d,d)", x, y, z);
        PyObject *res = PyObject_Call(v3_class, args, NULL);

        Py_DECREF(args);

        return res;
    }
    else {
        if (py_callback_v3 == NULL) {
            PyObject *v3_class = py_vector3_object();
            PyObject *args = PyTuple_New(0);
            py_callback_v3 = PyObject_Call(v3_class, args, NULL);

            Py_DECREF(args);
        }

        PyObject *pyx = PyFloat_FromDouble(x);
        PyObject *pyy = PyFloat_FromDouble(y);
        PyObject *pyz = PyFloat_FromDouble(z);

        PyObject_SetAttrString(py_callback_v3, "x", pyx);
        PyObject_SetAttrString(py_callback_v3, "y", pyy);
        PyObject_SetAttrString(py_callback_v3, "z", pyz);

        Py_DECREF(pyx);
        Py_DECREF(pyy);
        Py_DECREF(pyz);

        return py_callback_v3;
    }
}

static void py_user_material_func_wrap(vector3 x, void *user_data, medium_struct *medium) {
    PyObject *py_vec = vec2py(vector3_to_vec(x));

    PyObject *pyret = PyObject_CallFunctionObjArgs((PyObject *)user_data, py_vec, NULL);

    if (!pyret) {
        PyErr_PrintEx(0);
    }

    if (!pymedium_to_medium(pyret, medium)) {
        PyErr_PrintEx(0);
    }

    Py_DECREF(pyret);
}

static void py_epsilon_func_wrap(vector3 x, void *user_data, medium_struct *medium) {
    PyObject *py_vec = vec2py(vector3_to_vec(x));

    PyObject *pyret = PyObject_CallFunctionObjArgs((PyObject *)user_data, py_vec, NULL);

    if (!pyret) {
        PyErr_PrintEx(0);
    }

    double eps = PyFloat_AsDouble(pyret);

    medium->epsilon_diag.x = eps;
    medium->epsilon_diag.y = eps;
    medium->epsilon_diag.z = eps;

    Py_DECREF(pyret);
}

static std::string py_class_name_as_string(PyObject *po) {
    PyObject *py_type = PyObject_Type(po);
    PyObject *name = PyObject_GetAttrString(py_type, "__name__");

    const char *bytes = PyObject_ToCharPtr(name);

    std::string class_name(bytes);

    Py_XDECREF(py_type);
    Py_XDECREF(name);

    return class_name;
}

static int pyv3_to_v3(PyObject *po, vector3 *v) {

    PyObject *py_x = PyObject_GetAttrString(po, "x");
    PyObject *py_y = PyObject_GetAttrString(po, "y");
    PyObject *py_z = PyObject_GetAttrString(po, "z");

    if (!py_x || !py_y || !py_z) {
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

static int pyv3_to_cv3(PyObject *po, cvector3 *v) {
    PyObject *py_x = PyObject_GetAttrString(po, "x");
    PyObject *py_y = PyObject_GetAttrString(po, "y");
    PyObject *py_z = PyObject_GetAttrString(po, "z");

    if (!py_x || !py_y || !py_z) {
        PyErr_SetString(PyExc_ValueError, "Vector3 is not initialized");
        return 0;
    }

    std::complex<double> x = std::complex<double>(PyComplex_RealAsDouble(py_x),
                                                  PyComplex_ImagAsDouble(py_x));
    std::complex<double> y = std::complex<double>(PyComplex_RealAsDouble(py_y),
                                                  PyComplex_ImagAsDouble(py_y));
    std::complex<double> z = std::complex<double>(PyComplex_RealAsDouble(py_z),
                                                  PyComplex_ImagAsDouble(py_z));
    Py_DECREF(py_x);
    Py_DECREF(py_y);
    Py_DECREF(py_z);

    v->x.re = x.real();
    v->x.im = x.imag();
    v->y.re = y.real();
    v->y.im = y.imag();
    v->z.re = z.real();
    v->z.im = z.imag();

    return 1;
}

static PyObject* v3_to_pyv3(vector3 *v) {
    PyObject *v3_class = py_vector3_object();
    PyObject *args = Py_BuildValue("(ddd)", v->x, v->y, v->z);
    PyObject *py_v = PyObject_Call(v3_class, args, NULL);

    Py_DECREF(args);

    return py_v;
}

static int get_attr_v3(PyObject *py_obj, vector3 *v, const char *name) {
    PyObject *py_attr = PyObject_GetAttrString(py_obj, name);

    if (!py_attr) {
        PyErr_Format(PyExc_ValueError, "Class attribute '%s' is None\n", name);
        return 0;
    }

    if (!pyv3_to_v3(py_attr, v)) {
        return 0;
    }

    Py_XDECREF(py_attr);
    return 1;
}

static int get_attr_v3_cmplx(PyObject *py_obj, cvector3 *v, const char *name) {

    PyObject *py_attr = PyObject_GetAttrString(py_obj, name);

    if (!py_attr) {
        PyErr_Format(PyExc_ValueError, "Class attribute '%s' is None\n", name);
        return 0;
    }

    if (!pyv3_to_cv3(py_attr, v)) {
        return 0;
    }

    Py_XDECREF(py_attr);
    return 1;
}

static int get_attr_dbl(PyObject *py_obj, double *result, const char *name) {
    PyObject *py_attr = PyObject_GetAttrString(py_obj, name);

    if (!py_attr) {
        PyErr_Format(PyExc_ValueError, "Class attribute '%s' is None\n", name);
        return 0;
    }

    *result = PyFloat_AsDouble(py_attr);
    Py_XDECREF(py_attr);
    return 1;
}

static int get_attr_int(PyObject *py_obj, int *result, const char *name) {
    PyObject *py_attr = PyObject_GetAttrString(py_obj, name);

    if (!py_attr) {
        PyErr_Format(PyExc_ValueError, "Class attribute '%s' is None\n", name);
        return 0;
    }

    *result = PyInteger_AsLong(py_attr);
    Py_XDECREF(py_attr);
    return 1;
}

static int get_attr_material(PyObject *po, material_type *m) {
    PyObject *py_material = PyObject_GetAttrString(po, "material");

    if (!py_material) {
        PyErr_SetString(PyExc_ValueError, "Object's material is not set\n");
        return 0;
    }

    if (!pymaterial_to_material(py_material, m)) {
        return 0;
    }

    Py_XDECREF(py_material);

    return 1;
}

static int pytransition_to_transition(PyObject *py_trans, transition *trans) {

    int from, to;
    double trans_rate, freq, gamma, pump_rate;
    vector3 sigma_diag;

    if (!get_attr_int(py_trans, &from, "from_level")            ||
        !get_attr_int(py_trans, &to, "to_level")                ||
        !get_attr_dbl(py_trans, &trans_rate, "transition_rate") ||
        !get_attr_dbl(py_trans, &freq, "frequency")             ||
        !get_attr_dbl(py_trans, &gamma, "gamma")                ||
        !get_attr_dbl(py_trans, &pump_rate, "pumping_rate")     ||
        !get_attr_v3(py_trans, &sigma_diag, "sigma_diag")) {

        return 0;
    }

    trans->from_level = from;
    trans->to_level = to;
    trans->transition_rate = trans_rate;
    trans->frequency = freq;
    trans->gamma = gamma;
    trans->pumping_rate = pump_rate;
    trans->sigma_diag.x = sigma_diag.x;
    trans->sigma_diag.y = sigma_diag.y;
    trans->sigma_diag.z = sigma_diag.z;

    return 1;
}

static int py_susceptibility_to_susceptibility(PyObject *po, susceptibility_struct *s) {
    if (!get_attr_v3(po, &s->sigma_diag, "sigma_diag") ||
        !get_attr_v3(po, &s->sigma_offdiag, "sigma_offdiag")) {

        return 0;
    }

    s->frequency = 0;
    s->gamma = 0;
    s->noise_amp = 0;
    s->transitions.resize(0);
    s->initial_populations.resize(0);

    if (PyObject_HasAttrString(po, "frequency")) {
        if(!get_attr_dbl(po, &s->frequency, "frequency")) {
            return 0;
        }
    }

    if (PyObject_HasAttrString(po, "gamma")) {
        if(!get_attr_dbl(po, &s->gamma, "gamma")) {
            return 0;
        }
    }

    if (PyObject_HasAttrString(po, "noise_amp")) {
        if(!get_attr_dbl(po, &s->noise_amp, "noise_amp")) {
            return 0;
        }
    }

    if (PyObject_HasAttrString(po, "transitions")) {
        // MultilevelAtom
        PyObject *py_trans = PyObject_GetAttrString(po, "transitions");
        if (!py_trans) {
            return 0;
        }
        int length = PyList_Size(py_trans);
        s->transitions.resize(length);

        for (int i = 0; i < length; ++i) {
            if (!pytransition_to_transition(PyList_GetItem(py_trans, i), &s->transitions[i])) {
                return 0;
            }
        }
        Py_DECREF(py_trans);

        PyObject *py_pop = PyObject_GetAttrString(po, "initial_populations");
        if (!py_pop) {
            return 0;
        }
        length = PyList_Size(py_pop);
        s->initial_populations.resize(length);

        for (int i = 0; i < length; ++i) {
            s->initial_populations[i] = PyFloat_AsDouble(PyList_GetItem(py_pop, i));
        }
        Py_DECREF(py_pop);
    }

    std::string class_name = py_class_name_as_string(po);

    if (class_name.find(std::string("Drude")) != std::string::npos) {
        s->drude = true;
    } else {
        s->drude = false;
    }

    s->is_file = false;

    return 1;
}

static int py_list_to_susceptibility_list(PyObject *po, susceptibility_list *sl) {
    if (!PyList_Check(po)) {
        PyErr_SetString(PyExc_TypeError, "Expected a list\n");
        return 0;
    }

    int length = PyList_Size(po);
    sl->num_items = length;
    sl->items = new susceptibility_struct[length];

    for(int i = 0; i < length; i++) {
        if (!py_susceptibility_to_susceptibility(PyList_GetItem(po, i), &sl->items[i])) {
            return 0;
        }
    }

    return 1;
}

static int pymaterial_to_material(PyObject *po, material_type *mt) {
    material_data *md;

    if (PyObject_IsInstance(po, py_material_object())) {
        md = make_dielectric(1);
        if (!pymedium_to_medium(po, &md->medium)) {
            return 0;
        }
    } else if (PyFunction_Check(po)) {
        PyObject *eps = PyObject_GetAttrString(po, "eps");
        if (!eps) {
            return 0;
        }
        if (eps == Py_True) {
            md = make_user_material(py_epsilon_func_wrap, po);
        } else {
            md = make_user_material(py_user_material_func_wrap, po);
        }
        Py_DECREF(eps);
    } else if (IsPyString(po)) {
        const char *eps_input_file = PyObject_ToCharPtr(po);
        md = make_file_material(eps_input_file);
    } else if (PyArray_Check(po)) {
        PyArrayObject *pao = (PyArrayObject*)po;

        if (!PyArray_ISCARRAY(pao)) {
            PyErr_SetString(PyExc_TypeError, "Numpy array must be C-style contiguous.");
            return 0;
        }
        md = new material_data();
        md->which_subclass=material_data::MATERIAL_FILE;
        md->epsilon_dims[0] = md->epsilon_dims[1] = md->epsilon_dims[2] = 1;
        md->epsilon_data = new realnum[PyArray_SIZE(pao)];
        memcpy(md->epsilon_data, (realnum*)PyArray_DATA(pao), PyArray_SIZE(pao) * sizeof(realnum));

        for (int i = 0; i < PyArray_NDIM(pao); ++i) {
            md->epsilon_dims[i] = (size_t)PyArray_DIMS(pao)[i];
        }

        master_printf("read in %zdx%zdx%zd numpy array for epsilon\n",
                      md->epsilon_dims[0],
                      md->epsilon_dims[1],
                      md->epsilon_dims[2]);
    } else {
        PyErr_SetString(PyExc_TypeError, "Expected a Medium, a function, or a filename");
        return 0;
    }

    *mt = md;

    return 1;
}

template<class T>
static void set_v3_on_pyobj(PyObject *py_obj, T *v3, const char *attr) {
    PyObject *v3_class = py_vector3_object();
    PyObject *v3_args = Py_BuildValue("(d,d,d)", v3->x, v3->y, v3->z);
    PyObject *pyv3 = PyObject_Call(v3_class, v3_args, NULL);
    PyObject_SetAttrString(py_obj, attr, pyv3);

    Py_DECREF(v3_args);
    Py_DECREF(pyv3);
}

static PyObject *susceptibility_to_py_obj(susceptibility_struct *s) {
    PyObject *geom_mod = get_geom_mod();

    PyObject *res;
    PyObject *args = PyTuple_New(0);

    if (s->noise_amp == 0) {
        if (s->drude) {
            PyObject *py_drude_class = PyObject_GetAttrString(geom_mod, "DrudeSusceptibility");
            res = PyObject_Call(py_drude_class, args, NULL);
            Py_DECREF(py_drude_class);
        }
        else {
            PyObject *py_lorentz_class = PyObject_GetAttrString(geom_mod, "LorentzianSusceptibility");
            res = PyObject_Call(py_lorentz_class, args, NULL);
            Py_DECREF(py_lorentz_class);
        }
    }
    else {
        if (s->drude) {
            PyObject *py_noisy_drude_class = PyObject_GetAttrString(geom_mod, "NoisyDrudeSusceptibility");
            res = PyObject_Call(py_noisy_drude_class, args, NULL);
            Py_DECREF(py_noisy_drude_class);
        }
        else {
            PyObject *py_noisy_lorentz_class = PyObject_GetAttrString(geom_mod, "NoisyLorentzianSusceptibility");
            res = PyObject_Call(py_noisy_lorentz_class, args, NULL);
            Py_DECREF(py_noisy_lorentz_class);
        }
        PyObject *py_noise = PyFloat_FromDouble(s->noise_amp);
        PyObject_SetAttrString(res, "noise_amp", py_noise);
        Py_DECREF(py_noise);
    }

    set_v3_on_pyobj<vector3>(res, &s->sigma_diag, "sigma_diag");
    set_v3_on_pyobj<vector3>(res, &s->sigma_offdiag, "sigma_offdiag");

    PyObject *py_freq = PyFloat_FromDouble(s->frequency);
    PyObject *py_gamma = PyFloat_FromDouble(s->gamma);

    PyObject_SetAttrString(res, "frequency", py_freq);
    PyObject_SetAttrString(res, "gamma", py_gamma);

    Py_DECREF(args);
    Py_DECREF(py_freq);
    Py_DECREF(py_gamma);

    return res;
}

static PyObject *susceptibility_list_to_py_list(susceptibility_list *sl) {
    PyObject *res = PyList_New(sl->num_items);

    for (Py_ssize_t i = 0; i < sl->num_items; ++i) {
        PyList_SetItem(res, i, susceptibility_to_py_obj(&sl->items[i]));
    }

    return res;
}

static PyObject *material_to_py_material(material_type mat) {
    switch (mat->which_subclass) {
    case meep_geom::material_data::MEDIUM: {
        PyObject *geom_mod = get_geom_mod();
        PyObject *medium_class = PyObject_GetAttrString(geom_mod, "Medium");

        PyObject *medium_args = PyTuple_New(0);
        PyObject *py_mat = PyObject_Call(medium_class, medium_args, NULL);

        PyObject *py_E_sus = susceptibility_list_to_py_list(&mat->medium.E_susceptibilities);
        PyObject *py_H_sus = susceptibility_list_to_py_list(&mat->medium.H_susceptibilities);
        PyObject_SetAttrString(py_mat, "E_susceptibilities", py_E_sus);
        PyObject_SetAttrString(py_mat, "H_susceptibilities", py_H_sus);

        set_v3_on_pyobj<vector3>(py_mat, &mat->medium.epsilon_diag, "epsilon_diag");
        set_v3_on_pyobj<vector3>(py_mat, &mat->medium.mu_diag, "mu_diag");
        set_v3_on_pyobj<vector3>(py_mat, &mat->medium.E_chi2_diag, "E_chi2_diag");
        set_v3_on_pyobj<vector3>(py_mat, &mat->medium.E_chi3_diag, "E_chi3_diag");
        set_v3_on_pyobj<vector3>(py_mat, &mat->medium.H_chi2_diag, "H_chi2_diag");
        set_v3_on_pyobj<vector3>(py_mat, &mat->medium.H_chi3_diag, "H_chi3_diag");
        set_v3_on_pyobj<vector3>(py_mat, &mat->medium.D_conductivity_diag, "D_conductivity_diag");
        set_v3_on_pyobj<vector3>(py_mat, &mat->medium.B_conductivity_diag, "B_conductivity_diag");
        set_v3_on_pyobj<cvector3>(py_mat, &mat->medium.epsilon_offdiag, "epsilon_offdiag");
        set_v3_on_pyobj<cvector3>(py_mat, &mat->medium.mu_offdiag, "mu_offdiag");

        Py_DECREF(medium_args);
        Py_DECREF(medium_class);
        Py_DECREF(py_E_sus);
        Py_DECREF(py_H_sus);

        return py_mat;
    }
    default:
        // Only Medium is supported at this time.
        PyErr_SetString(PyExc_NotImplementedError, "Can only convert C++ medium_struct to Python");
        return NULL;
    }
}

static int pymedium_to_medium(PyObject *po, medium_struct *m) {
    if (!get_attr_v3(po, &m->epsilon_diag, "epsilon_diag") ||
        !get_attr_v3(po, &m->mu_diag, "mu_diag")) {

        return 0;
    }
    if (!get_attr_v3_cmplx(po, &m->mu_offdiag, "mu_offdiag") ||
        !get_attr_v3_cmplx(po, &m->epsilon_offdiag, "epsilon_offdiag")) {

        return 0;
    }

    PyObject *py_e_susceptibilities = PyObject_GetAttrString(po, "E_susceptibilities");
    PyObject *py_h_susceptibilities = PyObject_GetAttrString(po, "H_susceptibilities");

    if (!py_e_susceptibilities || !py_h_susceptibilities) {
        return 0;
    }

    if (!py_list_to_susceptibility_list(py_e_susceptibilities, &m->E_susceptibilities) ||
       !py_list_to_susceptibility_list(py_h_susceptibilities, &m->H_susceptibilities)) {

        return 0;
    }

    Py_XDECREF(py_e_susceptibilities);
    Py_XDECREF(py_h_susceptibilities);

    if (!get_attr_v3(po, &m->E_chi2_diag, "E_chi2_diag") ||
        !get_attr_v3(po, &m->E_chi3_diag, "E_chi3_diag") ||
        !get_attr_v3(po, &m->H_chi2_diag, "H_chi2_diag") ||
        !get_attr_v3(po, &m->H_chi3_diag, "H_chi3_diag") ||
        !get_attr_v3(po, &m->D_conductivity_diag, "D_conductivity_diag") ||
        !get_attr_v3(po, &m->B_conductivity_diag, "B_conductivity_diag")) {

        return 0;
    }

    return 1;
}

static int pysphere_to_sphere(PyObject *py_sphere, geometric_object *go) {

    material_type material;
    vector3 center;
    double radius;

    if (!get_attr_v3(py_sphere, &center, "center") ||
       !get_attr_dbl(py_sphere, &radius, "radius") ||
       !get_attr_material(py_sphere, &material)) {

        go->subclass.sphere_data = NULL;
        return 0;
    }

    *go = make_sphere(material, center, radius);

    return 1;
}

static int pycylinder_to_cylinder(PyObject *py_cyl, geometric_object *cyl) {
    material_type material;
    vector3 center, axis;
    double radius, height;

    if (!get_attr_v3(py_cyl, &center, "center") ||
       !get_attr_v3(py_cyl, &axis, "axis") ||
       !get_attr_dbl(py_cyl, &radius, "radius") ||
       !get_attr_dbl(py_cyl, &height, "height") ||
       !get_attr_material(py_cyl, &material)) {

        cyl->subclass.cylinder_data = NULL;
        return 0;
    }

    *cyl = make_cylinder(material, center, radius, height, axis);

    return 1;
}

static int pywedge_to_wedge(PyObject *py_wedge, geometric_object *wedge) {
    geometric_object cyl;
    if (!pycylinder_to_cylinder(py_wedge, &cyl)) {
        return 0;
    }

    double wedge_angle;
    vector3 wedge_start;

    if (!get_attr_dbl(py_wedge, &wedge_angle, "wedge_angle") ||
       !get_attr_v3(py_wedge, &wedge_start, "wedge_start")) {

        wedge->subclass.cylinder_data = NULL;
        geometric_object_destroy(cyl);
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
    if (!pycylinder_to_cylinder(py_cone, &cyl)) {
        return 0;
    }

    double radius2;
    if (!get_attr_dbl(py_cone, &radius2, "radius2")) {
        cone->subclass.cylinder_data = NULL;
        geometric_object_destroy(cyl);
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

    if (!get_attr_material(py_blk, &material) ||
       !get_attr_v3(py_blk, &center, "center") ||
       !get_attr_v3(py_blk, &e1, "e1") ||
       !get_attr_v3(py_blk, &e2, "e2") ||
       !get_attr_v3(py_blk, &e3, "e3") ||
       !get_attr_v3(py_blk, &size, "size")) {

        blk->subclass.block_data = NULL;
        return 0;
    }

    *blk = make_block(material, center, e1, e2, e3, size);
    return 1;
}

static int pyellipsoid_to_ellipsoid(PyObject *py_ell, geometric_object *e) {
    geometric_object blk;
    if (!pyblock_to_block(py_ell, &blk)) {
        return 0;
    }

    material_type material = (material_type) blk.material;
    vector3 center = blk.center;
    vector3 e1 = blk.subclass.block_data->e1;
    vector3 e2 = blk.subclass.block_data->e2;
    vector3 e3 = blk.subclass.block_data->e3;
    vector3 size = blk.subclass.block_data->size;

    *e = make_ellipsoid(material, center, e1, e2, e3, size);

    geometric_object_destroy(blk);

    return 1;
}

static int pyprism_to_prism(PyObject *py_prism, geometric_object *p) {
    material_type material;
    double height;
    vector3 axis, center;

    if (!get_attr_material(py_prism, &material)    ||
        !get_attr_dbl(py_prism, &height, "height") ||
        !get_attr_v3(py_prism, &center, "center")  ||
        !get_attr_v3(py_prism, &axis, "axis")) {

        return 0;
    }

    PyObject *py_vert_list = PyObject_GetAttrString(py_prism, "vertices");

    if (!py_vert_list) {
        PyErr_PrintEx(0);
        return 0;
    }

    if (!PyList_Check(py_vert_list)) {
        PyErr_SetString(PyExc_TypeError, "Expected Prism.vertices to be a list\n");
        return 0;
    }

    int num_vertices = PyList_Size(py_vert_list);
    vector3 *vertices = new vector3[num_vertices];

    for (Py_ssize_t i = 0; i < num_vertices; ++i) {
        vector3 v3;
        if (!pyv3_to_v3(PyList_GetItem(py_vert_list, i), &v3)) {
            return 0;
        }
        vertices[i] = v3;
    }

    *p = make_prism(material, vertices, num_vertices, height, axis);
    p->center = center;

    delete [] vertices;
    Py_DECREF(py_vert_list);

    return 1;
}

static int py_gobj_to_gobj(PyObject *po, geometric_object *o) {
    int success = 0;
    std::string go_type = py_class_name_as_string(po);

    if (go_type == "Sphere") {
        success = pysphere_to_sphere(po, o);
    }
    else if (go_type == "Cylinder") {
        success = pycylinder_to_cylinder(po, o);
    }
    else if (go_type == "Wedge") {
        success = pywedge_to_wedge(po, o);
    }
    else if (go_type == "Cone") {
        success = pycone_to_cone(po, o);
    }
    else if (go_type == "Block") {
        success = pyblock_to_block(po, o);
    }
    else if (go_type == "Ellipsoid") {
        success = pyellipsoid_to_ellipsoid(po, o);
    }
    else if (go_type == "Prism") {
        success = pyprism_to_prism(po, o);
    }
    else {
        PyErr_Format(PyExc_TypeError, "Error: %s is not a valid GeometricObject type\n", go_type.c_str());
        return 0;
    }

    return success;
}

static int py_list_to_gobj_list(PyObject *po, geometric_object_list *l) {
    if (!PyList_Check(po)) {
        PyErr_SetString(PyExc_TypeError, "Expected a list");
        return 0;
    }

    int length = PyList_Size(po);

    l->num_items = length;
    l->items = new geometric_object[length];

    for(int i = 0; i < length; i++) {
        PyObject *py_gobj = PyList_GetItem(po, i);
        geometric_object go;
        if (!py_gobj_to_gobj(py_gobj, &go)) {
            return 0;
        }
        geometric_object_copy(&go, &l->items[i]);
    }

    return 1;
}

static PyObject *gobj_to_py_obj(geometric_object *gobj) {
    switch (gobj->which_subclass) {
    case geometric_object::PRISM: {
        PyObject *geom_mod = get_geom_mod();
        PyObject *prism_class = PyObject_GetAttrString(geom_mod, "Prism");

        int num_verts = gobj->subclass.prism_data->vertices.num_items;
        prism *prsm = gobj->subclass.prism_data;

        PyObject *py_verts = PyList_New(num_verts);
        for (int i = 0; i < num_verts; ++i)
         PyList_SetItem(py_verts, i, v3_to_pyv3(prsm->vertices.items + i));

        double height = prsm->height;
        vector3 axis = prsm->axis;
        PyObject *py_axis = v3_to_pyv3(&axis);

        PyObject *py_mat = material_to_py_material((meep_geom::material_type)gobj->material);

        PyObject *args = Py_BuildValue("(OdO)", py_verts, height, py_axis);
        PyObject *kwargs = Py_BuildValue("{s:O}", "material", py_mat);
        PyObject *res = PyObject_Call(prism_class, args, kwargs);

        Py_DECREF(prism_class);
        Py_DECREF(args);
        Py_DECREF(kwargs);
        Py_DECREF(py_verts);
        Py_DECREF(py_axis);
        Py_DECREF(py_mat);

        return res;
    }
    case geometric_object::BLOCK:
    case geometric_object::SPHERE:
    case geometric_object::CYLINDER:
    default:
        // We currently only have the need to create python Prisms from C++.
        // Other geometry can be added as needed.
        PyErr_SetString(PyExc_RuntimeError, "Conversion of non-prism geometric_object to Python is not supported");
        return NULL;
    }
}

static PyObject* gobj_list_to_py_list(geometric_object_list *objs) {

    PyObject *py_res = PyList_New(objs->num_items);

    for (int i = 0; i < objs->num_items; ++i) {
        PyList_SetItem(py_res, i, gobj_to_py_obj(&objs->items[i]));
        geometric_object_destroy(objs->items[i]);
    }

    free(objs->items);

    return py_res;
}
