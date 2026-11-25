

/* Copyright (C) 2005-2026 Massachusetts Institute of Technology
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


// SWIG warnings filtered for the declarations in meep/vec.hpp:
// 322: Redundant redeclaration.  This is ok: the wrappers are created correctly.
// 503: Can't wrap a friend function unless renamed.  We don't need class-specific
//      wrappers for friend functions; the namespace-level declarations are wrapped.
// 509: Overloaded method shadowed by another one that Python can't tell apart.
%warnfilter(322,509) meep::component_direction;
%warnfilter(322,509) meep::direction_component;
%warnfilter(322,503) meep::zero_vec;
%warnfilter(322,503) meep::veccyl;
%warnfilter(322,503) meep::zero_ivec;
%warnfilter(322,503) meep::one_ivec;
%warnfilter(322,503) meep::iveccyl;
%warnfilter(322) meep::abort;
%warnfilter(503) meep::one_vec;
%warnfilter(503) meep::volcyl;
%warnfilter(503) meep::volone;
%warnfilter(503) meep::vol1d;
%warnfilter(503) meep::voltwo;
%warnfilter(503) meep::vol2d;
%warnfilter(503) meep::vol3d;
%warnfilter(503) meep::identity;
%warnfilter(503) meep::rotate4;
%warnfilter(503) meep::rotate2;
%warnfilter(503) meep::mirror;
%warnfilter(503) meep::r_to_minus_r_symmetry;
%warnfilter(509) meep::component_name;
%warnfilter(509) meep::coordinate_mismatch;
%warnfilter(509) meep::ivec::ivec;
%warnfilter(509) meep::symmetry::transform;
%warnfilter(509) meep::symmetry::phase_shift;

// Renaming python builtins
%rename(meep_type) meep::type;
%rename(vec_abs) meep::abs;
%rename(vec_max) meep::max;
%rename(vec_min) meep::min;
%rename(print_grid_volume) meep::grid_volume::print;
%rename(symmetry_reduce) meep::symmetry::reduce;

// Operator renaming
%rename(__contains__) meep::volume::operator&&;
%rename(__getitem__) meep::grid_volume::operator[];
%rename(symmetry_assign) meep::symmetry::operator=;

%rename(vec_dim_val) meep::vec::vec(ndim, double);
%rename(vec_dim) meep::vec::vec(ndim);

// ignore methods that conflict with python built-ins
%ignore meep::vec::str;
%ignore meep::ivec::str;
%ignore meep::volume::str;
%ignore meep::grid_volume::str;


%{
// Append Python's shortest round-trip representation of x (as in repr(float)).
// Returns false, with a Python exception set, on failure.
static bool meep_append_double_repr(std::string &s, double x) {
    char *p = PyOS_double_to_string(x, 'r', 0, Py_DTSF_ADD_DOT_0, NULL);
    if (!p) return false;
    s += p;
    PyMem_Free(p);
    return true;
}

// Same format as meep::vec::str(), but without rounding to six decimal places.
static bool meep_append_vec_repr(std::string &s, const meep::vec &v) {
    static const meep::direction cyl_dirs[] = {meep::R, meep::Z};
    static const meep::direction dirs[] = {meep::X, meep::Y, meep::Z};
    const meep::direction *ds = v.dim == meep::Dcyl ? cyl_dirs : dirs;
    const int n = v.dim == meep::Dcyl ? 2 : 3;
    s += '{';
    for (int i = 0; i < n; ++i) {
        if (i > 0) s += ',';
        if (!meep_append_double_repr(s, v.in_direction(ds[i]))) return false;
    }
    s += '}';
    return true;
}
%}

// The __repr__ methods keep the GIL even when SWIG threads are enabled: the str()
// methods format into static buffers, and the helpers above call the Python C API.
%nothread meep::vec::__repr__;
%nothread meep::ivec::__repr__;
%nothread meep::volume::__repr__;
%nothread meep::grid_volume::__repr__;

// Include the C++ class declarations
%include "meep/vec.hpp"

%extend meep::vec {
   PyObject *__repr__() {
       std::string s;
       if (!meep_append_vec_repr(s, *$self)) return NULL;
       return PyUnicode_FromString(s.c_str());
   }
}

%extend meep::ivec {
    const char *__repr__() { return $self->str(); }
}

%extend meep::volume {
    PyObject *__repr__() {
        std::string s = "min_corner:";
        if (!meep_append_vec_repr(s, $self->get_min_corner())) return NULL;
        s += ", max_corner:";
        if (!meep_append_vec_repr(s, $self->get_max_corner())) return NULL;
        return PyUnicode_FromString(s.c_str());
  }
}

%extend meep::grid_volume {
    const char *__repr__() { return $self->str(); }
}
