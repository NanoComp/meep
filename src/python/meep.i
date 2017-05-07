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

%module(directors="1") meep

%{
#define SWIG_FILE_WITH_INIT
#include "meep/vec.hpp"
#include "meep.hpp"
using namespace meep;
%}

%include "numpy.i"

%{
#if NPY_API_VERSION < 0x00000007
#define NPY_ARRAY_C_CONTIGUOUS NPY_C_CONTIGUOUS
#define NPY_ARRAY_ALIGNED  NPY_ALIGNED
#endif
%}

%init %{
  import_array();
%}

// TODO: apply necessary numpy typemaps

// Rename python builtins
%rename(br_apply) meep::boundary_region::apply;
%rename(_is) meep::dft_chunk::is;
%rename(Meep_None) meep::None;

// Operator renaming
%rename(boundary_region_assign) meep::boundary_region::operator=;

// TODO:  Fix these with a typemap when necesary
%feature("immutable") meep::fields_chunk::connections;
%feature("immutable") meep::fields_chunk::num_connections;

%feature("director") meep::Callback;

// TODO: Ignore with build flags once all have been analyzed and determined benign.
%warnfilter(509);
%warnfilter(322,509) meep::component_direction;
%warnfilter(322,509) meep::direction_component;
%warnfilter(322,503) meep::zero_vec;
%warnfilter(322,503) meep::veccyl;
%warnfilter(322,503) meep::zero_ivec;
%warnfilter(322,503) meep::one_ivec;
%warnfilter(322,503) meep::iveccyl;
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
%rename(ftype) meep::type;
%rename(vec_abs) meep::abs;
%rename(vec_max) meep::max;
%rename(vec_min) meep::min;
%rename(print_grid_volume) meep::grid_volume::print;
%rename(symmetry_reduce) meep::symmetry::reduce;

// Operator renaming
%rename(volume_and) meep::volume::operator&&;
%rename(grid_volume_getitem) meep::grid_volume::operator[];
%rename(symmetry_assign) meep::symmetry::operator=;


%include "meep/vec.hpp"
%include "meep.hpp"
