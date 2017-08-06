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


// 322:  Redundant declarations are ok. The wrappers are created correctly.
// 503:  We don't need to create class-specific wrappers for friend functions
// TODO: Check all 509's individually.
// TODO: Move warnfilters to separate 'warnings.i' file.
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
%rename(meep_type) meep::type;
%rename(vec_abs) meep::abs;
%rename(vec_max) meep::max;
%rename(vec_min) meep::min;
%rename(print_grid_volume) meep::grid_volume::print;
%rename(symmetry_reduce) meep::symmetry::reduce;

// Operator renaming
// TODO: Test these
%rename(volume_and) meep::volume::operator&&;
%rename(grid_volume_getitem) meep::grid_volume::operator[];
%rename(symmetry_assign) meep::symmetry::operator=;

%rename(vec_from_dim) meep::vec::vec(ndim, double);

%include "meep/vec.hpp"
