/* Copyright (C) 2005-2026 Massachusetts Institute of Technology
%
%  This program is free software; you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation; either version 2, or (at your option)
%  any later version.
%
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
*/

#include "meep/backend_hooks.hpp"

namespace meep {

/* Process-global table.  Zero-initialized: every hook starts as a null
 * function pointer, which the call sites read as "no backend installed,
 * fall through to the CPU path". */
backend_hooks meep_backend = {};

} /* namespace meep */
