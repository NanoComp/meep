/* Copyright (C) 2005-2015 Massachusetts Institute of Technology
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
%
%  You should have received a copy of the GNU General Public License
%  along with this program; if not, write to the Free Software Foundation,
%  Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
*/

#ifndef MEEPTYPES_HPP_
#define MEEPTYPES_HPP_

#include <stddef.h>

/*
 * Integer type used by meep for indexing arrays. It should be set to an SIGNED
 * INTEGER TYPE that can hold all possible memory addresses on the particular
 * system, e.g. a 32 bit int on a 32 bit system and a 64 bit int on a 64 bit
 * system. Usually ptrdiff_t should be a reasonable default as it fulfills those
 * requirements.
 */
namespace meep{
	typedef ptrdiff_t integer;
} /* namespace meep */

#endif /* MEEPTYPES_HPP_ */
