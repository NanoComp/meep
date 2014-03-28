/* Copyright (C) 2005-2014 Massachusetts Institute of Technology.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

#ifndef BICGSTAB_H
#define BICGSTAB_H

#include "meep.hpp"

namespace meep {

typedef void (*bicgstab_op)(const realnum *x, realnum *y, void *data);

int bicgstabL(const int L, 
	      const int n, realnum *x,
	      bicgstab_op A, void *Adata, const realnum *b,
	      const double tol, 
	      int *iters, // input *iters = max iters, output = actual iters
	      realnum *work, // if you pass work=NULL, bicgstab returns nwork
	      const bool quiet);

} // namespace meep

#endif /* BICGSTAB_H */
