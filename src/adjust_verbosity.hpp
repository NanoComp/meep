/* Copyright (C) 2005-2023 Massachusetts Institute of Technology
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
#ifndef ADJUST_VERBOSITY_H
#define ADJUST_VERBOSITY_H

// NOTE: This header assumes it has been #included *after* headers that declare
// verbosity, HAVE_MPB, and mpb_verbosity.

namespace meep {

// This is a RAII (Resource Allocation Is Initialization) class which adjusts
// the mpb_verbosity level when created, and restores it when deleted.
class adjust_mpb_verbosity {
public:
  adjust_mpb_verbosity() {
#if defined(HAVE_MPB) &&                                                                           \
    (MPB_VERSION_MAJOR > 1 || (MPB_VERSION_MAJOR == 1 && MPB_VERSION_MINOR >= 11))
    old_level = mpb_verbosity;
    mpb_verbosity = verbosity - 1;
    if (mpb_verbosity < 0) mpb_verbosity = 0;
    if (mpb_verbosity > 3) mpb_verbosity = 3;
#else
    // avoid warnings
    (void)old_level;
#endif
  }

  ~adjust_mpb_verbosity() {
#if defined(HAVE_MPB) &&                                                                           \
    (MPB_VERSION_MAJOR > 1 || (MPB_VERSION_MAJOR == 1 && MPB_VERSION_MINOR >= 11))
    mpb_verbosity = old_level;
#endif
  }

private:
  int old_level;
};

} // namespace meep

#endif // ADJUST_VERBOSITY_H
