#ifndef PYMPB_H
#define PYMPB_H

#include "mpb.h"
#include "mpb/maxwell.h"
// #include "mpb/scalar.h"
#include "meepgeom.hpp"

namespace py_mpb {
void add_eigenmode_source(int band_num, const vector3 &kpoint, bool match_frequency,
                          int parity, double resolution, double eigensolver_tol);
} // namespace py_mpb
#endif
