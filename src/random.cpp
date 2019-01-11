/* Copyright (C) 2005-2019 Massachusetts Institute of Technology.
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

#include "meep.hpp"
#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif

#ifdef HAVE_LIBGSL
#include <gsl/gsl_randist.h>
#else
#include <stdlib.h>
#endif
#include <time.h>

using namespace std;

namespace meep {

static bool rand_inited = false;
#ifdef HAVE_LIBGSL
static gsl_rng *rng = NULL;
#endif

static void init_rand(void) {
  if (!rand_inited) {
#ifdef HAVE_LIBGSL
    if (rng) gsl_rng_free(rng);
    rng = gsl_rng_alloc(gsl_rng_mt19937);
#endif
    rand_inited = true; // no infinite loop since rand_inited == true
    set_random_seed(time(NULL) * (1 + my_global_rank()));
  }
}

void set_random_seed(unsigned long seed) {
  init_rand();
  seed = ((unsigned long)broadcast(0, (int)seed)) + my_rank();
#ifdef HAVE_LIBGSL
  gsl_rng_set(rng, seed);
#else
  srand(seed);
#endif
}

int random_int(int a, int b) {
  init_rand();
#ifdef HAVE_LIBGSL
  return ((int)gsl_rng_uniform_int(rng, b - a + 1)) + a;
#else
  return a + rand() % (b - a + 1);
#endif
}

double uniform_random(double a, double b) {
  init_rand();
#ifdef HAVE_LIBGSL
  return a + gsl_rng_uniform(rng) * (b - a);
#else
  return a + rand() * (b - a) / RAND_MAX;
#endif
}

double gaussian_random(double mean, double stddev) {
  init_rand();
#ifdef HAVE_LIBGSL
  return mean + gsl_ran_gaussian(rng, stddev);
#else
  // Box-Muller algorithm to generate Gaussian from uniform
  // see Knuth vol II algorithm P, sec. 3.4.1
  double v1, v2, s;
  do {
    v1 = uniform_random(-1, 1);
    v2 = uniform_random(-1, 1);
    s = v1 * v1 + v2 * v2;
  } while (s >= 1.0);
  if (s == 0) {
    return mean;
  } else {
    return mean + v1 * sqrt(-2 * log(s) / s) * stddev;
  }
#endif
}

} // namespace meep
