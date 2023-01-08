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

#include "meep.hpp"
#include "config.h"

#include "support/meep_mt.h"
#include <time.h>

using namespace std;

namespace meep {

static bool rand_inited = false;

static void init_rand(void) {
  if (!rand_inited) {
    rand_inited = true; // no infinite loop since rand_inited == true
    set_random_seed(time(NULL) * (1 + my_global_rank()));
  }
}

void set_random_seed(unsigned long seed) {
  init_rand();
  meep_mt_init_genrand(seed);
}

void restore_random_seed() {
  init_rand();
  meep_mt_restore_genrand();
}

int random_int(int a, int b) {
  init_rand();
  return a + meep_mt_genrand_int32() % (b - a + 1);
}

double uniform_random(double a, double b) {
  init_rand();
  return a + meep_mt_genrand_res53() * (b - a);
}

double gaussian_random(double mean, double stddev) {
  init_rand();
  // Box-Muller algorithm to generate Gaussian from uniform
  // see Knuth vol II algorithm P, sec. 3.4.1
  double v1, v2, s;
  do {
    v1 = 2 * meep_mt_genrand_res53() - 1;
    v2 = 2 * meep_mt_genrand_res53() - 1;
    s = v1 * v1 + v2 * v2;
  } while (s >= 1.0);
  if (s == 0) { return mean; }
  else { return mean + v1 * sqrt(-2 * log(s) / s) * stddev; }
}

} // namespace meep
