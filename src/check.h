/* Copyright (C) 2000 Massachusetts Institute of Technology.
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

#ifndef CHECK_H
#define CHECK_H

/* some useful error-checking macros: */

#define CHECK(condition, message) do { \
     if (!(condition))  { \
          fprintf(stderr, "CHECK failure on line %d of " __FILE__ ": " \
		  message "\n", __LINE__);  exit(EXIT_FAILURE); \
     } \
} while (0)

#define CHK_MALLOC(p, t, n) do {                              \
     size_t CHK_MALLOC_n_tmp = (n);                           \
     (p) = (t *) malloc(sizeof(t) * CHK_MALLOC_n_tmp);        \
     CHECK((p) || CHK_MALLOC_n_tmp == 0, "out of memory!");   \
} while (0)

#endif /* CHECK_H */

