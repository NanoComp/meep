/* Copyright (C) 1999, 2000, 2001, 2002, Massachusetts Institute of Technology.
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

#ifndef SMOB_H
#define SMOB_H 1

/* use new smob functions from Guile 1.4+ */

#ifdef __cplusplus
extern "C" {
#endif

#include <libguile/smob.h>

#ifdef __cplusplus
} /* extern "C" */
#endif

#ifndef HAVE_SCM_MAKE_SMOB_TYPE
/* Thanks to Greg Badros for posting a Guile smob tutorial; see
   http://sources.redhat.com/ml/guile/1999-04/msg00107.html
   However, this way of creating smobs no longer works as of Guile 1.4. */
#define REGISTER_SMOBFUNS(T) \
  do { scm_tc16_smob_ ## T = scm_newsmob(& T ## _smobfuns); } while (0)
#define MAKE_SMOBFUNS(T) \
  static scm_smobfuns T ## _smobfuns = { \
    &mark_ ## T, \
    &free_ ## T, \
    &print_ ## T,  0 }
#endif /* ! HAVE_SCM_MAKE_SMOB_TYPE */

#ifdef HAVE_SCM_NEWSMOB
#  define NEWCELL_SMOB(ANSWER,T,PSMOB) \
       SCM_NEWSMOB(ANSWER, scm_tc16_smob_ ## T, PSMOB)
#else
#  define NEWCELL_SMOB(ANSWER,T,PSMOB) do { \
     SCM_NEWCELL((ANSWER)); \
     SCM_SETCAR((ANSWER),scm_tc16_smob_ ## T); \
     SCM_SETCDR((ANSWER),(SCM) (PSMOB)); \
   } while (0)
#endif

/* T_SMOB_P(T, X) is true iff X is an instance of the T SMOB type */
#ifdef HAVE_SCM_SMOB_PREDICATE
#  define T_SMOB_P(T, X) SCM_SMOB_PREDICATE(scm_tc16_smob_ ## T, X)
#else
#  define T_SMOB_P(T, X) (SCM_NIMP (X) \
                          && SCM_TYP16 (X) == (scm_tc16_smob_ ## T))
#endif

/* T_SMOB(T, X) returns the T * with the guts of the X instance; it
   assumes X is a T SMOB instance, and could crash if it is not */
#ifdef HAVE_SCM_SMOB_DATA
#  define T_SMOB(T, X)  ((T *) SCM_SMOB_DATA(X))
#else
#  define T_SMOB(T, X)  ((T *)SCM2PTR(SCM_CDR(X)))
#endif

/* Since T_SMOB(X) can be dangerous if X is not a T
   object, we also have a SAFE_T_SMOB macro: */
#define SAFE_T_SMOB(T, X) (T_SMOB_P(T,X) ? T_SMOB(T,X) : NULL)

#endif /* SMOB_H */
