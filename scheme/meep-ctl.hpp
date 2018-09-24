// -*- C++ -*-
#ifndef MEEP_CTL_H
#define MEEP_CTL_H 

#include "meep.hpp"

#include "meep-ctl-const.hpp"

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif
#include "ctl-io.h"

#include "meep-ctl-swig.hpp"

extern int verbose; // in main.c

/***************************************************************************/

#define CK(ex, msg) \
    (void)((ex) || (meep::abort(msg), 0))


#endif /* MEEP_CTL_H */
