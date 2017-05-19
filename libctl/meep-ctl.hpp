// -*- C++ -*-
#ifndef MEEP_CTL_H
#define MEEP_CTL_H 

#include "meep.hpp"

#include "meep-ctl-const.hpp"

#include "config.h"
#include "ctl-io.h"

#include "meep-ctl-swig.hpp"

extern int verbose;      // in meep.cpp
extern int libctl_quiet; // in meep.cpp

/***************************************************************************/

#define CK(ex, msg) \
    (void)((ex) || (meep::abort(msg), 0))


#endif /* MEEP_CTL_H */
