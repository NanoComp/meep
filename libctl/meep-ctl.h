#ifndef MEEP_CTL_H
#define MEEP_CTL_H 

#include "meep.h"

#include "meep-ctl-const.h"

#include "config.h"
#include "ctl-io.h"

#include "meep-ctl-swig.h"

extern int verbose; // in main.c

/***************************************************************************/

#define CK(ex, msg) \
    (void)((ex) || (meep::abort(msg), 0))


#endif /* MEEP_CTL_H */
