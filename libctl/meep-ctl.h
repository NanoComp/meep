#ifndef MEEP_CTL_H
#define MEEP_CTL_H 

#include "meep.h"

#include "meep-ctl-const.h"

#include "config.h"
#include "ctl-io.h"

namespace ctlio {
  vector3 vec2vector3(const meep::vec &v);
}

extern int verbose; // in main.c

/***************************************************************************/

#define CK(ex, msg) \
    (void)((ex) || (meep::abort(msg), 0))


#endif /* MEEP_CTL_H */
