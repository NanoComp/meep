#ifndef MEEP_CTL_H
#define MEEP_CTL_H 

#include "meep.h"

#include "my-smob.h"
#include "meep-ctl-const.h"

extern int verbose; // in main.c

/***************************************************************************/
typedef meep::structure structure_smob;

extern long scm_tc16_smob_structure_smob;
#define STRUCTURE_P(X) T_SMOB_P(structure_smob, X)
#define STRUCTURE(X) T_SMOB(structure_smob, X)
#define SAFE_STRUCTURE(X) SAFE_T_SMOB(structure_smob, X)

extern void register_structure_smobs(void);
extern structure_smob *assert_structure_smob(SCM fo);

/***************************************************************************/

#define CK(ex, msg) \
    (void)((ex) || (meep::abort(msg), 0))


#endif /* MEEP_CTL_H */
