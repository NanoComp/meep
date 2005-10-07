#include "meep-ctl.hpp"

using namespace meep;

/**************************************************************************/

/* The following are hook functions called from main() when
   starting the program and just before exiting.  */

static initialize *meep_init = 0;

void ctl_start_hook(int *argc, char ***argv)
{
  meep_init = new initialize(*argc, *argv);
}

void ctl_stop_hook(void)
{
  delete meep_init;
}

extern "C" void SWIG_init();

void ctl_export_hook(void)
{
  SWIG_init();
}

/**************************************************************************/
