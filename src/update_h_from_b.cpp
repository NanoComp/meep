#include <string.h>

#include "meep.hpp"
#include "meep_internals.hpp"

namespace meep {

void fields::update_h_from_b() {
  for (int i=0;i<num_chunks;i++)
    if (chunks[i]->is_mine())
      if (chunks[i]->update_h_from_b())
	chunk_connections_valid = false; // H allocated - reconnect chunks

  /* synchronize to avoid deadlocks if one process decides it needs
     to allocate H ... */
  chunk_connections_valid = and_to_all(chunk_connections_valid);
}

// set H = B/mu ... return true if any H field was allocated
bool fields_chunk::update_h_from_b() {
  bool allocated_h = false;

  // now, update H from B with the help of step_generic.cpp
  DOCMP FOR_H_AND_B(hc,bc) if (f[hc][cmp]) {
      const direction d_hc = component_direction(hc);
      const int s_hc = stride_any_direction[d_hc];
      const direction d_1 = cycle_direction(v.dim, d_hc, 1);
      const int s_1 = stride_any_direction[d_1];
      const direction d_2 = cycle_direction(v.dim, d_hc, 2);
      const int s_2 = stride_any_direction[d_2];
      
      component bc_1 = direction_component(bc,d_1);
      component bc_2 = direction_component(bc,d_2);
      
      direction dsig = d_2;
      direction dsigg = d_hc;
      direction dsig1 = d_1;
      direction dsig1inv = d_hc;
      direction dsig2 = d_2;
      direction dsig2inv = d_1;
      
      // lazily allocate any H fields that are needed (H==B initially)
      if (f[hc][cmp] == f[bc][cmp]
	  && (s->invmu[hc][d_hc]
	      || s->sigsize[dsig] > 1
	      || s->sigsize[dsigg] > 1
	      || (s->sigsize[dsig1] > 1
		  && (s->invmu[hc][d_1] || s->invmu[hc][d_2])))) {
	f[hc][cmp] = new double[v.ntot()];
	memcpy(f[hc][cmp], f[bc][cmp], v.ntot() * sizeof(double));
	allocated_h = true;
      }

      if (f[hc][cmp] != f[bc][cmp])
	step_update_EDHB(f[hc][cmp], hc, v, 
			 f[bc][cmp], f[bc_1][cmp], f[bc_2][cmp],
			 f_prev[bc][cmp], f_prev[bc_1][cmp], f_prev[bc_2][cmp],
			 s->invmu[hc][d_hc],
			 s->invmu[hc][d_1], s->invmu[hc][d_2],
			 s_hc, s_1, s_2, NULL, NULL,
			 dsig, s->sig[dsig], s->siginv[dsig],
			 dsigg, s->sig[dsigg],
			 dsig1, s->sig[dsig1],
			 dsig1inv, s->sig[dsig1inv],
			 dsig2, s->sig[dsig2],
			 dsig2inv, s->sig[dsig2inv],
			 s->sigsize[dsig],s->sigsize[dsigg],s->sigsize[dsig1]);
    }
  
  return allocated_h;
}

} // namespace meep
