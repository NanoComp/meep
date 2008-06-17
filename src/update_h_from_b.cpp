#include "meep.hpp"
#include "meep_internals.hpp"

namespace meep {
  
void fields::update_h_from_b() {
  for (int i=0;i<num_chunks;i++)
    if (chunks[i]->is_mine())
      chunks[i]->update_h_from_b();
}

void fields_chunk::update_h_from_b() {
  
    DOCMP FOR_H_AND_B(hc,bc) if (f[hc][cmp]) {
      const direction d_hc = component_direction(hc);
      const int s_hc = stride_any_direction[d_hc];
      const direction d_1 = direction((d_hc+1)%3);
      const component hc_1 = direction_component(hc,d_1);
      const int s_1 = stride_any_direction[d_1];
      const direction d_2 = direction((d_hc+2)%3);
      const component hc_2 = direction_component(hc,d_2);
      const int s_2 = stride_any_direction[d_2];
      
      component bc_1 = direction_component(bc,d_1);
      component bc_2 = direction_component(bc,d_2);
      
      direction dsig = (direction)((d_hc+2)%3);
      direction dsigg = (direction)(d_hc);
      direction dsig1 = (direction)((d_hc+1)%3);
      direction dsig1inv = (direction)(d_hc);
      direction dsig2 = (direction)((d_hc+2)%3);
      direction dsig2inv = (direction)((d_hc+1)%3);
      
      step_update_EDHB(f[hc][cmp], hc, v, 
		       f[bc][cmp], f[bc_1][cmp], f[bc_2][cmp],
		       f_prev[bc][cmp], f_prev[bc_1][cmp], f_prev[bc_2][cmp],
		       s->invmu[hc][d_hc], s->invmu[hc][d_1], s->invmu[hc][d_2],
		       s_hc, s_1, s_2, NULL, NULL,
		       dsig, s->sig[dsig], s->siginv[dsig],
		       dsigg, s->sig[dsigg],
		       dsig1, s->sig[dsig1],
		       dsig1inv, s->sig[dsig1inv],
		       dsig2, s->sig[dsig2],
		       dsig2inv, s->sig[dsig2inv],
		       s->sigsize[dsig],s->sigsize[dsigg],s->sigsize[dsig1]);
    }
  
}

} // namespace meep
