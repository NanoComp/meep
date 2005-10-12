// -*- C++ -*-
/* These are functions for the libctl front-end which are exported
   via SWIG. */

vector3 vec2vector3(const meep::vec &v);

meep::structure *make_structure(int dims, vector3 size, double resolution,
				ctlio::geometric_object_list geometry,
				ctlio::material_type default_mat,
				ctlio::pml_list pml_layers,
				ctlio::symmetry_list symmetries,
				int num_chunks, double Courant);

ctlio::cvector3_list do_harminv(ctlio::cnumber_list vals, double dt, 
				double fmin, double fmax, int maxbands);
