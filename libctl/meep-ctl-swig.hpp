// -*- C++ -*-
/* These are functions for the libctl front-end which are exported
   via SWIG. */

vector3 vec_to_vector3(const meep::vec &v);
meep::vec vector3_to_vec(const vector3 v3);
void set_dimensions(int dims);

meep::structure *make_structure(int dims, vector3 size, double resolution,
				bool ensure_periodicity_p,
				ctlio::geometric_object_list geometry,
				ctlio::material_type default_mat,
				ctlio::pml_list pml_layers,
				ctlio::symmetry_list symmetries,
				int num_chunks, double Courant);

ctlio::cvector3_list do_harminv(ctlio::cnumber_list vals, double dt, 
				double fmin, double fmax, int maxbands);

ctlio::number_list dft_flux_flux(meep::dft_flux *f);
