// -*- C++ -*-
/* These are functions for the libctl front-end which are exported
   via SWIG. */

#ifndef MEEP_CTL_SWIG_HPP
#define MEEP_CTL_SWIG_HPP 1

vector3 vec_to_vector3(const meep::vec &v);
meep::vec vector3_to_vec(const vector3 v3);
void set_dimensions(int dims);

meep::structure *make_structure(int dims, vector3 size, vector3 center,
				double resolution, bool enable_averaging,
				double subpixel_tol, int subpixel_maxeval,
				bool ensure_periodicity_p,
				ctlio::geometric_object_list geometry,
				ctlio::material_type_list extra_materials,
				ctlio::material_type default_mat,
				ctlio::pml_list pml_layers,
				ctlio::symmetry_list symmetries,
				int num_chunks, double Courant,
				double global_D_conductivity_diag_,
				double global_B_conductivity_diag_);

ctlio::cvector3_list do_harminv(ctlio::cnumber_list vals, double dt, 
				double fmin, double fmax, int maxbands);

ctlio::number_list dft_flux_flux(meep::dft_flux *f);

ctlio::cnumber_list make_casimir_g(double T, double dt, double sigma, meep::field_type ft,
				   complex<double> (*eps_func)(complex<double> omega) = 0,
				   double Tfft = 0);

// wrapper around constructor to fool SWIG
meep::geometric_volume_list
  *make_geometric_volume_list(const meep::geometric_volume &gv,
			      int c, complex<double> weight,
			      meep::geometric_volume_list *next);

#endif // MEEP_CTL_SWIG_HPP
