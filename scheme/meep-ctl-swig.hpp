// -*- C++ -*-
/* These are functions for the libctl front-end which are exported
   via SWIG. */

#ifndef MEEP_CTL_SWIG_HPP
#define MEEP_CTL_SWIG_HPP 1

vector3 vec_to_vector3(const meep::vec &);
meep::vec vector3_to_vec(const vector3 v3);
void set_dimensions(int dims);

meep::structure *make_structure(int dims, vector3 size, vector3 center,
				double resolution, bool enable_averaging,
				double subpixel_tol, int subpixel_maxeval,
				bool ensure_periodicity_p,
				ctlio::geometric_object_list geometry,
				ctlio::material_type_list extra_materials,
				ctlio::material_type default_mat,
				const char *eps_input_file,
				ctlio::pml_list pml_layers,
				ctlio::symmetry_list symmetries,
				int num_chunks, double Courant,
				double global_D_conductivity_diag_,
				double global_B_conductivity_diag_);

ctlio::cvector3_list do_harminv(ctlio::cnumber_list vals, double dt, 
				double fmin, double fmax, int maxbands,
                                double spectral_density, double Q_thresh,
                                double rel_err_thresh, double err_thresh,
                                double rel_amp_thresh, double amp_thresh);

ctlio::cnumber_list dft_flux_flux(meep::dft_flux *f);
ctlio::number_list dft_force_force(meep::dft_force *f);
ctlio::number_list dft_ldos_ldos(meep::dft_ldos *f);
ctlio::cnumber_list dft_ldos_F(meep::dft_ldos *f);
ctlio::cnumber_list dft_ldos_J(meep::dft_ldos *f);
ctlio::cnumber_list dft_near2far_farfield(meep::dft_near2far *f, const meep::vec &x);

ctlio::cnumber_list make_casimir_g(double T, double dt, double sigma, meep::field_type ft,
				   std::complex<double> (*eps_func)(std::complex<double> omega) = 0,
				   double Tfft = 0);

ctlio::cnumber_list make_casimir_g_kz(double T, double dt, double sigma, meep::field_type ft);

// wrapper around constructor to fool SWIG
meep::volume_list *make_volume_list(const meep::volume &v,
				    int c, std::complex<double> weight,
				    meep::volume_list *next);

#endif // MEEP_CTL_SWIG_HPP
