// operators must be explicitly renamed for SWIG to work

%rename(meep_symmetry_add) meep::symmetry::operator+;
%rename(meep_symmetry_mul) meep::symmetry::operator*;
%rename(meep_symmetry_sub) meep::symmetry::operator-;
%rename(meep_symmetry_negate) meep::symmetry::operator-();
%rename(meep_symmetry_eq) meep::symmetry::operator==;
%rename(meep_symmetry_neq) meep::symmetry::operator!=;

%rename(meep_boundary_region_add) meep::boundary_region::operator+;
%rename(meep_boundary_region_mul) meep::boundary_region::operator*;

%rename(meep_dft_chunk_subeq) meep::dft_chunk::operator-=;
%rename(meep_dft_flux_subeq) meep::dft_flux::operator-=;
