/* These are functions for the libctl front-end which are exported
   via SWIG. */

vector3 vec2vector3(const meep::vec &v);

meep::structure *make_structure(int dims, vector3 size, double resolution,
				ctlio::geometric_object_list geometry,
				int desired_num_chunks,
				const meep::symmetry &S);
