#include <cassert>
#include <meep.hpp>
using namespace meep;
using namespace meep_geom;

double eps(const vec pt) {
    return 1.0;
}

int main(int argc, char **argv) {
    initialize mpi(argc, argv);

    double resolution = 20;
    grid_volume v = vol2d(5, 10, resolution);
    structure s(v, eps, pml(1.0));

    geometric_object objects[1];
    objects[0] = make_sphere();
    geometric_object_list geom_list = {1, objects};
    set_materials_from_geometry(&s, geom_list);

    fields f(&s);

    f.output_hdf5(Dielectric, v.surroundings());

    double freq = 0.3;
    double fwidth = 0.1;
    gaussian_src_time src(freq, fwidth);
    f.add_point_source(Ey, src, vec(1.1, 2.3));
    while (f.time() < f.last_source_time()) {
        f.step();
    }

    return 0;
}
