#include <assert.h>
#include "meep.hpp"
#include "meepgeom.hpp"

using namespace meep;
using namespace meep_geom;

int main(int argc, char **argv) {

    grid_volume gv1 = vol1d(5, 11);
    gv1.center_origin();
    std::vector<geom_box> b1 = split_grid_volume_into_boxes(&gv1, 10);
    assert(b1.size() == 6);

    grid_volume gv2 = vol2d(5, 5, 11);
    gv2.center_origin();
    std::vector<geom_box> b2 = split_grid_volume_into_boxes(&gv2, 10);
    assert(b2.size() == 6 * 6);

    grid_volume gv3 = vol3d(5, 5, 5, 11);
    gv3.center_origin();
    std::vector<geom_box> b3 = split_grid_volume_into_boxes(&gv3, 10);
    assert(b3.size() == 6 * 6 * 6);

    grid_volume gvcyl = volcyl(5, 5, 11);
    gvcyl.center_origin();
    std::vector<geom_box> bcyl = split_grid_volume_into_boxes(&gvcyl, 10);
    assert(bcyl.size() == 6 * 6);

    return 0;
}
