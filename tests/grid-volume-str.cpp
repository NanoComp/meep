#include <cstdio>
#include <cstring>
#include <type_traits>

#include "meep.hpp"

using namespace meep;

static_assert(std::is_same<decltype(&ivec::str), const char *(ivec::*)()>::value,
              "ivec::str must not expose a caller-provided buffer");
static_assert(std::is_same<decltype(&vec::str), const char *(vec::*)()>::value,
              "vec::str must not expose a caller-provided buffer");
static_assert(std::is_same<decltype(&volume::str), const char *(volume::*)()>::value,
              "volume::str must not expose a caller-provided buffer");
static_assert(std::is_same<decltype(&grid_volume::str), const char *(grid_volume::*)()>::value,
              "grid_volume::str must not expose a caller-provided buffer");

int main() {
  ivec integer_vector(1, 2, 3);
  if (std::strcmp(integer_vector.str(), "{1,2,3}") != 0) return 1;

  vec minimum(1.0, 2.0, 3.0);
  vec maximum(4.0, 5.0, 6.0);
  volume bounds(minimum, maximum);
  if (std::strcmp(bounds.str(), "min_corner:{1.000000,2.000000,3.000000}, "
                                "max_corner:{4.000000,5.000000,6.000000}") != 0)
    return 1;

  grid_volume gv = vol3d(10.0, 11.0, 9.0, 1.0);
  const char *description = gv.str();
  if (std::strstr(description, "grid_volume {") == nullptr ||
      std::strstr(description, "dim:3D") == nullptr || std::strstr(description, "\n}") == nullptr) {
    std::fprintf(stderr, "str() did not produce the expected complete output\n");
    return 1;
  }
  return 0;
}
