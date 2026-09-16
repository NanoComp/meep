#include <algorithm>
#include <array>
#include <cstdio>
#include <cstring>

#include "meep.hpp"

using namespace meep;

static bool check_bounded_write(grid_volume &gv, size_t buflen) {
  constexpr unsigned char sentinel = 0xa5;
  std::array<unsigned char, 512> storage;
  storage.fill(sentinel);

  char *buffer = reinterpret_cast<char *>(storage.data());
  if (gv.str(buffer, buflen) != buffer) {
    std::fprintf(stderr, "str() did not return the caller-provided buffer\n");
    return false;
  }

  if (!std::all_of(storage.begin() + buflen, storage.end(),
                   [sentinel](unsigned char value) { return value == sentinel; })) {
    std::fprintf(stderr, "str() wrote beyond a %zu-byte buffer\n", buflen);
    return false;
  }

  if (buflen > 0 && std::memchr(buffer, '\0', buflen) == nullptr) {
    std::fprintf(stderr, "str() did not terminate a %zu-byte buffer\n", buflen);
    return false;
  }
  return true;
}

int main() {
  grid_volume gv = vol3d(10.0, 11.0, 9.0, 1.0);
  for (size_t buflen : {size_t(0), size_t(1), size_t(8), size_t(64), size_t(128), size_t(256)})
    if (!check_bounded_write(gv, buflen)) return 1;

  char buffer[512];
  gv.str(buffer, sizeof(buffer));
  if (std::strstr(buffer, "grid_volume {") == nullptr || std::strstr(buffer, "dim:3D") == nullptr ||
      std::strstr(buffer, "\n}") == nullptr) {
    std::fprintf(stderr, "str() did not produce the expected complete output\n");
    return 1;
  }
  return 0;
}
