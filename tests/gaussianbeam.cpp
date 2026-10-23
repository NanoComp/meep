// Gaussian-beam fields must satisfy Maxwell's equations and overwrite their output.
#include <algorithm>
#include <array>
#include <cmath>
#include <complex>

#include <meep.hpp>

using namespace meep;
typedef std::complex<double> cd;
typedef std::array<cd, 6> field_values;

static field_values field_at(const gaussianbeam &beam, const vec &p, cd initial = cd(0.0)) {
  field_values result;
  result.fill(initial);
  beam.get_fields(result.data(), p);
  return result;
}

static void check_component(cd actual, cd expected, const char *law, int component, int rotation,
                            double waist, int material, int polarization) {
  const double error = std::abs(actual - expected);
  const double scale = std::max(1.0, std::max(std::abs(actual), std::abs(expected)));
  if (!std::isfinite(error) || error > 1e-6 * scale)
    meep::abort("gaussianbeam %s component %d: error=%g scale=%g; rotation=%d waist=%g "
                "material=%d polarization=%d\n",
                law, component, error, scale, rotation, waist, material, polarization);
}

int main(int argc, char **argv) {
  initialize mpi(argc, argv);
  const double freq = 1.0;
  const double h = 1e-5;
  const direction axes[3] = {X, Y, Z};
  const double waists[2] = {0.8, 2.0};
  const cd polarizations[4][2] = {{cd(1.0), cd(0.0)},
                                  {cd(0.0, 1.0), cd(0.0)},
                                  {cd(1.0, 0.4), cd(0.2, 1.0)},
                                  {cd(1.0), cd(0.0, 1.0)}};

  for (int rotation = 0; rotation < 2; ++rotation) {
    // The rotated basis is orthonormal and right-handed, with ez x ex = ey.
    const vec ex = rotation ? vec(0.0, 0.0, 1.0) : vec(1.0, 0.0, 0.0);
    const vec ey = rotation ? vec(0.8, -0.6, 0.0) : vec(0.0, 1.0, 0.0);
    const vec ez = rotation ? vec(0.6, 0.8, 0.0) : vec(0.0, 0.0, 1.0);
    const vec x0(0.13, -0.17, 0.08);
    // Both transverse offsets must be nonzero to expose the missing local Hx.
    const vec p = x0 + ex * 0.3 + ey * 0.4 + ez * 0.5;

    for (int wi = 0; wi < 2; ++wi)
      for (int material = 0; material < 2; ++material)
        for (int polarization = 0; polarization < 4; ++polarization) {
          const double eps = material ? 2.25 : 1.0;
          const double mu = material ? 1.3 : 1.0;
          cd E0[3];
          for (int j = 0; j < 3; ++j)
            E0[j] = polarizations[polarization][0] * ex.in_direction(axes[j]) +
                    polarizations[polarization][1] * ey.in_direction(axes[j]);

          gaussianbeam beam(x0, ez, waists[wi], freq, eps, mu, E0);
          const field_values center = field_at(beam, p);
          cd derivative[3][6];
          for (int axis = 0; axis < 3; ++axis) {
            vec delta(0.0, 0.0, 0.0);
            delta.set_direction(axes[axis], h);
            const field_values hi = field_at(beam, p + delta);
            const field_values lo = field_at(beam, p - delta);
            for (int j = 0; j < 6; ++j)
              derivative[axis][j] = (hi[j] - lo[j]) / (2 * h);
          }

          for (int j = 0; j < 3; ++j) {
            const int a = (j + 1) % 3, b = (j + 2) % 3;
            const cd from_E = (derivative[a][b] - derivative[b][a]) / cd(0.0, 2 * pi * freq * mu);
            const cd from_H =
                -(derivative[a][b + 3] - derivative[b][a + 3]) / cd(0.0, 2 * pi * freq * eps);
            check_component(center[j + 3], from_E, "Faraday", j, rotation, waists[wi], material,
                            polarization);
            check_component(center[j], from_H, "Ampere", j, rotation, waists[wi], material,
                            polarization);
          }

          // get_fields must overwrite its output, independent of previous contents.
          const field_values reused = field_at(beam, p, cd(7.0, 3.0));
          for (int j = 0; j < 6; ++j)
            check_component(reused[j], center[j], "output reuse", j, rotation, waists[wi], material,
                            polarization);
        }
  }
  master_printf("gaussianbeam Maxwell and output-reuse checks passed (32 cases).\n");
  return 0;
}
