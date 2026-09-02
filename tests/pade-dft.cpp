/* Copyright (C) 2005-2026 Massachusetts Institute of Technology
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 */

#include <algorithm>
#include <cmath>
#include <complex>
#include <vector>

#include <meep.hpp>

using namespace meep;

double vacuum_eps(const vec &) { return 1.0; }

static void require(bool condition, const char *message) {
  if (!condition) meep::abort("Padé DFT test failed: %s", message);
}

static void require_close(std::complex<double> actual, std::complex<double> expected,
                          double relative_tolerance, const char *message) {
  const double scale = std::max(1.0, std::abs(expected));
  if (!(std::abs(actual - expected) <= relative_tolerance * scale))
    meep::abort("Padé DFT test failed: %s (actual %.16g%+.16gi, expected %.16g%+.16gi)",
                message, actual.real(), actual.imag(), expected.real(), expected.imag());
}

static size_t local_chunk_count(dft_chunk *chunks) {
  size_t count = 0;
  for (dft_chunk *cur = chunks; cur; cur = cur->next_in_dft)
    ++count;
  return count;
}

static void set_monitor_field(dft_chunk *chunks, double value) {
  for (dft_chunk *cur = chunks; cur; cur = cur->next_in_dft) {
    PLOOP_OVER_IVECS(cur->fc->gv, cur->is, cur->ie, idx) {
      cur->fc->f[cur->c][0][idx] = realnum(value);
      if (cur->fc->f[cur->c][1]) cur->fc->f[cur->c][1][idx] = 0.0;
    }
  }
}

static void append_sample(dft_chunk *chunks, double value, double time) {
  set_monitor_field(chunks, value);
  for (dft_chunk *cur = chunks; cur; cur = cur->next_in_dft)
    cur->update_dft(time);
}

static double modal_sample(const std::vector<double> &amplitudes,
                           const std::vector<double> &ratios, size_t n) {
  double value = 0.0;
  for (size_t i = 0; i < amplitudes.size(); ++i)
    value += amplitudes[i] * std::pow(ratios[i], double(n));
  return value;
}

static std::complex<double> modal_sum(const std::vector<double> &amplitudes,
                                      const std::vector<double> &ratios,
                                      std::complex<double> z, size_t count) {
  std::complex<double> sum = 0.0;
  for (size_t i = 0; i < amplitudes.size(); ++i) {
    const std::complex<double> rz = ratios[i] * z;
    sum += amplitudes[i] * (1.0 - std::pow(rz, double(count))) / (1.0 - rz);
  }
  return sum;
}

static std::complex<double> modal_infinite_sum(const std::vector<double> &amplitudes,
                                               const std::vector<double> &ratios,
                                               std::complex<double> z) {
  std::complex<double> sum = 0.0;
  for (size_t i = 0; i < amplitudes.size(); ++i)
    sum += amplitudes[i] / (1.0 - ratios[i] * z);
  return sum;
}

static void test_analytic_tail(component c, const std::vector<double> &amplitudes,
                               const std::vector<double> &ratios, size_t history_size,
                               size_t sample_count, double frequency, double first_time,
                               int decimation_factor) {
  grid_volume gv = volone(4.0, 10.0);
  gv.center_origin();
  structure s(gv, vacuum_eps);
  fields f(&s);
  f.use_real_fields();
  f.add_point_source(c, gaussian_src_time(0.2, 0.4), vec(0.0), 1.0);

  component components[] = {c};
  const volume point(vec(0.0), vec(0.0));
  dft_fields monitor = f.add_dft_fields(components, 1, point, &frequency, 1, false,
                                        decimation_factor, false, history_size);
  const double sample_period = decimation_factor * f.dt;
  for (size_t n = 0; n < sample_count; ++n)
    append_sample(monitor.chunks, modal_sample(amplitudes, ratios, n),
                  first_time + n * sample_period);

  require(sum_to_all(local_chunk_count(monitor.chunks)) > 0, "analytic monitor has no chunks");
  const double tolerance = sizeof(realnum) == sizeof(float) ? 2e-4 : 2e-10;
  for (dft_chunk *cur = monitor.chunks; cur; cur = cur->next_in_dft) {
    require(cur->get_pade_count() == history_size, "analytic history did not wrap to capacity");
    const std::complex<double> z = std::polar(1.0, cur->omega[0] * sample_period);
    const std::complex<double> phase = std::polar(1.0, cur->omega[0] * first_time);
    const std::complex<double> raw_expected =
        cur->scale * phase * modal_sum(amplitudes, ratios, z, sample_count);
    const std::complex<double> corrected_expected =
        cur->scale * phase * modal_infinite_sum(amplitudes, ratios, z);
    const std::complex<double> raw(cur->dft[0]);
    const std::complex<double> corrected(cur->dft_values()[0]);
    require_close(raw, raw_expected, tolerance, "ordinary DFT has the wrong sampled value");
    require_close(corrected, corrected_expected, tolerance,
                  "Padé correction does not equal the analytic infinite sum");
    require_close(corrected - raw, corrected_expected - raw_expected, tolerance,
                  "Padé correction double-counted already accumulated samples");
    const dft_pade_error error = cur->get_pade_error();
    require(error.ready, "analytic fit was not ready");
    require(error.invalid_fits == 0, "analytic damped-mode fit was rejected");
  }
}

static void test_rejected_tail(double ratio, const char *message) {
  const size_t history_size = 8;
  const double frequency = 0.0;
  grid_volume gv = volone(4.0, 10.0);
  gv.center_origin();
  structure s(gv, vacuum_eps);
  fields f(&s);
  f.use_real_fields();
  f.add_point_source(Ex, gaussian_src_time(0.2, 0.4), vec(0.0), 1.0);
  component components[] = {Ex};
  dft_fields monitor = f.add_dft_fields(components, 1, volume(vec(0.0), vec(0.0)), &frequency, 1,
                                        false, 1, false, history_size);
  for (size_t n = 0; n < history_size; ++n)
    append_sample(monitor.chunks, std::pow(ratio, double(n)), 0.75 + n * f.dt);

  for (dft_chunk *cur = monitor.chunks; cur; cur = cur->next_in_dft) {
    const std::complex<realnum> raw = cur->dft[0];
    const std::complex<realnum> corrected = cur->dft_values()[0];
    require(corrected == raw, message);
    require(cur->get_pade_error().invalid_fits > 0, "rejected Padé fit was not reported");
  }
}

static void test_diagnostic_checkpoint_isolation() {
  const size_t history_size = 8;
  const double frequency = 0.11;
  grid_volume gv = volone(4.0, 10.0);
  gv.center_origin();
  structure s(gv, vacuum_eps);
  fields f(&s);
  f.use_real_fields();
  f.add_point_source(Ex, gaussian_src_time(0.2, 0.4), vec(0.0), 1.0);
  component components[] = {Ex};
  dft_fields monitor = f.add_dft_fields(components, 1, volume(vec(0.0), vec(0.0)), &frequency, 1,
                                        false, 1, false, history_size);
  for (size_t n = 0; n < history_size; ++n)
    append_sample(monitor.chunks, std::pow(0.8, double(n)) + 0.02 * std::pow(0.3, double(n)),
                  0.5 + n * f.dt);

  for (dft_chunk *cur = monitor.chunks; cur; cur = cur->next_in_dft) {
    const std::complex<double> checkpoint(cur->dft_values()[0]);
    const dft_pade_error first = cur->get_pade_error();
    require(std::isinf(first.drift[0]), "first diagnostic unexpectedly had a drift baseline");

    append_sample(cur, 2.0 * std::pow(0.8, double(history_size)),
                  0.5 + history_size * f.dt);
    const std::complex<double> intervening(cur->dft_values()[0]);
    append_sample(cur, 0.25 * std::pow(0.8, double(history_size + 1)),
                  0.5 + (history_size + 1) * f.dt);
    const std::complex<double> current(cur->dft_values()[0]);
    const double expected = std::abs(current - checkpoint) /
                            std::max(std::abs(current), std::abs(checkpoint));
    const double adjacent = std::abs(current - intervening) /
                            std::max(std::abs(current), std::abs(intervening));
    const dft_pade_error second = cur->get_pade_error();
    require(std::abs(second.drift[0] - expected) <= 1e-10 * std::max(1.0, expected),
            "ordinary corrected-value reads advanced the diagnostic baseline");
    require(std::abs(expected - adjacent) > 1e-8,
            "diagnostic isolation fixture did not distinguish adjacent generations");
    require(cur->get_pade_error().drift[0] == second.drift[0],
            "repeated diagnostics changed the same-generation baseline");
    break;
  }
}

static void test_lifecycle_and_bounded_storage() {
  const size_t history_size = 8;
  const double frequency = 0.2;
  grid_volume gv = volone(4.0, 10.0);
  gv.center_origin();
  structure s(gv, vacuum_eps);
  fields f(&s);
  f.use_real_fields();
  f.add_point_source(Ex, gaussian_src_time(frequency, 0.4), vec(0.0), 1.0);

  component components[] = {Ex};
  dft_fields enabled_a = f.add_dft_fields(components, 1, f.v, &frequency, 1, true, 1, false,
                                          history_size);
  dft_fields enabled_b = f.add_dft_fields(components, 1, f.v, &frequency, 1, true, 1, false,
                                          history_size);
  dft_fields disabled = f.add_dft_fields(components, 1, f.v, &frequency, 1, true, 1, false, 0);

  for (size_t step = 0; step < 3 * history_size; ++step)
    f.step();

  require(sum_to_all(local_chunk_count(enabled_a.chunks)) > 0, "monitor has no chunks");
  for (dft_chunk *cur = enabled_a.chunks; cur; cur = cur->next_in_dft) {
    require(cur->pade_enabled(), "enabled monitor did not allocate Padé state");
    require(cur->get_pade_samples() == history_size, "history capacity changed");
    require(cur->get_pade_count() == history_size, "history did not remain bounded");
    const std::complex<realnum> *values = cur->dft_values();
    require(values == cur->dft_values(), "corrected-value cache was not reused");
  }
  for (dft_chunk *cur = disabled.chunks; cur; cur = cur->next_in_dft) {
    require(!cur->pade_enabled(), "zero history unexpectedly enabled Padé");
    require(cur->dft_values() == cur->dft, "disabled accessor did not return raw storage");
  }

  std::vector<std::vector<std::complex<realnum> > > raw_before_scale;
  for (dft_chunk *cur = enabled_a.chunks; cur; cur = cur->next_in_dft)
    raw_before_scale.emplace_back(cur->dft, cur->dft + cur->N * cur->omega.size());
  enabled_a.scale_dfts(2.0);
  size_t chunk_index = 0;
  for (dft_chunk *cur = enabled_a.chunks; cur; cur = cur->next_in_dft, ++chunk_index) {
    require(cur->get_pade_count() == 0, "scaling did not clear history");
    for (size_t i = 0; i < cur->N * cur->omega.size(); ++i)
      require(cur->dft[i] == realnum(2.0) * raw_before_scale[chunk_index][i],
              "scaling materialized a corrected prediction");
  }

  dft_fields enabled_c = f.add_dft_fields(components, 1, f.v, &frequency, 1, true, 1, false,
                                          history_size);
  dft_fields enabled_d = f.add_dft_fields(components, 1, f.v, &frequency, 1, true, 1, false,
                                          history_size);
  for (size_t step = 0; step < history_size; ++step)
    f.step();
  if (enabled_c.chunks && enabled_d.chunks) {
    *enabled_c.chunks -= *enabled_d.chunks;
    for (dft_chunk *cur = enabled_c.chunks; cur; cur = cur->next_in_dft) {
      require(cur->get_pade_count() == 0, "subtraction did not clear history");
      for (size_t i = 0; i < cur->N * cur->omega.size(); ++i)
        require(cur->dft[i] == std::complex<realnum>(0.0, 0.0),
                "identical corrected snapshots did not subtract to zero");
    }
  }

  f.t = 0;
  f.zero_fields();
  f.step();
  for (dft_chunk *cur = enabled_b.chunks; cur; cur = cur->next_in_dft)
    require(cur->get_pade_count() == 1, "time discontinuity did not start a fresh history");
}

static void test_exported_overloads_link() {
  typedef dft_chunk *(fields::*old_add_dft)(component, const volume &, const double *, size_t, bool,
                                            std::complex<double>, dft_chunk *, bool,
                                            std::complex<double>, bool, int, int, bool);
  typedef dft_chunk *(fields::*new_add_dft)(component, const volume &, const double *, size_t, bool,
                                            std::complex<double>, dft_chunk *, bool,
                                            std::complex<double>, bool, int, int, bool, size_t);
  typedef dft_flux (fields::*old_flux)(const volume_list *, const double *, size_t, bool, bool, int);
  typedef dft_flux (fields::*new_flux)(const volume_list *, const double *, size_t, bool, bool, int,
                                       size_t);
  typedef dft_flux (fields::*old_mode)(direction, const volume &, const double *, size_t, bool, int);
  typedef dft_flux (fields::*new_mode)(direction, const volume &, const double *, size_t, bool, int,
                                       size_t);
  typedef dft_fields (fields::*old_fields)(component *, int, const volume, const double *, size_t,
                                           bool, int, bool);
  typedef dft_fields (fields::*new_fields)(component *, int, const volume, const double *, size_t,
                                           bool, int, bool, size_t);
  typedef dft_energy (fields::*old_energy)(const volume_list *, const double *, size_t, int);
  typedef dft_energy (fields::*new_energy)(const volume_list *, const double *, size_t, int, size_t);
  typedef dft_force (fields::*old_force)(const volume_list *, const double *, size_t, int);
  typedef dft_force (fields::*new_force)(const volume_list *, const double *, size_t, int, size_t);
  typedef dft_near2far (fields::*old_n2f)(const volume_list *, const double *, size_t, int, int);
  typedef dft_near2far (fields::*new_n2f)(const volume_list *, const double *, size_t, int, int,
                                          size_t);

  old_add_dft a0 = static_cast<old_add_dft>(&fields::add_dft);
  new_add_dft a1 = static_cast<new_add_dft>(&fields::add_dft);
  old_flux f0 = static_cast<old_flux>(&fields::add_dft_flux);
  new_flux f1 = static_cast<new_flux>(&fields::add_dft_flux);
  old_mode m0 = static_cast<old_mode>(&fields::add_mode_monitor);
  new_mode m1 = static_cast<new_mode>(&fields::add_mode_monitor);
  old_fields d0 = static_cast<old_fields>(&fields::add_dft_fields);
  new_fields d1 = static_cast<new_fields>(&fields::add_dft_fields);
  old_energy e0 = static_cast<old_energy>(&fields::add_dft_energy);
  new_energy e1 = static_cast<new_energy>(&fields::add_dft_energy);
  old_force s0 = static_cast<old_force>(&fields::add_dft_force);
  new_force s1 = static_cast<new_force>(&fields::add_dft_force);
  old_n2f n0 = static_cast<old_n2f>(&fields::add_dft_near2far);
  new_n2f n1 = static_cast<new_n2f>(&fields::add_dft_near2far);
  require(a0 && a1 && f0 && f1 && m0 && m1 && d0 && d1 && e0 && e1 && s0 && s1 && n0 && n1,
          "old/new exported overloads were not linkable");
}

int main(int argc, char **argv) {
  initialize mpi(argc, argv);
  verbosity = 0;

  test_exported_overloads_link();
  test_analytic_tail(Ex, {1.0}, {0.72}, 8, 13, 0.17, 1.25, 3);
  test_analytic_tail(Ex, {1.2, -0.35}, {0.61, -0.28}, 10, 15, 0.09, 0.8, 2);
  // Magnetic monitors are sampled half a Yee timestep before electric ones.
  test_analytic_tail(Hy, {0.9}, {0.68}, 8, 12, 0.13, 1.5 - 0.5 / 20.0, 1);
  test_rejected_tail(1.0 - 1.0e-10, "near-pole fit did not fail closed to the raw DFT");
  test_rejected_tail(1.0 - 1.0e-7, "unstable amplified tail did not fail closed to the raw DFT");
  test_diagnostic_checkpoint_isolation();
  test_lifecycle_and_bounded_storage();
  return 0;
}
