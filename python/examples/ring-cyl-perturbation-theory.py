import meep as mp
import numpy as np
import argparse


def main(args):
    n = 3.4                 # index of waveguide
    r = 1
    a = r                   # inner radius of ring
    w = 1                   # width of waveguide
    b = a + w               # outer radius of ring
    pad = 4                 # padding between waveguide and edge of PML

    dpml = 2                # thickness of PML
    pml_layers = [mp.PML(dpml)]

    resolution = args.res

    sr = b + pad + dpml            # radial size (cell is from 0 to sr)
    dimensions = mp.CYLINDRICAL    # coordinate system is (r,phi,z) instead of (x,y,z)
    cell = mp.Vector3(sr, 0, 0)

    geometry = [mp.Block(center=mp.Vector3(a + (w / 2)),
                         size=mp.Vector3(w, 1e20, 1e20),
                         material=mp.Medium(index=n))]

    # Finding a resonance mode with a high Q-value (calculated with Harminv)

    fcen = 0.15         # pulse center frequency
    df = 0.1            # pulse width (in frequency)

    m = 5

    if args.perp:
        component = mp.Hz
        component_name = 'meep.Hz'
    else:
        component = mp.Ez
        component_name = 'meep.Ez'

    sources = [mp.Source(mp.GaussianSource(fcen, fwidth=df),
                         component=component,
                         center=mp.Vector3(r+0.1),
                         amplitude=1)]

    sim = mp.Simulation(cell_size=cell,
                        geometry=geometry,
                        boundary_layers=pml_layers,
                        resolution=resolution,
                        sources=sources,
                        dimensions=dimensions,
                        m=m)

    h = mp.Harminv(component, mp.Vector3(r+0.1), fcen, df)
    sim.run(mp.after_sources(h),
            until_after_sources=200)

    harminv_freq_at_r = h.modes[0].freq

    sim.reset_meep()

    # unperturbed geometry with narrowband source centered at resonant frequency

    fcen = harminv_freq_at_r
    df = 0.05 * fcen

    sources = [mp.Source(mp.GaussianSource(fcen, fwidth=df),
                         component=component,
                         center=mp.Vector3(r + 0.1),
                         amplitude=1)]

    sim = mp.Simulation(cell_size=cell,
                        geometry=geometry,
                        boundary_layers=pml_layers,
                        resolution=resolution,
                        sources=sources,
                        dimensions=dimensions,
                        m=m)

    sim.run(until_after_sources=200)

    # now need to calculate the surface integrals that go into dw/dR. Fields parallel and perpendicular to the interface
    # AND at the inner and outer surfaces are treated differently, so each will be calculated separately.

    npts = 100
    angles = 2 * np.pi / npts * np.arange(npts)
    deps_inner = 1 - n ** 2
    deps_inv_inner = 1 - 1/(n**2)
    deps_outer = n ** 2 - 1
    deps_inv_outer = -1 + 1 / (n ** 2)

    parallel_fields_inner = []
    perpendicular_fields_inner = []
    parallel_fields_outer = []
    perpendicular_fields_outer = []

    for angle in angles:
        # section for fields at inner surface
        point = mp.Vector3(a, angle)

        # section for fields parallel to interface (Ez and Ep)
        e_z_field = abs(sim.get_field_point(mp.Ez, point))**2
        e_p_field = abs(sim.get_field_point(mp.Ep, point))**2
        parallel_field = e_z_field + e_p_field
        # fields have to be multiplied by Δε to get integrand
        parallel_field_integrand = deps_inner * parallel_field
        parallel_fields_inner.append(parallel_field_integrand)

        # section for fields perpendicular to interface (Er)
        e_r_field = abs(sim.get_field_point(mp.Er, point))**2
        perpendicular_field = e_r_field
        # fields have to be multiplied by Δ(1/ε) and ε**2 to get integrand
        perpendicular_field_integrand = deps_inv_inner * (abs(sim.get_epsilon_point(point, harminv_freq_at_r))**2) * perpendicular_field
        perpendicular_fields_inner.append(perpendicular_field_integrand)

        # section for fields at outer surface
        point = mp.Vector3(b, angle)

        # section for fields parallel to interface (Ez and Ep)
        e_z_field = abs(sim.get_field_point(mp.Ez, point))**2
        e_p_field = abs(sim.get_field_point(mp.Ep, point))**2
        parallel_field = e_z_field + e_p_field
        # fields have to be multiplied by Δε to get integrand
        parallel_field_integrand = deps_outer * parallel_field
        parallel_fields_outer.append(parallel_field_integrand)

        # section for fields perpendicular to interface (Er)
        e_r_field = abs(sim.get_field_point(mp.Er, point))
        perpendicular_field = e_r_field**2
        # fields have to be multiplied by Δ(1/ε) and ε**2 to get integrand
        perpendicular_field_integrand = deps_inv_outer * (abs(sim.get_epsilon_point(point, harminv_freq_at_r))**2) * perpendicular_field
        perpendicular_fields_outer.append(perpendicular_field_integrand)

    # numerator_surface_integral = 2 * np.pi * b * (mean([mean(parallel_fields_inner), mean(parallel_fields_outer)]) - mean([mean(perpendicular_fields_inner), mean(perpendicular_fields_outer)]))
    numerator_surface_integral = (np.sum(parallel_fields_outer)-np.sum(perpendicular_fields_outer))*2*np.pi*b/npts \
                                 + (np.sum(parallel_fields_inner)-np.sum(perpendicular_fields_inner))*2*np.pi*a/npts
    denominator_surface_integral = sim.electric_energy_in_box(center=mp.Vector3(), size=mp.Vector3(2*sr))
    perturb_theory_dw_dr = -harminv_freq_at_r * numerator_surface_integral / (2 * denominator_surface_integral)

    # perturbed geometry with narrowband source

    dr = 0.02

    sim.reset_meep()
    a = r + dr  # inner radius of ring
    b = a + w  # outer radius of ring

    sources = [mp.Source(mp.GaussianSource(fcen, fwidth=df),
                         component=component,
                         center=mp.Vector3(r + dr + 0.1),
                         amplitude=1)]

    geometry = [mp.Block(center=mp.Vector3(a + (w / 2)),
                         size=mp.Vector3(w, mp.inf, mp.inf),
                         material=mp.Medium(index=n))]

    sim = mp.Simulation(cell_size=cell,
                        geometry=geometry,
                        boundary_layers=pml_layers,
                        resolution=resolution,
                        sources=sources,
                        dimensions=dimensions,
                        m=m)

    h = mp.Harminv(mp.Ez, mp.Vector3(r + dr + 0.1), fcen, df)
    sim.run(mp.after_sources(h),
            until_after_sources=100)

    harminv_freq_at_r_plus_dr = h.modes[0].freq

    finite_diff_dw_dr = (harminv_freq_at_r_plus_dr - harminv_freq_at_r) / dr

    print("res: {}, dr: {}\ndwdR:, {} (pert. theory), {} (finite diff.)".format(resolution, dr, perturb_theory_dw_dr, finite_diff_dw_dr))

    # The rest of main() is not explicitly included in the tutorial, but this shows you how we did
    # some error calculations
    # center_diff_dw_dr = []
    # harminv_freqs_at_r_plus_dr = []
    #
    # drs = np.logspace(start=-3, stop=-1, num=10)
    #
    # for dr in drs:
    #     sim.reset_meep()
    #     w = 1 + dr  # width of waveguide
    #     b = a + w
    #     print(f'The current dr is dr={dr}')
    #     if len(harminv_freqs_at_r_plus_dr) == 0:
    #         fcen = harminv_freq_at_r
    #     else:
    #         fcen = harminv_freqs_at_r_plus_dr[-1]
    #     df = 0.01
    #
    #     sources = [mp.Source(mp.GaussianSource(fcen, fwidth=df), component, mp.Vector3(r + 0.1))]
    #
    #     geometry = [mp.Block(center=mp.Vector3(a + (w / 2)),
    #                          size=mp.Vector3(w, 1e20, 1e20),
    #                          material=mp.Medium(index=n))]
    #
    #     sim = mp.Simulation(cell_size=cell,
    #                         geometry=geometry,
    #                         boundary_layers=pml_layers,
    #                         resolution=resolution,
    #                         sources=sources,
    #                         dimensions=dimensions,
    #                         m=m)
    #
    #     h = mp.Harminv(component, mp.Vector3(r + 0.1), fcen, df)
    #     sim.run(mp.after_sources(h), until_after_sources=200)
    #
    #     harminv_freq_at_r_plus_dr = h.modes[0].freq
    #     harminv_freqs_at_r_plus_dr.append(harminv_freq_at_r_plus_dr)
    #
    #     dw_dr = (harminv_freq_at_r_plus_dr - harminv_freq_at_r) / dr
    #     center_diff_dw_dr.append(dw_dr)
    #
    # relative_errors_dw_dr = [abs((dw_dR - perturb_theory_dw_dr) / dw_dR) for dw_dR in center_diff_dw_dr]
    #
    # perturb_predicted_freqs_at_r_plus_dr = [dr * perturb_theory_dw_dr + harminv_freq_at_r for dr in drs]
    # relative_errors_freqs_at_r_plus_dr = [
    #     abs((perturb_predicted_freqs_at_r_plus_dr[i] - harminv_freqs_at_r_plus_dr[i]) / harminv_freqs_at_r_plus_dr[i])
    #     for i in range(len(harminv_freqs_at_r_plus_dr))]
    #
    # results_string  = '\ncomponent={}\n'.format(component_name)
    # results_string += 'perturb_theory_dw_dr={}\n'.format(perturb_theory_dw_dr)
    # results_string += 'drs={}\n'.format(drs)
    # results_string += 'center_diff_dw_dr={}\n'.format(center_diff_dw_dr)
    # results_string += 'harminv_freqs_at_r_plus_dr={}\n'.format(harminv_freqs_at_r_plus_dr)
    # results_string += 'relative_errors_dw_dr={}\n'.format(relative_errors_dw_dr)
    # results_string += 'relative_errors_freqs_at_r_plus_dr={}'.format(relative_errors_freqs_at_r_plus_dr)
    #
    # if mp.am_master():
    #     f = open('ring_cyl_perturbation_theory.dat', 'a')
    #     f.write(results_string)
    #     f.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-perp', type=bool, default=False, help='controls if in-plane perpendicular fields are excited or not. (default: False, for Ez source)')
    parser.add_argument('-res', type=int, default=100, help='resolution (default: 100 pixels/um)')
    args = parser.parse_args()
    main(args)
