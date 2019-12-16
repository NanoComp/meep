from __future__ import division

import meep as mp
import numpy as np
from statistics import mean
from math import floor
import matplotlib.pyplot as plt
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

    resolution = 10

    sr = b + pad + dpml            # radial size (cell is from 0 to sr)
    dimensions = mp.CYLINDRICAL    # coordinate system is (r,phi,z) instead of (x,y,z)
    cell = mp.Vector3(sr, 0, 0)

    m = 4

    geometry = [mp.Block(center=mp.Vector3(a + (w / 2)),
                         size=mp.Vector3(w, 1e20, 1e20),
                         material=mp.Medium(index=n))]

    # Finding a resonance mode with a high Q-value (calculated with Harminv)

    fcen = 0.18         # pulse center frequency
    df = 0.1            # pulse width (in frequency)

    if (args.perpendicular):
        component = mp.Hz
        component_name = 'meep.Hz'
    else:
        component = mp.Ez
        component_name = 'meep.Ez'

    sources = [mp.Source(mp.GaussianSource(fcen, fwidth=df), component, mp.Vector3(r+0.1), amplitude=1)]

    sim = mp.Simulation(cell_size=cell,
                        geometry=geometry,
                        boundary_layers=pml_layers,
                        resolution=resolution,
                        sources=sources,
                        dimensions=dimensions,
                        m=m)

    h = mp.Harminv(component, mp.Vector3(r+0.1), fcen, df)
    sim.run(mp.after_sources(h), until_after_sources=200)

    Harminv_freq_at_R = h.modes[0].freq

    sim.reset_meep()

    # now running the simulation that will be used with perturbation theory to calculate dw/dR

    fcen = Harminv_freq_at_R
    df = 0.01

    sources = [mp.Source(mp.GaussianSource(fcen, fwidth=df), component, mp.Vector3(r+0.1), amplitude=1)]

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

    # section for fields at inner surface
    npts_inner = 10
    angles_inner = 2 * np.pi / npts_inner * np.arange(npts_inner)
    deps_inner = 1 - n ** 2
    deps_inv_inner = 1 - 1/(n**2)

    # section for fields parallel to interface (Ez and Ep)
    parallel_fields_inner = []
    for angle in angles_inner:
        point = mp.Vector3(a, angle)
        e_z_field = abs(sim.get_field_point(mp.Ez, point))**2
        e_p_field = abs(sim.get_field_point(mp.Ep, point))**2
        e_parallel_field = e_z_field + e_p_field
        # fields have to be multiplied by Δε
        e_parallel_field = deps_inner * e_parallel_field
        parallel_fields_inner.append(e_parallel_field)

    # section for fields perpendicular to interface (Er)
    perpendicular_fields_inner = []
    for angle in angles_inner:
        point = mp.Vector3(a, angle)
        e_r_field = abs(sim.get_field_point(mp.Er, point))**2
        e_perpendicular_field = e_r_field
        # fields have to be multiplied by Δ(1/ε) and ε**2
        e_perpendicular_field = deps_inv_inner * (abs(sim.get_epsilon_point(point, Harminv_freq_at_R))**2) * e_perpendicular_field
        perpendicular_fields_inner.append(e_perpendicular_field)

    # section for fields at outer surface
    npts_outer = npts_inner
    angles_outer = 2 * np.pi / npts_outer * np.arange(npts_outer)
    deps_outer = n ** 2 - 1
    deps_inv_outer = -1 + 1/(n**2)

    # section for fields parallel to interface (Ez and Ep)
    parallel_fields_outer = []
    for angle in angles_outer:
        point = mp.Vector3(b, angle)
        e_z_field = abs(sim.get_field_point(mp.Ez, point))**2
        e_p_field = abs(sim.get_field_point(mp.Ep, point))**2
        e_parallel_field = e_z_field + e_p_field
        # fields have to be multiplied by Δε
        e_parallel_field = deps_outer * e_parallel_field
        parallel_fields_outer.append(e_parallel_field)

    # section for fields perpendicular to interface (Er)
    perpendicular_fields_outer = []
    for angle in angles_inner:
        point = mp.Vector3(b, angle)
        e_r_field = abs(sim.get_field_point(mp.Er, point))
        e_perpendicular_field = e_r_field**2
        # fields have to be multiplied by Δ(1/ε) and ε**2
        e_perpendicular_field = deps_inv_outer * (abs(sim.get_epsilon_point(point, Harminv_freq_at_R))**2) * e_perpendicular_field
        perpendicular_fields_outer.append(e_perpendicular_field)

    numerator_surface_integral = 2 * np.pi * b * (mean([mean(parallel_fields_inner), mean(parallel_fields_outer)]) - mean([mean(perpendicular_fields_inner), mean(perpendicular_fields_outer)]))
    denominator_surface_integral = sim.electric_energy_in_box(center=mp.Vector3((b + pad/2) / 2), size=mp.Vector3(b + pad/2))
    perturb_theory_dw_dR_at_R = -Harminv_freq_at_R * numerator_surface_integral / (4 * denominator_surface_integral)

    center_diff_dw_dR_at_R_plus_dR = []
    perturb_theory_dw_dR_at_R_plus_dR = []
    Harminv_freqs_at_R_plus_dR = []

    resolutions = np.logspace(start=2, stop=2.5, num=1)
    resolutions = [floor(resolution) for resolution in resolutions]

    for resolution in resolutions:
        sim.reset_meep()
        dr = 1e-2
        a = a + dr
        w = 1  # width of waveguide
        b = a + w
        print(f'The current resolution is resolution={resolution}')
        if len(Harminv_freqs_at_R_plus_dR) == 0:
            fcen = Harminv_freq_at_R
        else:
            fcen = Harminv_freqs_at_R_plus_dR[-1]
        df = 0.01

        sources = [mp.Source(mp.GaussianSource(fcen, fwidth=df),
                             component,
                             mp.Vector3(r + dr + 0.1))]

        geometry = [mp.Block(center=mp.Vector3(a + (w / 2)),
                             size=mp.Vector3(w, 1e20, 1e20),
                             material=mp.Medium(index=n))]

        sim = mp.Simulation(cell_size=cell,
                            geometry=geometry,
                            boundary_layers=pml_layers,
                            resolution=resolution,
                            sources=sources,
                            dimensions=dimensions,
                            m=m)

        h = mp.Harminv(component, mp.Vector3(r + dr + 0.1), fcen, df)
        sim.run(mp.after_sources(h), until_after_sources=200)

        Harminv_freq_at_R_plus_dR = h.modes[0].freq
        Harminv_freqs_at_R_plus_dR.append(Harminv_freq_at_R_plus_dR)

        dw_dR = (Harminv_freq_at_R_plus_dR - Harminv_freq_at_R) / dr
        center_diff_dw_dR_at_R_plus_dR.append(dw_dR)

        # now need to calculate the surface integrals that go into dw/dR. Fields parallel and perpendicular to the interface
        # AND at the inner and outer surfaces are treated differently, so each will be calculated separately.

        # section for fields at inner surface
        npts_inner = resolution
        angles_inner = 2 * np.pi / npts_inner * np.arange(npts_inner)
        deps_inner = 1 - n ** 2
        deps_inv_inner = 1 - 1 / (n ** 2)

        # section for fields parallel to interface (Ez and Ep)
        parallel_fields_inner = []
        for angle in angles_inner:
            point = mp.Vector3(a, angle)
            e_z_field = abs(sim.get_field_point(mp.Ez, point)) ** 2
            e_p_field = abs(sim.get_field_point(mp.Ep, point)) ** 2
            e_parallel_field = e_z_field + e_p_field
            # fields have to be multiplied by Δε
            e_parallel_field = deps_inner * e_parallel_field
            parallel_fields_inner.append(e_parallel_field)

        # section for fields perpendicular to interface (Er)
        perpendicular_fields_inner = []
        for angle in angles_inner:
            point = mp.Vector3(a, angle)
            e_r_field = abs(sim.get_field_point(mp.Er, point)) ** 2
            e_perpendicular_field = e_r_field
            # fields have to be multiplied by Δ(1/ε) and ε**2
            e_perpendicular_field = deps_inv_inner * (
                        abs(sim.get_epsilon_point(point, Harminv_freq_at_R)) ** 2) * e_perpendicular_field
            perpendicular_fields_inner.append(e_perpendicular_field)

        # section for fields at outer surface
        npts_outer = npts_inner
        angles_outer = 2 * np.pi / npts_outer * np.arange(npts_outer)
        deps_outer = n ** 2 - 1
        deps_inv_outer = -1 + 1 / (n ** 2)

        # section for fields parallel to interface (Ez and Ep)
        parallel_fields_outer = []
        for angle in angles_outer:
            point = mp.Vector3(b, angle)
            e_z_field = abs(sim.get_field_point(mp.Ez, point)) ** 2
            e_p_field = abs(sim.get_field_point(mp.Ep, point)) ** 2
            e_parallel_field = e_z_field + e_p_field
            # fields have to be multiplied by Δε
            e_parallel_field = deps_outer * e_parallel_field
            parallel_fields_outer.append(e_parallel_field)

        # section for fields perpendicular to interface (Er)
        perpendicular_fields_outer = []
        for angle in angles_inner:
            point = mp.Vector3(b, angle)
            e_r_field = abs(sim.get_field_point(mp.Er, point))
            e_perpendicular_field = e_r_field ** 2
            # fields have to be multiplied by Δ(1/ε) and ε**2
            e_perpendicular_field = deps_inv_outer * (
                        abs(sim.get_epsilon_point(point, Harminv_freq_at_R)) ** 2) * e_perpendicular_field
            perpendicular_fields_outer.append(e_perpendicular_field)

        numerator_surface_integral = 2 * np.pi * b * (
                    mean([mean(parallel_fields_inner), mean(parallel_fields_outer)]) - mean(
                [mean(perpendicular_fields_inner), mean(perpendicular_fields_outer)]))
        denominator_surface_integral = sim.electric_energy_in_box(center=mp.Vector3((b + pad / 2) / 2),
                                                                  size=mp.Vector3(b + pad / 2))
        temp_perturb_theory_dw_dR_at_R = -Harminv_freq_at_R * numerator_surface_integral / (4 * denominator_surface_integral)
        perturb_theory_dw_dR_at_R_plus_dR.append(temp_perturb_theory_dw_dR_at_R)

    relative_errors_dw_dR = [abs((center_diff_dw_dR_at_R_plus_dR[i] - perturb_theory_dw_dR_at_R_plus_dR[i]) / center_diff_dw_dR_at_R_plus_dR[i]) for i in range(len(resolutions))]

    results_string = '\ncomponent={}\nperturb_theory_dw_dR_at_R={}\nresolutions={}\ndr={}\nperturb_theory_dw_dR_at_R_plus_dR={}\ncenter_diff_dw_dR={}\nHarminv_freqs_at_R_plus_dR={}\nrelative_errors_dw_dR={}\n'.format(component_name,perturb_theory_dw_dR_at_R,resolutions,dr,perturb_theory_dw_dR_at_R_plus_dR,center_diff_dw_dR_at_R_plus_dR,Harminv_freqs_at_R_plus_dR,relative_errors_dw_dR)

    if mp.am_master():
        f = open('ring_cyl_perturbation_theory_covergence_test.dat', 'a')
        f.write(results_string)
        f.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-perpendicular', type=bool, default=False, help='True for perpendicular field source, False for parallel field source. False default.' )
    args = parser.parse_args()
    main(args)