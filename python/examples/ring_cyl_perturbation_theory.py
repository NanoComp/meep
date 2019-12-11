from __future__ import division

import meep as mp
import numpy as np
from statistics import mean
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

    resolution = 100

    sr = b + pad + dpml            # radial size (cell is from 0 to sr)
    dimensions = mp.CYLINDRICAL    # coordinate system is (r,phi,z) instead of (x,y,z)
    cell = mp.Vector3(sr, 0, 0)

    m = 4

    geometry = [mp.Block(center=mp.Vector3(a + (w / 2)),
                         size=mp.Vector3(w, 1e20, 1e20),
                         material=mp.Medium(index=n))]

    # Finding a resonance mode with a high Q-value (calculated with Harminv)

    fcen = 0.15         # pulse center frequency
    df = 0.1            # pulse width (in frequency)

    if (args.perpendicular):
        component = mp.Hz
    else:
        component = mp.Ez

    sources = [mp.Source(mp.GaussianSource(fcen, fwidth=df), component=component, mp.Vector3(r+0.1), amplitude=1)]

    sim = mp.Simulation(cell_size=cell,
                        geometry=geometry,
                        boundary_layers=pml_layers,
                        resolution=resolution,
                        sources=sources,
                        dimensions=dimensions,
                        m=m)

    h = mp.Harminv(c=component, mp.Vector3(r+0.1), fcen, df)
    sim.run(mp.after_sources(h), until_after_sources=200)

    Harminv_freq_at_R = h.modes[0].freq

    sim.reset_meep()

    # now running the simulation that will be used with perturbation theory to calculate dw/dR

    fcen = Harminv_freq_at_R
    df = 0.01

    sources = [mp.Source(mp.GaussianSource(fcen, fwidth=df), component=component, mp.Vector3(r+0.1), amplitude=1)]

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
    npts_inner = 100
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
    perturb_theory_dw_dR = -Harminv_freq_at_R * numerator_surface_integral / (4 * denominator_surface_integral)

    center_diff_dw_dR = []
    Harminv_freqs_at_R_plus_dR = []

    drs = np.logspace(start=-3, stop=-1, num=10)

    for dr in drs:
        sim.reset_meep()
        w = 1 + dr  # width of waveguide
        b = a + w
        print(f'The current dr is dr={dr}')
        if len(Harminv_freqs_at_R_plus_dR) == 0:
            fcen = Harminv_freq_at_R
        else:
            fcen = Harminv_freqs_at_R_plus_dR[-1]
        df = 0.01

        sources = [mp.Source(mp.GaussianSource(fcen, fwidth=df), component=component, mp.Vector3(r + 0.1))]

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

        h = mp.Harminv(c=component, mp.Vector3(r + 0.1), fcen, df)
        sim.run(mp.after_sources(h), until_after_sources=200)

        Harminv_freq_at_R_plus_dR = h.modes[0].freq
        Harminv_freqs_at_R_plus_dR.append(Harminv_freq_at_R_plus_dR)

        dw_dR = (Harminv_freq_at_R_plus_dR - Harminv_freq_at_R) / dr
        center_diff_dw_dR.append(dw_dR)

    relative_errors_dw_dR = [abs((dw_dR - perturb_theory_dw_dR) / dw_dR) for dw_dR in center_diff_dw_dR]

    perturb_predicted_freqs_at_R_plus_dR = [dr * perturb_theory_dw_dR + Harminv_freq_at_R for dr in drs]
    relative_errors_freqs_at_R_plus_dR = [
        abs((perturb_predicted_freqs_at_R_plus_dR[i] - Harminv_freqs_at_R_plus_dR[i]) / Harminv_freqs_at_R_plus_dR[i])
        for i in range(len(Harminv_freqs_at_R_plus_dR))]

    results_string = 'component={}\ndrs={}\nrelative_errors_dw_dR={}\nrelative_errors_freqs_at_R_plus_dR={}'.format(component,drs,relative_errors_dw_dR,relative_errors_freqs_at_R_plus_dR)

    if mp.am_master():
        f = open('ring_cyl_perturbation_theory.dat', 'a')
        f.write(results_string)
        f.close

        #plt.figure(dpi=150)
        #plt.loglog(drs, relative_errors_dw_dR, 'bo-', label='relative error')
        #plt.grid(True, which='both', ls='-')
        #plt.xlabel('perturbation amount $dR$')
        #plt.ylabel('relative error between $dω/dR$')
        #plt.legend(loc='upper right')
        #plt.title('Comparison of Perturbation Theory and \nCenter-Difference Calculations in Finding $dω/dR$')
        #plt.tight_layout()
        # plt.show()
        #plt.savefig('ring_Hz_perturbation_theory.dw_dR_error.png')
        #plt.clf()

        #plt.figure(dpi=150)
        #plt.loglog(drs, relative_errors_freqs_at_R_plus_dR, 'bo-', label='relative error')
        #plt.grid(True, which='both', ls='-')
        #plt.xlabel('perturbation amount $dR$')
        #plt.ylabel('relative error between $ω(R+dR)$')
        #plt.legend(loc='upper left')
        #plt.title('Comparison of resonance frequencies at $R+dR$ predicted by\nperturbation theory and found with Harminv')
        #plt.tight_layout()
        # plt.show()
        #plt.savefig('ring_Hz_perturbation_theory.freqs_error.png')
        #plt.clf()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-perpendicular', type=bool, default=False, help='True for perpendicular field source, False for parallel field source. False default.' )
    args = parser.parse_args()
    main(args)