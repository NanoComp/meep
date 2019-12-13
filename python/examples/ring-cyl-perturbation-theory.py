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

    fcen = 0.18         # pulse center frequency
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
    df = 0.01

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

    npts = 50
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

        # section for fields perpendicular to interface (Dr)
        d_r_field = abs(sim.get_field_point(mp.Dr, point))**2
        perpendicular_field = d_r_field
        # fields have to be multiplied by Δ(1/ε) and ε**2 to get integrand
        perpendicular_field_integrand = deps_inv_inner * perpendicular_field
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

        # section for fields perpendicular to interface (Dr)
        d_r_field = abs(sim.get_field_point(mp.Dr, point))
        perpendicular_field = d_r_field**2
        # fields have to be multiplied by Δ(1/ε) and ε**2 to get integrand
        perpendicular_field_integrand = deps_inv_outer * perpendicular_field
        perpendicular_fields_outer.append(perpendicular_field_integrand)

    numerator_surface_integral = (np.sum(parallel_fields_outer)-np.sum(perpendicular_fields_outer))*2*np.pi*b/npts \
                                 + (np.sum(parallel_fields_inner)-np.sum(perpendicular_fields_inner))*2*np.pi*a/npts
    denominator_surface_integral = sim.electric_energy_in_box(center=mp.Vector3(), size=mp.Vector3(2*sr))
    perturb_theory_dw_dr = -harminv_freq_at_r * numerator_surface_integral / (2 * denominator_surface_integral)

    # perturbed geometry with narrowband source

    dr = args.dr

    sim.reset_meep()
    a = r + dr  # inner radius of ring
    b = a + w  # outer radius of ring

    fcen = perturb_theory_dw_dr * dr + harminv_freq_at_r
    df = 0.01

    sources = [mp.Source(mp.GaussianSource(fcen, fwidth=df),
                         component=component,
                         center=mp.Vector3(r + dr + 0.1))]

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
    sim.run(mp.after_sources(h),
            until_after_sources=200)

    harminv_freq_at_r_plus_dr = h.modes[0].freq

    finite_diff_dw_dr = (harminv_freq_at_r_plus_dr - harminv_freq_at_r) / dr

    print("component:, {}".format(component_name))
    print("res:, {}".format(resolution))
    print("dr:, {}".format(dr))
    print("w:, {} (unperturbed), {} (perturbed)".format(harminv_freq_at_r, harminv_freq_at_r_plus_dr))
    print("dwdR:, {} (pert. theory), {} (finite diff.)".format(perturb_theory_dw_dr, finite_diff_dw_dr))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-perp', type=bool, default=False, help='controls if in-plane perpendicular fields are excited or not. (default: False, for Ez source)')
    parser.add_argument('-res', type=int, default=100, help='resolution (default: 100 pixels/um)')
    parser.add_argument('-dr', type=float, default=0.02, help='dr (default: 0.02)')
    args = parser.parse_args()
    main(args)
