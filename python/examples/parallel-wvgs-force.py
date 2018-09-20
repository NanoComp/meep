import meep as mp
import argparse

def main(args):

    resolution = 30
    nSi = 3.45
    Si = mp.Medium(index=nSi)
    dpml = 1.0
    sx = 5
    sy = 3

    cell = mp.Vector3(sx + 2 * dpml, sy + 2 * dpml)
    pml_layers = mp.PML(dpml)

    a = 1.0     # waveguide width
    s = args.s  # waveguide separation distance

    geometry = [mp.Block(center=mp.Vector3(-0.5 * (s + a)),
                         size=mp.Vector3(a, a, mp.inf),
                         material=Si),
                mp.Block(center=mp.Vector3(0.5 * (s + a)),
                         size=mp.Vector3(a, a, mp.inf),
                         material=Si)]

    xodd = args.xodd
    symmetries = [mp.Mirror(mp.X, phase=-1.0 if xodd else 1.0),
                  mp.Mirror(mp.Y, phase=-1.0)]

    beta = 0.5
    k_point = mp.Vector3(z=beta)

    fcen = 0.22
    df = 0.06
    sources = [mp.Source(src=mp.GaussianSource(fcen, fwidth=df), component=mp.Ey,
                         center=mp.Vector3(-0.5 * (s + a)), size=mp.Vector3(a, a)),
               mp.Source(src=mp.GaussianSource(fcen, fwidth=df), component=mp.Ey,
                         center=mp.Vector3(0.5 * (s + a)), size=mp.Vector3(a, a),
                         amplitude=-1.0 if xodd else 1.0)]

    sim = mp.Simulation(resolution=resolution,
                        cell_size=cell,
                        boundary_layers=[pml_layers],
                        geometry=geometry,
                        symmetries=symmetries,
                        k_point=k_point,
                        sources=sources)

    h = mp.Harminv(mp.Ey, mp.Vector3(0.5 * (s + a)), fcen, df)

    sim.run(mp.after_sources(h), until_after_sources=200)

    f = h.modes[0].freq
    print("freq:, {}, {}".format(s, f))

    sim.reset_meep()

    new_sources = [
        mp.EigenModeSource(src=mp.GaussianSource(f, fwidth=df), component=mp.Ey,
                           size=mp.Vector3(a, a), center=mp.Vector3(-0.5 * (s + a)),
                           eig_kpoint=k_point, eig_match_freq=True, eig_parity=mp.ODD_Y),
        mp.EigenModeSource(src=mp.GaussianSource(f, fwidth=df), component=mp.Ey,
                           size=mp.Vector3(a, a), center=mp.Vector3(0.5 * (s + a)),
                           eig_kpoint=k_point, eig_match_freq=True, eig_parity=mp.ODD_Y,
                           amplitude=-1.0 if xodd else 1.0)
    ]

    sim.change_sources(new_sources)

    flx_reg = mp.FluxRegion(direction=mp.Z, center=mp.Vector3(), size=mp.Vector3(1.2 * (2 * a + s), 1.2 * a))
    wvg_pwr = sim.add_flux(f, 0, 1, flx_reg)

    frc_reg1 = mp.ForceRegion(mp.Vector3(0.5 * s), direction=mp.X, weight=1.0, size=mp.Vector3(y=a))
    frc_reg2 = mp.ForceRegion(mp.Vector3(0.5 * s + a), direction=mp.X, weight=-1.0, size=mp.Vector3(y=a))
    wvg_force = sim.add_force(f, 0, 1, frc_reg1, frc_reg2)

    runtime = 5000
    sim.run(until_after_sources=runtime)
    sim.display_fluxes(wvg_pwr)
    sim.display_forces(wvg_force)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-xodd', type=bool, default=True, help='odd mirror-symmetry plane in the X direction?')
    parser.add_argument('-s', type=float, default=1.0, help='waveguide separation distance')
    args = parser.parse_args()
    main(args)
