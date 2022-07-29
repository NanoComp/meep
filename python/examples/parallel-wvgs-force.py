import matplotlib.pyplot as plt
import numpy as np

import meep as mp

resolution = 40  # pixels/μm

Si = mp.Medium(index=3.45)

dpml = 1.0
pml_layers = [mp.PML(dpml)]

sx = 5
sy = 3
cell = mp.Vector3(sx + 2 * dpml, sy + 2 * dpml, 0)

a = 1.0  # waveguide width/height

k_point = mp.Vector3(z=0.5)


def parallel_waveguide(s, xodd):
    geometry = [
        mp.Block(
            center=mp.Vector3(-0.5 * (s + a)),
            size=mp.Vector3(a, a, mp.inf),
            material=Si,
        ),
        mp.Block(
            center=mp.Vector3(0.5 * (s + a)), size=mp.Vector3(a, a, mp.inf), material=Si
        ),
    ]

    symmetries = [mp.Mirror(mp.X, phase=-1 if xodd else 1), mp.Mirror(mp.Y, phase=-1)]

    sim = mp.Simulation(
        resolution=resolution,
        cell_size=cell,
        geometry=geometry,
        boundary_layers=pml_layers,
        symmetries=symmetries,
        k_point=k_point,
    )

    sim.init_sim()
    EigenmodeData = sim.get_eigenmode(
        0.22,
        mp.Z,
        mp.Volume(center=mp.Vector3(), size=mp.Vector3(sx, sy)),
        2 if xodd else 1,
        k_point,
        match_frequency=False,
        parity=mp.ODD_Y,
    )

    fcen = EigenmodeData.freq
    print(f'freq:, {"xodd" if xodd else "xeven"}, {s}, {fcen}')

    sim.reset_meep()

    eig_sources = [
        mp.EigenModeSource(
            src=mp.GaussianSource(fcen, fwidth=0.1 * fcen),
            size=mp.Vector3(sx, sy),
            center=mp.Vector3(),
            eig_band=2 if xodd else 1,
            eig_kpoint=k_point,
            eig_match_freq=False,
            eig_parity=mp.ODD_Y,
        )
    ]

    sim.change_sources(eig_sources)

    flux_reg = mp.FluxRegion(
        direction=mp.Z, center=mp.Vector3(), size=mp.Vector3(sx, sy)
    )
    wvg_flux = sim.add_flux(fcen, 0, 1, flux_reg)

    force_reg1 = mp.ForceRegion(
        mp.Vector3(0.49 * s), direction=mp.X, weight=1, size=mp.Vector3(y=sy)
    )
    force_reg2 = mp.ForceRegion(
        mp.Vector3(0.5 * s + 1.01 * a), direction=mp.X, weight=-1, size=mp.Vector3(y=sy)
    )
    wvg_force = sim.add_force(fcen, 0, 1, force_reg1, force_reg2)

    sim.run(until_after_sources=1500)

    flux = mp.get_fluxes(wvg_flux)[0]
    force = mp.get_forces(wvg_force)[0]
    print(
        f'data:, {"xodd" if xodd else "xeven"}, {s}, {flux}, {force}, {-force / flux}'
    )

    sim.reset_meep()
    return flux, force


s = np.arange(0.05, 1.05, 0.05)
fluxes_odd = np.zeros(s.size)
forces_odd = np.zeros(s.size)
fluxes_even = np.zeros(s.size)
forces_even = np.zeros(s.size)

for k in range(len(s)):
    fluxes_odd[k], forces_odd[k] = parallel_waveguide(s[k], True)
    fluxes_even[k], forces_even[k] = parallel_waveguide(s[k], False)

plt.figure(dpi=150)
plt.plot(s, -forces_odd / fluxes_odd, "rs", label="anti-symmetric")
plt.plot(s, -forces_even / fluxes_even, "bo", label="symmetric")
plt.grid(True)
plt.xlabel("waveguide separation s/a")
plt.ylabel("optical force (F/L)(ac/P)")
plt.legend(loc="upper right")
plt.show()
