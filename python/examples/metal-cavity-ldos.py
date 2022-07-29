import math

import matplotlib.pyplot as plt
import numpy as np

import meep as mp


def metal_cavity(w):
    resolution = 50
    sxy = 2
    dpml = 1
    sxy += 2 * dpml
    cell = mp.Vector3(sxy, sxy)

    pml_layers = [mp.PML(dpml)]
    a = 1
    t = 0.1
    geometry = [
        mp.Block(mp.Vector3(a + 2 * t, a + 2 * t, mp.inf), material=mp.metal),
        mp.Block(mp.Vector3(a, a, mp.inf), material=mp.air),
    ]

    geometry.append(
        mp.Block(
            center=mp.Vector3(a / 2), size=mp.Vector3(2 * t, w, mp.inf), material=mp.air
        )
    )

    fcen = math.sqrt(0.5) / a
    df = 0.2
    sources = [
        mp.Source(
            src=mp.GaussianSource(fcen, fwidth=df), component=mp.Ez, center=mp.Vector3()
        )
    ]

    symmetries = [mp.Mirror(mp.Y)]

    sim = mp.Simulation(
        cell_size=cell,
        geometry=geometry,
        boundary_layers=pml_layers,
        sources=sources,
        symmetries=symmetries,
        resolution=resolution,
    )

    h = mp.Harminv(mp.Ez, mp.Vector3(), fcen, df)
    sim.run(mp.after_sources(h), until_after_sources=500)

    m = h.modes[0]
    f = m.freq
    Q = m.Q
    Vmode = 0.25 * a * a
    ldos_1 = Q / Vmode / (2 * math.pi * f * math.pi * 0.5)

    sim.reset_meep()

    T = 2 * Q * (1 / f)
    sim.run(mp.dft_ldos(f, 0, 1), until_after_sources=T)
    ldos_2 = sim.ldos_data[0]

    return ldos_1, ldos_2


ws = np.arange(0.2, 0.5, 0.1)
ldos_1 = np.zeros(len(ws))
ldos_2 = np.zeros(len(ws))

for j in range(len(ws)):
    ldos_1[j], ldos_2[j] = metal_cavity(ws[j])
    print(f"ldos:, {ldos_1[j]}, {ldos_2[2]}")

plt.figure(dpi=150)
plt.semilogy(1 / ws, ldos_1, "bo-", label="2Q/(πωV)")
plt.semilogy(1 / ws, ldos_2, "rs-", label="LDOS")
plt.xlabel("a/w")
plt.ylabel("2Q/(πωW) or LDOS")
plt.show()
