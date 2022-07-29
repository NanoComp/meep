# Material dispersion example, from the Meep tutorial.  Here, we simply
# simulate homogenous space filled with a dispersive material, and compute
# its modes as a function of wavevector k.  Since omega/c = k/n, we can
# extract the dielectric function epsilon(omega) = (ck/omega)^2.
import meep as mp

cell = mp.Vector3()
resolution = 20

# We'll use a dispersive material with two polarization terms, just for
# illustration.  The first one is a strong resonance at omega=1.1,
# which leads to a polaritonic gap in the dispersion relation.  The second
# one is a weak resonance at omega=0.5, whose main effect is to add a
# small absorption loss around that frequency.

susceptibilities = [
    mp.LorentzianSusceptibility(frequency=1.1, gamma=1e-5, sigma=0.5),
    mp.LorentzianSusceptibility(frequency=0.5, gamma=0.1, sigma=2e-5),
]

default_material = mp.Medium(epsilon=2.25, E_susceptibilities=susceptibilities)

fcen = 1.0
df = 2.0

sources = [
    mp.Source(mp.GaussianSource(fcen, fwidth=df), component=mp.Ez, center=mp.Vector3())
]

kmin = 0.3
kmax = 2.2
k_interp = 99

kpts = mp.interpolate(k_interp, [mp.Vector3(kmin), mp.Vector3(kmax)])

sim = mp.Simulation(
    cell_size=cell,
    geometry=[],
    sources=sources,
    default_material=default_material,
    resolution=resolution,
)

all_freqs = sim.run_k_points(200, kpts)  # a list of lists of frequencies

for fs, kx in zip(all_freqs, [v.x for v in kpts]):
    for f in fs:
        print(f"eps:, {f.real:.6g}, {f.imag:.6g}, {(kx / f) ** 2:.6g}")
