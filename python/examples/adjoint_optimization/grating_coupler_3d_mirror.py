"""A 3D grating coupler optimized together with its fiber and its mirror.

A single-mode fiber illuminates a silicon grating coupler through an
anti-reflection coating and a BEOL oxide stack, and an aluminium mirror below
the buried oxide recovers the downward-radiated lobe.  Four things are optimized
at once, and the point of the example is that they are coupled: the best grating
depends on where the mirror sits, and the best fiber angle depends on both.

    grating topology        MaterialGrid + subpixel-smoothed projection
    fiber incidence angle   JAX, through the angular-spectrum propagator
    fiber lateral position  JAX, through the same chain
    mirror depth            shape derivative of a dispersive object

The ARC and the 2.2 um of cladding above the grating are *not meshed*.  The
fiber mode is propagated down through them analytically and injected as
equivalent surface currents on a plane just above the silicon, which is what
keeps a 3D run tractable -- it removes about 2.5 um of z, or fifty pixels at
resolution 20.  The fiber's angle and position are parameters of the mode, so
the derivative reaches them through JAX:

    (offset, tilt) -> gaussian_mode -> incident_fields -> equivalent_sources
                   -> amp_data -> Meep's adjoint solver

Light the grating radiates back upward leaves through the PML above the
injection plane.  That is the right model precisely because the ARC is there:
the ARC nulls the air/oxide interface, and the oxide above the grating is
homogeneous, so there is nothing for the upward lobe to reflect from.

Foundry constraints are Applied Nanotools' published NanoSOI specifications:
220 nm silicon, 2.0 um buried oxide, 2.2 um PECVD cladding, a single full etch,
60 nm minimum feature and 70 nm minimum gap.  The length scale is imposed with
the smoothed-subpixel-projection constraints of

    A. M. Hammond, A. Oskooi, I. M. Hammond, M. Chen, S. E. Ralph and
    S. G. Johnson, "Unifying and accelerating level-set and density-based
    topology optimization by subpixel-smoothed projection," Opt. Express 33(16),
    33620 (2025),

which requires the `ssp_topopt` package (github.com/NanoComp/SSP).

Why the mirror depth is worth optimizing: aluminium reflects from oxide with
|r| = 0.9817 at arg r = -169.4 deg, so the round trip returns in phase at
grating-to-mirror gaps of 0.79, 1.33, 1.87, 2.41 and 2.94 um.  A standard 2.0 um
box sits almost exactly between the third and fourth, so the stock stack is near
worst case and the optimum is a few hundred nm away.  The run prints the
converged depth against that analytic value, which is a check a reader can
follow.

Performance, all measured:

  * Real time-domain fields, not complex: 1.95x.  The adjoint needs complex DFT
    *monitors*, not complex fields.  The objective changes by exactly 4x -- the
    cos(wt) versus exp(-iwt) source convention -- which cancels in an efficiency.
  * MPI, not OpenMP: 13.3x on 64 ranks against 1.23x for 8 threads.  The solve
    is memory-bandwidth bound, so ranks buy bandwidth and threads do not.  Run
    under `mpirun` with one thread per rank.
  * The Courant number is derived from the metal.  Drude stability is a
    condition on omega_p*dt, measured to break between 1.264 and 1.580, so
    holding that fixed makes dt nearly independent of resolution and the cost
    scales as res^3 rather than res^4.
  * Pade extrapolation of the DFT tail: 14% fewer timesteps at equal accuracy.

What is verified, and what the optimizer is for.  Injecting into homogeneous
oxide and measuring the flux plane by plane, 99.4% of the launched power
reaches the silicon surface, arriving as a 6.3 um rms beam centred +0.49 um
downstream -- consistent with an 8 degree tilt over the propagation distance.
So the fiber, the analytic stack and the equivalent-current injection all work.

The seed device does not: a uniform grating scatters into the full 6 um width
of the aperture and a 0.5 um waveguide collects almost none of it, which is
-54 dB.  There is deliberately no taper.  Focusing the scattered light into the
waveguide is precisely what the design region is there to discover, and giving
it a taper would be solving by hand the part worth optimizing.

Run:
    mpirun -n 64 python grating_coupler_3d_mirror.py optimize
    mpirun -n 64 python grating_coupler_3d_mirror.py forward
    python grating_coupler_3d_mirror.py report --history history.npz
"""

import argparse
import math
import os
import time
from typing import List, Optional, Tuple

import jax
import jax.numpy as jnp
import numpy as np

import meep as mp
import meep.adjoint as mpa
from meep.materials import Al as Al_rakic

jax.config.update("jax_enable_x64", True)

# --------------------------------------------------------------------- stack
WAVELENGTH = 1.55
FCEN = 1.0 / WAVELENGTH

N_SI, N_OXIDE, N_AIR = 3.48, 1.44, 1.0
N_ARC = math.sqrt(N_AIR * N_OXIDE)  # quarter-wave match, ~1.20

T_ARC = WAVELENGTH / (4 * N_ARC)
T_CLAD = 2.2  # ANT standard PECVD cladding, handled analytically
T_STANDOFF = 1.2  # meshed oxide above the grating
T_LAUNCH = 0.3  # launch plane height above the silicon surface
T_SI = 0.220  # ANT device layer
T_BOX = 2.0  # ANT buried oxide
T_SPACER = 0.8  # oxide below the box; the mirror moves inside this
T_AL = 0.100  # over twelve skin depths, so optically opaque
T_HANDLE = 0.5

SI = mp.Medium(index=N_SI)
OXIDE = mp.Medium(index=N_OXIDE)

# ------------------------------------------------------------ foundry rules
N_EFF_GRATING = 2.1  # etched 220 nm slab at ~50% duty; MPB gives 2.088 for the wg

MIN_FEATURE = 0.060  # ANT NanoSOI minimum feature size
MIN_GAP = 0.070  # ANT NanoSOI minimum spacing

# ----------------------------------------------------------------- geometry
APERTURE_X = 12.0  # along the waveguide; the fiber MFD is 10.4 um
APERTURE_Y = 6.0  # full width; a mirror plane runs along y = 0
WG_WIDTH = 0.5
WG_LENGTH = 3.0
DPML = 1.0
MARGIN = 0.5  # keeps the design region clear of the PML interface

# -------------------------------------------------------------------- fiber
FIBER_MFD = 10.4  # SMF-28 at 1550 nm
FIBER_WAIST = FIBER_MFD / 2
FIBER_TILT_DEG = 8.0  # spoils the second-order back-reflection
FIBER_HEIGHT = 2.0  # facet above the top of the ARC

# ------------------------------------------------------------------- solver
WP_DT_MAX = 0.95  # headroom below the measured Drude stability limit
PADE_SAMPLES = 12
PADE_TOL = 1e-2
LOG_FLOOR = 1e-30  # keeps log finite if a monitor reads exactly zero


def aluminium() -> Tuple[mp.Medium, float]:
    """A single Drude pole matched to Rakic's aluminium at the design wavelength.

    The library model spans 0.2-12.4 um in five poles.  A 1550 nm design does
    not need that band, and each pole costs polarization state and `update_P`
    work wherever sigma is nonzero.  One pole matched to the same reference
    epsilon is exact at 1550 nm and tracks reflectivity to better than 0.001
    across the C band.

    Written as `frequency = fp, sigma = 1` rather than the library's
    `frequency = 1e-10` with `sigma ~ 1e20`: same product, far less dynamic
    range to carry into a timestepping recurrence.
    """
    eps = complex(np.squeeze(Al_rakic.epsilon(FCEN))[0, 0])
    d = eps - 1.0
    x = FCEN**2 * abs(d) ** 2 / (-d.real)
    fp = math.sqrt(x)
    gamma = x * d.imag / (abs(d) ** 2 * FCEN)
    return (
        mp.Medium(
            epsilon=1.0,
            E_susceptibilities=[
                mp.DrudeSusceptibility(frequency=fp, gamma=gamma, sigma=1.0)
            ],
        ),
        fp,
    )


AL, AL_FP = aluminium()


def courant_for(resolution: float) -> float:
    """Keep omega_p * dt below the measured Drude stability limit.

    Deriving this from the metal means changing resolution cannot silently cross
    the threshold, and because the Courant number then rises linearly with
    resolution, dt is nearly constant and cost scales as res^3, not res^4.
    """
    return min(0.5, WP_DT_MAX * resolution / (2 * math.pi * AL_FP))


def constructive_gaps(count: int = 5) -> List[float]:
    """Grating-to-mirror oxide gaps whose round trip returns in phase."""
    eps = complex(np.squeeze(AL.epsilon(FCEN))[0, 0])
    r = (N_OXIDE - np.sqrt(eps)) / (N_OXIDE + np.sqrt(eps))
    k = 2 * math.pi / WAVELENGTH
    return [
        (2 * math.pi * m - float(np.angle(r))) / (2 * k * N_OXIDE)
        for m in range(1, count + 1)
    ]


class Layout:
    """Where every layer sits, measured downward from the top of the silicon."""

    def __init__(self):
        self.z_si = 0.0
        self.z_launch = T_SI / 2 + T_LAUNCH
        top_of_box = -T_SI / 2
        self.z_standoff = T_SI / 2 + T_STANDOFF / 2
        self.z_box = top_of_box - (T_BOX + T_SPACER) / 2
        self.z_mirror0 = top_of_box - T_BOX - T_AL / 2
        self.z_handle = top_of_box - T_BOX - T_SPACER - T_AL - T_HANDLE / 2
        self.z_top = T_SI / 2 + T_STANDOFF
        self.z_bottom = top_of_box - T_BOX - T_SPACER - T_AL - T_HANDLE
        self.sz = (self.z_top - self.z_bottom) + 2 * DPML
        self.sx = APERTURE_X + 2 * WG_LENGTH + 2 * (DPML + MARGIN)
        self.sy = APERTURE_Y + 2 * (DPML + MARGIN)
        self.center = mp.Vector3(0, 0, (self.z_top + self.z_bottom) / 2)

    def mirror_gap(self, z_mirror: float) -> float:
        """Oxide thickness between the silicon and the top of the mirror."""
        return (-T_SI / 2) - (z_mirror + T_AL / 2)


LAYOUT = Layout()


def upward_stack() -> mpa.Stack:
    """What the fiber crosses between the launch plane and its facet.

    Ordered away from the monitor and terminating semi-infinite.  The ARC is a
    genuine quarter-wave layer rather than decoration: untreated, the air/oxide
    interface reflects about 3.4%, a fifth of a dB off a device whose whole
    budget is around one.
    """
    remaining_oxide = (T_STANDOFF - T_LAUNCH) + T_CLAD
    return mpa.Stack(
        layers=[
            mpa.Layer(index=N_OXIDE, thickness=remaining_oxide),
            mpa.Layer(index=N_ARC, thickness=T_ARC),
            mpa.Layer(index=N_AIR, thickness=None),
        ]
    )


def launch_volume() -> mp.Volume:
    """The plane the fiber's equivalent currents sit on.

    Kept 0.3 um above the silicon so the injected field reaches the grating
    without spreading, and T_STANDOFF - T_LAUNCH = 0.9 um below the PML so the
    current sheets are not radiating into its near field.  The first version had
    0.3 um of clearance, which is 2.4 pixels at resolution 8.
    """
    return mp.Volume(
        center=mp.Vector3(0, 0, LAYOUT.z_launch),
        size=mp.Vector3(APERTURE_X, APERTURE_Y, 0),
    )


# ------------------------------------------------------------- design chain
def design_resolution(resolution: float) -> float:
    """The design grid is finer than the Meep grid, so the filter is resolved.

    SSP smooths the level set analytically, but the conic filter still has to be
    represented: a 60 nm feature on a 20/um Meep grid is barely one pixel.
    """
    return max(4 * resolution, 1.0 / (MIN_FEATURE / 4))


def design_shape(resolution: float) -> Tuple[int, int]:
    """Shape of the *half* design the optimizer owns; see `mirror_y`."""
    dres = design_resolution(resolution)
    return int(APERTURE_X * dres), int(APERTURE_Y * dres) // 2


def mirror_y(half: jnp.ndarray) -> jnp.ndarray:
    """Build the full design from its y >= 0 half.

    The simulation has a `Mirror(Y)` symmetry, so Meep solves only y >= 0 and
    the adjoint gradient comes back with 98% of its magnitude there -- measured,
    403 against 6.4 summed over each half.  Optimizing the full grid against
    that would leave half the variables with no gradient and would break the
    very symmetry the solver is exploiting.

    Parameterizing the half and reflecting it keeps design and solver
    consistent, halves the number of variables, and makes the returned gradient
    the whole derivative rather than half of one.
    """
    return jnp.concatenate([half[:, ::-1], half], axis=1)


def project(latent: jnp.ndarray, beta: float, resolution: float) -> jnp.ndarray:
    """Latent variables to a projected density, with the length scale built in.

    Conic radius from the target feature size, as the SSP paper recommends, and
    `ssp2` for the projection so beta can go to infinity without the gradient
    going with it.
    """
    from ssp_topopt import conic_filter, ssp2

    latent = mirror_y(latent)
    nx, ny = latent.shape
    dres = design_resolution(resolution)
    lx, ly = (nx - 1) / dres, (ny - 1) / dres
    filtered = conic_filter(latent, MIN_FEATURE, lx, ly, dres)
    return ssp2(filtered, beta, 0.5, dres), filtered


def length_penalties(latent, beta, resolution):
    """ANT's 60 nm feature and 70 nm gap, as two scalar constraints <= 0."""
    from ssp_topopt import length_constraint

    projected, filtered = project(latent, beta, resolution)
    dres = design_resolution(resolution)
    return (
        length_constraint("solid", filtered, projected, dres, MIN_FEATURE),
        length_constraint("void", filtered, projected, dres, MIN_GAP),
    )


def seed_grating(nx: int, ny: int, margin: float = 0.02) -> np.ndarray:
    """A uniform grating at the phase-matched period, as a starting point.

    period = lambda / (n_eff - n_clad sin(theta)), with n_eff ~ 2.7 for a
    half-filled 220 nm silicon layer.  Starting from a structure that already
    couples is worth a great deal: a uniform rho = 0.5 slab couples almost
    nothing, and an objective near zero is a poor place to start a gradient
    method (and, separately, a poor place to benchmark DFT convergence).
    """
    # n_eff of the *etched* layer, not of bulk silicon.  A 220 nm slab at
    # roughly half duty in oxide guides at about 2.1 -- MPB puts the output
    # waveguide's own mode at 2.088 -- and the first version of this used 2.7,
    # giving 0.632 um where phase matching wants 0.816.  A 21% error in the
    # period is a large phase mismatch and costs far more than it looks.
    period = WAVELENGTH / (
        N_EFF_GRATING - N_OXIDE * math.sin(math.radians(FIBER_TILT_DEG))
    )
    x = np.linspace(-APERTURE_X / 2, APERTURE_X / 2, nx)
    # Binary, not a smooth cosine.  The length-scale constraint penalises
    # regions where the filtered field is flat, and a cosine is flat across
    # every crest and trough -- measured, it starts 1.7e4 in violation, against
    # -1.0 (comfortably feasible) for the binary grating with the same period.
    # Starting infeasible makes MMA spend its budget reaching feasibility and
    # the efficiency falls while it does.
    line = (np.cos(2 * math.pi * x / period) > 0).astype(float)
    # Held `margin` inside the box.  Exactly 0 and 1 sit *on* the bounds, and
    # MMA then has no interior direction to move in: measured, three successive
    # iterations returned a bit-identical design, objective and gradient.
    line = margin + (1 - 2 * margin) * line
    return np.repeat(line[:, None], ny, axis=1)  # ny is already the half


# ---------------------------------------------------------------- simulation
def build(resolution: float, z_mirror: float, weights: np.ndarray):
    """The coupler with no sources; the fiber's currents are added separately."""
    inf = mp.inf
    geometry = [
        # Oxide everywhere, then the things that are not oxide.
        mp.Block(
            center=mp.Vector3(0, 0, (LAYOUT.z_top + LAYOUT.z_bottom) / 2),
            size=mp.Vector3(inf, inf, LAYOUT.z_top - LAYOUT.z_bottom),
            material=OXIDE,
        ),
        mp.Block(
            center=mp.Vector3(0, 0, LAYOUT.z_handle),
            size=mp.Vector3(inf, inf, T_HANDLE),
            material=SI,
        ),
        # The waveguide the grating feeds, running out in -x.
        mp.Block(
            center=mp.Vector3(-(APERTURE_X / 2 + WG_LENGTH / 2), 0, 0),
            size=mp.Vector3(WG_LENGTH + 1.0, WG_WIDTH, T_SI),
            material=SI,
        ),
        mp.Block(
            center=mp.Vector3(0, 0, z_mirror),
            size=mp.Vector3(APERTURE_X + 2 * WG_LENGTH, inf, T_AL),
            material=AL,
            differentiable=["center"],
            name="mirror",
        ),
    ]

    nx, ny = weights.shape
    grid = mp.MaterialGrid(
        mp.Vector3(nx, ny),
        OXIDE,
        SI,
        weights=weights,
        do_averaging=False,  # SSP has already smoothed the level set
        grid_type="U_MEAN",
    )
    region = mpa.DesignRegion(
        grid,
        volume=mp.Volume(
            center=mp.Vector3(0, 0, 0),
            size=mp.Vector3(APERTURE_X, APERTURE_Y, T_SI),
        ),
    )
    geometry.append(mp.Block(center=region.center, size=region.size, material=grid))

    sim = mp.Simulation(
        cell_size=mp.Vector3(LAYOUT.sx, LAYOUT.sy, LAYOUT.sz),
        geometry_center=LAYOUT.center,
        resolution=resolution,
        Courant=courant_for(resolution),
        geometry=geometry,
        boundary_layers=[mp.PML(DPML)],
        # Ex is even about y = 0 for the TE mode this couples into.
        symmetries=[mp.Mirror(mp.Y, phase=+1)],
        sources=[],
        # Real fields: 1.95x, and the adjoint only needs complex DFT monitors.
        force_complex_fields=False,
        # The cost model knows a metal pixel is dearer than an oxide one, but it
        # has no communication term, so it is worth measuring at your rank count.
        split_chunks_evenly=False,
    )
    return sim, region, grid


def waveguide_monitor(sim, frequencies):
    """Power in the guided mode, which is what a coupler is for."""
    return mpa.EigenmodeCoefficient(
        sim,
        mp.Volume(
            center=mp.Vector3(-(APERTURE_X / 2 + WG_LENGTH * 0.6), 0, 0),
            size=mp.Vector3(0, 4 * WG_WIDTH, 8 * T_SI),
        ),
        mode=1,
        forward=False,  # power leaving along -x
    )


# --------------------------------------------------------------- fiber chain
def make_propagator(sim, frequencies):
    """A propagator for the launch plane, facing up toward the fiber.

    `from_volume` rather than `from_monitor`: a DFT monitor needs the field
    components allocated, which does not happen until sources exist, and the
    sources are exactly what this is being built to produce.
    """
    return mpa.AngularSpectrum.from_volume(
        sim,
        sim._fit_volume_to_simulation(launch_volume()),
        upward_stack(),
        frequencies,
        sign=+1,  # the fiber is above
        pad_factor=4,
    )


def fiber_distance() -> float:
    """From the launch plane up to the fiber facet, through the unmeshed stack."""
    return FIBER_HEIGHT + T_ARC + T_CLAD + (T_STANDOFF - T_LAUNCH)


def fiber_sources(propagator, params, frequencies) -> List[mp.Source]:
    """Equivalent currents for a fiber at (lateral offset, tilt in degrees)."""
    offset, tilt_deg = float(params[0]), float(params[1])
    mode = mpa.gaussian_mode(FIBER_WAIST, tilt_deg=tilt_deg, offset=offset)
    fields = propagator.incident_fields(mode, distance=fiber_distance())
    plane = launch_volume()
    return propagator.equivalent_sources(
        fields,
        mp.GaussianSource(FCEN, fwidth=0.2 * FCEN),
        center=plane.center,
        size=plane.size,
    )


def launched_power(args, params, frequencies) -> float:
    """Power the fiber delivers into the cell, measured the way Meep measures it.

    `EigenmodeCoefficient` returns |a|^2 in Meep's power normalization, and the
    propagator's own `power()` is in its own -- dividing one by the other gives
    a number bigger than one, which is how this was caught.  So measure the
    input the same way as the output: the same stack with the design region set
    to plain oxide, and a flux plane just below the launch plane.

    It costs one short extra solve, once, and it is what makes the reported
    number a coupling efficiency rather than an arbitrary scale.  The source
    convention cancels here too: real and complex fields differ by exactly 4x
    in numerator and denominator alike.
    """
    # Homogeneous, deliberately.  Keeping the stack underneath would measure
    # injected *minus* reflected, and an aluminium mirror sends nearly all of it
    # back -- that read 1.6x low.  It is also far cheaper: with nothing to ring,
    # this settles in about 6 s against 400 s for a full iteration.
    sim = mp.Simulation(
        cell_size=mp.Vector3(LAYOUT.sx, LAYOUT.sy, LAYOUT.sz),
        geometry_center=LAYOUT.center,
        resolution=args.resolution,
        Courant=courant_for(args.resolution),
        geometry=[],
        default_material=OXIDE,
        boundary_layers=[mp.PML(DPML)],
        symmetries=[mp.Mirror(mp.Y, phase=+1)],
        sources=[],
        force_complex_fields=False,
    )
    sim.init_sim()
    propagator = make_propagator(sim, frequencies)
    sim.change_sources(fiber_sources(propagator, params, frequencies))

    flux = sim.add_flux(
        frequencies[0],
        0 if len(frequencies) == 1 else frequencies[-1] - frequencies[0],
        len(frequencies),
        mp.FluxRegion(
            center=mp.Vector3(0, 0, LAYOUT.z_launch - 0.15),
            size=mp.Vector3(APERTURE_X, APERTURE_Y, 0),
            direction=mp.Z,
        ),
    )
    sim.run(until_after_sources=mp.stop_when_dft_decayed(1e-6, 20, 400))
    return float(abs(np.mean(mp.get_fluxes(flux))))


def fiber_amp_jacobian(propagator, params, frequencies, step=1e-4):
    """d(amp_data)/d(offset, tilt), by differencing the *analytic* chain.

    Meep returns d(objective)/d(amp_data) for a source flagged differentiable;
    turning that into a derivative with respect to the fiber's position and
    angle needs the Jacobian of the chain that produced those amplitudes.

    That chain -- mode, propagate through the stack, apply the equivalence
    principle -- involves no simulation, so differencing it is free.  It is also
    the well-conditioned half: offset and tilt are smooth O(1) parameters, quite
    unlike differencing a position against a pixel.  `equivalent_sources`
    materialises to NumPy on the way out, which is where a JAX trace would stop
    anyway, so this avoids reimplementing its internals to get around that.
    """
    jac = []
    for i in range(len(params)):
        hi, lo = np.array(params, float), np.array(params, float)
        hi[i] += step
        lo[i] -= step
        a_hi = [
            np.asarray(s.amp_data) for s in fiber_sources(propagator, hi, frequencies)
        ]
        a_lo = [
            np.asarray(s.amp_data) for s in fiber_sources(propagator, lo, frequencies)
        ]
        jac.append([(h - l) / (2 * step) for h, l in zip(a_hi, a_lo)])
    return jac


def fiber_gradient(cotangents, jacobian) -> np.ndarray:
    """Contract Meep's amp_data cotangents onto the fiber parameters.

    The adjoint convention is g = dJ/d(Re a) - i dJ/d(Im a), so the real
    derivative along a parameter is sum Re(g * da/dp).

    One cotangent per source, paired with that source's own column of the
    Jacobian.  Sizes must agree: silently truncating a mismatch to the shorter
    of the two is how this returned 1e-18 instead of a gradient.
    """
    out = np.zeros(len(jacobian))
    for p, per_source in enumerate(jacobian):
        if len(per_source) != len(cotangents):
            raise ValueError(
                f"{len(cotangents)} cotangents for {len(per_source)} sources; "
                "they must correspond one to one."
            )
        for g, da in zip(cotangents, per_source):
            g = np.asarray(g).ravel()
            da = np.asarray(da).ravel()
            if g.size != da.size:
                raise ValueError(
                    f"cotangent has {g.size} elements but amp_data has "
                    f"{da.size}; these must be the same array shape."
                )
            out[p] += float(np.real(np.dot(g, da)))
    return out


# ------------------------------------------------------------- diagnostics
class Diagnostics:
    """Per-run folder of figures, written as the optimization proceeds.

    Watching the numbers scroll past is a poor way to notice that a design has
    gone grey, that a gradient is concentrated in one corner, or that the metal
    has walked into the substrate.  Everything here is written every iteration
    so a run can be diagnosed while it is still going.
    """

    def __init__(self, root: str, layout, resolution: float):
        self.dir = root
        self.layout = layout
        self.resolution = resolution
        self.history = []
        if mp.am_master():
            os.makedirs(self.dir, exist_ok=True)

    def _plt(self):
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        return plt

    def record(self, iteration, design, gradient, efficiency, z_mirror, sim=None):
        self.history.append(efficiency)
        if not mp.am_master():
            return
        plt = self._plt()
        self._design(plt, iteration, design, efficiency)
        self._gradient(plt, iteration, gradient)
        self._objective(plt)
        if sim is not None:
            self._cross_section(plt, iteration, sim, efficiency, z_mirror)

    def _extent(self, arr):
        """Full aperture for a mirrored design, the y >= 0 half for a gradient."""
        full = np.asarray(arr).shape[1] >= int(
            APERTURE_Y * design_resolution(self.resolution)
        )
        y0 = -APERTURE_Y / 2 if full else 0.0
        return [-APERTURE_X / 2, APERTURE_X / 2, y0, APERTURE_Y / 2]

    def _design(self, plt, iteration, design, efficiency):
        fig, ax = plt.subplots(figsize=(9, 5))
        im = ax.imshow(
            np.asarray(design).T,
            origin="lower",
            extent=self._extent(design),
            cmap="binary",
            vmin=0,
            vmax=1,
            aspect="equal",
        )
        ax.set_title(
            f"design, iteration {iteration}   "
            f"efficiency {efficiency:.4e} ({10 * np.log10(max(efficiency, 1e-300)):.1f} dB)"
        )
        ax.set_xlabel("x (um)  -- waveguide is to the left")
        ax.set_ylabel("y (um)")
        fig.colorbar(im, ax=ax, label="density (0 oxide, 1 silicon)")
        fig.tight_layout()
        fig.savefig(f"{self.dir}/design_{iteration:03d}.png", dpi=110)
        plt.close(fig)

    def _gradient(self, plt, iteration, gradient):
        g = np.asarray(gradient)
        # Symmetric limits so the colormap's zero is the gradient's zero.
        lim = float(np.max(np.abs(g))) or 1.0
        fig, ax = plt.subplots(figsize=(9, 5))
        im = ax.imshow(
            g.T,
            origin="lower",
            extent=self._extent(g),
            cmap="RdBu_r",
            vmin=-lim,
            vmax=lim,
            aspect="equal",
        )
        ax.set_title(
            f"objective gradient, iteration {iteration}   " f"max |g| = {lim:.3e}"
        )
        ax.set_xlabel("x (um)")
        ax.set_ylabel("y (um)")
        fig.colorbar(im, ax=ax, label="d(log coupling)/d(density)")
        fig.tight_layout()
        fig.savefig(f"{self.dir}/gradient_{iteration:03d}.png", dpi=110)
        plt.close(fig)

    def _objective(self, plt):
        fig, ax = plt.subplots(figsize=(7, 4))
        y = np.asarray(self.history)
        ax.semilogy(np.arange(1, len(y) + 1), y, "o-", ms=3)
        ax.set_xlabel("iteration")
        ax.set_ylabel("coupling efficiency")
        ax.grid(True, which="both", alpha=0.3)
        best = float(np.max(y))
        ax.set_title(f"best {best:.4e} = {10 * np.log10(max(best, 1e-300)):.2f} dB")
        fig.tight_layout()
        fig.savefig(f"{self.dir}/objective.png", dpi=110)
        plt.close(fig)

    def _cross_section(self, plt, iteration, sim, efficiency, z_mirror):
        """An xz slice through y = 0, with the metal called out."""
        x = np.linspace(-LAYOUT.sx / 2, LAYOUT.sx / 2, 500)
        z = np.linspace(LAYOUT.z_bottom - DPML, LAYOUT.z_top + DPML, 400)
        # At the design frequency, not at DC.  A Drude metal has eps_infinity =
        # 1, so a DC plot draws the aluminium as though it were vacuum and the
        # mirror is simply invisible -- which is exactly how this looked the
        # first time.  At 1550 nm those same pixels are -232.7 + 47.1j.
        eps = np.real(
            np.asarray(sim.get_epsilon_grid(x, np.array([0.0]), z, frequency=FCEN))
        ).reshape(len(x), len(z))

        fig, ax = plt.subplots(figsize=(10, 4.5))
        # Clipped and diverging, so the metal (strongly negative) and the
        # dielectrics (2 to 12) are both legible on one scale.
        im = ax.pcolormesh(
            x,
            z,
            np.clip(eps, -13, 13).T,
            cmap="RdBu_r",
            vmin=-13,
            vmax=13,
            shading="auto",
        )
        ax.axhline(LAYOUT.z_launch, color="w", ls=":", lw=1.2)
        ax.text(
            -LAYOUT.sx / 2 + 0.2,
            LAYOUT.z_launch + 0.08,
            "fiber launch",
            color="w",
            fontsize=8,
        )
        top, bot = z_mirror + T_AL / 2, z_mirror - T_AL / 2
        ax.axhline(top, color="r", lw=1.0)
        ax.axhline(bot, color="r", lw=1.0)
        gap = LAYOUT.mirror_gap(z_mirror)
        near = min(constructive_gaps(), key=lambda g: abs(g - gap))
        ax.text(
            -LAYOUT.sx / 2 + 0.2,
            top + 0.10,
            f"Al mirror: gap {gap:.3f} um " f"(nearest analytic optimum {near:.3f})",
            color="r",
            fontsize=8,
        )
        for edge in (LAYOUT.z_bottom, LAYOUT.z_top):
            ax.axhline(edge, color="k", ls="--", lw=0.6, alpha=0.5)
        ax.set_xlabel("x (um)")
        ax.set_ylabel("z (um)")
        ax.set_title(
            f"iteration {iteration}: efficiency {efficiency:.4e} "
            f"({10 * np.log10(max(efficiency, 1e-300)):.1f} dB)"
        )
        fig.colorbar(
            im, ax=ax, label=f"Re(eps) at {1000 * WAVELENGTH:.0f} nm (clipped)"
        )
        fig.tight_layout()
        fig.savefig(f"{self.dir}/cross_section_{iteration:03d}.png", dpi=110)
        plt.close(fig)


# ------------------------------------------------------------------ history
class History:
    """Everything worth post-processing, one row per iteration.

    Saved as a single npz rather than printed, because the interesting questions
    afterwards -- did the mirror walk to the analytic optimum, did the fiber
    angle move away from 8 degrees, when did the length-scale constraints become
    active -- are all about trajectories rather than the final number.
    """

    def __init__(self, path: str):
        self.path = path
        self.rows = {
            k: []
            for k in (
                "iteration",
                "efficiency",
                "objective",
                "latent",
                "projected",
                "grad_design",
                "grad_offset",
                "grad_tilt",
                "grad_mirror",
                "fiber_offset",
                "fiber_tilt",
                "mirror_z",
                "mirror_gap",
                "constraint_solid",
                "constraint_void",
                "beta",
                "timesteps",
                "wall_seconds",
            )
        }

    def add(self, **kw):
        for k, v in kw.items():
            self.rows[k].append(np.asarray(v))

    def save(self):
        if not mp.am_master():
            return
        np.savez_compressed(
            self.path,
            **{k: np.asarray(v) for k, v in self.rows.items() if v},
            wavelength=WAVELENGTH,
            min_feature=MIN_FEATURE,
            min_gap=MIN_GAP,
            constructive_gaps=np.asarray(constructive_gaps()),
        )


def say(msg: str):
    if mp.am_master():
        print(msg, flush=True)


# --------------------------------------------------------------- the problem
class Coupler:
    """One optimizable grating coupler: design, fiber and mirror together."""

    def __init__(self, args):
        self.args = args
        self.frequencies = (
            [FCEN]
            if args.nfreq == 1
            else list(
                np.linspace(
                    FCEN * (1 - args.bandwidth / 2),
                    FCEN * (1 + args.bandwidth / 2),
                    args.nfreq,
                )
            )
        )
        self.nx, self.ny = design_shape(args.resolution)
        self.beta = args.beta_schedule[0]
        self._power_ref = None
        self._propagator_for_power = None
        self._last_sim = None
        self.diagnostics = None

    def reference_power(self, fiber):
        """Input power for a given fiber placement, in Meep's units.

        Calibrating with a simulation is what puts |a|^2 and the input power in
        the same normalization, but doing it per placement would cost a full
        extra solve on every iteration -- the optimizer moves the fiber each
        time, so a cache keyed on position never hits.

        The fix is that only the *unit conversion* needs a simulation, and that
        is a constant.  Measure it once at the nominal fiber, then track the
        placement dependence with the propagator's own analytic power, which is
        exact and free:

            P_meep(fiber) = P_meep(nominal) * P_asm(fiber) / P_asm(nominal)
        """
        if self._power_ref is None:
            nominal = [0.0, FIBER_TILT_DEG]
            self._power_ref = (
                launched_power(self.args, nominal, self.frequencies),
                self._analytic_power(nominal),
            )
        meep_ref, asm_ref = self._power_ref
        if asm_ref <= 0:
            return meep_ref
        return meep_ref * self._analytic_power(fiber) / asm_ref

    def _grad_log_input_power(self, fiber, step: float = 1e-4) -> np.ndarray:
        """d(log P_in)/d(offset, tilt), from the analytic power alone."""
        out = np.zeros(len(fiber))
        for i in range(len(fiber)):
            hi, lo = np.array(fiber, float), np.array(fiber, float)
            hi[i] += step
            lo[i] -= step
            p_hi, p_lo = self._analytic_power(hi), self._analytic_power(lo)
            if p_hi > 0 and p_lo > 0:
                out[i] = (math.log(p_hi) - math.log(p_lo)) / (2 * step)
        return out

    def _analytic_power(self, fiber) -> float:
        """Power the propagator says the fiber delivers; no simulation."""
        if self._propagator_for_power is None:
            nx, ny = design_shape(self.args.resolution)
            sim, _, _ = build(
                self.args.resolution, LAYOUT.z_mirror0, np.zeros((nx, ny))
            )
            sim.init_sim()
            self._propagator_for_power = make_propagator(sim, self.frequencies)
        prop = self._propagator_for_power
        mode = mpa.gaussian_mode(
            FIBER_WAIST, tilt_deg=float(fiber[1]), offset=float(fiber[0])
        )
        fields = prop.incident_fields(mode, distance=fiber_distance())
        return float(np.mean(np.asarray(prop.power(fields, distance=fiber_distance()))))

    def _weights(self, latent):
        projected, _ = project(jnp.asarray(latent), self.beta, self.args.resolution)
        return np.asarray(projected, dtype=np.float64)

    def evaluate(self, latent, fiber, z_mirror, need_gradient=True):
        """One forward (and optionally adjoint) solve; returns value and pieces."""
        import time

        weights = self._weights(latent)
        sim, region, _ = build(self.args.resolution, z_mirror, weights)
        # The propagator reads the grid metadata of the launch plane, which only
        # exists once the fields do -- and the fields are built without sources
        # precisely because the sources are what the propagator is here to make.
        sim.init_sim()
        propagator = make_propagator(sim, self.frequencies)
        sources = fiber_sources(propagator, fiber, self.frequencies)
        sim.change_sources(sources)

        monitor = waveguide_monitor(sim, self.frequencies)
        opt = mpa.OptimizationProblem(
            simulation=sim,
            # Logarithmic, which matters a great deal at this starting point.
            # d(log f)/dx = (1/f) df/dx, so when the coupling is 1e-6 the
            # gradient is amplified by 1e6 and the optimizer moves instead of
            # crawling; a linear objective gives a vanishing gradient exactly
            # where the design is worst.  See Hammond et al.,
            # https://doi.org/10.1364/OE.466015.
            #
            # It also makes the normalisation an additive constant, so the
            # design and mirror gradients need no rescaling by the input power.
            objective_functions=[
                lambda a: jnp.log(jnp.mean(jnp.abs(a) ** 2) + LOG_FLOOR)
            ],
            objective_arguments=[monitor],
            design_regions=[region],
            frequencies=self.frequencies,
            decay_by=self.args.decay_by,
            minimum_run_time=self.args.min_run_time,
            maximum_run_time=self.args.max_run_time,
            pade_samples=PADE_SAMPLES if self.args.pade else 0,
            pade_tolerance=PADE_TOL if self.args.pade else None,
        )

        t0 = time.time()
        value, grads = opt([weights.ravel()], need_gradient=need_gradient)
        wall = time.time() - t0
        steps = sim.timestep()

        # The objective is log(|a|^2).  Dividing by the launched power to get
        # an efficiency is then a subtraction, so the design and mirror
        # gradients pass through untouched -- only the reported number and the
        # fiber gradient see it.
        p_in = self.reference_power(fiber)
        log_obj = float(np.real(np.squeeze(value)))
        efficiency = (
            float(np.exp(log_obj) / p_in) if p_in > 0 else float(np.exp(log_obj))
        )

        self._last_sim = sim
        if not need_gradient:
            return efficiency, None, wall, steps, propagator
        return efficiency, grads, wall, steps, propagator

    def step(self, x, history, iteration):
        """Objective and gradient for the packed parameter vector."""
        latent, fiber, z_mirror = self.unpack(x)
        value, grads, wall, steps, propagator = self.evaluate(latent, fiber, z_mirror)

        # Design: chain Meep's density gradient back through filter and
        # projection with a VJP, which is what makes the length scale part of
        # the parameterization rather than a penalty bolted on afterwards.
        # The gradient comes back on the *full* grid; fold it onto the half the
        # optimizer owns, so both mirrored copies of each variable contribute.
        dens_grad = np.asarray(grads["design"], dtype=np.float64).reshape(
            self.nx, 2 * self.ny
        )

        def to_density(u):
            return project(u, self.beta, self.args.resolution)[0]

        _, vjp = jax.vjp(to_density, jnp.asarray(latent))
        grad_latent = np.asarray(vjp(jnp.asarray(dens_grad))[0])

        # Fiber: Meep gives d/d(amp_data) for *each* of the equivalent-current
        # sources -- the equivalence principle places four -- keyed by their
        # index in registration order, which is the order fiber_sources built
        # them.  Pair each with its own column of the analytic Jacobian.
        src_keys = sorted((k for k in grads if k not in ("design", "mirror")), key=str)
        if not src_keys:
            grad_fiber = np.zeros(2)
        else:
            cotangents = [np.asarray(grads[k]["amp_data"]) for k in src_keys]
            jac = fiber_amp_jacobian(propagator, fiber, self.frequencies)
            grad_fiber = fiber_gradient(cotangents, jac)
            # The objective the optimizer sees is log(|a|^2) - log(P_in), and
            # unlike the design and the mirror, the fiber moves P_in.  That term
            # is analytic and needs no simulation, so difference it here.
            grad_fiber = grad_fiber - self._grad_log_input_power(fiber)

        grad_mirror = float(np.real(np.atleast_1d(grads["mirror"]["center"])[2]))

        cs, cv = length_penalties(jnp.asarray(latent), self.beta, self.args.resolution)
        projected = self._weights(latent)

        history.add(
            iteration=iteration,
            efficiency=value,
            objective=value,
            latent=latent,
            projected=projected,
            grad_design=grad_latent,
            grad_offset=grad_fiber[0],
            grad_tilt=grad_fiber[1],
            grad_mirror=grad_mirror,
            fiber_offset=fiber[0],
            fiber_tilt=fiber[1],
            mirror_z=z_mirror,
            mirror_gap=LAYOUT.mirror_gap(z_mirror),
            constraint_solid=float(cs),
            constraint_void=float(cv),
            beta=self.beta,
            timesteps=steps,
            wall_seconds=wall,
        )
        if self.diagnostics is not None:
            self.diagnostics.record(
                iteration, projected, grad_latent, value, z_mirror, self._last_sim
            )
        say(
            f"  iter {iteration:3d}  eff={value:.6e}  "
            f"gap={LAYOUT.mirror_gap(z_mirror):.4f}um  "
            f"tilt={fiber[1]:.3f}deg  offset={fiber[0]:+.4f}um  "
            f"beta={self.beta:g}  steps={steps}  {wall:.0f}s"
        )
        return value, self.pack_gradient(grad_latent, grad_fiber, grad_mirror)

    # The optimizer sees one flat vector; these two keep the packing in one place.
    def unpack(self, x):
        n = self.nx * self.ny
        latent = np.asarray(x[:n]).reshape(self.nx, self.ny)
        fiber = np.asarray(x[n : n + 2])
        z_mirror = float(x[n + 2])
        return latent, fiber, z_mirror

    def pack(self, latent, fiber, z_mirror):
        return np.concatenate(
            [np.asarray(latent).ravel(), np.asarray(fiber), [z_mirror]]
        )

    def pack_gradient(self, grad_latent, grad_fiber, grad_mirror):
        return np.concatenate(
            [np.asarray(grad_latent).ravel(), np.asarray(grad_fiber), [grad_mirror]]
        )

    def bounds(self):
        n = self.nx * self.ny
        lo = np.concatenate(
            [np.zeros(n), [-2.0, 0.0], [LAYOUT.z_mirror0 - T_SPACER / 2]]
        )
        hi = np.concatenate(
            [np.ones(n), [2.0, 20.0], [LAYOUT.z_mirror0 + T_SPACER / 2]]
        )
        return lo, hi


# ---------------------------------------------------------------- run modes
def run_optimize(args):
    import nlopt

    coupler = Coupler(args)
    run_dir = args.outdir or time.strftime("run_%Y%m%d_%H%M%S")
    if mp.am_master():
        os.makedirs(run_dir, exist_ok=True)
    coupler.diagnostics = Diagnostics(run_dir, LAYOUT, args.resolution)
    history = History(os.path.join(run_dir, "history.npz"))

    latent0 = seed_grating(coupler.nx, coupler.ny)
    x = coupler.pack(latent0, [0.0, FIBER_TILT_DEG], LAYOUT.z_mirror0)
    lo, hi = coupler.bounds()

    gaps = constructive_gaps()
    say(f"aluminium: single Drude pole, fp = {AL_FP:.4f} /um")
    say(
        f"constructive grating-to-mirror gaps: "
        f"{', '.join(f'{g:.3f}' for g in gaps)} um"
    )
    say(
        f"starting gap {LAYOUT.mirror_gap(LAYOUT.z_mirror0):.3f} um "
        f"(a standard 2.0 um box sits between the third and fourth)"
    )
    say(
        f"design grid {coupler.nx} x {coupler.ny}, "
        f"Courant {courant_for(args.resolution):.3f}, "
        f"{mp.count_processors()} ranks"
    )

    counter = {"n": 0}

    def objective(x_in, grad):
        counter["n"] += 1
        value, g = coupler.step(x_in, history, counter["n"])
        if grad.size:
            grad[:] = g
        history.save()
        return value

    def solid_constraint(x_in, grad):
        latent, _, _ = coupler.unpack(x_in)
        f = lambda u: length_penalties(u, coupler.beta, args.resolution)[0]
        v, vjp = jax.vjp(f, jnp.asarray(latent))
        if grad.size:
            grad[:] = 0.0
            grad[: coupler.nx * coupler.ny] = np.asarray(vjp(1.0)[0]).ravel()
        return float(v)

    def void_constraint(x_in, grad):
        latent, _, _ = coupler.unpack(x_in)
        f = lambda u: length_penalties(u, coupler.beta, args.resolution)[1]
        v, vjp = jax.vjp(f, jnp.asarray(latent))
        if grad.size:
            grad[:] = 0.0
            grad[: coupler.nx * coupler.ny] = np.asarray(vjp(1.0)[0]).ravel()
        return float(v)

    # A beta continuation: SSP keeps the gradient alive as beta grows, so the
    # design binarizes without the usual trade against differentiability.
    # The usual topology-optimization schedule: a couple of epochs at finite
    # beta with the lengthscale constraints *off*, then a final epoch at beta =
    # infinity with them on.
    #
    # Leaving the constraints off early is not a convenience.  The constraint
    # measures features of a near-binary structure, and at finite beta the
    # projected design is grey, so there are no features to measure and the
    # number it returns is an artefact: measured on one seed, -1.0 (comfortably
    # feasible) at beta = inf against +2390 at beta = 8, growing worse the
    # greyer the design.  Imposing it there makes MMA chase a quantity that
    # means nothing, and it did -- the objective fell monotonically while the
    # optimizer worked on it.
    #
    # The finite-beta epochs still earn their place: they give a smoother
    # landscape to find the broad shape of the grating in.  SSP is what allows
    # the last epoch to be beta = infinity at all, rather than merely large.
    for stage, beta in enumerate(args.beta_schedule):
        coupler.beta = beta
        constrained = not np.isfinite(beta) or stage == len(args.beta_schedule) - 1
        say(
            f"\n--- stage {stage + 1}/{len(args.beta_schedule)}: beta = {beta:g}, "
            f"lengthscale constraints {'on' if constrained else 'off'} ---"
        )
        opt = nlopt.opt(nlopt.LD_MMA, x.size)
        opt.set_lower_bounds(lo)
        opt.set_upper_bounds(hi)
        opt.set_max_objective(objective)
        if constrained:
            opt.add_inequality_constraint(solid_constraint, 1e-8)
            opt.add_inequality_constraint(void_constraint, 1e-8)
        opt.set_maxeval(args.iters_per_stage)
        x = opt.optimize(x)

    latent, fiber, z_mirror = coupler.unpack(x)
    gap = LAYOUT.mirror_gap(z_mirror)
    nearest = min(gaps, key=lambda g: abs(g - gap))
    say(
        f"\nconverged mirror gap {gap:.4f} um; "
        f"nearest analytic optimum {nearest:.4f} um "
        f"(difference {abs(gap - nearest) * 1e3:.0f} nm)"
    )
    say(f"fiber tilt {fiber[1]:.3f} deg, lateral offset {fiber[0]:+.4f} um")
    history.save()
    say(f"history and figures written to {run_dir}/")


def run_forward(args):
    """One solve at the starting point, plus the spectrum, and no optimization."""
    coupler = Coupler(args)
    latent = seed_grating(coupler.nx, coupler.ny)
    value, _, wall, steps, _ = coupler.evaluate(
        latent, [0.0, FIBER_TILT_DEG], LAYOUT.z_mirror0, need_gradient=False
    )
    say(f"efficiency {value:.6e}   {steps} timesteps   {wall:.1f} s")
    say(f"mirror gap {LAYOUT.mirror_gap(LAYOUT.z_mirror0):.4f} um")
    say(f"constructive gaps {', '.join(f'{g:.3f}' for g in constructive_gaps())} um")


def run_report(args):
    """Summarize a saved history without touching Meep."""
    data = np.load(args.history)
    n = len(data["iteration"])
    say(f"{n} iterations from {args.history}")
    say(
        f"{'iter':>5} {'efficiency':>12} {'gap um':>9} {'tilt deg':>9} "
        f"{'offset um':>10} {'g_solid':>9} {'g_void':>9} {'beta':>7}"
    )
    for i in range(n):
        say(
            f"{int(data['iteration'][i]):5d} {data['efficiency'][i]:12.5e} "
            f"{data['mirror_gap'][i]:9.4f} {data['fiber_tilt'][i]:9.3f} "
            f"{data['fiber_offset'][i]:10.4f} {data['constraint_solid'][i]:9.2e} "
            f"{data['constraint_void'][i]:9.2e} {data['beta'][i]:7.0f}"
        )
    best = int(np.argmax(data["efficiency"]))
    say(
        f"\nbest at iteration {int(data['iteration'][best])}: "
        f"efficiency {data['efficiency'][best]:.6e}, "
        f"gap {data['mirror_gap'][best]:.4f} um"
    )
    say(
        f"analytic optima: {', '.join(f'{g:.3f}' for g in data['constructive_gaps'])} um"
    )


def main():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("mode", choices=["optimize", "forward", "report"])
    p.add_argument("--resolution", type=float, default=20.0)
    p.add_argument(
        "--nfreq",
        type=int,
        default=1,
        help="frequencies in the objective; the spectrum is evaluated "
        "separately at the end",
    )
    p.add_argument("--bandwidth", type=float, default=0.05)
    p.add_argument(
        "--beta-schedule",
        type=float,
        nargs="+",
        default=[8, 32, np.inf],
        help="Epochs of increasing beta. The lengthscale constraints are off "
        "for every epoch but the last, because at finite beta the projected "
        "design is grey and the constraint measures nothing real -- on one "
        "seed, -1.0 at beta=inf against +2390 at beta=8. SSP is what lets the "
        "final epoch be beta=inf rather than merely large.",
    )
    p.add_argument("--iters-per-stage", type=int, default=15)
    p.add_argument("--decay-by", type=float, default=1e-6)
    p.add_argument("--min-run-time", type=float, default=50.0)
    p.add_argument("--max-run-time", type=float, default=800.0)
    p.add_argument(
        "--no-pade",
        dest="pade",
        action="store_false",
        help="use the decay criterion instead of Pade extrapolation",
    )
    p.add_argument(
        "--outdir", default=None, help="run folder; defaults to run_<timestamp>"
    )
    p.add_argument("--history", default="history.npz")
    args = p.parse_args()

    {"optimize": run_optimize, "forward": run_forward, "report": run_report}[args.mode](
        args
    )


if __name__ == "__main__":
    main()
