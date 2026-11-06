"""Grating coupler radiating through a thick glass superstrate into a fiber.

A two-etch silicon grating coupler radiates upward, propagates a few hundred
microns through silica, crosses the silica/air interface, and is collected by a
fiber. The propagation and the interface are handled analytically by
`meep.adjoint.AngularSpectrum` rather than being meshed, which is what makes the
problem tractable: at resolution 20 a 300 um tall cell would be several hundred
million pixels, and the interface would put Meep's near-to-far transformation out
of reach anyway, since it requires a homogeneous medium.

The FDTD cell therefore stops less than a micron above the device layer. Only the
grating is simulated; everything above the monitor is a layer stack.

Two modes:

    python grating_coupler_asm.py forward
        One simulation, then the far field. Prints the coupling efficiency and
        the monitor diagnostics, and sweeps fiber tilt, working distance and
        superstrate thickness -- none of which re-runs the simulation, because
        they live entirely in the analytic propagator.

    python grating_coupler_asm.py optimize
        Topology optimization of the two etch layers against the fiber overlap.

Run with --help for the geometry and solver knobs.
"""

import argparse
import math

import jax
import jax.numpy as jnp
import numpy as np

import meep as mp
import meep.adjoint as mpa

jax.config.update("jax_enable_x64", True)

# ------------------------------------------------------------------ geometry
N_SI, N_OXIDE, N_AIR = 3.48, 1.444, 1.0
T_DEVICE = 0.22  # silicon device layer
T_DEEP = 0.15  # lower etch level: the deep grating teeth
T_SHALLOW = T_DEVICE - T_DEEP  # upper etch level, for directionality
T_BOX = 2.0  # buried oxide
T_HANDLE = 1.0  # how much silicon handle to include below the box

WAVEGUIDE_LENGTH = 4.0
MONITOR_STANDOFF = 1.0  # monitor height above the device layer
DPML = 1.0

# The monitor plane must be clear of the near field, and the transform is
# periodic, so the monitor also has to be wide enough that the radiated field
# has decayed by its ends -- otherwise it wraps around. `report()` measures both,
# and `check_monitor` below refuses to report numbers built on a wrapped field.
# The cell is therefore padded well beyond the grating aperture.
PAD = 12.0
MONITOR_MARGIN = 0.5
EDGE_AMPLITUDE_LIMIT = 1e-2


def build_arguments():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("mode", choices=["forward", "optimize"])
    parser.add_argument(
        "--resolution", type=float, default=20.0, help="pixels per micron"
    )
    parser.add_argument(
        "--aperture", type=float, default=20.0, help="width of the grating, in microns"
    )
    parser.add_argument(
        "--superstrate",
        type=float,
        default=300.0,
        help="silica thickness above the device layer, in microns",
    )
    parser.add_argument(
        "--working-distance",
        type=float,
        default=5.0,
        help="fiber facet height above the silica/air interface",
    )
    parser.add_argument(
        "--fiber-waist",
        type=float,
        default=None,
        help="fiber 1/e field radius; defaults to the diffracted "
        "beam size at the facet",
    )
    parser.add_argument(
        "--tilt", type=float, default=8.0, help="fiber tilt from normal, in degrees"
    )
    parser.add_argument(
        "--wavelengths", type=float, nargs="+", default=[1.52, 1.55, 1.58]
    )
    parser.add_argument(
        "--pad-factor",
        type=int,
        default=16,
        help="zero padding before the transform; a beam that "
        "spreads over hundreds of microns needs a wide window",
    )
    parser.add_argument("--iterations", type=int, default=15)
    parser.add_argument("--design-resolution", type=float, default=None)
    return parser.parse_args()


def cell_geometry(args):
    """Returns the cell size, center, and the x extent of the design regions."""
    x_min = -(args.aperture / 2 + WAVEGUIDE_LENGTH + PAD + DPML)
    x_max = args.aperture / 2 + PAD + DPML
    y_min = -(T_BOX + T_HANDLE + DPML)
    # Leave a wavelength of oxide above the monitor before the PML starts;
    # a plane pressed up against the absorber picks up its residual
    # reflection, and near-grazing rays have nowhere to go.
    y_max = T_DEVICE + MONITOR_STANDOFF + 1.2 + DPML
    size = mp.Vector3(x_max - x_min, y_max - y_min)
    center = mp.Vector3(0.5 * (x_min + x_max), 0.5 * (y_min + y_max))
    return size, center


def design_regions(args):
    """The two etch levels, stacked, spanning the grating aperture.

    Splitting the device layer into a deep and a shallow level is what lets the
    grating radiate preferentially upward instead of symmetrically; a single
    fully etched layer is limited to roughly half its power going the wrong way.
    """
    resolution = args.design_resolution or (2 * args.resolution)
    nx = int(round(args.aperture * resolution)) + 1
    oxide, silicon = mp.Medium(index=N_OXIDE), mp.Medium(index=N_SI)

    grids, volumes = [], []
    for y_lo, y_hi in ((0.0, T_DEEP), (T_DEEP, T_DEVICE)):
        grids.append(
            mp.MaterialGrid(
                mp.Vector3(nx, 1, 1),
                oxide,
                silicon,
                weights=np.ones((nx,)),
                do_averaging=False,
                beta=0,
            )
        )
        volumes.append(
            mp.Volume(
                center=mp.Vector3(0, 0.5 * (y_lo + y_hi)),
                size=mp.Vector3(args.aperture, y_hi - y_lo),
            )
        )
    return grids, volumes, nx


def build_simulation(args, weights=None):
    """Assembles the simulation. The cell stops just above the device layer."""
    size, center = cell_geometry(args)
    grids, volumes, nx = design_regions(args)
    if weights is not None:
        for grid, w in zip(grids, weights):
            grid.update_weights(np.asarray(w).ravel())

    silicon, oxide = mp.Medium(index=N_SI), mp.Medium(index=N_OXIDE)
    geometry = [
        # silicon handle below the buried oxide
        mp.Block(
            center=mp.Vector3(center.x, -(T_BOX + T_HANDLE + DPML) / 2 - T_BOX / 2),
            size=mp.Vector3(mp.inf, T_HANDLE + DPML),
            material=silicon,
        ),
        # the input waveguide, running in from the left
        mp.Block(
            center=mp.Vector3(
                -(args.aperture / 2 + WAVEGUIDE_LENGTH + PAD + DPML) / 2
                - args.aperture / 4,
                T_DEVICE / 2,
            ),
            size=mp.Vector3(2 * (WAVEGUIDE_LENGTH + PAD + DPML), T_DEVICE),
            material=silicon,
        ),
    ] + [
        mp.Block(center=volume.center, size=volume.size, material=grid)
        for grid, volume in zip(grids, volumes)
    ]

    frequencies = [1.0 / w for w in args.wavelengths]
    center_frequency = float(np.mean(frequencies))
    source = [
        mp.EigenModeSource(
            mp.GaussianSource(center_frequency, fwidth=0.2 * center_frequency),
            center=mp.Vector3(
                -(args.aperture / 2 + WAVEGUIDE_LENGTH * 0.6), T_DEVICE / 2
            ),
            size=mp.Vector3(0, 6 * T_DEVICE),
            eig_band=1,
            eig_parity=mp.ODD_Z,
            eig_match_freq=True,
        )
    ]

    simulation = mp.Simulation(
        cell_size=size,
        geometry_center=center,
        resolution=args.resolution,
        boundary_layers=[mp.PML(DPML)],
        geometry=geometry,
        sources=source,
        default_material=oxide,
        dimensions=2,
    )
    return simulation, grids, volumes, nx, frequencies


def monitor_volume(args):
    """The plane the radiated field is read on, inside the homogeneous oxide."""
    size, center = cell_geometry(args)
    width = size.x - 2 * DPML - 2 * MONITOR_MARGIN
    return mp.Volume(
        center=mp.Vector3(center.x, T_DEVICE + MONITOR_STANDOFF),
        size=mp.Vector3(width, 0),
    )


def build_stack(args):
    """Everything above the monitor: the rest of the silica, then air."""
    remaining_oxide = args.superstrate - MONITOR_STANDOFF
    if remaining_oxide <= 0:
        raise ValueError(
            f"--superstrate {args.superstrate} must exceed the monitor standoff "
            f"{MONITOR_STANDOFF}."
        )
    return (
        mpa.Stack([mpa.Layer(N_OXIDE, remaining_oxide), mpa.Layer(N_AIR)]),
        remaining_oxide,
    )


def diffracted_waist(args):
    """A rough estimate of the beam radius at the fiber facet.

    Used only to pick a sensible default fiber mode: a beam launched from an
    aperture this size and allowed to diffract for this far is not going to be
    collected by a standard single-mode fiber.
    """
    waist = args.aperture / 4
    wavelength = float(np.mean(args.wavelengths))
    rayleigh = math.pi * waist**2 * N_OXIDE / wavelength
    distance = args.superstrate - MONITOR_STANDOFF + args.working_distance
    return waist * math.sqrt(1 + (distance / rayleigh) ** 2)


def initial_weights(args, nx):
    """A uniform grating, as a starting point for the optimizer."""
    x = np.linspace(-args.aperture / 2, args.aperture / 2, nx)
    period = float(np.mean(args.wavelengths)) / (
        N_OXIDE * math.sin(math.radians(args.tilt)) + 2.6
    )
    # A binary grating, not a smoothly graded one: a sinusoidal index
    # modulation is a much weaker scatterer and radiates very little.
    deep = (np.cos(2 * math.pi * x / period) > 0).astype(float)
    shallow = (np.cos(2 * math.pi * x / period - 0.6) > -0.3).astype(float)
    return [deep, shallow]


def check_monitor(propagator, simulation, monitor):
    """Prints the monitor diagnostics, and complains if they are bad.

    These are worth reading before any efficiency is believed. The transform is
    periodic, so a field that has not decayed by the ends of the monitor wraps
    around and contaminates everything downstream -- and it does so quietly,
    producing a plausible-looking number that does not improve with resolution.
    """
    report = propagator.report_monitor(simulation, monitor)
    print("\nmonitor diagnostics (per wavelength)")
    for key, value in report.items():
        print(f"  {key:22s} " + "  ".join(f"{v:.2e}" for v in np.atleast_1d(value)))

    edge = float(np.max(report["edge_amplitude"]))
    if edge > 0.1:
        print(f"\n  WARNING: the field is still {edge:.1%} of its peak at the ends of")
        print("  the monitor, so the transform is wrapping badly and the numbers")
        print("  below are not trustworthy. Widen the cell by increasing PAD.")
    elif edge > EDGE_AMPLITUDE_LIMIT:
        print(f"\n  Note: {edge:.1%} of the peak field remains at the ends of the")
        print("  monitor. A few percent is normal here and is not the beam tail --")
        print("  it is near-grazing radiation, which travels sideways rather than")
        print("  decaying, so widening the cell barely helps. It bounds the")
        print("  accuracy of the efficiencies below at a similar level. Power at")
        print("  those angles never reaches the fiber in any case.")
    else:
        print(f"  edge amplitude {edge:.1e} is small; the field has decayed.")
    return report


# --------------------------------------------------------------------- modes
def run_forward(args):
    """One simulation, then everything the analytic propagator gives for free."""
    simulation, grids, volumes, nx, frequencies = build_simulation(
        args, weights=initial_weights(args, design_regions(args)[2])
    )
    volume = monitor_volume(args)
    monitor = simulation.add_dft_fields([mp.Ez, mp.Hx], frequencies, where=volume)
    simulation.run(
        until_after_sources=mp.stop_when_dft_decayed(1e-9, minimum_run_time=50)
    )

    stack, remaining_oxide = build_stack(args)
    propagator = mpa.AngularSpectrum.from_monitor(
        simulation, monitor, stack, volume, pad_factor=args.pad_factor
    )

    check_monitor(propagator, simulation, monitor)

    waist = args.fiber_waist or diffracted_waist(args)
    distance = remaining_oxide + args.working_distance
    fiber = mpa.gaussian_mode(waist, tilt_deg=args.tilt)
    efficiency = propagator.overlap_monitor(simulation, monitor, fiber, distance)

    print(f"\nfiber waist {waist:.2f} um at {distance:.1f} um, tilt {args.tilt}deg")
    for wavelength, value in zip(args.wavelengths, np.atleast_1d(efficiency)):
        print(
            f"  {wavelength*1e3:.0f} nm   modal purity = {value:.4f}"
            f"   ({10*np.log10(max(value,1e-12)):+.2f} dB)"
        )

    # None of the sweeps below re-runs the simulation: the fiber, the working
    # distance and the superstrate are all parameters of the analytic stack.
    print("\nvs fiber tilt (no further simulation)")
    for tilt in np.arange(args.tilt - 4, args.tilt + 4.1, 2.0):
        value = propagator.overlap_monitor(
            simulation,
            monitor,
            mpa.gaussian_mode(waist, tilt_deg=float(tilt)),
            distance,
        )
        print(f"  {tilt:5.1f} deg   {np.mean(value):.4f}")

    print("\nvs working distance (no further simulation)")
    for extra in (0.0, 5.0, 20.0, 50.0):
        value = propagator.overlap_monitor(
            simulation, monitor, fiber, remaining_oxide + extra
        )
        print(f"  {extra:5.1f} um    {np.mean(value):.4f}")

    print("\nvs superstrate thickness (rebuilds only the stack)")
    for thickness in (100.0, 200.0, 300.0, 500.0):
        alternative = mpa.Stack(
            [mpa.Layer(N_OXIDE, thickness - MONITOR_STANDOFF), mpa.Layer(N_AIR)]
        )
        other = mpa.AngularSpectrum.from_monitor(
            simulation, monitor, alternative, volume, pad_factor=args.pad_factor
        )
        value = other.overlap_monitor(
            simulation,
            monitor,
            fiber,
            thickness - MONITOR_STANDOFF + args.working_distance,
        )
        print(f"  {thickness:5.0f} um   {np.mean(value):.4f}")


def run_optimize(args):
    """Topology optimization of both etch layers against the fiber overlap."""
    import nlopt

    simulation, grids, volumes, nx, frequencies = build_simulation(args)
    volume = monitor_volume(args)
    stack, remaining_oxide = build_stack(args)
    waist = args.fiber_waist or diffracted_waist(args)
    distance = remaining_oxide + args.working_distance

    # get_array_metadata needs the fields allocated, and it is also the only
    # way to learn how many samples Meep will actually put on the plane after
    # snapping it to the grid.
    simulation.init_sim()
    propagator = mpa.AngularSpectrum(
        stack,
        frequencies,
        pitch=1.0 / args.resolution,
        num_points=len(simulation.get_array_metadata(vol=volume)[0]),
        normal=mp.Y,
        pad_factor=args.pad_factor,
    )
    fiber = mpa.gaussian_mode(waist, tilt_deg=args.tilt)

    # The objective is written with jax.numpy; OptimizationProblem recognizes
    # that from the array type it returns and differentiates it with jax.vjp.
    def objective(*monitor_values):
        fields = propagator.take(monitor_values)
        return jnp.mean(propagator.overlap(fields, fiber, distance))

    optimization = mpa.OptimizationProblem(
        simulation=simulation,
        objective_functions=objective,
        objective_arguments=propagator.objective_arguments(simulation, volume),
        design_regions=[
            mpa.DesignRegion(grid, volume=v) for grid, v in zip(grids, volumes)
        ],
        frequencies=frequencies,
        decay_by=1e-7,
    )

    weights = initial_weights(args, nx)
    history = []

    def evaluate(x, gradient):
        split = np.split(x, 2)
        value, gradients = optimization([split[0], split[1]])
        if gradient.size > 0:
            # One gradient per design region, each (num weights, num frequencies)
            # -- or (num weights,) when Meep squeezes a single frequency away.
            # Reshaping rather than atleast_2d keeps the weights on the first
            # axis in both cases.
            gradient[:] = np.concatenate(
                [np.asarray(g).reshape(nx, -1).sum(axis=1) for g in gradients]
            )
        history.append(float(value))
        print(f"  iteration {len(history):3d}   mean efficiency = {float(value):.5f}")
        return float(value)

    solver = nlopt.opt(nlopt.LD_MMA, 2 * nx)
    solver.set_lower_bounds(0.0)
    solver.set_upper_bounds(1.0)
    solver.set_max_objective(evaluate)
    solver.set_maxeval(args.iterations)
    solver.optimize(np.concatenate(weights))

    print(f"\nstarted at {history[0]:.5f}, finished at {max(history):.5f}")


if __name__ == "__main__":
    arguments = build_arguments()
    if arguments.mode == "forward":
        run_forward(arguments)
    else:
        run_optimize(arguments)
