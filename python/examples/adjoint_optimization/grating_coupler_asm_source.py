"""Grating coupler driven *from the fiber*, with the fiber outside the cell.

The companion example, `grating_coupler_asm.py`, runs the coupler in the
transmit direction: light comes in along the waveguide, radiates up, and the
angular-spectrum propagator carries the radiated field hundreds of microns to a
fiber. This one runs it in the receive direction. A fiber mode is specified
where it is physically meaningful -- at the facet, far above the chip and
outside the FDTD cell -- carried *down* through the silica and across the
silica/air interface analytically, and injected at a plane just above the
device layer as equivalent surface currents.

Nothing about the cell changes: the superstrate and the working distance stay
out of the simulation in both directions. What changes is which end of the
chain is the source.

The whole path is differentiable. The fiber's waist, tilt and lateral offset
are ordinary JAX values that reach Meep only through the injected currents, and
Meep returns the cotangent with respect to those currents; JAX carries it back
the rest of the way. So the grating and the fiber alignment can be optimized
together, which is the thing the transmit-direction example cannot do.

    # what the fiber launches, and where it is pointing
    python grating_coupler_asm_source.py forward

    # co-optimize the grating and the fiber alignment
    python grating_coupler_asm_source.py optimize --iterations 20

The geometry helpers are imported from the transmit-direction example rather
than duplicated; the two describe the same device.
"""

import argparse
import numpy as np

import meep as mp
import meep.adjoint as mpa

import grating_coupler_asm as tx

WAVEGUIDE_MODE_STANDOFF = 0.6  # where the waveguide mode is measured, as a
# fraction of the waveguide run-in


def build_arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("forward", "optimize"))
    parser.add_argument("--resolution", type=int, default=30)
    parser.add_argument("--aperture", type=float, default=12.0)
    parser.add_argument("--wavelengths", type=float, nargs="+", default=[1.55])
    parser.add_argument(
        "--superstrate",
        type=float,
        default=300.0,
        help="silica thickness above the device layer, in microns",
    )
    parser.add_argument(
        "--working-distance",
        type=float,
        default=50.0,
        help="air gap between the silica surface and the fiber facet",
    )
    parser.add_argument("--waist", type=float, default=5.2, help="fiber MFD / 2")
    parser.add_argument("--tilt", type=float, default=8.0, help="fiber tilt, degrees")
    parser.add_argument("--offset", type=float, default=0.0)
    parser.add_argument("--iterations", type=int, default=20)
    parser.add_argument("--learning-rate", type=float, default=0.02)
    # consumed by the geometry helpers imported from the transmit example
    parser.add_argument("--design-resolution", type=float, default=None)
    parser.add_argument("--pad-factor", type=int, default=16)
    parser.add_argument("--fiber-waist", type=float, default=None)
    return parser.parse_args()


def source_plane(args):
    """The plane the fiber's field is injected on.

    The same plane the transmit example reads the radiated field from: inside
    the homogeneous oxide, clear of the near field, and wide enough that the
    beam has decayed by its ends.
    """
    return tx.monitor_volume(args)


def injection_propagator(simulation, args, frequency):
    """A propagator for the source plane, facing up toward the fiber.

    `from_volume` rather than `from_monitor`: a DFT monitor needs the field
    components to be allocated, which does not happen until sources are added,
    and the sources are exactly what this is being built to produce.
    """
    stack, _ = tx.build_stack(args)
    return mpa.AngularSpectrum.from_volume(
        simulation,
        simulation._fit_volume_to_simulation(source_plane(args)),
        stack,
        [frequency],
        sign=+1,  # the fiber is in the +y direction
        pad_factor=4,
    )


def fiber_distance(args):
    """From the source plane to the fiber facet."""
    return (args.superstrate - tx.MONITOR_STANDOFF) + args.working_distance


def waveguide_monitor(args):
    """Where the coupled power is measured, in the input waveguide."""
    return mp.Volume(
        center=mp.Vector3(
            -(args.aperture / 2 + tx.WAVEGUIDE_LENGTH * WAVEGUIDE_MODE_STANDOFF),
            tx.T_DEVICE / 2,
        ),
        size=mp.Vector3(0, 6 * tx.T_DEVICE),
    )


def build_simulation(args, weights=None):
    """The coupler with no sources; the fiber's currents are added separately."""
    simulation, grids, volumes, nx, frequencies = tx.build_simulation(args, weights)
    # The transmit example drives the waveguide. Here the fiber drives the
    # grating instead, so that source goes away.
    simulation.sources = []
    return simulation, grids, volumes, nx, frequencies


def fiber_sources(simulation, args, frequency, mode=None):
    """Equivalent surface currents for the fiber mode, at the source plane."""
    propagator = injection_propagator(simulation, args, frequency)
    if mode is None:
        mode = mpa.gaussian_mode(
            waist=args.waist, tilt_deg=args.tilt, offset=args.offset
        )
    fields = propagator.incident_fields(mode, distance=fiber_distance(args))
    plane = source_plane(args)
    return (
        propagator,
        fields,
        propagator.equivalent_sources(
            fields,
            mp.GaussianSource(frequency, fwidth=0.2 * frequency),
            center=plane.center,
            size=plane.size,
        ),
    )


def run_forward(args):
    """Inject the fiber mode and report how much reaches the waveguide."""
    simulation, _, _, _, frequencies = build_simulation(
        args, weights=tx.initial_weights(args, tx.design_regions(args)[2])
    )
    frequency = float(np.mean(frequencies))
    propagator, fields, sources = fiber_sources(simulation, args, frequency)

    report = propagator.report(fields)
    print("\nthe field the fiber puts on the source plane")
    print(
        f"  travelling toward the chip   {float(np.ravel(report['downgoing_fraction'])[0]):.4f}"
    )
    print(
        f"  evanescent                   {float(np.ravel(report['evanescent_fraction'])[0]):.3e}"
    )
    print(
        f"  amplitude at the plane edges {float(np.ravel(report['edge_amplitude'])[0]):.3e}"
    )
    if float(np.ravel(report["edge_amplitude"])[0]) > tx.EDGE_AMPLITUDE_LIMIT:
        print(
            "  ^ the beam has not decayed by the ends of the plane, so the "
            "injected field wraps around. Widen the cell or narrow the beam."
        )

    simulation.change_sources(sources)
    monitor = simulation.add_mode_monitor(
        [frequency],
        mp.ModeRegion(volume=waveguide_monitor(args)),
    )
    simulation.run(until_after_sources=mp.stop_when_dft_decayed(1e-9))

    coefficients = simulation.get_eigenmode_coefficients(
        monitor, [1], eig_parity=mp.ODD_Z
    )
    # the backward-going coefficient: the waveguide runs off to the left
    coupled = abs(coefficients.alpha[0, 0, 1]) ** 2
    print(f"\npower coupled into the waveguide mode: {coupled:.6e}")
    return coupled


def run_optimize(args):
    """Co-optimize the grating and the fiber alignment."""
    import jax

    jax.config.update("jax_enable_x64", True)
    import jax.numpy as jnp
    import optax

    simulation, grids, volumes, nx, frequencies = build_simulation(args)
    frequency = float(np.mean(frequencies))
    plane = source_plane(args)
    propagator = injection_propagator(simulation, args, frequency)
    distance = fiber_distance(args)

    design_regions = [
        mpa.DesignRegion(grid, volume=volume) for grid, volume in zip(grids, volumes)
    ]
    monitor = mpa.EigenmodeCoefficient(
        simulation, waveguide_monitor(args), mode=1, forward=False
    )

    # A placeholder source, replaced on every iteration by the currents JAX
    # computes. Its amplitudes are what the wrapper differentiates with respect
    # to; the fiber parameters live upstream and Meep never sees them.
    _, _, sources = fiber_sources(simulation, args, frequency)

    wrapper = mpa.MeepJaxWrapper(
        simulation,
        sources,
        [monitor],
        design_regions,
        [frequency],
        minimum_run_time=200.0,
    )

    def currents(fiber):
        """The amp_data array for each equivalent-current sheet, in JAX."""
        mode = mpa.gaussian_mode(
            waist=fiber["waist"], tilt_deg=fiber["tilt"], offset=fiber["offset"]
        )
        fields = propagator.incident_fields(mode, distance=distance)
        # get_equiv_sources puts n x H on the electric components and -n x E on
        # the magnetic ones, so the two are swapped relative to each other.
        swap = {
            mp.Ex: mp.Hx,
            mp.Ey: mp.Hy,
            mp.Ez: mp.Hz,
            mp.Hx: mp.Ex,
            mp.Hy: mp.Ey,
            mp.Hz: mp.Ez,
        }
        out = []
        for s in sources:
            paired = swap[s.component]
            store = fields.H if mp.is_electric(s.component) else fields.E
            out.append(jnp.ravel(store[paired][0]))
        return out

    def loss(params):
        values = wrapper(
            [jnp.asarray(w) for w in params["weights"]], currents(params["fiber"])
        )
        return -jnp.sum(jnp.abs(jnp.asarray(values)) ** 2)

    params = {
        "weights": [jnp.asarray(w) for w in tx.initial_weights(args, nx)],
        "fiber": {
            "waist": jnp.asarray(args.waist),
            "tilt": jnp.asarray(args.tilt),
            "offset": jnp.asarray(args.offset),
        },
    }
    optimizer = optax.adam(args.learning_rate)
    state = optimizer.init(params)

    for iteration in range(args.iterations):
        value, gradient = jax.value_and_grad(loss)(params)
        updates, state = optimizer.update(gradient, state, params)
        params = optax.apply_updates(params, updates)
        params["weights"] = [jnp.clip(w, 0.0, 1.0) for w in params["weights"]]
        fiber = params["fiber"]
        print(
            f"{iteration:3d}  coupling {-float(value):.6e}  "
            f"waist {float(fiber['waist']):.3f}  "
            f"tilt {float(fiber['tilt']):.3f}  "
            f"offset {float(fiber['offset']):+.3f}"
        )
    return params


if __name__ == "__main__":
    arguments = build_arguments()
    if arguments.mode == "forward":
        run_forward(arguments)
    else:
        run_optimize(arguments)
