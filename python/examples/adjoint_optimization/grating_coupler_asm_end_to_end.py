"""End-to-end check of angular-spectrum injection, in a small 2D grating coupler.

A fiber mode is specified in air above a silica superstrate, outside the FDTD
cell, carried down through the stack by `AngularSpectrum.incident_fields`, and
injected as equivalent surface currents just above a silicon grating. Three
independent checks, each against something the injection path does not use:

  1. The field Meep launches matches the propagator's prediction a few microns
     below the injection plane, in a uniform cell, with no fitted parameters:
     only the source's own spectrum, which Meep applies to every source.
  2. Reciprocity. Coupling fiber -> waveguide (inject the fiber mode, measure
     the waveguide mode) equals coupling waveguide -> fiber (launch the
     waveguide mode, project the radiated field onto the fiber mode through
     the same stack), at every fiber tilt. The two agree to discretization
     error, which shrinks at second order with resolution.
  3. The gradient of the coupled power with respect to the fiber's waist, tilt
     and offset -- through `incident_fields`, `equivalent_currents`,
     `MeepJaxWrapper` and Meep's source adjoint -- matches finite differences
     of the full simulation.

    python grating_coupler_asm_end_to_end.py            # ~20 minutes, one core
    python grating_coupler_asm_end_to_end.py --quick    # coarser, ~2 minutes
    python grating_coupler_asm_end_to_end.py --replot   # redraw from saved results
"""

import argparse

import jax
import jax.numpy as jnp
import numpy as np

import meep as mp
import meep.adjoint as mpa
from meep.adjoint import angular_spectrum as asm
from meep.adjoint import source_gradient

jax.config.update("jax_enable_x64", True)
mp.verbosity(0)

WAVELENGTH = 1.55
FCEN = 1 / WAVELENGTH
N_OXIDE, N_AIR, N_SILICON = 1.444, 1.0, 3.48

CELL = mp.Vector3(26, 9)
PML = 1.0
PLANE = mp.Volume(center=mp.Vector3(0, 2.5), size=mp.Vector3(12, 0))  # injection
STACK = mpa.Stack([mpa.Layer(N_OXIDE, 2.0), mpa.Layer(N_AIR)])  # above PLANE
FIBER_DISTANCE = 12.0  # from PLANE: 2 um of silica, then 10 um of air
WAIST = 2.5

WG_Y, WG_T = -1.5, 0.22  # silicon waveguide, fully etched grating
PERIOD, NUM_TEETH, GRATING_X0 = 0.62, 10, -3.0
SOURCE_X, MONITOR_X = -11.0, -10.0  # far from the beam, so free-space light
WG_LINE = mp.Volume(center=mp.Vector3(MONITOR_X, WG_Y), size=mp.Vector3(0, 2.5))

SRC_TIME = mp.GaussianSource(FCEN, fwidth=0.2 * FCEN)
DECAY = mp.stop_when_dft_decayed(1e-9)


def device():
    silicon = mp.Medium(index=N_SILICON)
    geometry = [
        mp.Block(
            center=mp.Vector3(GRATING_X0 - 50, WG_Y),
            size=mp.Vector3(100, WG_T),
            material=silicon,
        )
    ]
    for i in range(NUM_TEETH):
        geometry.append(
            mp.Block(
                center=mp.Vector3(GRATING_X0 + (i + 0.25) * PERIOD, WG_Y),
                size=mp.Vector3(0.5 * PERIOD, WG_T),
                material=silicon,
            )
        )
    return geometry


def simulation(resolution, geometry=(), sources=()):
    return mp.Simulation(
        cell_size=CELL,
        resolution=resolution,
        boundary_layers=[mp.PML(PML)],
        default_material=mp.Medium(index=N_OXIDE),
        geometry=list(geometry),
        sources=list(sources),
        force_complex_fields=True,
    )


def fiber_mode(tilt, offset=0.0, waist=WAIST):
    return mpa.gaussian_mode(waist, tilt_deg=tilt, offset=offset)


def injection(resolution, mode):
    """The propagator for PLANE and the fiber's incident fields there."""
    sim = simulation(resolution)
    propagator = mpa.AngularSpectrum.from_volume(
        sim, sim._fit_volume_to_simulation(PLANE), STACK, [FCEN], pad_factor=8
    )
    return propagator, propagator.incident_fields(mode, distance=FIBER_DISTANCE)


def fiber_sources(propagator, fields):
    return propagator.equivalent_sources(
        fields, SRC_TIME, center=PLANE.center, size=PLANE.size
    )


def homogeneous(propagator, sign):
    """A propagator through uniform silica on the same samples."""
    return mpa.AngularSpectrum(
        mpa.Stack([mpa.Layer(N_OXIDE)]),
        [FCEN],
        propagator.pitch,
        propagator.num_points,
        normal=mp.Y,
        sign=sign,
        pad_factor=8,
    )


def waveguide_power(sim, flux, direction):
    result = sim.get_eigenmode_coefficients(flux, [1], eig_parity=mp.ODD_Z)
    return abs(result.alpha[0, 0, 0 if direction > 0 else 1]) ** 2


# ------------------------------------------------------------ 1. the field


def check_field(resolution, tilt=10.0, below=3.5):
    """Meep's launched field against the propagator's prediction."""
    propagator, fields = injection(resolution, fiber_mode(tilt))
    sim = simulation(resolution, sources=fiber_sources(propagator, fields))
    line = mp.Volume(
        center=PLANE.center - mp.Vector3(0, below), size=mp.Vector3(PLANE.size.x, 0)
    )
    region = mp.Volume(center=mp.Vector3(), size=CELL - mp.Vector3(2 * PML, 2 * PML))
    cut = sim.add_dft_fields([mp.Ez], [FCEN], where=line, yee_grid=False)
    whole = sim.add_dft_fields([mp.Ez], [FCEN], where=region, yee_grid=False)
    sim.run(until_after_sources=DECAY)

    x = np.asarray(sim.get_array_metadata(vol=line)[0])
    measured = np.asarray(sim.get_dft_array(cut, mp.Ez, 0))
    # Meep scales every source by the DTFT of its time profile; that is the
    # only factor applied to the prediction.
    spectrum = complex(source_gradient.time_profile_dtft(sim, SRC_TIME, [FCEN])[0])
    down = homogeneous(propagator, sign=-1)
    outgoing = mpa.TangentialFields(E=fields.E, H=fields.H, normal=mp.Y, sign=-1)
    predicted = spectrum * np.asarray(
        down.propagate(outgoing, below, coordinates=x)[mp.Ez][0]
    )
    meta = sim.get_array_metadata(vol=region)
    return dict(
        x=x,
        measured=measured,
        predicted=predicted,
        map=np.abs(np.asarray(sim.get_dft_array(whole, mp.Ez, 0))),
        map_x=np.asarray(meta[0]),
        map_y=np.asarray(meta[1]),
        line_y=PLANE.center.y - below,
    )


# ------------------------------------------------------- 2. reciprocity


def transmit_efficiency(resolution, tilts):
    """Waveguide -> fiber, for every tilt from one simulation."""
    sim = simulation(
        resolution,
        device(),
        [
            mp.EigenModeSource(
                SRC_TIME,
                center=mp.Vector3(SOURCE_X, WG_Y),
                size=mp.Vector3(0, 2.5),
                eig_band=1,
                eig_parity=mp.ODD_Z,
                direction=mp.X,
            )
        ],
    )
    flux_in = sim.add_flux(FCEN, 0, 1, mp.FluxRegion(volume=WG_LINE))
    monitor = sim.add_dft_fields([mp.Ez, mp.Hx], [FCEN], where=PLANE, yee_grid=False)
    flux_up = sim.add_flux(FCEN, 0, 1, mp.FluxRegion(volume=PLANE))
    sim.run(until_after_sources=DECAY)

    power_in = waveguide_power(sim, flux_in, +1)
    propagator = mpa.AngularSpectrum.from_monitor(
        sim, monitor, STACK, PLANE, pad_factor=8
    )
    fields = propagator.fields_from_monitor(sim, monitor)
    # Meep's flux per unit of the propagator's power, from the same plane
    units = mp.get_fluxes(flux_up)[0] / float(
        homogeneous(propagator, sign=+1).power(fields, 0.0)[0]
    )
    return np.array(
        [
            units
            * float(
                propagator.overlap(
                    fields, fiber_mode(t), FIBER_DISTANCE, incident_power=jnp.ones(1)
                )[0]
            )
            / power_in
            for t in tilts
        ]
    )


def receive_efficiency(resolution, tilt):
    """Fiber -> waveguide: coupled power over the power the fiber emits."""
    mode = fiber_mode(tilt)
    propagator, fields = injection(resolution, mode)

    # the power launched into the cell, measured by Meep without the device
    reference = simulation(resolution, sources=fiber_sources(propagator, fields))
    below = mp.Volume(
        center=PLANE.center - mp.Vector3(0, 0.5),
        size=mp.Vector3(CELL.x - 2 * PML, 0),
    )
    launched = reference.add_flux(FCEN, 0, 1, mp.FluxRegion(volume=below, weight=-1))
    reference.run(until_after_sources=DECAY)

    sim = simulation(resolution, device(), fiber_sources(propagator, fields))
    flux = sim.add_flux(FCEN, 0, 1, mp.FluxRegion(volume=WG_LINE))
    sim.run(until_after_sources=DECAY)

    # The fiber emits more than reaches the cell: the stack reflects some.
    # The ratio comes from the propagator, in its own units.
    k0 = 2 * np.pi * np.array([FCEN])
    a = np.squeeze(
        np.asarray(
            mode.spectrum(
                propagator, propagator.kt, k0[:, None], propagator._indices[-1]
            )
        )
    )
    admittance = [
        np.real(np.asarray(propagator._admittances[-1][p])[0])
        for p in (asm.S_POLARIZATION, asm.P_POLARIZATION)
    ]
    emitted = (
        0.5
        * propagator._spectral_measure()
        * np.sum(
            admittance[0] * abs(a[..., 0]) ** 2 + admittance[1] * abs(a[..., 1]) ** 2
        )
    )
    outgoing = mpa.TangentialFields(E=fields.E, H=fields.H, normal=mp.Y, sign=-1)
    arriving = float(homogeneous(propagator, sign=-1).power(outgoing, 0.0)[0])

    return waveguide_power(sim, flux, -1) / (
        mp.get_fluxes(launched)[0] * emitted / arriving
    )


# ------------------------------------------------------------ 3. gradients


def coupled_power_loss(resolution, run_time=400.0):
    """Waveguide power as a JAX function of (waist, tilt, offset)."""
    sim = simulation(resolution, device())
    inert = mp.MaterialGrid(
        mp.Vector3(4, 4), mp.Medium(index=N_OXIDE), mp.Medium(index=N_OXIDE)
    )
    # MeepJaxWrapper needs a design region; this one changes nothing
    region = mpa.DesignRegion(
        inert, volume=mp.Volume(center=mp.Vector3(8, -2.5), size=mp.Vector3(0.4, 0.4))
    )
    sim.geometry = sim.geometry + [
        mp.Block(center=region.center, size=region.size, material=inert)
    ]
    propagator, seed = injection(resolution, fiber_mode(10.0))
    sources = fiber_sources(propagator, seed)
    monitor = mpa.EigenmodeCoefficient(
        sim, WG_LINE, mode=1, forward=False, eig_parity=mp.ODD_Z
    )
    wrapper = mpa.MeepJaxWrapper(
        sim,
        sources,
        [monitor],
        [region],
        [FCEN],
        minimum_run_time=run_time,
        maximum_run_time=run_time,
    )

    def loss(fiber):
        mode = fiber_mode(fiber[1], offset=fiber[2], waist=fiber[0])
        fields = propagator.incident_fields(mode, distance=FIBER_DISTANCE)
        sheet = propagator.equivalent_currents(fields)
        values = wrapper(
            [jnp.full((4, 4, 1), 0.5)], [jnp.ravel(sheet[s.component]) for s in sources]
        )
        return jnp.sum(jnp.abs(jnp.asarray(values)) ** 2)

    return loss


def check_gradient(resolution):
    x0 = jnp.array([WAIST, 10.0, 0.3])  # waist (um), tilt (deg), offset (um)
    _, adjoint = jax.value_and_grad(coupled_power_loss(resolution))(x0)
    reference = []
    for i, step in enumerate((0.02, 0.1, 0.02)):
        e = jnp.zeros(3).at[i].set(step)
        plus = float(coupled_power_loss(resolution)(x0 + e))
        minus = float(coupled_power_loss(resolution)(x0 - e))
        reference.append((plus - minus) / (2 * step))
    return np.asarray(adjoint), np.array(reference)


# ------------------------------------------------------------------ figure

INK, INK_2, GRID, SURFACE = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"
BLUE, ORANGE = "#2a78d6", "#eb6834"


def figure(field, tilts, eta_tx, eta_rx, resolutions, gaps, adjoint, reference, path):
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap

    plt.rcParams.update(
        {
            "font.size": 9,
            "axes.edgecolor": INK_2,
            "axes.labelcolor": INK,
            "axes.titlecolor": INK,
            "xtick.color": INK_2,
            "ytick.color": INK_2,
            "axes.grid": True,
            "grid.color": GRID,
            "grid.linewidth": 0.6,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "figure.facecolor": SURFACE,
            "axes.facecolor": SURFACE,
            "legend.frameon": False,
        }
    )
    fig = plt.figure(figsize=(12, 7))
    grid = fig.add_gridspec(2, 3, hspace=0.45, wspace=0.32)

    # (a) the launched beam
    ax = fig.add_subplot(grid[0, 0])
    blues = LinearSegmentedColormap.from_list(
        "blues", [SURFACE, "#9cc1ee", BLUE, "#123a6b"]
    )
    ax.pcolormesh(
        field["map_x"], field["map_y"], field["map"].T, cmap=blues, shading="auto"
    )
    ax.axhline(PLANE.center.y, color=INK_2, lw=1, ls="--")
    ax.axhline(field["line_y"], color=ORANGE, lw=1.5)
    for y, label in (
        (PLANE.center.y, "injection plane"),
        (field["line_y"], "comparison line"),
    ):
        ax.text(6.8, y + 0.2, label, color=INK_2, fontsize=8, ha="right")
    ax.set(xlabel="x (µm)", ylabel="y (µm)", xlim=(-8, 7), aspect="equal")
    ax.grid(False)
    ax.set_title("(a) Launched beam |Ez|, 10° fiber tilt", loc="left")

    # (b) launched field against the propagator
    ax = fig.add_subplot(grid[0, 1])
    ax.plot(field["x"], field["measured"].real, color=BLUE, lw=2, label="Meep (FDTD)")
    ax.plot(
        field["x"][:: max(1, len(field["x"]) // 25)],
        field["predicted"].real[:: max(1, len(field["x"]) // 25)],
        "o",
        color=ORANGE,
        ms=4,
        label="angular spectrum",
    )
    error = np.linalg.norm(field["measured"] - field["predicted"]) / np.linalg.norm(
        field["measured"]
    )
    ax.text(
        0.98,
        0.06,
        f"relative difference\n{error:.1%}",
        transform=ax.transAxes,
        color=INK_2,
        ha="right",
    )
    ax.set(xlabel="x (µm)", ylabel="Re Ez (a.u.)")
    ax.legend(loc="upper right")
    ax.set_title("(b) Field below the plane, no fitted scale", loc="left")

    # (c) gradients
    ax = fig.add_subplot(grid[0, 2])
    names = ["waist", "tilt", "offset"]
    position = np.arange(3)
    scale = np.abs(reference)
    ax.bar(
        position - 0.18, adjoint / scale, 0.34, color=BLUE, label="adjoint (JAX + Meep)"
    )
    ax.bar(
        position + 0.18,
        reference / scale,
        0.34,
        color=ORANGE,
        label="finite difference",
    )
    for i in range(3):
        rel = abs(adjoint[i] - reference[i]) / abs(reference[i])
        y = np.sign(reference[i]) * 1.08
        ax.text(
            i,
            y,
            f"{rel:.0e}",
            ha="center",
            va="bottom" if y > 0 else "top",
            color=INK_2,
            fontsize=8,
        )
    ax.axhline(0, color=INK_2, lw=0.8)
    ax.set_xticks(position, names)
    ax.set_ylim(-1.45, 1.9)
    ax.set_ylabel("dP / dparameter, normalized")
    ax.legend(loc="upper left", fontsize=8)
    ax.set_title("(c) Gradient w.r.t. fiber parameters", loc="left")

    # (d) reciprocity
    ax = fig.add_subplot(grid[1, 0:2])
    ax.plot(tilts, 100 * eta_tx, color=BLUE, lw=2, label="waveguide → fiber (transmit)")
    ax.plot(
        tilts,
        100 * eta_rx,
        "o",
        color=ORANGE,
        ms=6,
        label="fiber → waveguide (receive)",
    )
    ax.set(xlabel="fiber tilt (degrees)", ylabel="coupling efficiency (%)")
    ax.set_ylim(0, 1.35 * 100 * max(np.max(eta_tx), np.max(eta_rx)))
    ax.legend(loc="upper center", ncol=2)
    ax.set_title(
        f"(d) Reciprocity: both directions, resolution {resolutions[-1]}", loc="left"
    )

    # (e) convergence of the remaining gap
    ax = fig.add_subplot(grid[1, 2])
    ax.loglog(resolutions, 100 * np.asarray(gaps), "o-", color=BLUE, lw=2, ms=6)
    guide = 100 * gaps[0] * (np.asarray(resolutions) / resolutions[0]) ** -2.0
    ax.loglog(resolutions, guide, "--", color=INK_2, lw=1)
    ax.text(resolutions[1], guide[1] * 0.7, "second order", color=INK_2, fontsize=8)
    ax.set_xticks(resolutions, [str(r) for r in resolutions])
    ticks = 100 * np.asarray(gaps)
    ax.set_yticks(ticks, [f"{t:.0f}%" if t >= 10 else f"{t:.1f}%" for t in ticks])
    ax.minorticks_off()
    ax.set(xlabel="resolution (pixels/µm)", ylabel="transmit vs receive gap")
    ax.set_title("(e) The gap is discretization error", loc="left")

    fig.savefig(path, dpi=200, bbox_inches="tight")
    print(f"wrote {path}")


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--quick", action="store_true", help="coarser grids")
    parser.add_argument("--output", default="grating_coupler_asm_end_to_end.png")
    parser.add_argument(
        "--replot",
        action="store_true",
        help="redraw the figure from the results saved by a previous run",
    )
    args = parser.parse_args()
    saved = args.output.rsplit(".", 1)[0] + ".npz"

    if args.replot:
        r = dict(np.load(saved, allow_pickle=True))
        field = r.pop("field").item()
        figure(field, **r, path=args.output)
        return

    resolutions = [10, 15, 20] if args.quick else [20, 30, 40]
    tilts = np.arange(4.0, 18.1, 2.0)

    field = check_field(resolutions[-1])
    error = np.linalg.norm(field["measured"] - field["predicted"])
    print(f"field: relative difference {error / np.linalg.norm(field['measured']):.2%}")

    eta_tx = transmit_efficiency(resolutions[-1], tilts)
    eta_rx = np.array([receive_efficiency(resolutions[-1], t) for t in tilts])
    for t, a, b in zip(tilts, eta_tx, eta_rx):
        print(f"tilt {t:4.1f}  transmit {a:.4%}  receive {b:.4%}")

    gaps = []
    for r in resolutions:
        if r == resolutions[-1]:
            i = int(np.argmin(np.abs(tilts - 10.0)))
            tx, rx = eta_tx[i], eta_rx[i]
        else:
            tx, rx = transmit_efficiency(r, [10.0])[0], receive_efficiency(r, 10.0)
        gaps.append(abs(tx - rx) / max(tx, rx))
        print(f"resolution {r}: transmit/receive gap {gaps[-1]:.1%}")

    adjoint, reference = check_gradient(resolutions[0])
    for name, a, b in zip(("waist", "tilt", "offset"), adjoint, reference):
        print(f"d/d{name}: adjoint {a:+.6g}  finite difference {b:+.6g}")

    if mp.am_master():
        results = dict(
            tilts=tilts,
            eta_tx=eta_tx,
            eta_rx=eta_rx,
            resolutions=np.array(resolutions),
            gaps=np.array(gaps),
            adjoint=adjoint,
            reference=reference,
        )
        np.savez(saved, field=np.array(field, dtype=object), **results)
        figure(field, **results, path=args.output)


if __name__ == "__main__":
    main()
