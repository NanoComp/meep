import warnings
import functools

import numpy as np

from meep.geom import Vector3, check_nonnegative

import meep as mp


def check_positive(prop, val):
    if val > 0:
        return val
    else:
        raise ValueError(f"{prop} must be positive. Got {val}")


class Source:
    """
    The `Source` class is used to specify the current sources via the `Simulation.sources`
    attribute. Note that all sources in Meep are separable in time and space, i.e. of the
    form $\\mathbf{J}(\\mathbf{x},t) = \\mathbf{A}(\\mathbf{x}) \\cdot f(t)$ for some functions
    $\\mathbf{A}$ and $f$. Non-separable sources can be simulated, however, by modifying
    the sources after each time step. When real fields are being used (which is the
    default in many cases; see `Simulation.force_complex_fields`), only the real part of
    the current source is used.

    **Important note**: These are *current* sources (**J** terms in Maxwell's equations),
    even though they are labelled by electric/magnetic field components. They do *not*
    specify a particular electric/magnetic field which would be what is called a "hard"
    source in the FDTD literature. There is no fixed relationship between the current
    source and the resulting field amplitudes; it depends on the surrounding geometry, as
    described in the
    [FAQ](FAQ.md#how-does-the-current-amplitude-relate-to-the-resulting-field-amplitude)
    and in Section 4.4 ("Currents and Fields: The Local Density of States") in [Chapter
    4](http://arxiv.org/abs/arXiv:1301.5366) ("Electromagnetic Wave Source Conditions") of
    the book [Advances in FDTD Computational Electrodynamics: Photonics and
    Nanotechnology](https://www.amazon.com/Advances-FDTD-Computational-Electrodynamics-Nanotechnology/dp/1608071707).
    """

    def __init__(
        self,
        src,
        component,
        center=None,
        volume=None,
        size=Vector3(),
        amplitude=1.0,
        amp_func=None,
        amp_func_file="",
        amp_data=None,
    ):
        """
        Construct a `Source`.

        + **`src` [`SourceTime` class ]** — Specify the time-dependence of the source (see
          below). No default.

        + **`component` [`component` constant ]** — Specify the direction and type of the
          current component: e.g. `mp.Ex`, `mp.Ey`, etcetera for an electric-charge
          current, and `mp.Hx`, `mp.Hy`, etcetera for a magnetic-charge current. Note that
          currents pointing in an arbitrary direction are specified simply as multiple
          current sources with the appropriate amplitudes for each component. No default.

        + **`center` [`Vector3`]** — The location of the center of the current source in
          the cell. No default.

        + **`size` [`Vector3`]** — The size of the current distribution along each
          direction of the cell. Default is `(0,0,0)`: a point-dipole source.

        + **`volume` [`Volume`]** — A `meep.Volume` can be used to specify the source
          region instead of a `center` and a `size`.

        + **`amplitude` [`complex`]** — An overall complex amplitude multiplying the
          current source. Default is 1.0. Note that specifying a complex `amplitude`
          imparts a phase shift to the real part of the overall current and thus
          does *not* require using complex fields for the entire simulation
          (via `force_complex_fields=True`).

        + **`amp_func` [`function`]** — A Python function of a single argument, that takes
          a `Vector3` giving a position and returns a complex current amplitude for that
          point. The position argument is *relative* to the `center` of the current
          source, so that you can move your current around without changing your function.
          Default is `None`, meaning that a constant amplitude of 1.0 is used. Note that
          your amplitude function (if any) is *multiplied* by the `amplitude` property, so
          both properties can be used simultaneously.

        + **`amp_func_file` [`string`]** — String of the form
          `path_to_h5_file.h5:dataset`. The `.h5` extension is optional. Meep will read
          the HDF5 file and create an amplitude function that interpolates into the grid
          specified by the file. Meep expects the data to be split into real and imaginary
          parts, so in the above example it will look for `dataset.re` and `dataset.im` in
          the file `path_to_h5_file.h5`. Defaults to the empty string.

        + **`amp_data` [`numpy.ndarray with dtype=numpy.complex128`]** — Like
          `amp_func_file` above, but instead of interpolating into an HDF5 file,
          interpolates into a complex NumPy array. The array should be three dimensions.
          For a 2d simulation, just pass 1 for the third dimension, e.g., `arr =
          np.zeros((N, M, 1), dtype=np.complex128)`. Defaults to `None`.

        As described in Section 4.2 ("Incident Fields and Equivalent Currents") in
        [Chapter 4](http://arxiv.org/abs/arXiv:1301.5366) ("Electromagnetic Wave Source
        Conditions") of the book [Advances in FDTD Computational Electrodynamics:
        Photonics and Nanotechnology](https://www.amazon.com/Advances-FDTD-Computational-Electrodynamics-Nanotechnology/dp/1608071707),
        it is also possible to supply a source that is designed to couple exclusively into
        a single waveguide mode (or other mode of some cross section or periodic region)
        at a single frequency, and which couples primarily into that mode as long as the
        bandwidth is not too broad. This is possible if you have
        [MPB](https://mpb.readthedocs.io) installed: Meep will call MPB to compute the
        field profile of the desired mode, and use the field profile to produce an
        equivalent current source. The mode-launcher feature does *not* work in
        cylindrical coordinates. To use the mode launcher, instead of a `source` you
        should use an `EigenModeSource`.
        """
        if center is None and volume is None:
            raise ValueError("Source requires either center or volume")

        if volume:
            self.center = volume.center
            self.size = volume.size
        else:
            self.center = Vector3(*center)
            self.size = Vector3(*size)

        self.src = src
        self.component = component
        self.amplitude = complex(amplitude)
        self.amp_func = amp_func
        self.amp_func_file = amp_func_file
        self.amp_data = amp_data

    def add_source(self, sim):
        where = mp.Volume(
            self.center,
            self.size,
            dims=sim.dimensions,
            is_cylindrical=sim.is_cylindrical,
        ).swigobj
        add_vol_src_args = [self.component, self.src.swigobj, where]
        add_vol_src = functools.partial(sim.fields.add_volume_source, *add_vol_src_args)

        if self.amp_func_file:
            fname_dset = self.amp_func_file.rsplit(":", 1)
            if len(fname_dset) != 2:
                err_msg = "Expected a string of the form 'h5filename:dataset'. Got '{}'"
                raise ValueError(err_msg.format(self.amp_func_file))

            fname, dset = fname_dset
            if not fname.endswith(".h5"):
                fname += ".h5"

            add_vol_src(fname, dset, self.amplitude * 1.0)
        elif self.amp_func:
            add_vol_src(self.amp_func, self.amplitude * 1.0)
        elif self.amp_data is not None:
            add_vol_src(self.amp_data, self.amplitude * 1.0)
        else:
            add_vol_src(self.amplitude * 1.0)


class SourceTime:
    """
    This is the parent for classes describing the time dependence of sources; it should
    not be instantiated directly.
    """

    def __init__(self, is_integrated=False):
        self.is_integrated = is_integrated


class ContinuousSource(SourceTime):
    """
    A continuous-wave (CW) source is proportional to $\\exp(-i\\omega t)$, possibly with a
    smooth (exponential/tanh) turn-on/turn-off. In practice, the CW source [never produces
    an exact single-frequency
    response](FAQ.md#why-doesnt-the-continuous-wave-cw-source-produce-an-exact-single-frequency-response).
    """

    def __init__(
        self,
        frequency=None,
        start_time=0,
        end_time=1.0e20,
        width=0,
        fwidth=float("inf"),
        cutoff=3.0,
        wavelength=None,
        is_integrated=False,
        **kwargs,
    ):
        """
        Construct a `ContinuousSource`.

        + **`frequency` [`number`]** — The frequency *f* in units of $c$/distance or ω in
          units of 2π$c$/distance. See [Units](Introduction.md#units-in-meep). No default
          value. You can instead specify `wavelength=x` or `period=x`, which are both a
          synonym for `frequency=1/x`; i.e. 1/ω in these units is the vacuum wavelength or
          the temporal period.

        + **`start_time` [`number`]** — The starting time for the source. Default is 0
          (turn on at $t=0$).

        + **`end_time` [`number`]** — The end time for the source. Default is
          10<sup>20</sup> (never turn off).

        + **`width` [`number`]** — Roughly, the temporal width of the smoothing
          (technically, the inverse of the exponential rate at which the current turns off
          and on). Default is 0 (no smoothing). You can instead specify `fwidth=x`, which
          is a synonym for `width=1/x` (i.e. the frequency width is proportional to the
          inverse of the temporal width).

        + **`slowness` [`number`]** — Controls how far into the exponential tail of the
          tanh function the source turns on. Default is 3.0. A larger value means that the
          source turns on more gradually at the beginning.

        + **`is_integrated` [`boolean`]** — If `True`, the source is the integral of the
          current (the [dipole
          moment](https://en.wikipedia.org/wiki/Electric_dipole_moment)) which oscillates
          but does not increase for a sinusoidal current. In practice, there is little
          difference between integrated and non-integrated sources *except* for
          [planewaves extending into
          PML](Perfectly_Matched_Layer.md#planewave-sources-extending-into-pml). Default
          is `False`.
        """

        if frequency is None and wavelength is None:
            raise ValueError(
                f"Must set either frequency or wavelength in {self.__class__.__name__}."
            )

        super().__init__(is_integrated=is_integrated, **kwargs)
        self.frequency = 1 / wavelength if wavelength else float(frequency)
        self.start_time = start_time
        self.end_time = end_time
        self.width = max(width, 1 / fwidth)
        self.cutoff = cutoff
        self.swigobj = mp.continuous_src_time(
            self.frequency, self.width, self.start_time, self.end_time, self.cutoff
        )
        self.swigobj.is_integrated = self.is_integrated


class GaussianSource(SourceTime):
    """
    A Gaussian-pulse source roughly proportional to $\\exp(-i\\omega t - (t-t_0)^2/2w^2)$.
    Technically, the "Gaussian" sources in Meep are the (discrete-time) derivative of a
    Gaussian, i.e. they are $(-i\\omega)^{-1} \\frac{\\partial}{\\partial t} \\exp(-i\\omega t -
    (t-t_0)^2/2w^2)$, but the difference between this and a true Gaussian is usually
    irrelevant.
    """

    def __init__(
        self,
        frequency=None,
        width=0,
        fwidth=float("inf"),
        start_time=0,
        cutoff=5.0,
        is_integrated=False,
        wavelength=None,
        **kwargs,
    ):
        """
        Construct a `GaussianSource`.

        + **`frequency` [`number`]** — The center frequency $f$ in units of $c$/distance
          (or $\\omega$ in units of $2\\pi c$/distance). See [Units](Introduction.md#units-in-meep).
          No default value. You can instead specify `wavelength=x` or `period=x`, which
          are both a synonym for `frequency=1/x`; i.e. $1/\\omega$ in these units is the vacuum
          wavelength or the temporal period.

        + **`width` [`number`]** — The width $w$ used in the Gaussian. No default value.
          You can instead specify `fwidth=x`, which is a synonym for `width=1/x` (i.e. the
          frequency width is proportional to the inverse of the temporal width).

        + **`start_time` [`number`]** — The starting time for the source; default is 0
          (turn on at $t=0$). This is not the time of the peak. See below.

        + **`cutoff` [`number`]** — How many `width`s the current decays for before it is
          cut off and set to zero &mdash; this applies for both turn-on and turn-off of
          the pulse. Default is 5.0. A larger value of `cutoff` will reduce the amount of
          high-frequency components that are introduced by the start/stop of the source,
          but will of course lead to longer simulation times. The peak of the Gaussian is
          reached at the time $t_0$=`start_time + cutoff*width`.

        + **`is_integrated` [`boolean`]** — If `True`, the source is the integral of the
          current (the [dipole moment](https://en.wikipedia.org/wiki/Electric_dipole_moment))
          which is guaranteed to be zero after the current turns off. In practice, there
          is little difference between integrated and non-integrated sources *except* for
          [planewaves extending into PML](Perfectly_Matched_Layer.md#planewave-sources-extending-into-pml).
          Default is `False`.

        + **`fourier_transform(f)`** — Returns the Fourier transform of the current
          evaluated at frequency $f$ ($\\omega=2\\pi f$) given by:
          $$
          \\widetilde G(\\omega) \\equiv \\frac{1}{\\sqrt{2\\pi}}
          \\int e^{i\\omega t}G(t)\\,dt \\equiv
          \\frac{1}{\\Delta f}
          e^{i\\omega t_0 -\\frac{(\\omega-\\omega_0)^2}{2\\Delta f^2}}
          $$
          where $G(t)$ is the current (not the dipole moment). In this formula, $\\Delta f$
          is the `fwidth` of the source, $\\omega_0$ is $2\\pi$ times its `frequency,` and
          $t_0$ is the peak time discussed above. Note that this does not include any
          `amplitude` or `amp_func` factor that you specified for the source.
        """
        if frequency is None and wavelength is None:
            raise ValueError(
                f"Must set either frequency or wavelength in {self.__class__.__name__}."
            )

        super().__init__(is_integrated=is_integrated, **kwargs)
        self.frequency = 1 / wavelength if wavelength else float(frequency)
        self.width = max(width, 1 / fwidth)
        self.start_time = start_time
        self.cutoff = cutoff

        self.swigobj = mp.gaussian_src_time(
            self.frequency,
            self.width,
            self.start_time,
            self.start_time + 2 * self.width * self.cutoff,
        )
        self.swigobj.is_integrated = self.is_integrated

    def fourier_transform(self, freq):
        return self.swigobj.fourier_transform(freq)


class CustomSource(SourceTime):
    """
    A user-specified source function $f(t)$. You can also specify start/end times at which
    point your current is set to zero whether or not your function is actually zero. These
    are optional, but you must specify an `end_time` explicitly if you want `run`
    functions like `until_after_sources` to work, since they need to know when your source
    turns off. To use a custom source within an `EigenModeSource`, you must specify the
    `center_frequency` parameter, since Meep does not know the frequency content of the
    `CustomSource`. The resultant eigenmode is calculated at this frequency only. For a
    demonstration of a [linear-chirped pulse](FAQ.md#how-do-i-create-a-chirped-pulse), see
    [`examples/chirped_pulse.py`](https://github.com/NanoComp/meep/blob/master/python/examples/chirped_pulse.py).
    """

    def __init__(
        self,
        src_func,
        start_time=-1.0e20,
        end_time=1.0e20,
        is_integrated=False,
        center_frequency=0,
        fwidth=0,
        **kwargs,
    ):
        """
        Construct a `CustomSource`.

        + **`src_func` [`function`]** — The function $f(t)$ specifying the time-dependence
          of the source. It should take one argument (the time in Meep units) and return a
          complex number.

        + **`start_time` [`number`]** — The starting time for the source. Default is
          -10<sup>20</sup>: turn on at $t=-\\infty$. Note, however, that the simulation
          normally starts at $t=0$ with zero fields as the initial condition, so there is
          implicitly a sharp turn-on at $t=0$ whether you specify it or not.

        + **`end_time` [`number`]** — The end time for the source. Default is
          10<sup>20</sup> (never turn off).

        + **`is_integrated` [`boolean`]** — If `True`, the source is the integral of the
          current (the [dipole
          moment](https://en.wikipedia.org/wiki/Electric_dipole_moment)) which is
          guaranteed to be zero after the current turns off. In practice, there is little
          difference between integrated and non-integrated sources *except* for
          [planewaves extending into
          PML](Perfectly_Matched_Layer.md#planewave-sources-extending-into-pml). Default
          is `False`.

        + **`center_frequency` [`number`]** — Optional center frequency so that the
          `CustomSource` can be used within an `EigenModeSource`. Defaults to 0.

        + **`fwidth` [`number`]** — Optional bandwidth in frequency units.
          Default is 0. For bandwidth-limited sources, this parameter is used to
          automatically determine the decimation factor of the time-series updates
          of the DFT fields monitors (if any).
        """
        super().__init__(is_integrated=is_integrated, **kwargs)
        self.src_func = src_func
        self.start_time = start_time
        self.end_time = end_time
        self.fwidth = fwidth
        self.center_frequency = center_frequency
        self.swigobj = mp.custom_py_src_time(
            src_func, start_time, end_time, center_frequency, fwidth
        )
        self.swigobj.is_integrated = self.is_integrated


class EigenModeSource(Source):
    """
    This is a subclass of `Source` and has **all of the properties** of `Source` above.
    However, you normally do not specify a `component`. Instead of `component`, the
    current source components and amplitude profile are computed by calling MPB to compute
    the modes, $\\mathbf{u}_{n,\\mathbf{k}}(\\mathbf{r}) e^{i \\mathbf{k} \\cdot \\mathbf{r}}$,
    of the dielectric profile in the region given by the `size` and `center` of the
    source, with the modes computed as if the *source region were repeated periodically in
    all directions*. If an `amplitude` and/or `amp_func` are supplied, they are
    *multiplied* by this current profile. The desired eigenmode and other features are
    specified by the properties shown in `__init__`.

    Eigenmode sources are normalized so that in the case of a time-harmonic simulation
    with all sources and fields having monochromatic time dependence $e^{-i 2\\pi f_m t}$
    where $f_m$ is the frequency of the eigenmode, the total time-average power of the
    fields — the integral of the normal Poynting vector over the entire cross-sectional
    line or plane — is equal to 1. This convention has two use cases:

    + For [frequency-domain
      calculations](Python_User_Interface.md#frequency-domain-solver) involving a
      `ContinuousSource` time dependence, the time-average power of the fields is 1.

    + For time-domain calculations involving a time dependence $W(t)$ which is typically a
      [Gaussian](#gaussiansource), the amplitude of the fields at frequency $f$ will be
      multiplied by $\\widetilde W(f)$, the Fourier transform of $W(t)$, while
      field-bilinear quantities like the [Poynting flux](#flux-spectra) and [energy
      density](#energy-density-spectra) are multiplied by $|\\widetilde W(f)|^2$. For the
      particular case of a Gaussian time dependence, the Fourier transform at $f$ can be
      obtained via the `fourier_transform` class method.

    In either case, the `eig_power` method returns the total power at frequency `f`.
    However, for a user-defined [`CustomSource`](#customsource), `eig_power` will *not*
    include the $|\\widetilde W(f)|^2$ factor since Meep does not know the Fourier
    transform of your source function $W(t)$. You will have to multiply by this yourself
    if you need it.

    **Note:** Due to discretization effects, the normalization of eigenmode sources to
    yield unit power transmission is only approximate: at any finite resolution, the power
    of the fields as measured using [DFT flux](#flux-spectra) monitors will not precisely
    match that of calling `eig_power` but will rather include discretization errors that
    decrease with resolution.  Generally, the most reliable procedure is to normalize your
    calculations by the power computed in a separate normalization run at the same
    resolution, as shown in several of the tutorial examples.

    Note that Meep's MPB interface only supports dispersionless non-magnetic materials but
    it does support anisotropic $\\varepsilon$. Any nonlinearities, magnetic responses $\\mu$,
    conductivities $\\sigma$, or dispersive polarizations in your materials will be *ignored* when
    computing the eigenmode source. PML will also be ignored.

    The `SourceTime` object (`Source.src`), which specifies the time dependence of the
    source, can be one of `ContinuousSource`, `GaussianSource` or `CustomSource`.
    """

    def __init__(
        self,
        src,
        center=None,
        volume=None,
        eig_lattice_size=None,
        eig_lattice_center=None,
        component=mp.ALL_COMPONENTS,
        direction=mp.AUTOMATIC,
        eig_band=1,
        eig_kpoint=Vector3(),
        eig_match_freq=True,
        eig_parity=mp.NO_PARITY,
        eig_resolution=0,
        eig_tolerance=1e-12,
        **kwargs,
    ):
        """
        Construct an `EigenModeSource`.

        + **`eig_band` [`integer` or `DiffractedPlanewave` class]** — Either the index $n$
          (1,2,3,...) of the desired band $\\omega_n(\\mathbf{k})$ to compute in MPB where
          1 denotes the lowest-frequency band at a given $\\mathbf{k}$ point, and so on,
          or alternatively a diffracted planewave in homogeneous media.

        + **`direction` [`mp.X`, `mp.Y`, or `mp.Z;` default `mp.AUTOMATIC`],
          `eig_match_freq` [`boolean;` default `True`], `eig_kpoint` [`Vector3`]** — By
          default (if `eig_match_freq` is `True`), Meep tries to find a mode with the same
          frequency $\\omega_n(\\mathbf{k})$ as the `src` property (above), by scanning
          $\\mathbf{k}$ vectors in the given `direction` using MPB's `find_k` functionality.
          Alternatively, if `eig_kpoint` is supplied, it is used as an initial guess for
          $\\mathbf{k}$. By default, `direction` is the direction normal to the source region,
          assuming `size` is $d$–1 dimensional in a $d$-dimensional simulation (e.g. a
          plane in 3d). If `direction` is set to `mp.NO_DIRECTION`, then `eig_kpoint` is
          not only the initial guess and the search direction of the $\\mathbf{k}$ vectors, but is
          also taken to be the direction of the waveguide, allowing you to [launch modes
          in oblique ridge waveguides](Python_Tutorials/Eigenmode_Source.md#oblique-waveguides)
          (not perpendicular to the source plane).  If `eig_match_freq` is `False`, then the
          $\\mathbf{k}$ vector of the desired mode is specified with  `eig_kpoint` (in Meep units
          of 2π/(unit length)). Also, the eigenmode frequency computed by MPB overwrites
          the `frequency` parameter of the `src` property for a `GaussianSource` and
          `ContinuousSource` but not `CustomSource` (the `width` or any other parameter of
          `src` is unchanged). By default, the $\\mathbf{k}$ components in the plane of the source
          region are zero.  However, if the source region spans the *entire* cell in some
          directions, and the cell has Bloch-periodic boundary conditions via the
          `k_point` parameter, then the mode's $\\mathbf{k}$ components in those directions will
          match `k_point` so that the mode satisfies the Meep boundary conditions,
          regardless of `eig_kpoint`. Note that once $\\mathbf{k}$ is either found by MPB, or
          specified by `eig_kpoint`, the field profile used to create the current sources
          corresponds to the [Bloch mode](https://en.wikipedia.org/wiki/Bloch_wave),
          $\\mathbf{u}_{n,\\mathbf{k}}(\\mathbf{r})$, multiplied by the appropriate
          exponential factor, $e^{i \\mathbf{k} \\cdot \\mathbf{r}}$.

        + **`eig_parity` [`mp.NO_PARITY` (default), `mp.EVEN_Z`, `mp.ODD_Z`, `mp.EVEN_Y`,
          `mp.ODD_Y`]** — The parity (= polarization in 2d) of the mode to calculate,
          assuming the structure has $z$ and/or $y$ mirror symmetry *in the source
          region*, with respect to the `center` of the source region.  (In particular, it
          does not matter if your simulation as a whole has that symmetry, only the cross
          section where you are introducing the source.) If the structure has both $y$ and
          $z$ mirror symmetry, you can combine more than one of these, e.g. `EVEN_Z + ODD_Y`.
          Default is `NO_PARITY`, in which case MPB computes all of the bands
          which will still be even or odd if the structure has mirror symmetry, of course.
          This is especially useful in 2d simulations to restrict yourself to a desired
          polarization.

        + **`eig_resolution` [`integer`, defaults to `2*resolution` ]** — The spatial
          resolution to use in MPB for the eigenmode calculations. This defaults to twice
          the Meep `resolution` in which case the structure is linearly interpolated from
          the Meep pixels.

        + **`eig_tolerance` [`number`, defaults to 10<sup>–12</sup> ]** — The tolerance to
          use in the MPB eigensolver. MPB terminates when the eigenvalues stop changing to
          less than this fractional tolerance.

        + **`component` [as above, but defaults to `ALL_COMPONENTS`]** — Once the MPB
          modes are computed, equivalent electric and magnetic sources are created within
          Meep. By default, these sources include magnetic and electric currents in *all*
          transverse directions within the source region, corresponding to the mode fields
          as described in Section 4.2 ("Incident Fields and Equivalent Currents") in
          [Chapter 4](http://arxiv.org/abs/arXiv:1301.5366) ("Electromagnetic Wave Source
          Conditions") of the book [Advances in FDTD Computational Electrodynamics:
          Photonics and Nanotechnology](https://www.amazon.com/Advances-FDTD-Computational-Electrodynamics-Nanotechnology/dp/1608071707).
          If you specify a `component` property, however, you can include only one
          component of these currents if you wish. Most users won't need this feature.

        + **`eig_lattice_size` [`Vector3`], `eig_lattice_center` [`Vector3`]** — Normally,
          the MPB computational unit cell is the same as the source volume given by the
          `size` and `center` parameters. However, occasionally you want the unit cell to
          be larger than the source volume. For example, to create an eigenmode source in
          a periodic medium, you need to pass MPB the entire unit cell of the periodic
          medium, but once the mode is computed then the actual current sources need only
          lie on a cross section of that medium. To accomplish this, you can specify the
          optional `eig_lattice_size` and `eig_lattice_center`, which define a volume
          (which must enclose `size` and `center`) that is used for the unit cell in MPB
          with the dielectric function ε taken from the corresponding region in the Meep
          simulation.
        """

        super().__init__(src, component, center, volume, **kwargs)
        self.eig_lattice_size = eig_lattice_size
        self.eig_lattice_center = eig_lattice_center
        self.component = component
        self.direction = direction
        self.eig_band = eig_band
        self.eig_kpoint = mp.Vector3(*eig_kpoint)
        self.eig_match_freq = eig_match_freq
        self.eig_parity = eig_parity
        self.eig_resolution = eig_resolution
        self.eig_tolerance = eig_tolerance

    @property
    def eig_lattice_size(self):
        return self._eig_lattice_size

    @eig_lattice_size.setter
    def eig_lattice_size(self, val):
        self._eig_lattice_size = self.size if val is None else val

    @property
    def eig_lattice_center(self):
        return self._eig_lattice_center

    @eig_lattice_center.setter
    def eig_lattice_center(self, val):
        self._eig_lattice_center = self.center if val is None else val

    @property
    def component(self):
        return self._component

    @component.setter
    def component(self, val):
        if val != mp.ALL_COMPONENTS:
            warnings.warn(
                "EigenModeSource component is not ALL_COMPONENTS (the default), which makes it non-unidirectional.",
                RuntimeWarning,
            )
        self._component = val

    @property
    def eig_band(self):
        return self._eig_band

    @eig_band.setter
    def eig_band(self, val):
        if isinstance(val, int):
            self._eig_band = check_positive("EigenModeSource.eig_band", val)
        else:
            self._eig_band = val

    @property
    def eig_resolution(self):
        return self._eig_resolution

    @eig_resolution.setter
    def eig_resolution(self, val):
        self._eig_resolution = check_nonnegative("EigenModeSource.eig_resolution", val)

    @property
    def eig_tolerance(self):
        return self._eig_tolerance

    @eig_tolerance.setter
    def eig_tolerance(self, val):
        self._eig_tolerance = check_positive("EigenModeSource.eig_tolerance", val)

    def eig_power(self, freq):
        """
        Returns the total power of the fields from the eigenmode source at frequency `freq`.
        """
        amp = self.amplitude
        if callable(getattr(self.src, "fourier_transform", None)):
            amp *= self.src.fourier_transform(freq)
        return abs(amp) ** 2

    def add_source(self, sim):
        where = mp.Volume(
            self.center,
            self.size,
            dims=sim.dimensions,
            is_cylindrical=sim.is_cylindrical,
        ).swigobj
        if self.direction < 0:
            direction = sim.fields.normal_direction(where)
        else:
            direction = self.direction

        eig_vol = mp.Volume(
            self.eig_lattice_center,
            self.eig_lattice_size,
            sim.dimensions,
            is_cylindrical=sim.is_cylindrical,
        ).swigobj

        if isinstance(self.eig_band, mp.DiffractedPlanewave):
            eig_band = 1
            diffractedplanewave = mp.simulation.bands_to_diffractedplanewave(
                where, self.eig_band
            )
        elif isinstance(self.eig_band, int):
            eig_band = self.eig_band

        add_eig_src_args = [
            self.component,
            self.src.swigobj,
            direction,
            where,
            eig_vol,
            eig_band,
            mp.py_v3_to_vec(
                sim.dimensions, self.eig_kpoint, is_cylindrical=sim.is_cylindrical
            ),
            self.eig_match_freq,
            self.eig_parity,
            self.eig_resolution,
            self.eig_tolerance,
            self.amplitude,
        ]
        add_eig_src = functools.partial(
            sim.fields.add_eigenmode_source, *add_eig_src_args
        )

        if isinstance(self.eig_band, mp.DiffractedPlanewave):
            add_eig_src(self.amp_func, diffractedplanewave)
        else:
            add_eig_src(self.amp_func)


class GaussianBeam3DSource(Source):
    """
    This is a subclass of `Source` and has **all of the properties** of `Source` above. However, the `component` parameter of the `Source` object is ignored. The [Gaussian beam](https://en.wikipedia.org/wiki/Gaussian_beam) is a transverse electromagnetic mode for which the source region must be a *line* (in 2d) or *plane* (in 3d). For a beam polarized in the $x$ direction with propagation along $+z$, the electric field is defined by $\\mathbf{E}(r,z)=E_0\\hat{x}\\frac{w_0}{w(z)}\\exp\\left(\\frac{-r^2}{w(z)^2}\\right)\\exp\\left(-i\\left(kz + k\\frac{r^2}{2R(z)}\\right)\\right)$ where $r$ is the radial distance from the center axis of the beam, $z$ is the axial distance from the beam's focus (or "waist"), $k=2\\pi n/\\lambda$ is the wavenumber (for a free-space wavelength $\\lambda$ and refractive index $n$ of the homogeneous, lossless medium in which the beam propagates), $E_0$ is the electric-field amplitude at the origin, $w(z)$ is the radius at which the field amplitude decays by $1/e$ of its axial values, $w_0$ is the beam waist radius, and $R(z)$ is the radius of curvature of the beam's wavefront at $z$. The only independent parameters that need to be specified are $w_0$, $E_0$, $k$, and the location of the beam focus (i.e., the origin: $r=z=0$).

    In 3d, we use a ["complex point-source" method](https://doi.org/10.1364/JOSAA.16.001381) to define a source that generates an exact Gaussian-beam solution.  In 2d, we currently use the simple approximation of taking a cross-section of the 3d beam.  In both cases, the beam is most accurate near the source's center frequency.) To use the true solution for a 2d Gaussian beam, use the `GaussianBeam2DSource` class instead.

    The `SourceTime` object (`Source.src`), which specifies the time dependence of the source, should normally be a narrow-band `ContinuousSource` or `GaussianSource`.  (For a `CustomSource`, the beam frequency is determined by the source's `center_frequency` parameter.
    """

    def __init__(
        self,
        src,
        center=None,
        volume=None,
        component=mp.ALL_COMPONENTS,
        beam_x0=Vector3(),
        beam_kdir=Vector3(),
        beam_w0=None,
        beam_E0=Vector3(),
        **kwargs,
    ):
        """
        Construct a `GaussianBeamSource`.

        + **`beam_x0` [`Vector3`]** — The location of the beam focus *relative* to the center of the source. The beam focus does *not* need to lie within the source region (i.e., the beam focus can be anywhere, inside or outside the cell, independent of the position of the source).

        + **`beam_kdir` [`Vector3`]** — The propagation direction of the beam. The length is *ignored*. The wavelength of the beam is determined by the center frequency of the `Source.src` object and the refractive index (real part only) at the center of the source region.

        + **`beam_w0` [`number`]** — The beam waist radius.

        + **`beam_E0` [`Vector3`]** — The polarization vector of the beam. Elements can be complex valued (i.e., for circular polarization). The polarization vector must be *parallel* to the source region in order to generate a transverse mode.
        """

        super().__init__(src, component, center, volume, **kwargs)
        self._beam_x0 = beam_x0
        self._beam_kdir = beam_kdir
        self._beam_w0 = beam_w0
        self._beam_E0 = beam_E0

    @property
    def beam_x0(self):
        return self._beam_x0

    @property
    def beam_kdir(self):
        return self._beam_kdir

    @property
    def beam_w0(self):
        return self._beam_w0

    @property
    def beam_E0(self):
        return self._beam_E0

    def add_source(self, sim):
        where = mp.Volume(
            self.center,
            self.size,
            dims=sim.dimensions,
            is_cylindrical=sim.is_cylindrical,
        ).swigobj
        gaussianbeam_args = [
            mp.py_v3_to_vec(
                sim.dimensions, self.beam_x0, is_cylindrical=sim.is_cylindrical
            ),
            mp.py_v3_to_vec(
                sim.dimensions, self.beam_kdir, is_cylindrical=sim.is_cylindrical
            ),
            self.beam_w0,
            self.src.swigobj.frequency().real,
            sim.fields.get_eps(
                mp.py_v3_to_vec(sim.dimensions, self.center, sim.is_cylindrical)
            ).real,
            sim.fields.get_mu(
                mp.py_v3_to_vec(sim.dimensions, self.center, sim.is_cylindrical)
            ).real,
            np.array(
                [self.beam_E0.x, self.beam_E0.y, self.beam_E0.z], dtype=np.complex128
            ),
        ]
        gaussianbeam = mp.gaussianbeam(*gaussianbeam_args)
        add_vol_src_args = [self.src.swigobj, where, gaussianbeam]
        add_vol_src = functools.partial(sim.fields.add_volume_source, *add_vol_src_args)
        add_vol_src()


def get_equiv_sources(field, normal_vec, time_src, center, size):
    """Given the fields along a slice, returns the equivalent sources as meep source objects for the beam."""
    # Get fields
    Ex, Ey, Ez, Hx, Hy, Hz = field
    nHat = normal_vec

    # Electric current K = nHat x H
    Kx = nHat[1] * Hz - nHat[2] * Hy
    Ky = nHat[2] * Hx - nHat[0] * Hz
    Kz = nHat[0] * Hy - nHat[1] * Hx

    # Mangnetic current N = - nHat x E
    Nx = nHat[2] * Ey - nHat[1] * Ez
    Ny = nHat[0] * Ez - nHat[2] * Ex
    Nz = nHat[1] * Ex - nHat[0] * Ey

    # Source components
    components = {mp.Ex: Kx, mp.Ey: Ky, mp.Ez: Kz, mp.Hx: Nx, mp.Hy: Ny, mp.Hz: Nz}

    # Make sources
    sources = [
        mp.Source(
            time_src,
            field_comp,
            center=center,
            size=size,
            amp_data=source_comp,
        )
        for field_comp, source_comp in components.items()
        if np.sum(np.abs(source_comp))
    ]

    return sources


class GaussianBeam2DSource(GaussianBeam3DSource):
    """
    Identical to `GaussianBeam3DSource` except that the beam is defined in 2d.
    This is useful for 2d simulations where the 3d beam is not exact.
    """

    def get_fields(self, sim):
        """Calls green2d under various conditions (incoming vs outgoing) providing the correct Hankel functions and returns the fields at the slice provided by the source."""
        from scipy.special import hankel1e, hankel2e, jve
        from sys import float_info

        # Beam parameters
        freq = self.src.swigobj.frequency().real
        eps = sim.fields.get_eps(
            mp.py_v3_to_vec(sim.dimensions, self.center, sim.is_cylindrical)
        ).real
        mu = sim.fields.get_mu(
            mp.py_v3_to_vec(sim.dimensions, self.center, sim.is_cylindrical)
        ).real
        k = 2 * np.pi * freq * np.sqrt(eps * mu)

        # Get this coordinate system
        center = mp.py_v3_to_vec(sim.dimensions, self.center, sim.is_cylindrical)
        size = mp.py_v3_to_vec(sim.dimensions, self.size, sim.is_cylindrical)
        beam_x0 = mp.py_v3_to_vec(sim.dimensions, self.beam_x0, sim.is_cylindrical)
        beam_kdir = mp.py_v3_to_vec(sim.dimensions, self.beam_kdir, sim.is_cylindrical)
        beam_E0 = self.beam_E0

        # Check for errors
        if size.x() and size.y():
            raise Exception(
                "GaussianBeam2DSource should be a line source, not a plane. Either set size.x or size.y to zero."
            )

        # Complex point source
        x0 = center.x() + beam_x0.x()
        y0 = center.y() + beam_x0.y()
        z0 = k * self.beam_w0**2 / 2
        kdir = beam_kdir.x() + 1j * beam_kdir.y()
        kdir = kdir / np.abs(kdir)
        jx0 = 1j * z0 * np.real(kdir)
        jy0 = 1j * z0 * np.imag(kdir)
        X0 = np.array([x0 + jx0, y0 + jy0]).astype(complex)

        # Create grid
        x = (
            np.linspace(
                center.x() - size.x() / 2,
                center.x() + size.x() / 2,
                int(2 * sim.resolution * size.x()),
            )
            if size.x() > 0
            else np.array([center.x()])
        )
        y = (
            np.linspace(
                center.y() - size.y() / 2,
                center.y() + size.y() / 2,
                int(2 * sim.resolution * size.y()),
            )
            if size.y() > 0
            else np.array([center.y()])
        )
        xx, yy = np.meshgrid(x, y)
        X = np.transpose(np.array([xx, yy])[:, :, :], axes=(0, 2, 1))

        # Find which points are incoming vs outgoing
        incoming_arg, outgoing_arg, waist_points = self.incoming_mask(
            x0, y0, beam_kdir.x(), beam_kdir.y(), X
        )
        waist_xy = np.array([np.real(X0[0]), np.real(X0[1])])[:, np.newaxis, np.newaxis]

        # Find large imaginary argument hankel shift
        r, _ = self.get_r_rhat(X, X0)
        kr = (k * r).flatten()
        max_imag = np.abs(kr.imag[np.argmax(np.abs(kr.imag))])
        max_float = int(np.log(float_info.max) / 2)
        shift = 0 if np.abs(max_imag) < max_float else max_imag - max_float

        # Fix hankel functions for large imaginary arguments and remove nans and infs when non-indexed
        scaled_hankel1 = lambda o, kr: hankel1e(o, kr) * np.exp(1j * kr - shift)
        scaled_hankel2 = lambda o, kr: hankel2e(o, kr) * np.exp(-1j * kr - shift)
        scaled_jv = lambda o, kr: jve(o, kr) * np.exp(np.abs(kr.imag) - shift)

        # Get field for outgoing points (hankel2)
        o_fields2D = self.green2d(
            X[:, outgoing_arg][..., np.newaxis],
            freq,
            eps,
            mu,
            X0,
            kdir,
            scaled_hankel2,
            beam_E0,
        )[:, :, :, np.newaxis]

        # Get field for incoming points (hankel1))
        i_fields2D = self.green2d(
            X[:, incoming_arg][..., np.newaxis],
            freq,
            eps,
            mu,
            X0,
            kdir,
            scaled_hankel1,
            beam_E0,
        )[:, :, :, np.newaxis]

        # Get field for waist points (jv)
        w_fields2D = self.green2d(
            X[:, waist_points][..., np.newaxis],
            freq,
            eps,
            mu,
            X0,
            kdir,
            scaled_jv,
            beam_E0,
        )[:, :, :, np.newaxis]
        w_fields2D_norm = self.green2d(
            waist_xy, freq, eps, mu, X0, kdir, scaled_jv, beam_E0
        )[:, :, :, np.newaxis]

        # Test for overflow and get normalizations properly
        E_fields_norm = np.sqrt(
            np.sum(np.abs(w_fields2D_norm[:3]) ** 2, axis=0).item(0)
        )

        # Normalize fields
        i_fields2D = i_fields2D / E_fields_norm  # hankel1
        o_fields2D = o_fields2D / E_fields_norm  # hankel2
        w_fields2D = w_fields2D / E_fields_norm  # jv

        # Combine fields
        fields2D = np.zeros(
            (6, incoming_arg.shape[0], incoming_arg.shape[1], 1), dtype=np.complex128
        )
        fields2D[:, incoming_arg] += i_fields2D[..., 0]
        fields2D[:, outgoing_arg] += o_fields2D[..., 0]
        fields2D[:, waist_points] += w_fields2D[..., 0]

        return fields2D

    def incoming_mask(self, x0, y0, kx, ky, X):
        """Given a beam with a waist at (x0, y0) and a direction of propagation (kx, ky) returns the boolean masks along the meshgrid X for incoming waves, outgoing waves, and waist points."""

        # Create the plane where the beam is incident
        kx, ky = np.round(kx, 8), np.round(ky, 8)
        plane = lambda x, y: ky * (y - y0) + kx * (x - x0)

        # Find the sign of our points of interest
        point_sign = np.sign(plane(X[0], X[1]))[:, :]

        # Incoming and outgoing points
        incoming_points = point_sign == 1
        outgoing_points = point_sign == -1
        waist_points = point_sign == 0

        return incoming_points, outgoing_points, waist_points

    def get_r_rhat(self, X, X0):
        """Returns r and rhat before normalizing rhat to be used in green2d and get_fields for overflow prediction"""
        rhat = X - np.repeat(
            np.repeat(X0[:, np.newaxis, np.newaxis], X.shape[1], axis=1),
            X.shape[2],
            axis=2,
        )
        r = np.sqrt(np.sum(rhat * rhat, axis=0))
        return r, rhat

    def green2d(self, X, freq, eps, mu, X0, kdir, hankel, beam_E0):
        """Produces the 2D Green's function for an arbitrary complex point source at X0 along the meshgrid X."""

        # Position variables
        EH = np.zeros([6] + list(X.shape)[1:], dtype=complex)
        r, rhat = self.get_r_rhat(X, X0)
        rhat = rhat / r[np.newaxis, :, :]

        # Frequency variables
        omega = 2 * np.pi * freq
        k = omega * np.sqrt(eps * mu)
        ik = 1j * k
        kr = k * r
        Z = np.sqrt(mu / eps)
        H0_kr = hankel(0, kr)
        H1_kr = hankel(1, kr)
        ikH1 = 0.25 * ik * H1_kr

        # E and H source hankel and vector components
        H2_kr = hankel(2, kr)
        p = np.zeros(3, dtype=complex)
        p[0] = beam_E0.x
        p[1] = beam_E0.y
        p[2] = beam_E0.z
        pdotrhat = p[0] * rhat[0] + p[1] * rhat[1]
        rhatcrossp = rhat[0] * p[1] - rhat[1] * p[0]

        # First fill Electric Source fields
        eHx = -rhat[1] * ikH1 * p[2]
        eHy = rhat[0] * ikH1 * p[2]
        eHz = -rhatcrossp * ikH1
        eEx = -(rhat[0] * (pdotrhat / r * 0.25 * Z)) * H1_kr + (
            rhat[1] * (rhatcrossp * omega * mu * 0.125)
        ) * (H0_kr - H2_kr)
        eEy = -(rhat[1] * (pdotrhat / r * 0.25 * Z)) * H1_kr - (
            rhat[0] * (rhatcrossp * omega * mu * 0.125)
        ) * (H0_kr - H2_kr)
        eEz = (-0.25 * omega * mu) * H0_kr * p[2]

        # Create new p vector for magnetic source
        p = -np.array(
            [
                -p[2] * np.imag(kdir),
                p[2] * np.real(kdir),
                p[0] * np.imag(kdir) - p[1] * np.real(kdir),
            ]
        )
        pdotrhat = p[0] * rhat[0] + p[1] * rhat[1]
        rhatcrossp = rhat[0] * p[1] - rhat[1] * p[0]

        # H sources
        hEx = rhat[1] * ikH1 * p[2]
        hEy = -rhat[0] * ikH1 * p[2]
        hEz = rhatcrossp * ikH1
        hHx = -(rhat[0] * (pdotrhat / r * 0.25 / Z)) * H1_kr + (
            rhat[1] * (rhatcrossp * omega * eps * 0.125)
        ) * (H0_kr - H2_kr)
        hHy = -(rhat[1] * (pdotrhat / r * 0.25 / Z)) * H1_kr - (
            rhat[0] * (rhatcrossp * omega * eps * 0.125)
        ) * (H0_kr - H2_kr)
        hHz = (-0.25 * omega * eps) * H0_kr * p[2]

        # Fill arrays
        EH[:] = np.array([eEx, eEy, eEz, eHx, eHy, eHz]) + np.array(
            [hEx, hEy, hEz, hHx, hHy, hHz]
        )  # Ex cross Hz = kdir

        return EH

    def add_source(self, sim):
        """Calls the add_source method for each equivalent source."""
        fields = self.get_fields(sim)

        # Get normal vector nHat
        size = mp.py_v3_to_vec(sim.dimensions, self.size, sim.is_cylindrical)
        if size.x():
            nHat = mp.Vector3(0, 1) * np.sign(self._beam_kdir.y)
        elif size.y():
            nHat = mp.Vector3(1, 0) * np.sign(self._beam_kdir.x)
        sources = get_equiv_sources(fields, nHat, self.src, self.center, self.size)

        for source in sources:
            source.add_source(sim)


class GaussianBeamSource(GaussianBeam3DSource):
    """
    Wrapper for GaussianBeam3DSource to warn the user when running 2D simulations that this behavior is deprecated and they should be using GaussianBeam2DSource.
    """

    def add_source(self, sim):
        if sim.dimensions == 2:
            warnings.warn(
                "GaussianBeamSource is deprecated for 2D simulations. For more accurate results, use GaussianBeam2DSource instead. In the future, this will be the default behavior.",
                DeprecationWarning,
            )
        super().add_source(sim)


class IndexedSource(Source):
    """
    created a source object using (SWIG-wrapped mp::srcdata*) srcdata.
    """

    def __init__(self, src, srcdata, amp_arr, needs_boundary_fix=False):
        self.src = src
        self.num_pts = len(amp_arr)
        self.srcdata = srcdata
        self.amp_arr = amp_arr
        self.needs_boundary_fix = needs_boundary_fix

    def add_source(self, sim):
        sim.fields.register_src_time(self.src.swigobj)
        sim.fields.add_srcdata(
            self.srcdata,
            self.src.swigobj,
            self.num_pts,
            self.amp_arr,
            self.needs_boundary_fix,
        )
        return
