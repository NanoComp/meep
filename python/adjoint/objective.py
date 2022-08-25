"""Handling of objective functions and objective quantities."""
import abc
from collections import namedtuple
from autograd import numpy as npa
import numpy as np
from meep.simulation import py_v3_to_vec

import meep as mp

from .filter_source import FilteredSource

Grid = namedtuple("Grid", ["x", "y", "z", "w"])

# 3 possible components for E x n and H x n
# signs are handled in code
EH_TRANSVERSE = [
    [mp.Hz, mp.Hy, mp.Ez, mp.Ey],
    [mp.Hx, mp.Hz, mp.Ex, mp.Ez],
    [mp.Hy, mp.Hx, mp.Ey, mp.Ex],
]

# Holds the components for each current source
# for the cases of x,y, and z normal vectors.
# This is the same as swapping H and E in the above list
JK_TRANSVERSE = [
    [mp.Ey, mp.Ez, mp.Hy, mp.Hz],
    [mp.Ez, mp.Ex, mp.Hz, mp.Hx],
    [mp.Ex, mp.Ey, mp.Hx, mp.Hy],
]

# Holds the amplitudes used in Poynting Flux adjoint sources
FLUX_AMPLITUDES = np.array([1 / 4, -1 / 4, -1 / 4, 1 / 4], dtype=np.complex128)


class ObjectiveQuantity(abc.ABC):
    """A differentiable objective quantity.

    Attributes:
        sim: the Meep simulation object with which the objective quantity is registered.
        frequencies: the frequencies at which the objective quantity is evaluated.
        num_freq: the number of frequencies at which the objective quantity is evaluated.
    """

    def __init__(self, sim):
        self.sim = sim
        self._eval = None
        self._frequencies = None

    @property
    def frequencies(self):
        return self._frequencies

    @property
    def num_freq(self):
        return len(self.frequencies)

    @abc.abstractmethod
    def __call__(self):
        """Evaluates the objective quantity."""

    @abc.abstractmethod
    def register_monitors(self, frequencies):
        """Registers monitors in the forward simulation."""

    @abc.abstractmethod
    def place_adjoint_source(self, dJ):
        """Places appropriate sources for the adjoint simulation."""

    def get_evaluation(self):
        """Evaluates the objective quantity."""
        if self._eval is not None:
            return self._eval
        else:
            raise RuntimeError(
                "You must first run a forward simulation before requesting the evaluation of an objective quantity."
            )

    def _adj_src_scale(self, include_resolution=True):
        """Calculates the scale for the adjoint sources."""
        T = self.sim.meep_time()
        dt = self.sim.fields.dt
        src = self._create_time_profile()

        if include_resolution:
            num_dims = self.sim._infer_dimensions(self.sim.k_point)
            dV = 1 / self.sim.resolution**num_dims
        else:
            dV = 1

        iomega = (1.0 - np.exp(-1j * (2 * np.pi * self._frequencies) * dt)) * (
            1.0 / dt
        )  # scaled frequency factor with discrete time derivative fix

        # an ugly way to calcuate the scaled dtft of the forward source
        y = np.array(
            [src.swigobj.current(t, dt) for t in np.arange(0, T, dt)]
        )  # time domain signal
        fwd_dtft = (
            np.matmul(
                np.exp(
                    1j
                    * 2
                    * np.pi
                    * self._frequencies[:, np.newaxis]
                    * np.arange(y.size)
                    * dt
                ),
                y,
            )
            * dt
            / np.sqrt(2 * np.pi)
        )  # dtft

        # Interestingly, the real parts of the DTFT and fourier transform match, but the imaginary parts are very different...
        # fwd_dtft = src.fourier_transform(src.frequency)
        """
        For some reason, there seems to be an additional phase
        factor at the center frequency that needs to be applied
        to *all* frequencies...
        """
        src_center_dtft = (
            np.matmul(
                np.exp(
                    1j
                    * 2
                    * np.pi
                    * np.array([src.frequency])[:, np.newaxis]
                    * np.arange(y.size)
                    * dt
                ),
                y,
            )
            * dt
            / np.sqrt(2 * np.pi)
        )
        adj_src_phase = np.exp(1j * np.angle(src_center_dtft)) * self.fwidth_scale

        if self._frequencies.size == 1:
            # Single frequency simulations. We need to drive it with a time profile.
            scale = dV * iomega / fwd_dtft / adj_src_phase  # final scale factor
        else:
            # multi frequency simulations
            scale = dV * iomega / adj_src_phase
        # compensate for the fact that real fields take the real part of the current,
        # which halves the Fourier amplitude at the positive frequency (Re[J] = (J + J*)/2)
        if self.sim.using_real_fields():
            scale *= 2
        return scale

    def _create_time_profile(self, fwidth_frac=0.1, adj_cutoff=5):
        """Creates a time domain waveform for normalizing the adjoint source(s).

        For single frequency objective functions, we should generate a guassian pulse with a reasonable
        bandwidth centered at said frequency.

        TODO:
        The user may specify a scalar valued objective function across multiple frequencies (e.g. MSE) in
        which case we should check that all the frequencies fit in the specified bandwidth.
        """
        self.fwidth_scale = np.exp(-2j * np.pi * adj_cutoff / fwidth_frac)
        return mp.GaussianSource(
            np.mean(self._frequencies),
            fwidth=fwidth_frac * np.mean(self._frequencies),
            cutoff=adj_cutoff,
        )


class EigenmodeCoefficient(ObjectiveQuantity):
    """A frequency-dependent eigenmode coefficient.
    Attributes:
        volume: the volume over which the eigenmode coefficient is calculated.
        mode: the eigenmode number.
        forward: whether the forward or backward mode coefficient is returned as
          the result of the evaluation.
        kpoint_func: an optional k-point function to use when evaluating the eigenmode
          coefficient. When specified, this overrides the effect of `forward`.
        kpoint_func_overlap_idx: the index of the mode coefficient to return when
          specifying `kpoint_func`. When specified, this overrides the effect of
          `forward` and should have a value of either 0 or 1.
    """

    def __init__(
        self,
        sim,
        volume,
        mode,
        forward=True,
        kpoint_func=None,
        kpoint_func_overlap_idx=0,
        decimation_factor=0,
        **kwargs
    ):
        super().__init__(sim)
        if kpoint_func_overlap_idx not in [0, 1]:
            raise ValueError(
                "`kpoint_func_overlap_idx` should be either 0 or 1, but got %d"
                % (kpoint_func_overlap_idx,)
            )
        self.volume = volume
        self.mode = mode
        self.forward = forward
        self.kpoint_func = kpoint_func
        self.kpoint_func_overlap_idx = kpoint_func_overlap_idx
        self.eigenmode_kwargs = kwargs
        self._monitor = None
        self._cscale = None
        self.decimation_factor = decimation_factor

    def register_monitors(self, frequencies):
        self._frequencies = np.asarray(frequencies)
        self._monitor = self.sim.add_mode_monitor(
            frequencies,
            mp.ModeRegion(center=self.volume.center, size=self.volume.size),
            yee_grid=True,
            decimation_factor=self.decimation_factor,
        )
        return self._monitor

    def place_adjoint_source(self, dJ):
        dJ = np.atleast_1d(dJ)
        if dJ.ndim == 2:
            dJ = np.sum(dJ, axis=1)
        time_src = self._create_time_profile()
        da_dE = 0.5 * self._cscale
        scale = self._adj_src_scale()

        if self.kpoint_func:
            eig_kpoint = -1 * self.kpoint_func(time_src.frequency, self.mode)
        else:
            center_frequency = 0.5 * (
                np.min(self.frequencies) + np.max(self.frequencies)
            )
            direction = mp.Vector3(
                *(np.eye(3)[self._monitor.normal_direction] * np.abs(center_frequency))
            )
            eig_kpoint = -1 * direction if self.forward else direction

        if self._frequencies.size == 1:
            amp = da_dE * dJ * scale
            src = time_src
        else:
            scale = da_dE * dJ * scale
            src = FilteredSource(
                time_src.frequency,
                self._frequencies,
                scale,
                self.sim.fields.dt,
            )
            amp = 1
        source = mp.EigenModeSource(
            src,
            eig_band=self.mode,
            direction=mp.NO_DIRECTION,
            eig_kpoint=eig_kpoint,
            amplitude=amp,
            eig_match_freq=True,
            size=self.volume.size,
            center=self.volume.center,
            **self.eigenmode_kwargs,
        )
        return [source]

    def __call__(self):
        if self.kpoint_func:
            kpoint_func = self.kpoint_func
            overlap_idx = self.kpoint_func_overlap_idx
        else:
            center_frequency = 0.5 * (
                np.min(self.frequencies) + np.max(self.frequencies)
            )
            kpoint = mp.Vector3(
                *(np.eye(3)[self._monitor.normal_direction] * np.abs(center_frequency))
            )
            kpoint_func = lambda *not_used: kpoint if self.forward else -1 * kpoint
            overlap_idx = 0
        ob = self.sim.get_eigenmode_coefficients(
            self._monitor,
            [self.mode],
            direction=mp.NO_DIRECTION,
            kpoint_func=kpoint_func,
            **self.eigenmode_kwargs,
        )
        overlaps = ob.alpha.squeeze(axis=0)
        assert overlaps.ndim == 2
        self._eval = overlaps[:, overlap_idx]
        self._cscale = ob.cscale
        return self._eval


class FourierFields(ObjectiveQuantity):
    def __init__(self, sim, volume, component, yee_grid=False, decimation_factor=0):
        super().__init__(sim)
        self.volume = sim._fit_volume_to_simulation(volume)
        self.component = component
        self.yee_grid = yee_grid
        self.decimation_factor = decimation_factor

    def register_monitors(self, frequencies):
        self._frequencies = np.asarray(frequencies)
        self._monitor = self.sim.add_dft_fields(
            [self.component],
            self._frequencies,
            where=self.volume,
            yee_grid=self.yee_grid,
            decimation_factor=self.decimation_factor,
        )
        return self._monitor

    def place_adjoint_source(self, dJ):
        time_src = self._create_time_profile()
        sources = []

        mon_size = self.sim.fields.dft_monitor_size(
            self._monitor.swigobj, self.volume.swigobj, self.component
        )
        dJ = dJ.astype(np.complex128)
        if (
            np.prod(mon_size) * self.num_freq != dJ.size
            and np.prod(mon_size) * self.num_freq**2 != dJ.size
        ):
            raise ValueError("The format of J is incorrect!")

        # The objective function J is a vector. Each component corresponds to a frequency.
        if np.prod(mon_size) * self.num_freq**2 == dJ.size and self.num_freq > 1:
            dJ = np.sum(dJ, axis=1)
        """The adjoint solver requires the objective function
        to be scalar valued with regard to objective arguments
        and position, but the function may be vector valued
        with regard to frequency. In this case, the Jacobian
        will be of the form [F,F,...] where F is the number of
        frequencies. Because of linearity, we can sum across the
        second frequency dimension to calculate a frequency
        scale factor for each point (rather than a scale vector).
        """

        self.all_fouriersrcdata = self._monitor.swigobj.fourier_sourcedata(
            self.volume.swigobj, self.component, self.sim.fields, dJ
        )

        for fourier_data in self.all_fouriersrcdata:
            amp_arr = np.array(fourier_data.amp_arr).reshape(-1, self.num_freq)
            scale = amp_arr * self._adj_src_scale(include_resolution=False)

            if self.num_freq == 1:
                sources += [
                    mp.IndexedSource(
                        time_src, fourier_data, scale[:, 0], not self.yee_grid
                    )
                ]
            else:
                src = FilteredSource(
                    time_src.frequency, self._frequencies, scale, self.sim.fields.dt
                )
                (num_basis, num_pts) = src.nodes.shape
                for basis_i in range(num_basis):
                    sources += [
                        mp.IndexedSource(
                            src.time_src_bf[basis_i],
                            fourier_data,
                            src.nodes[basis_i],
                            not self.yee_grid,
                        )
                    ]
        return sources

    def __call__(self):
        self._eval = np.array(
            [
                self.sim.get_dft_array(self._monitor, self.component, i)
                for i in range(self.num_freq)
            ]
        )
        return self._eval


class Near2FarFields(ObjectiveQuantity):
    def __init__(self, sim, Near2FarRegions, far_pts, decimation_factor=0):
        super().__init__(sim)
        self.Near2FarRegions = Near2FarRegions
        self.far_pts = far_pts  # list of far pts
        self._nfar_pts = len(far_pts)
        self.decimation_factor = decimation_factor

    def register_monitors(self, frequencies):
        self._frequencies = np.asarray(frequencies)
        self._monitor = self.sim.add_near2far(
            self._frequencies,
            *self.Near2FarRegions,
            decimation_factor=self.decimation_factor,
        )
        return self._monitor

    def place_adjoint_source(self, dJ):
        time_src = self._create_time_profile()
        sources = []
        if dJ.ndim == 4:
            dJ = np.sum(dJ, axis=0)

        farpt_list = np.array([list(pi) for pi in self.far_pts]).flatten()
        far_pt0 = self.far_pts[0]
        far_pt_vec = py_v3_to_vec(
            self.sim.dimensions,
            far_pt0,
            self.sim.is_cylindrical,
        )

        all_nearsrcdata = self._monitor.swigobj.near_sourcedata(
            far_pt_vec, farpt_list, self._nfar_pts, dJ
        )
        for near_data in all_nearsrcdata:
            cur_comp = near_data.near_fd_comp
            amp_arr = np.array(near_data.amp_arr).reshape(-1, self.num_freq)
            scale = amp_arr * self._adj_src_scale(include_resolution=False)

            if self.num_freq == 1:
                sources += [mp.IndexedSource(time_src, near_data, scale[:, 0])]
            else:
                src = FilteredSource(
                    time_src.frequency,
                    self._frequencies,
                    scale,
                    self.sim.fields.dt,
                )
                (num_basis, num_pts) = src.nodes.shape
                for basis_i in range(num_basis):
                    sources += [
                        mp.IndexedSource(
                            src.time_src_bf[basis_i],
                            near_data,
                            src.nodes[basis_i],
                        )
                    ]
        return sources

    def __call__(self):
        self._eval = np.array(
            [self.sim.get_farfield(self._monitor, far_pt) for far_pt in self.far_pts]
        ).reshape((self._nfar_pts, self.num_freq, 6))
        return self._eval


class LDOS(ObjectiveQuantity):
    def __init__(self, sim, decimation_factor=0, **kwargs):
        super().__init__(sim)
        self.decimation_factor = decimation_factor
        self.srckwarg = kwargs

    def register_monitors(self, frequencies):
        self._frequencies = np.asarray(frequencies)
        self.forward_src = self.sim.sources
        return

    def place_adjoint_source(self, dJ):
        time_src = self._create_time_profile()
        if dJ.ndim == 2:
            dJ = np.sum(dJ, axis=1)
        dJ = dJ.flatten()
        sources = []
        forward_f_scale = np.array(
            [self.ldos_scale / self.ldos_Jdata[k] for k in range(self.num_freq)]
        )
        if self._frequencies.size == 1:
            amp = (dJ * self._adj_src_scale(False) * forward_f_scale)[0]
            src = time_src
        else:
            scale = dJ * self._adj_src_scale(False) * forward_f_scale
            src = FilteredSource(
                time_src.frequency,
                self._frequencies,
                scale,
                self.sim.fields.dt,
            )
            amp = 1
        for forward_src_i in self.forward_src:
            if isinstance(forward_src_i, mp.EigenModeSource):
                src_i = mp.EigenModeSource(
                    src,
                    component=forward_src_i.component,
                    eig_kpoint=forward_src_i.eig_kpoint,
                    amplitude=amp,
                    eig_band=forward_src_i.eig_band,
                    size=forward_src_i.size,
                    center=forward_src_i.center,
                    **self.srckwarg,
                )
            else:
                src_i = mp.Source(
                    src,
                    component=forward_src_i.component,
                    amplitude=amp,
                    size=forward_src_i.size,
                    center=forward_src_i.center,
                    **self.srckwarg,
                )
            if mp.is_electric(src_i.component):
                src_i.amplitude *= -1
            sources += [src_i]

        return sources

    def __call__(self):
        self._eval = self.sim.ldos_data
        self.ldos_scale = self.sim.ldos_scale
        self.ldos_Jdata = self.sim.ldos_Jdata
        return np.array(self._eval)


class PoyntingFlux(ObjectiveQuantity):
    """A frequency-dependent Poynting Flux adjoint source.
    Attributes:
        volume: The volume over which the Poynting Flux is calculated.
        This function currently only works for monitors with a defined
        normal vector (e.g. planes in 3d or lines in 2d). User supplied
        normal vectors may be implemented in the future. It also only
        works with monitors aligned to a coordinate direction.
        decimation_factor: Whether to skip points in the time series every
        decimation_factor timesteps. See "add_dft_fields" documentation.
        The biggest warning there is to be careful to avoid aliasing if
        the fields vary quickly in time.
        Note on yee_grid: For the Poynting Flux to work, H and E components
        must lie at the same points. Therefore, the Yee grid will always be false.
    """

    def __init__(self, sim, volume, scale=1, decimation_factor=0):
        super().__init__(sim)
        # _fit_volume_to_simulation increases the dimensionality of
        # the volume, so we'll use the user input volume
        self.volume = sim._fit_volume_to_simulation(volume)
        self.decimation_factor = decimation_factor
        self.scale = scale
        # get_normal returns an index for the two
        # dictionaries of cross products
        self.normal = self.get_normal(volume)

    def register_monitors(self, frequencies):
        self._frequencies = np.asarray(frequencies)
        self._monitor = []
        # List to hold FourierFields objects
        self.F_fields_list = []
        for comp in EH_TRANSVERSE[self.normal]:
            # instantiate the FourierFields monitors
            F_field = FourierFields(self.sim, self.volume, comp)
            self.F_fields_list.append(F_field)
            self._monitor.append(F_field.register_monitors(self._frequencies))
        return self._monitor

    def place_adjoint_source(self, dJ):
        source = []
        for pos, field in enumerate(self.F_fields_list):
            source.append(field.place_adjoint_source(np.flipud(dJ[pos].flatten()))[0])
        return source

    def __call__(self):
        self.field_component_evaluations = []
        for field in self.F_fields_list:
            # Get the dft evaluation from a call to the underlying
            # FourierFields object
            field_here = field()
            self.field_component_evaluations.append(field_here)
        # Get integration weights Meep uses
        self.metadata = self.sim.get_array_metadata(vol=self.volume)
        self.field_component_evaluations.append(
            np.array([self.metadata[3]]).astype(complex)
        )
        [H1, H2, E1, E2, meta] = self.field_component_evaluations
        self._eval = self.field_component_evaluations
        print("This is meta*E2")
        print(meta * E2)
        return np.array([H1, H2, E1, E2, meta])

    # takes in a 1x5xNxM vector where the size five array corresponds to
    # [H1,H2,E1,E1,meta]
    # multiple frequencies will be tested later
    @staticmethod
    def compute_flux(*inputs):
        flux = npa.sum(
            npa.real(
                inputs[0][4]
                * (
                    npa.conj(inputs[0][0]) * inputs[0][3]
                    - npa.conj(inputs[0][1]) * inputs[0][2]
                )
            )
        )
        return flux

    # returns 0,1, or 2 corresponding to x, y, or z normal vectors
    # TODO: Handle user-specified normal vectors and cases when 2d
    # has a zero-size dimension other than z
    def get_normal(self, volume):
        # I'll add cylindrical later (since the normal vector gets a little different)
        if volume.size.x == 0:
            return 0
        elif volume.size.y == 0:
            return 1
        else:
            return 2
