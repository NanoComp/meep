"""Handling of objective functions and objective quantities."""

from abc import ABC, abstractmethod
import numpy as np
import meep as mp
from .filter_source import FilteredSource
from .optimization_problem import Grid
from meep.simulation import py_v3_to_vec


class ObjectiveQuantitiy(ABC):
    @abstractmethod
    def __init__(self):
        return

    @abstractmethod
    def register_monitors(self):
        return

    @abstractmethod
    def place_adjoint_source(self):
        return

    @abstractmethod
    def __call__(self):
        return

    @abstractmethod
    def get_evaluation(self):
        return


class EigenmodeCoefficient(ObjectiveQuantitiy):
    def __init__(self,
                 sim,
                 volume,
                 mode,
                 forward=True,
                 kpoint_func=None,
                 **kwargs):
        '''
        '''
        self.sim = sim
        self.volume = volume
        self.mode = mode
        self.forward = 0 if forward else 1
        self.normal_direction = None
        self.kpoint_func = kpoint_func
        self.eval = None
        self.EigenMode_kwargs = kwargs
        return

    def register_monitors(self, frequencies):
        self.frequencies = np.asarray(frequencies)
        self.monitor = self.sim.add_mode_monitor(frequencies,
                                                 mp.ModeRegion(
                                                     center=self.volume.center,
                                                     size=self.volume.size),
                                                 yee_grid=True)
        self.normal_direction = self.monitor.normal_direction
        return self.monitor

    def place_adjoint_source(self, dJ):
        '''Places an equivalent eigenmode monitor facing the opposite direction. Calculates the
        correct scaling/time profile.
        dJ ........ the user needs to pass the dJ/dMonitor evaluation
        '''
        dJ = np.atleast_1d(dJ)
        dt = self.sim.fields.dt  # the timestep size from sim.fields.dt of the forward sim
        # determine starting kpoint for reverse mode eigenmode source
        direction_scalar = 1 if self.forward else -1
        if self.kpoint_func is None:
            if self.normal_direction == 0:
                k0 = direction_scalar * mp.Vector3(x=1)
            elif self.normal_direction == 1:
                k0 = direction_scalar * mp.Vector3(y=1)
            elif self.normal_direction == 2:
                k0 == direction_scalar * mp.Vector3(z=1)
        else:
            k0 = direction_scalar * self.kpoint_func(self.time_src.frequency,
                                                     1)
        if dJ.ndim == 2:
            dJ = np.sum(dJ, axis=1)
        da_dE = 0.5 * self.cscale  # scalar popping out of derivative

        scale = adj_src_scale(self, dt)

        if self.frequencies.size == 1:
            # Single frequency simulations. We need to drive it with a time profile.
            amp = da_dE * dJ * scale  # final scale factor
            src = self.time_src
        else:
            # multi frequency simulations
            scale = da_dE * dJ * scale
            src = FilteredSource(
                self.time_src.frequency,
                self.frequencies,
                scale,
                dt,
            )  # generate source from broadband response
            amp = 1

        # generate source object
        self.source = [
            mp.EigenModeSource(
                src,
                eig_band=self.mode,
                direction=mp.NO_DIRECTION,
                eig_kpoint=k0,
                amplitude=amp,
                eig_match_freq=True,
                size=self.volume.size,
                center=self.volume.center,
                **self.EigenMode_kwargs,
            )
        ]

        return self.source

    def __call__(self):
        # Eigenmode data
        self.time_src = create_time_profile(self)
        direction = mp.NO_DIRECTION if self.kpoint_func else mp.AUTOMATIC
        ob = self.sim.get_eigenmode_coefficients(
            self.monitor,
            [self.mode],
            direction=direction,
            kpoint_func=self.kpoint_func,
            **self.EigenMode_kwargs,
        )
        self.eval = np.squeeze(ob.alpha[:, :, self.forward]
                               )  # record eigenmode coefficients for scaling
        self.cscale = ob.cscale  # pull scaling factor

        return self.eval

    def get_evaluation(self):
        '''Returns the requested eigenmode coefficient.
        '''
        try:
            return self.eval
        except AttributeError:
            raise RuntimeError(
                "You must first run a forward simulation before resquesting an eigenmode coefficient."
            )


class FourierFields(ObjectiveQuantitiy):
    def __init__(self, sim, volume, component):
        self.sim = sim
        self.volume = volume
        self.eval = None
        self.component = component
        return

    def register_monitors(self, frequencies):
        self.frequencies = np.asarray(frequencies)
        self.num_freq = len(self.frequencies)
        self.monitor = self.sim.add_dft_fields(
            [self.component],
            self.frequencies,
            where=self.volume,
            yee_grid=False,
        )
        return self.monitor

    def place_adjoint_source(self, dJ):
        dt = self.sim.fields.dt  # the timestep size from sim.fields.dt of the forward sim
        self.sources = []
        scale = adj_src_scale(self, dt)

        x_dim, y_dim, z_dim = len(self.dg.x), len(self.dg.y), len(self.dg.z)

        if self.num_freq == 1:
            amp = -dJ[0].copy().reshape(x_dim, y_dim, z_dim) * scale
            src = self.time_src
            if self.component in [mp.Hx, mp.Hy, mp.Hz]:
                amp = -amp
            for zi in range(z_dim):
                for yi in range(y_dim):
                    for xi in range(x_dim):
                        if amp[xi, yi, zi] != 0:
                            self.sources += [
                                mp.Source(
                                    src,
                                    component=self.component,
                                    amplitude=amp[xi, yi, zi],
                                    center=mp.Vector3(self.dg.x[xi],
                                                      self.dg.y[yi],
                                                      self.dg.z[zi]),
                                )
                            ]
        else:
            '''The adjoint solver requires the objective function
            to be scalar valued with regard to objective arguments
            and position, but the function may be vector valued
            with regard to frequency. In this case, the Jacobian
            will be of the form [F,F,...] where F is the number of
            frequencies. Because of linearity, we can sum across the
            second frequency dimension to calculate a frequency
            scale factor for each point (rather than a scale vector).
            '''
            dJ = np.sum(
                dJ,
                axis=1)  # sum along first dimension bc Jacobian is always diag
            dJ_4d = np.array([
                dJ[f].copy().reshape(x_dim, y_dim, z_dim)
                for f in range(self.num_freq)
            ])
            if self.component in [mp.Hx, mp.Hy, mp.Hz]:
                dJ_4d = -dJ_4d
            for zi in range(z_dim):
                for yi in range(y_dim):
                    for xi in range(x_dim):
                        '''We only need to add a current source if the
                        jacobian is nonzero for all frequencies at 
                        that particular point. Otherwise, the fitting
                        algorithm is going to fail.
                        '''
                        if not np.all((dJ_4d[:, xi, yi, zi] == 0)):
                            final_scale = -dJ_4d[:, xi, yi, zi] * scale
                            src = FilteredSource(
                                self.time_src.frequency,
                                self.frequencies,
                                final_scale,
                                dt,
                            )
                            self.sources += [
                                mp.Source(
                                    src,
                                    component=self.component,
                                    amplitude=1,
                                    center=mp.Vector3(self.dg.x[xi],
                                                      self.dg.y[yi],
                                                      self.dg.z[zi]),
                                )
                            ]

        return self.sources

    def __call__(self):
        self.dg = Grid(*self.sim.get_array_metadata(dft_cell=self.monitor))
        self.eval = np.array([
            self.sim.get_dft_array(self.monitor, self.component, i)
            for i in range(self.num_freq)
        ])
        self.time_src = create_time_profile(self)
        return self.eval

    def get_evaluation(self):
        try:
            return self.eval
        except AttributeError:
            raise RuntimeError("You must first run a forward simulation.")


class Near2FarFields(ObjectiveQuantitiy):
    def __init__(self, sim, Near2FarRegions, far_pts):
        self.sim = sim
        self.Near2FarRegions = Near2FarRegions
        self.eval = None
        self.far_pts = far_pts  #list of far pts
        self.nfar_pts = len(far_pts)
        return

    def register_monitors(self, frequencies):
        self.frequencies = np.asarray(frequencies)
        self.num_freq = len(self.frequencies)
        self.monitor = self.sim.add_near2far(
            self.frequencies,
            *self.Near2FarRegions,
            yee_grid=True,
        )
        return self.monitor

    def place_adjoint_source(self, dJ):
        dt = self.sim.fields.dt  # the timestep size from sim.fields.dt of the forward sim
        self.sources = []
        if dJ.ndim == 4:
            dJ = np.sum(dJ, axis=0)
        dJ = dJ.flatten()
        farpt_list = np.array([list(pi) for pi in self.far_pts]).flatten()
        far_pt0 = self.far_pts[0]
        far_pt_vec = py_v3_to_vec(
            self.sim.dimensions,
            far_pt0,
            self.sim.is_cylindrical,
        )

        self.all_nearsrcdata = self.monitor.swigobj.near_sourcedata(
            far_pt_vec, farpt_list, self.nfar_pts, dJ)
        for near_data in self.all_nearsrcdata:
            cur_comp = near_data.near_fd_comp
            amp_arr = np.array(near_data.amp_arr).reshape(-1, self.num_freq)
            scale = amp_arr * adj_src_scale(self, dt, include_resolution=False)

            if self.num_freq == 1:
                self.sources += [
                    mp.IndexedSource(self.time_src, near_data, scale[:, 0])
                ]
            else:
                src = FilteredSource(
                    self.time_src.frequency,
                    self.frequencies,
                    scale,
                    dt,
                )
                (num_basis, num_pts) = src.nodes.shape
                for basis_i in range(num_basis):
                    self.sources += [
                        mp.IndexedSource(
                            src.time_src_bf[basis_i],
                            near_data,
                            src.nodes[basis_i],
                        )
                    ]

        return self.sources

    def __call__(self):
        self.time_src = create_time_profile(self)
        self.eval = np.array([
            self.sim.get_farfield(self.monitor, far_pt)
            for far_pt in self.far_pts
        ]).reshape((self.nfar_pts, self.num_freq, 6))
        return self.eval

    def get_evaluation(self):
        try:
            return self.eval
        except AttributeError:
            raise RuntimeError("You must first run a forward simulation.")


def create_time_profile(obj_quantity, fwidth_frac=0.1):
    '''
    For single frequency objective functions, we should
    generate a guassian pulse with a reasonable bandwidth
    centered at said frequency.
    
    TODO:
    The user may specify a scalar valued objective
    function across multiple frequencies (e.g. MSE)
    in which case we should check that all the frequencies
    fit in the specified bandwidth.
    '''
    return mp.GaussianSource(
        np.mean(obj_quantity.frequencies),
        fwidth=fwidth_frac * np.mean(obj_quantity.frequencies),
    )


def adj_src_scale(obj_quantity, dt, include_resolution=True):
    # -------------------------------------- #
    # Get scaling factor
    # -------------------------------------- #
    # leverage linearity and combine source for multiple frequencies
    T = obj_quantity.sim.meep_time()

    if not include_resolution:
        dV = 1
    elif obj_quantity.sim.cell_size.y == 0:
        dV = 1 / obj_quantity.sim.resolution
    elif obj_quantity.sim.cell_size.z == 0:
        dV = 1 / obj_quantity.sim.resolution * 1 / obj_quantity.sim.resolution
    else:
        dV = 1 / obj_quantity.sim.resolution * 1 / obj_quantity.sim.resolution * 1 / obj_quantity.sim.resolution

    iomega = (1.0 -
              np.exp(-1j * (2 * np.pi * obj_quantity.frequencies) * dt)) * (
                  1.0 / dt
              )  # scaled frequency factor with discrete time derivative fix

    src = obj_quantity.time_src

    # an ugly way to calcuate the scaled dtft of the forward source
    y = np.array([src.swigobj.current(t, dt)
                  for t in np.arange(0, T, dt)])  # time domain signal
    fwd_dtft = np.matmul(
        np.exp(1j * 2 * np.pi * obj_quantity.frequencies[:, np.newaxis] *
               np.arange(y.size) * dt), y) * dt / np.sqrt(2 * np.pi)  # dtft

    # Interestingly, the real parts of the DTFT and fourier transform match, but the imaginary parts are very different...
    #fwd_dtft = src.fourier_transform(src.frequency)
    '''
    For some reason, there seems to be an additional phase
    factor at the center frequency that needs to be applied
    to *all* frequencies...
    '''
    src_center_dtft = np.matmul(
        np.exp(1j * 2 * np.pi * np.array([src.frequency])[:, np.newaxis] *
               np.arange(y.size) * dt), y) * dt / np.sqrt(2 * np.pi)
    adj_src_phase = np.exp(1j * np.angle(src_center_dtft))

    if obj_quantity.frequencies.size == 1:
        # Single frequency simulations. We need to drive it with a time profile.
        scale = dV * iomega / fwd_dtft / adj_src_phase  # final scale factor
    else:
        # multi frequency simulations
        scale = dV * iomega / adj_src_phase
    return scale
