from __future__ import division, print_function

import functools
import math
import numbers
import os
import re
import subprocess
import sys
import warnings
from collections import namedtuple
from collections import Sequence

import numpy as np

import meep as mp
from meep.geom import Vector3
from meep.source import EigenModeSource, check_positive


try:
    basestring
except NameError:
    basestring = str


FluxData = namedtuple('FluxData', ['E', 'H'])
ForceData = namedtuple('ForceData', ['offdiag1', 'offdiag2', 'diag'])
NearToFarData = namedtuple('NearToFarData', ['F'])


def get_num_args(func):
    if isinstance(func, Harminv):
        return 2
    return func.__code__.co_argcount


def vec(*args):
    try:
        # Check for vec(x, [y, [z]])
        return mp._vec(*args)
    except (TypeError, NotImplementedError):
        try:
            # Check for vec(iterable)
            if len(args) != 1:
                raise TypeError

            return mp._vec(*args[0])
        except (TypeError, NotImplementedError):
            print("Expected an iterable with three or fewer floating point values")
            print("    or something of the form vec(x, [y, [z]])")
            raise


def py_v3_to_vec(dims, v3, is_cylindrical=False):
    if dims == 1:
        return mp.vec(v3.z)
    elif dims == 2:
        if is_cylindrical:
            return mp.veccyl(v3.x, v3.z)
        else:
            return mp.vec(v3.x, v3.y)
    elif dims == 3:
        return mp.vec(v3.x, v3.y, v3.z)
    else:
        raise ValueError("Invalid dimensions in Volume: {}".format(dims))


class PML(object):

    def __init__(self, thickness,
                 direction=mp.ALL,
                 side=mp.ALL,
                 R_asymptotic=1e-15,
                 mean_stretch=1.0,
                 pml_profile=lambda u: u * u):

        self.thickness = thickness
        self.direction = direction
        self.side = side
        self.R_asymptotic = R_asymptotic
        self.mean_stretch = mean_stretch
        self.pml_profile = pml_profile

        if direction == mp.ALL and side == mp.ALL:
            self.swigobj = mp.pml(thickness, R_asymptotic, mean_stretch)
        elif direction == mp.ALL:
            self.swigobj = mp.pml(thickness, side, R_asymptotic, mean_stretch)
        else:
            self.swigobj = mp.pml(thickness, direction, side, R_asymptotic, mean_stretch)

    @property
    def R_asymptotic(self):
        return self._R_asymptotic

    @R_asymptotic.setter
    def R_asymptotic(self, val):
        self._R_asymptotic = check_positive('PML.R_asymptotic', val)

    @property
    def mean_stretch(self):
        return self._mean_stretch

    @mean_stretch.setter
    def mean_stretch(self, val):
        if val >= 1:
            self._mean_stretch = val
        else:
            raise ValueError("PML.mean_stretch must be >= 1. Got {}".format(val))


class Absorber(PML):
    pass


class Symmetry(object):

    def __init__(self, direction, phase=1):
        self.direction = direction
        self.phase = complex(phase)
        self.swigobj = None


class Rotate2(Symmetry):
    pass


class Rotate4(Symmetry):
    pass


class Mirror(Symmetry):
    pass


class Identity(Symmetry):
    pass


class Volume(object):

    def __init__(self, center, size=Vector3(), dims=2, is_cylindrical=False):
        self.center = center
        self.size = size
        self.dims = dims

        v1 = center - size.scale(0.5)
        v2 = center + size.scale(0.5)

        vec1 = py_v3_to_vec(self.dims, v1, is_cylindrical)
        vec2 = py_v3_to_vec(self.dims, v2, is_cylindrical)

        self.swigobj = mp.volume(vec1, vec2)


class FluxRegion(object):

    def __init__(self, center, size=Vector3(), direction=mp.AUTOMATIC, weight=1.0):
        self.center = center
        self.size = size
        self.direction = direction
        self.weight = complex(weight)


ModeRegion = FluxRegion


class ForceRegion(object):

    def __init__(self, center, direction, size=mp.Vector3(), weight=1.0):
        self.center = center
        self.direction = direction
        self.size = size
        self.weight = complex(weight)


class Near2FarRegion(object):

    def __init__(self, center, size=mp.Vector3(), direction=mp.AUTOMATIC, weight=1.0):
        self.center = center
        self.size = size
        self.direction = direction
        self.weight = complex(weight)


class FieldsRegion(object):

    def __init__(self, where=None, center=None, size=None):
        self.where = where
        self.center = center
        self.size = size


class DftObj(object):

    def __init__(self, func, args):
        self.func = func
        self.args = args
        self.swigobj = None

    def swigobj_attr(self, attr):
        if self.swigobj is None:
            self.swigobj = self.func(*self.args)
        return getattr(self.swigobj, attr)

    @property
    def save_hdf5(self):
        return self.swigobj_attr('save_hdf5')

    @property
    def load_hdf5(self):
        return self.swigobj_attr('load_hdf5')

    @property
    def scale_dfts(self):
        return self.swigobj_attr('scale_dfts')

    @property
    def remove(self):
        return self.swigobj_attr('remove')

    @property
    def freq_min(self):
        return self.swigobj_attr('freq_min')

    @property
    def dfreq(self):
        return self.swigobj_attr('dfreq')

    @property
    def Nfreq(self):
        return self.swigobj_attr('Nfreq')

    @property
    def where(self):
        return self.swigobj_attr('where')


class DftFlux(DftObj):

    def __init__(self, func, args):
        super(DftFlux, self).__init__(func, args)
        self.nfreqs = args[2]
        self.regions = args[3]
        self.num_components = 4

    @property
    def flux(self):
        return self.swigobj_attr('flux')

    @property
    def E(self):
        return self.swigobj_attr('E')

    @property
    def H(self):
        return self.swigobj_attr('H')

    @property
    def cE(self):
        return self.swigobj_attr('cE')

    @property
    def cH(self):
        return self.swigobj_attr('cH')

    @property
    def normal_direction(self):
        return self.swigobj_attr('normal_direction')


class DftForce(DftObj):

    def __init__(self, func, args):
        super(DftForce, self).__init__(func, args)
        self.nfreqs = args[2]
        self.regions = args[3]
        self.num_components = 6

    @property
    def force(self):
        return self.swigobj_attr('force')

    @property
    def offdiag1(self):
        return self.swigobj_attr('offdiag1')

    @property
    def offdiag2(self):
        return self.swigobj_attr('offdiag2')

    @property
    def diag(self):
        return self.swigobj_attr('diag')


class DftNear2Far(DftObj):

    def __init__(self, func, args):
        super(DftNear2Far, self).__init__(func, args)
        self.nfreqs = args[2]
        self.regions = args[3]
        self.num_components = 4

    @property
    def farfield(self):
        return self.swigobj_attr('farfield')

    @property
    def save_farfields(self):
        return self.swigobj_attr('save_farfields')

    @property
    def F(self):
        return self.swigobj_attr('F')

    @property
    def eps(self):
        return self.swigobj_attr('eps')

    @property
    def mu(self):
        return self.swigobj_attr('mu')


class DftFields(DftObj):

    def __init__(self, func, args):
        super(DftFields, self).__init__(func, args)
        self.nfreqs = args[6]
        self.regions = [FieldsRegion(where=args[1], center=args[2], size=args[3])]
        self.num_components = len(args[0])

    @property
    def chunks(self):
        return self.swigobj_attr('chunks')


Mode = namedtuple('Mode', ['freq', 'decay', 'Q', 'amp', 'err'])


class EigenmodeData(object):

    def __init__(self, band_num, omega, group_velocity, k, swigobj):
        self.band_num = band_num
        self.omega = omega
        self.group_velocity = group_velocity
        self.k = k
        self.swigobj = swigobj

    def amplitude(self, point, component):
        swig_point = mp.vec(point.x, point.y, point.z)
        return mp.eigenmode_amplitude(self.swigobj, swig_point, component)


class Harminv(object):

    def __init__(self, c, pt, fcen, df, mxbands=None):
        self.c = c
        self.pt = pt
        self.fcen = fcen
        self.df = df
        self.mxbands = mxbands
        self.data = []
        self.data_dt = 0
        self.modes = []
        self.spectral_density = 1.1
        self.Q_thresh = 50.0
        self.rel_err_thresh = mp.inf
        self.err_thresh = 0.01
        self.rel_amp_thresh = -1.0
        self.amp_thresh = -1.0
        self.step_func = self._harminv()

    def __call__(self, sim, todo):
        self.step_func(sim, todo)

    def _collect_harminv(self):

        def _collect1(c, pt):
            self.t0 = 0

            def _collect2(sim):
                self.data_dt = sim.meep_time() - self.t0
                self.t0 = sim.meep_time()
                self.data.append(sim.get_field_point(c, pt))
            return _collect2
        return _collect1

    def _check_freqs(self, sim):
        source_freqs = [(s.src.frequency, 0 if s.src.width == 0 else 1 / s.src.width)
                        for s in sim.sources
                        if hasattr(s.src, 'frequency')]

        harminv_max = self.fcen + 0.5 * self.df
        harminv_min = self.fcen - 0.5 * self.df

        for sf in source_freqs:
            sf_max = sf[0] + 0.5 * sf[1]
            sf_min = sf[0] - 0.5 * sf[1]
            if harminv_max > sf_max:
                warn_fmt = "Harminv frequency {} is outside maximum Source frequency {}"
                warnings.warn(warn_fmt.format(harminv_max, sf_max), RuntimeWarning)
            if harminv_min < sf_min:
                warn_fmt = "Harminv frequency {} is outside minimum Source frequency {}"
                warnings.warn(warn_fmt.format(harminv_min, sf_min), RuntimeWarning)

    def _analyze_harminv(self, sim, maxbands):
        harminv_cols = ['frequency', 'imag. freq.', 'Q', '|amp|', 'amplitude', 'error']
        display_run_data(sim, 'harminv', harminv_cols)
        self._check_freqs(sim)

        dt = self.data_dt if self.data_dt is not None else sim.fields.dt

        bands = mp.py_do_harminv(self.data, dt, self.fcen - self.df / 2, self.fcen + self.df / 2, maxbands,
                                 self.spectral_density, self.Q_thresh, self.rel_err_thresh, self.err_thresh,
                                 self.rel_amp_thresh, self.amp_thresh)

        modes = []
        for freq, amp, err in bands:
            Q = freq.real / (-2 * freq.imag)
            modes.append(Mode(freq.real, freq.imag, Q, amp, err))
            display_run_data(sim, 'harminv', [freq.real, freq.imag, Q, abs(amp), amp, err])

        return modes

    def _harminv(self):

        def _harm(sim):

            if self.mxbands is None or self.mxbands == 0:
                mb = 100
            else:
                mb = self.mxbands
            self.modes = self._analyze_harminv(sim, mb)

        f1 = self._collect_harminv()

        return _combine_step_funcs(at_end(_harm), f1(self.c, self.pt))


class Simulation(object):

    def __init__(self,
                 cell_size,
                 resolution,
                 geometry=[],
                 sources=[],
                 eps_averaging=True,
                 dimensions=3,
                 boundary_layers=[],
                 symmetries=[],
                 verbose=False,
                 force_complex_fields=False,
                 default_material=mp.Medium(),
                 m=0,
                 k_point=False,
                 extra_materials=[],
                 material_function=None,
                 epsilon_func=None,
                 epsilon_input_file='',
                 progress_interval=4,
                 subpixel_tol=1e-4,
                 subpixel_maxeval=100000,
                 ensure_periodicity=True,
                 num_chunks=0,
                 Courant=0.5,
                 accurate_fields_near_cylorigin=False,
                 filename_prefix='',
                 output_volume=None,
                 output_single_precision=False,
                 load_structure=''):

        self.cell_size = cell_size
        self.geometry = geometry
        self.sources = sources
        self.resolution = resolution
        self.dimensions = dimensions
        self.boundary_layers = boundary_layers
        self.symmetries = symmetries
        self.geometry_center = Vector3()
        self.eps_averaging = eps_averaging
        self.subpixel_tol = subpixel_tol
        self.subpixel_maxeval = subpixel_maxeval
        self.ensure_periodicity = ensure_periodicity
        self.extra_materials = extra_materials
        self.default_material = default_material
        self.epsilon_input_file = epsilon_input_file
        self.num_chunks = num_chunks
        self.Courant = Courant
        self.global_d_conductivity = 0
        self.global_b_conductivity = 0
        self.special_kz = False
        self.k_point = k_point
        self.fields = None
        self.structure = None
        self.accurate_fields_near_cylorigin = accurate_fields_near_cylorigin
        self.m = m
        self.force_complex_fields = force_complex_fields
        self.verbose = verbose
        self.progress_interval = progress_interval
        self.init_sim_hooks = []
        self.run_index = 0
        self.filename_prefix = filename_prefix
        self.output_append_h5 = None
        self.output_single_precision = output_single_precision
        self.output_volume = output_volume
        self.last_eps_filename = ''
        self.output_h5_hook = lambda fname: False
        self.interactive = False
        self.is_cylindrical = False
        self.material_function = material_function
        self.epsilon_func = epsilon_func
        self.load_structure_file = load_structure
        self.dft_objects = []

    # To prevent the user from having to specify `dims` and `is_cylindrical`
    # to Volumes they create, the library will adjust them appropriately based
    # on the settings in the Simulation instance. This method must be called on
    # any user-defined Volume before passing it to meep via its `swigobj`.
    def _fit_volume_to_simulation(self, vol):
        return Volume(vol.center, vol.size, dims=self.dimensions, is_cylindrical=self.is_cylindrical)

    # Every function that takes a user volume can be specified either by a volume
    # (a Python Volume or a SWIG-wrapped meep::volume), or a center and a size
    def _volume_from_kwargs(self, vol=None, center=None, size=None):
        if vol:
            if isinstance(vol, Volume):
                # A pure Python Volume
                return self._fit_volume_to_simulation(vol).swigobj
            else:
                # A SWIG-wrapped meep::volume
                return vol
        elif size and center:
            return Volume(center, size=size, dims=self.dimensions, is_cylindrical=self.is_cylindrical).swigobj
        else:
            raise ValueError("Need either a Volume, or a size and center")

    def _infer_dimensions(self, k):
        if self.dimensions == 3:

            def use_2d(self, k):
                zero_z = self.cell_size.z == 0
                return zero_z and (not k or self.special_kz or k.z == 0)

            if use_2d(self, k):
                return 2
            else:
                return 3
        return self.dimensions

    def _get_valid_material_frequencies(self):
        fmin = float('-inf')
        fmax = float('inf')

        for mat in [go.material for go in self.geometry] + self.extra_materials:
            if isinstance(mat, mp.Medium) and mat.valid_freq_range:
                if mat.valid_freq_range.min > fmin:
                    fmin = mat.valid_freq_range.min
                if mat.valid_freq_range.max < fmax:
                    fmax = mat.valid_freq_range.max

        return fmin, fmax

    def _check_material_frequencies(self):

        min_freq, max_freq = self._get_valid_material_frequencies()

        source_freqs = [(s.src.frequency, 0 if s.src.width == 0 else 1 / s.src.width)
                        for s in self.sources
                        if hasattr(s.src, 'frequency')]

        dft_freqs = []
        for dftf in self.dft_objects:
            dft_freqs.append(dftf.freq_min)
            dft_freqs.append(dftf.freq_min + dftf.Nfreq * dftf.dfreq)

        warn_src = ('Note: your sources include frequencies outside the range of validity of the ' +
                    'material models. This is fine as long as you eventually only look at outputs ' +
                    '(fluxes, resonant modes, etc.) at valid frequencies.')

        warn_dft_fmt = "DFT frequency {} is out of material's range of {}-{}"

        for sf in source_freqs:
            if sf[0] + 0.5 * sf[1] > max_freq or sf[0] - 0.5 * sf[1] < min_freq:
                warnings.warn(warn_src, RuntimeWarning)

        for dftf in dft_freqs:
            if dftf > max_freq or dftf < min_freq:
                warnings.warn(warn_dft_fmt.format(dftf, min_freq, max_freq), RuntimeWarning)

    def _create_grid_volume(self, k):
        dims = self._infer_dimensions(k)

        if dims == 0 or dims == 1:
            gv = mp.vol1d(self.cell_size.z, self.resolution)
        elif dims == 2:
            self.dimensions = 2
            gv = mp.vol2d(self.cell_size.x, self.cell_size.y, self.resolution)
        elif dims == 3:
            gv = mp.vol3d(self.cell_size.x, self.cell_size.y, self.cell_size.z, self.resolution)
        elif dims == mp.CYLINDRICAL:
            gv = mp.volcyl(self.cell_size.x, self.cell_size.z, self.resolution)
            self.dimensions = 2
            self.is_cylindrical = True
        else:
            raise ValueError("Unsupported dimentionality: {}".format(dims))

        gv.center_origin()
        return gv

    def _create_symmetries(self, gv):
        sym = mp.symmetry()

        # Initialize swig objects for each symmetry and combine them into one
        for s in self.symmetries:
            if isinstance(s, Identity):
                s.swigobj = mp.identity()
            elif isinstance(s, Rotate2):
                s.swigobj = mp.rotate2(s.direction, gv)
                sym += s.swigobj * complex(s.phase.real, s.phase.imag)
            elif isinstance(s, Rotate4):
                s.swigobj = mp.rotate4(s.direction, gv)
                sym += s.swigobj * complex(s.phase.real, s.phase.imag)
            elif isinstance(s, Mirror):
                s.swigobj = mp.mirror(s.direction, gv)
                sym += s.swigobj * complex(s.phase.real, s.phase.imag)
            else:
                s.swigobj = mp.symmetry()

        return sym

    def _get_dft_volumes(self):
        volumes = [self._volume_from_kwargs(vol=r.where if hasattr(r, 'where') else None,
                                            center=r.center, size=r.size)
                   for dft in self.dft_objects
                   for r in dft.regions]

        return volumes

    def _compute_fragment_stats(self, gv):

        def convert_volumes(dft_obj):
            volumes = []
            for r in dft_obj.regions:
                volumes.append(self._volume_from_kwargs(vol=r.where if hasattr(r, 'where') else None,
                                                        center=r.center, size=r.size))
            return volumes

        dft_data_list = [mp.dft_data(o.nfreqs, o.num_components, convert_volumes(o))
                         for o in self.dft_objects]

        stats = mp.compute_fragment_stats(
            self.geometry,
            gv,
            self.cell_size,
            mp.Vector3(),
            self.default_material,
            dft_data_list,
            self.subpixel_tol,
            self.subpixel_maxeval,
            self.ensure_periodicity
        )

        return stats

    def _init_structure(self, k=False):
        print('-' * 11)
        print('Initializing structure...')

        gv = self._create_grid_volume(k)
        sym = self._create_symmetries(gv)
        br = _create_boundary_region_from_boundary_layers(self.boundary_layers, gv)

        if self.boundary_layers and type(self.boundary_layers[0]) is Absorber:
            absorbers = self.boundary_layers
        else:
            absorbers = None

        if self.material_function:
            self.material_function.eps = False
            self.default_material = self.material_function
        elif self.epsilon_func:
            self.epsilon_func.eps = True
            self.default_material = self.epsilon_func
        elif self.epsilon_input_file:
            self.default_material = self.epsilon_input_file

        self.fragment_stats = self._compute_fragment_stats(gv)

        self.structure = mp.structure(gv, None, br, sym, self.num_chunks, self.Courant,
                                      self.eps_averaging, self.subpixel_tol, self.subpixel_maxeval)
        self.structure.shared_chunks = True

        mp.set_materials_from_geometry(self.structure, self.geometry, self.eps_averaging, self.subpixel_tol,
                                       self.subpixel_maxeval, self.ensure_periodicity and not not self.k_point,
                                       False, self.default_material, absorbers, self.extra_materials)
        if self.load_structure_file:
            self.load_structure(self.load_structure_file)

    def set_materials(self, geometry=None, default_material=None):
        if self.fields:
            self.fields.remove_susceptibilities()

        have_absorbers = self.boundary_layers and type(self.boundary_layers[0]) is Absorber

        mp.set_materials_from_geometry(
            self.structure,
            geometry if geometry is not None else self.geometry,
            self.eps_averaging,
            self.subpixel_tol,
            self.subpixel_maxeval,
            self.ensure_periodicity,
            False,
            default_material if default_material else self.default_material,
            self.boundary_layers if have_absorbers else None,
            self.extra_materials
        )

    def load_structure(self, fname):
        if self.structure is None:
            raise ValueError("Fields must be initialized before calling load_structure")
        self.structure.load(fname)

    def dump_structure(self, fname):
        if self.structure is None:
            raise ValueError("Fields must be initialized before calling dump_structure")
        self.structure.dump(fname)

    def init_sim(self):

        materials = [g.material for g in self.geometry if isinstance(g.material, mp.Medium)]
        if isinstance(self.default_material, mp.Medium):
            materials.append(self.default_material)
        for med in materials:
            if ((med.epsilon_diag.x < 1 and med.epsilon_diag.x > -mp.inf) or
                (med.epsilon_diag.y < 1 and med.epsilon_diag.y > -mp.inf) or
                (med.epsilon_diag.z < 1 and med.epsilon_diag.z > -mp.inf)):

                eps_warning = ("Epsilon < 1 may require adjusting the Courant parameter. " +
                               "See the 'Numerical Stability' entry under the 'Materials' " +
                               "section of the documentation")
                warnings.warn(eps_warning, RuntimeWarning)

        if self.structure is None:
            self._init_structure(self.k_point)

        self.fields = mp.fields(
            self.structure,
            self.m if self.is_cylindrical else 0,
            self.k_point.z if self.special_kz and self.k_point else 0,
            not self.accurate_fields_near_cylorigin
        )

        if self.verbose:
            self.fields.verbose()

        def use_real(self):
            cond1 = self.is_cylindrical and self.m != 0
            cond2 = any([s.phase.imag for s in self.symmetries])
            cond3 = not self.k_point
            cond4 = self.special_kz and self.k_point.x == 0 and self.k_point.y == 0
            cond5 = not (cond3 or cond4 or self.k_point == Vector3())
            return not (self.force_complex_fields or cond1 or cond2 or cond5)

        if use_real(self):
            self.fields.use_real_fields()
        else:
            print("Meep: using complex fields.")

        if self.k_point:
            v = Vector3(self.k_point.x, self.k_point.y) if self.special_kz else self.k_point
            self.fields.use_bloch(py_v3_to_vec(self.dimensions, v, self.is_cylindrical))

        for s in self.sources:
            self.add_source(s)

        for hook in self.init_sim_hooks:
            hook()

    def init_fields(self):
        warnings.warn('init_fields is deprecated. Please use init_sim instead', DeprecationWarning)
        self.init_sim()

    def require_dimensions(self):
        if self.structure is None:
            mp.set_dimensions(self._infer_dimensions(self.k_point))

    def meep_time(self):
        if self.fields is None:
            self.init_sim()
        return self.fields.time()

    def round_time(self):
        if self.fields is None:
            self.init_sim()

        return self.fields.round_time()

    def get_field_point(self, c, pt):
        v3 = py_v3_to_vec(self.dimensions, pt, self.is_cylindrical)
        return self.fields.get_field_from_comp(c, v3)

    def get_epsilon_point(self, pt):
        v3 = py_v3_to_vec(self.dimensions, pt, self.is_cylindrical)
        return self.fields.get_eps(v3)

    def get_filename_prefix(self):
        if self.filename_prefix:
            return self.filename_prefix
        else:
            _, filename = os.path.split(sys.argv[0])

            if filename == 'ipykernel_launcher.py' or filename == '__main__.py':
                return ''
            else:
                return re.sub(r'\.py$', '', filename)

    def use_output_directory(self, dname=''):
        if not dname:
            dname = self.get_filename_prefix() + '-out'

        closure = {'trashed': False}

        def hook():
            print("Meep: using output directory '{}'".format(dname))
            self.fields.set_output_directory(dname)
            if not closure['trashed']:
                mp.trash_output_directory(dname)
            closure['trashed'] = True

        self.init_sim_hooks.append(hook)

        if self.fields is not None:
            hook()
        self.filename_prefix = ''

        return dname

    def _run_until(self, cond, step_funcs):
        self.interactive = False
        if self.fields is None:
            self.init_sim()

        if isinstance(cond, numbers.Number):
            stop_time = cond
            t0 = self.round_time()

            def stop_cond(sim):
                return sim.round_time() >= t0 + stop_time

            cond = stop_cond

            step_funcs = list(step_funcs)
            step_funcs.append(display_progress(t0, t0 + stop_time, self.progress_interval))

        while not cond(self):
            for func in step_funcs:
                _eval_step_func(self, func, 'step')
            self.fields.step()

        # Translating the recursive scheme version of run-until into an iterative version
        # (because python isn't tail-call-optimized) means we need one extra iteration to
        # be the same as scheme.
        for func in step_funcs:
            _eval_step_func(self, func, 'step')

        for func in step_funcs:
            _eval_step_func(self, func, 'finish')

        print("run {} finished at t = {} ({} timesteps)".format(self.run_index, self.meep_time(), self.fields.t))
        self.run_index += 1

    def _run_sources_until(self, cond, step_funcs):
        if self.fields is None:
            self.init_sim()

        ts = self.fields.last_source_time()

        if isinstance(cond, numbers.Number):
            new_cond = (ts - self.round_time()) + cond
        else:
            def f(sim):
                return cond(sim) and sim.round_time() >= ts
            new_cond = f

        self._run_until(new_cond, step_funcs)

    def _run_sources(self, step_funcs):
        self._run_sources_until(self, 0, step_funcs)

    def _run_k_point(self, t, k):
        components = [s.component for s in self.sources]
        pts = [s.center for s in self.sources]

        src_freqs_min = min([s.src.frequency - 1 / s.src.width / 2 if isinstance(s.src, mp.GaussianSource) else mp.inf
                             for s in self.sources])
        fmin = max(0, src_freqs_min)

        fmax = max([s.src.frequency + 1 / s.src.width / 2 if isinstance(s.src, mp.GaussianSource) else 0
                    for s in self.sources])

        if not components or fmin > fmax:
            raise ValueError("Running with k_points requires a 'GaussianSource' source")

        self.change_k_point(k)
        self.restart_fields()

        h = Harminv(components[0], pts[0], 0.5 * (fmin + fmax), fmax - fmin)
        self.run(after_sources(h), until_after_sources=t)

        return [complex(m.freq, m.decay) for m in h.modes]

    def run_k_points(self, t, k_points):
        k_index = 0
        all_freqs = []

        for k in k_points:
            k_index += 1

            if k_index == 1:
                self.init_sim()
                output_epsilon(self)

            freqs = self._run_k_point(t, k)

            print("freqs:, {}, {}, {}, {}, ".format(k_index, k.x, k.y, k.z), end='')
            print(', '.join([str(f.real) for f in freqs]))
            print("freqs-im:, {}, {}, {}, {}, ".format(k_index, k.x, k.y, k.z), end='')
            print(', '.join([str(f.imag) for f in freqs]))

            all_freqs.append(freqs)

        return all_freqs

    def set_epsilon(self, eps):
        if self.fields is None:
            self.init_sim()

        self.structure.set_epsilon(eps, self.eps_averaging, self.subpixel_tol, self.subpixel_maxeval)

    def add_source(self, src):
        if self.fields is None:
            self.init_sim()

        where = Volume(src.center, src.size, dims=self.dimensions,
                       is_cylindrical=self.is_cylindrical).swigobj

        if isinstance(src, EigenModeSource):
            if src.direction < 0:
                direction = self.fields.normal_direction(where)
            else:
                direction = src.direction

            eig_vol = Volume(src.eig_lattice_center, src.eig_lattice_size, self.dimensions,
                             is_cylindrical=self.is_cylindrical).swigobj

            add_eig_src_args = [
                src.component,
                src.src.swigobj,
                direction,
                where,
                eig_vol,
                src.eig_band,
                py_v3_to_vec(self.dimensions, src.eig_kpoint, is_cylindrical=self.is_cylindrical),
                src.eig_match_freq,
                src.eig_parity,
                src.eig_resolution,
                src.eig_tolerance,
                src.amplitude
            ]
            add_eig_src = functools.partial(self.fields.add_eigenmode_source, *add_eig_src_args)

            if src.amp_func is None:
                add_eig_src()
            else:
                add_eig_src(src.amp_func)
        else:
            add_vol_src_args = [src.component, src.src.swigobj, where]
            add_vol_src = functools.partial(self.fields.add_volume_source, *add_vol_src_args)

            if src.amp_func_file:
                fname_dset = src.amp_func_file.rsplit(':', 1)
                if len(fname_dset) != 2:
                    err_msg = "Expected a string of the form 'h5filename:dataset'. Got '{}'"
                    raise ValueError(err_msg.format(src.amp_func_file))

                fname, dset = fname_dset
                if not fname.endswith('.h5'):
                    fname += '.h5'

                add_vol_src(fname, dset, src.amplitude * 1.0,)
            elif src.amp_func:
                add_vol_src(src.amp_func, src.amplitude * 1.0)
            elif src.amp_data is not None:
                add_vol_src(src.amp_data, src.amplitude * 1.0,)
            else:
                add_vol_src(src.amplitude * 1.0)

    def _evaluate_dft_objects(self):
        for dft in self.dft_objects:
            if dft.swigobj is None:
                dft.swigobj = dft.func(*dft.args)

    def add_dft_fields(self, components, freq_min, freq_max, nfreq, where=None, center=None, size=None):
        dftf = DftFields(self._add_dft_fields, [components, where, center, size, freq_min, freq_max, nfreq])
        self.dft_objects.append(dftf)
        return dftf

    def _add_dft_fields(self, components, where, center, size, freq_min, freq_max, nfreq):
        if self.fields is None:
            self.init_sim()
        try:
            where = self._volume_from_kwargs(where, center, size)
        except ValueError:
            where = self.fields.total_volume()

        return self.fields.add_dft_fields(components, where, freq_min, freq_max, nfreq)

    def output_dft(self, dft_fields, fname):
        if self.fields is None:
            self.init_sim()

        if hasattr(dft_fields, 'swigobj'):
            dft_fields_swigobj = dft_fields.swigobj
        else:
            dft_fields_swigobj = dft_fields

        self.fields.output_dft(dft_fields_swigobj, fname)

    def get_dft_data(self, dft_chunk):
        n = mp._get_dft_data_size(dft_chunk)
        arr = np.zeros(n, np.complex128)
        mp._get_dft_data(dft_chunk, arr)
        return arr

    def add_near2far(self, fcen, df, nfreq, *near2fars):
        n2f = DftNear2Far(self._add_near2far, [fcen, df, nfreq, near2fars])
        self.dft_objects.append(n2f)
        return n2f

    def _add_near2far(self, fcen, df, nfreq, near2fars):
        if self.fields is None:
            self.init_sim()
        return self._add_fluxish_stuff(self.fields.add_dft_near2far, fcen, df, nfreq, near2fars)

    def get_farfield(self, f, v):
        return mp._get_farfield(f.swigobj, py_v3_to_vec(self.dimensions, v, is_cylindrical=self.is_cylindrical))

    def output_farfields(self, near2far, fname, resolution, where=None, center=None, size=None):
        vol = self._volume_from_kwargs(where, center, size)
        near2far.save_farfields(fname, self.get_filename_prefix(), vol, resolution)

    def load_near2far(self, fname, n2f):
        if self.fields is None:
            self.init_sim()
        n2f.load_hdf5(self.fields, fname, '', self.get_filename_prefix())

    def save_near2far(self, fname, n2f):
        if self.fields is None:
            self.init_sim()
        n2f.save_hdf5(self.fields, fname, '', self.get_filename_prefix())

    def load_minus_near2far(self, fname, n2f):
        self.load_near2far(fname, n2f)
        n2f.scale_dfts(-1.0)

    def get_near2far_data(self, n2f):
        return NearToFarData(F=self.get_dft_data(n2f.F))

    def load_near2far_data(self, n2f, n2fdata):
        mp._load_dft_data(n2f.F, n2fdata.F)

    def load_minus_near2far_data(self, n2f, n2fdata):
        self.load_near2far_data(n2f, n2fdata)
        n2f.scale_dfts(complex(1.0))

    def add_force(self, fcen, df, nfreq, *forces):
        force = DftForce(self._add_force, [fcen, df, nfreq, forces])
        self.dft_objects.append(force)
        return force

    def _add_force(self, fcen, df, nfreq, forces):
        if self.fields is None:
            self.init_sim()
        return self._add_fluxish_stuff(self.fields.add_dft_force, fcen, df, nfreq, forces)

    def display_forces(self, *forces):
        force_freqs = get_force_freqs(forces[0])
        display_csv(self, 'force', zip(force_freqs, *[get_forces(f) for f in forces]))

    def load_force(self, fname, force):
        if self.fields is None:
            self.init_sim()
        force.load_hdf5(self.fields, fname, '', self.get_filename_prefix())

    def save_force(self, fname, force):
        if self.fields is None:
            self.init_sim()
        force.save_hdf5(self.fields, fname, '', self.get_filename_prefix())

    def load_minus_force(self, fname, force):
        self.load_force(fname, force)
        force.scale_dfts(-1.0)

    def get_force_data(self, force):
        return ForceData(offdiag1=self.get_dft_data(force.offdiag1),
                         offdiag2=self.get_dft_data(force.offdiag2),
                         diag=self.get_dft_data(force.diag))

    def load_force_data(self, force, fdata):
        mp._load_dft_data(force.offdiag1, fdata.offdiag1)
        mp._load_dft_data(force.offdiag2, fdata.offdiag2)
        mp._load_dft_data(force.diag, fdata.diag)

    def load_minus_force_data(self, force, fdata):
        self.load_force_data(force, fdata)
        force.scale_dfts(complex(-1.0))

    def add_flux(self, fcen, df, nfreq, *fluxes):
        flux = DftFlux(self._add_flux, [fcen, df, nfreq, fluxes])
        self.dft_objects.append(flux)
        return flux

    def _add_flux(self, fcen, df, nfreq, fluxes):
        if self.fields is None:
            self.init_sim()
        return self._add_fluxish_stuff(self.fields.add_dft_flux, fcen, df, nfreq, fluxes)

    def add_mode_monitor(self, fcen, df, nfreq, *fluxes):
        flux = DftFlux(self._add_mode_monitor, [fcen, df, nfreq, fluxes])
        self.dft_objects.append(flux)
        return flux

    def _add_mode_monitor(self, fcen, df, nfreq, fluxes):
        if self.fields is None:
            self.init_sim()

        if len(fluxes) != 1:
            raise ValueError("add_mode_monitor expected just one ModeRegion. Got {}".format(len(fluxes)))

        region = fluxes[0]
        v = mp.Volume(region.center, region.size, dims=self.dimensions, is_cylindrical=self.is_cylindrical)
        d0 = region.direction
        d = self.fields.normal_direction(v.swigobj) if d0 < 0 else d0

        return self.fields.add_mode_monitor(d, v.swigobj, fcen - df / 2, fcen + df / 2, nfreq)

    def add_eigenmode(self, fcen, df, nfreq, *fluxes):
        warnings.warn('add_eigenmode is deprecated. Please use add_mode_monitor instead.', DeprecationWarning)
        return self.add_mode_monitor(fcen, df, nfreq, *fluxes)

    def display_fluxes(self, *fluxes):
        display_csv(self, 'flux', zip(get_flux_freqs(fluxes[0]), *[get_fluxes(f) for f in fluxes]))

    def load_flux(self, fname, flux):
        if self.fields is None:
            self.init_sim()

        flux.load_hdf5(self.fields, fname, '', self.get_filename_prefix())

    def save_flux(self, fname, flux):
        if self.fields is None:
            self.init_sim()

        flux.save_hdf5(self.fields, fname, '', self.get_filename_prefix())

    def load_minus_flux(self, fname, flux):
        self.load_flux(fname, flux)
        flux.scale_dfts(complex(-1.0))

    def get_flux_data(self, flux):
        return FluxData(E=self.get_dft_data(flux.E), H=self.get_dft_data(flux.H))

    def load_flux_data(self, flux, fdata):
        mp._load_dft_data(flux.E, fdata.E)
        mp._load_dft_data(flux.H, fdata.H)

    def load_minus_flux_data(self, flux, fdata):
        self.load_flux_data(flux, fdata)
        flux.scale_dfts(complex(-1.0))

    def flux_in_box(self, d, box=None, center=None, size=None):
        if self.fields is None:
            raise RuntimeError('Fields must be initialized before using flux_in_box')

        box = self._volume_from_kwargs(box, center, size)

        return self.fields.flux_in_box(d, box)

    def electric_energy_in_box(self, box=None, center=None, size=None):
        if self.fields is None:
            raise RuntimeError('Fields must be initialized before using electric_energy_in_box')

        box = self._volume_from_kwargs(box, center, size)

        return self.fields.electric_energy_in_box(box)

    def magnetic_energy_in_box(self, box=None, center=None, size=None):
        if self.fields is None:
            raise RuntimeError('Fields must be initialized before using magnetic_energy_in_box')

        box = self._volume_from_kwargs(box, center, size)

        return self.fields.magnetic_energy_in_box(box)

    def field_energy_in_box(self, box=None, center=None, size=None):
        if self.fields is None:
            raise RuntimeError('Fields must be initialized before using field_energy_in_box')

        box = self._volume_from_kwargs(box, center, size)

        return self.fields.field_energy_in_box(box)

    def modal_volume_in_box(self, box=None, center=None, size=None):
        if self.fields is None:
            raise RuntimeError('Fields must be initialized before using modal_volume_in_box')

        try:
            box = self._volume_from_kwargs(box, center, size)
        except ValueError:
            box = self.fields.total_volume()

        return self.fields.modal_volume_in_box(box)

    def solve_cw(self, tol=1e-8, maxiters=10000, L=2):
        if self.fields is None:
            raise RuntimeError('Fields must be initialized before using solve_cw')
        return self.fields.solve_cw(tol, maxiters, L)

    def _add_fluxish_stuff(self, add_dft_stuff, fcen, df, nfreq, stufflist):
        vol_list = None

        for s in stufflist:
            v = Volume(center=s.center, size=s.size, dims=self.dimensions,
                       is_cylindrical=self.is_cylindrical)
            d0 = s.direction
            d = self.fields.normal_direction(v.swigobj) if d0 < 0 else d0
            c = mp.direction_component(mp.Sx, d)
            v2 = Volume(center=s.center, size=s.size, dims=self.dimensions,
                        is_cylindrical=self.is_cylindrical).swigobj
            vol_list = mp.make_volume_list(v2, c, s.weight, vol_list)

        stuff = add_dft_stuff(vol_list, fcen - df / 2, fcen + df / 2, nfreq)
        vol_list.__swig_destroy__(vol_list)

        return stuff

    def output_component(self, c, h5file=None):
        if self.fields is None:
            raise RuntimeError("Fields must be initialized before calling output_component")

        vol = self.fields.total_volume() if self.output_volume is None else self.output_volume
        h5 = self.output_append_h5 if h5file is None else h5file
        append = h5file is None and self.output_append_h5 is not None

        self.fields.output_hdf5(c, vol, h5, append, self.output_single_precision, self.get_filename_prefix())

        if h5file is None:
            nm = self.fields.h5file_name(mp.component_name(c), self.get_filename_prefix(), True)
            if c == mp.Dielectric:
                self.last_eps_filename = nm
            self.output_h5_hook(nm)

    def output_components(self, fname, *components):
        if self.fields is None:
            raise RuntimeError("Fields must be initialized before calling output_component")

        if self.output_append_h5 is None:
            f = self.fields.open_h5file(fname, mp.h5file.WRITE, self.get_filename_prefix(), True)
        else:
            f = None

        for c in components:
            self.output_component(c, h5file=f)
            if self.output_append_h5 is None:
                f.prevent_deadlock()

        if self.output_append_h5 is None:
            self.output_h5_hook(self.fields.h5file_name(fname, self.get_filename_prefix(), True))

    def h5topng(self, rm_h5, option, *step_funcs):
        opts = "h5topng {}".format(option)
        cmd = re.sub(r'\$EPS', self.last_eps_filename, opts)
        return convert_h5(rm_h5, cmd, *step_funcs)

    def get_array(self, vol=None, center=None, size=None, component=mp.Ez, cmplx=None, arr=None):
        dim_sizes = np.zeros(3, dtype=np.uintp)

        if vol is None and center is None and size is None:
            v = self.fields.total_volume()
        else:
            v = self._volume_from_kwargs(vol, center, size)

        self.fields.get_array_slice_dimensions(v, dim_sizes)

        dims = [s for s in dim_sizes if s != 0]

        if cmplx is None:
            cmplx = component < mp.Dielectric and not self.fields.is_real

        if arr is not None:
            if cmplx and not np.iscomplexobj(arr):
                raise ValueError("Requested a complex slice, but provided array of type {}.".format(arr.dtype))

            for a, b in zip(arr.shape, dims):
                if a != b:
                    fmt = "Expected dimensions {}, but got {}"
                    raise ValueError(fmt.format(dims, arr.shape))

            arr = np.require(arr, requirements=['C', 'W'])

        else:
            arr = np.zeros(dims, dtype=np.complex128 if cmplx else np.float64)

        if np.iscomplexobj(arr):
            self.fields.get_complex_array_slice(v, component, arr)
        else:
            self.fields.get_array_slice(v, component, arr)

        return arr

    def get_dft_array(self, dft_obj, component, num_freq):
        if hasattr(dft_obj, 'swigobj'):
            dft_swigobj = dft_obj.swigobj
        else:
            dft_swigobj = dft_obj

        if type(dft_swigobj) is mp.dft_fields:
            return mp.get_dft_fields_array(self.fields, dft_swigobj, component, num_freq)
        elif type(dft_swigobj) is mp.dft_flux:
            return mp.get_dft_flux_array(self.fields, dft_swigobj, component, num_freq)
        elif type(dft_swigobj) is mp.dft_force:
            return mp.get_dft_force_array(self.fields, dft_swigobj, component, num_freq)
        elif type(dft_swigobj) is mp.dft_near2far:
            return mp.get_dft_near2far_array(self.fields, dft_swigobj, component, num_freq)
        else:
            raise ValueError("Invalid type of dft object: {}".format(dft_swigobj))

    def get_eigenmode_coefficients(self, flux, bands, eig_parity=mp.NO_PARITY,
                                   eig_vol=None, eig_resolution=0, eig_tolerance=1e-7, kpoint_func=None):
        if self.fields is None:
            raise ValueError("Fields must be initialized before calling get_eigenmode_coefficients")
        if eig_vol is None:
            eig_vol = flux.where
        else:
            eig_vol = self._volume_from_kwargs(vol=eig_vol)

        num_bands = len(bands)
        coeffs = np.zeros(2 * num_bands * flux.Nfreq, dtype=np.complex128)
        vgrp = np.zeros(num_bands * flux.Nfreq)

        kpoints = mp.get_eigenmode_coefficients_and_kpoints(
            self.fields,
            flux.swigobj,
            eig_vol,
            np.array(bands, dtype=np.intc),
            eig_parity,
            eig_resolution,
            eig_tolerance,
            coeffs,
            vgrp,
            kpoint_func
        )

        return np.reshape(coeffs, (num_bands, flux.Nfreq, 2)), vgrp, kpoints

    def get_eigenmode(self, omega_src, direction, where, band_num, kpoint, eig_vol=None, match_frequency=True,
                      parity=mp.NO_PARITY, resolution=0, eigensolver_tol=1e-7, verbose=False):

        if self.fields is None:
            raise ValueError("Fields must be initialized before calling get_eigenmode")

        where = self._volume_from_kwargs(vol=where)
        if eig_vol is None:
            eig_vol = where
        else:
            eig_vol = self._volume_from_kwargs(vol=eig_vol)

        swig_kpoint = mp.vec(kpoint.x, kpoint.y, kpoint.z)
        emdata = mp._get_eigenmode(self.fields, omega_src, direction, where, eig_vol, band_num, swig_kpoint,
                                   match_frequency, parity, resolution, eigensolver_tol, verbose)

        return EigenmodeData(emdata.band_num, emdata.omega, emdata.group_velocity, emdata.Gk, emdata.data)

    def output_field_function(self, name, cs, func, real_only=False, h5file=None):
        if self.fields is None:
            raise RuntimeError("Fields must be initialized before calling output_field_function")

        ov = self.output_volume if self.output_volume else self.fields.total_volume()
        h5 = self.output_append_h5 if h5file is None else h5file
        append = h5file is None and self.output_append_h5 is not None

        self.fields.output_hdf5(name, [cs, func], ov, h5, append, self.output_single_precision,
                                self.get_filename_prefix(), real_only)
        if h5file is None:
            self.output_h5_hook(self.fields.h5file_name(name, self.get_filename_prefix(), True))

    def _get_field_function_volume(self, where=None, center=None, size=None):
        try:
            where = self._volume_from_kwargs(where, center, size)
        except ValueError:
            where = self.fields.total_volume()

        return where

    def integrate_field_function(self, cs, func, where=None, center=None, size=None):
        where = self._get_field_function_volume(where, center, size)
        return self.fields.integrate([cs, func], where)

    def integrate2_field_function(self, fields2, cs1, cs2, func, where=None, center=None, size=None):
        where = self._get_field_function_volume(where, center, size)
        return self.fields.integrate2(fields2, [cs1, cs2, func], where)

    def max_abs_field_function(self, cs, func, where=None, center=None, size=None):
        where = self._get_field_function_volume(where, center, size)
        return self.fields.max_abs([cs, func], where)

    def change_k_point(self, k):
        self.k_point = k

        if self.fields:
            needs_complex_fields = not (not self.k_point or self.k_point == mp.Vector3())

            if needs_complex_fields and self.fields.is_real:
                self.fields = None
                self.init_sim()
            else:
                if self.k_point:
                    self.fields.use_bloch(py_v3_to_vec(self.dimensions, self.k_point, self.is_cylindrical))

    def change_sources(self, new_sources):
        self.sources = new_sources
        if self.fields:
            self.fields.remove_sources()
            for s in self.sources:
                self.add_source(s)

    def reset_meep(self):
        self.fields = None
        self.structure = None
        self.dft_objects = []

    def restart_fields(self):
        if self.fields is not None:
            self.fields.t = 0
            self.fields.zero_fields()
        else:
            self.init_sim()

    def run(self, *step_funcs, **kwargs):
        until = kwargs.pop('until', None)
        until_after_sources = kwargs.pop('until_after_sources', None)

        if self.fields is None:
            self.init_sim()

        self._evaluate_dft_objects()
        self._check_material_frequencies()

        if kwargs:
            raise ValueError("Unrecognized keyword arguments: {}".format(kwargs.keys()))

        if until_after_sources is not None:
            self._run_sources_until(until_after_sources, step_funcs)
        elif until is not None:
            self._run_until(until, step_funcs)
        else:
            raise ValueError("Invalid run configuration")

    def get_epsilon(self):
        return self.get_array(component=mp.Dielectric)

    def get_mu(self):
        return self.get_array(component=mp.Permeability)

    def get_hpwr(self):
        return self.get_array(component=mp.H_EnergyDensity)

    def get_dpwr(self):
        return self.get_array(component=mp.D_EnergyDensity)

    def get_tot_pwr(self):
        return self.get_array(component=mp.EnergyDensity)

    def get_hfield(self):
        if self.is_cylindrical:
            r = self.get_array(cmplx=True, component=mp.Hr)
            p = self.get_array(cmplx=True, component=mp.Hp)
            return np.stack([r, p], axis=-1)
        else:
            x = self.get_array(cmplx=True, component=mp.Hx)
            y = self.get_array(cmplx=True, component=mp.Hy)
            z = self.get_array(cmplx=True, component=mp.Hz)
            return np.stack([x, y, z], axis=-1)

    def get_hfield_x(self):
        return self.get_array(cmplx=True, component=mp.Hx)

    def get_hfield_y(self):
        return self.get_array(cmplx=True, component=mp.Hy)

    def get_hfield_z(self):
        return self.get_array(cmplx=True, component=mp.Hz)

    def get_hfield_r(self):
        return self.get_array(cmplx=True, component=mp.Hr)

    def get_hfield_p(self):
        return self.get_array(cmplx=True, component=mp.Hp)

    def get_bfield(self):
        if self.is_cylindrical:
            r = self.get_array(cmplx=True, component=mp.Br)
            p = self.get_array(cmplx=True, component=mp.Bp)
            return np.stack([r, p], axis=-1)
        else:
            x = self.get_array(cmplx=True, component=mp.Bx)
            y = self.get_array(cmplx=True, component=mp.By)
            z = self.get_array(cmplx=True, component=mp.Bz)
            return np.stack([x, y, z], axis=-1)

    def get_bfield_x(self):
        return self.get_array(cmplx=True, component=mp.Bx)

    def get_bfield_y(self):
        return self.get_array(cmplx=True, component=mp.By)

    def get_bfield_z(self):
        return self.get_array(cmplx=True, component=mp.Bz)

    def get_bfield_r(self):
        return self.get_array(cmplx=True, component=mp.Br)

    def get_bfield_p(self):
        return self.get_array(cmplx=True, component=mp.Bp)

    def get_efield(self):
        if self.is_cylindrical:
            r = self.get_array(cmplx=True, component=mp.Er)
            p = self.get_array(cmplx=True, component=mp.Ep)
            return np.stack([r, p], axis=-1)
        else:
            x = self.get_array(cmplx=True, component=mp.Ex)
            y = self.get_array(cmplx=True, component=mp.Ey)
            z = self.get_array(cmplx=True, component=mp.Ez)
            return np.stack([x, y, z], axis=-1)

    def get_efield_x(self):
        return self.get_array(cmplx=True, component=mp.Ex)

    def get_efield_y(self):
        return self.get_array(cmplx=True, component=mp.Ey)

    def get_efield_z(self):
        return self.get_array(cmplx=True, component=mp.Ez)

    def get_efield_r(self):
        return self.get_array(cmplx=True, component=mp.Er)

    def get_efield_p(self):
        return self.get_array(cmplx=True, component=mp.Ep)

    def get_dfield(self):
        if self.is_cylindrical:
            r = self.get_array(cmplx=True, component=mp.Dr)
            p = self.get_array(cmplx=True, component=mp.Dp)
            return np.stack([r, p], axis=-1)
        else:
            x = self.get_array(cmplx=True, component=mp.Dx)
            y = self.get_array(cmplx=True, component=mp.Dy)
            z = self.get_array(cmplx=True, component=mp.Dz)
            return np.stack([x, y, z], axis=-1)

    def get_dfield_x(self):
        return self.get_array(cmplx=True, component=mp.Dx)

    def get_dfield_y(self):
        return self.get_array(cmplx=True, component=mp.Dy)

    def get_dfield_z(self):
        return self.get_array(cmplx=True, component=mp.Dz)

    def get_dfield_r(self):
        return self.get_array(cmplx=True, component=mp.Dr)

    def get_dfield_p(self):
        return self.get_array(cmplx=True, component=mp.Dp)

    def get_sfield(self):
        if self.is_cylindrical:
            r = self.get_array(cmplx=True, component=mp.Sr)
            p = self.get_array(cmplx=True, component=mp.Sp)
            return np.stack([r, p], axis=-1)
        else:
            x = self.get_array(cmplx=True, component=mp.Sx)
            y = self.get_array(cmplx=True, component=mp.Sy)
            z = self.get_array(cmplx=True, component=mp.Sz)
            return np.stack([x, y, z], axis=-1)

    def get_sfield_x(self):
        return self.get_array(cmplx=True, component=mp.Sx)

    def get_sfield_y(self):
        return self.get_array(cmplx=True, component=mp.Sy)

    def get_sfield_z(self):
        return self.get_array(cmplx=True, component=mp.Sz)

    def get_sfield_r(self):
        return self.get_array(cmplx=True, component=mp.Sr)

    def get_sfield_p(self):
        return self.get_array(cmplx=True, component=mp.Sp)


def _create_boundary_region_from_boundary_layers(boundary_layers, gv):
    br = mp.boundary_region()

    for layer in boundary_layers:

        if isinstance(layer, Absorber):
            continue

        boundary_region_args = [
            mp.boundary_region.PML,
            layer.thickness,
            layer.R_asymptotic,
            layer.mean_stretch,
            mp.py_pml_profile,
            layer.pml_profile,
            1 / 3,
            1 / 4,
        ]

        if layer.direction == mp.ALL:
            d = mp.start_at_direction(gv.dim)
            loop_stop_directi = mp.stop_at_direction(gv.dim)

            while d < loop_stop_directi:
                if layer.side == mp.ALL:
                    b = mp.High
                    loop_stop_bi = mp.Low

                    while b != loop_stop_bi:
                        br += mp.boundary_region(*(boundary_region_args + [d, b]))
                        b = (b + 1) % 2
                        loop_stop_bi = mp.High
                else:
                    br += mp.boundary_region(*(boundary_region_args + [d, layer.side]))
                d += 1
        else:
            if layer.side == mp.ALL:
                b = mp.High
                loop_stop_bi = mp.Low

                while b != loop_stop_bi:
                    br += mp.boundary_region(*(boundary_region_args + [layer.direction, b]))
                    b = (b + 1) % 2
                    loop_stop_bi = mp.High
            else:
                br += mp.boundary_region(*(boundary_region_args + [layer.direction, layer.side]))
    return br


# Private step functions

def _combine_step_funcs(*step_funcs):
    def _combine(sim, todo):
        for func in step_funcs:
            _eval_step_func(sim, func, todo)
    return _combine


def _eval_step_func(sim, func, todo):
    num_args = get_num_args(func)

    if num_args != 1 and num_args != 2:
        raise ValueError("Step function '{}'' requires 1 or 2 arguments".format(func.__name__))
    elif num_args == 1:
        if todo == 'step':
            func(sim)
    elif num_args == 2:
        func(sim, todo)


def _when_true_funcs(cond, *step_funcs):
    def _true(sim, todo):
        if todo == 'finish' or cond(sim):
            for f in step_funcs:
                _eval_step_func(sim, f, todo)
    return _true


# Public step functions

def after_sources(*step_funcs):
    def _after_sources(sim, todo):
        time = sim.fields.last_source_time()
        if sim.round_time() >= time:
            for func in step_funcs:
                _eval_step_func(sim, func, todo)
    return _after_sources


def after_sources_and_time(t, *step_funcs):
    def _after_s_and_t(sim, todo):
        time = sim.fields.last_source_time() + t - sim.round_time()
        if sim.round_time() >= time:
            for func in step_funcs:
                _eval_step_func(sim, func, todo)
    return _after_s_and_t


def after_time(t, *step_funcs):
    def _after_t(sim):
        return sim.round_time() >= t
    return _when_true_funcs(_after_t, *step_funcs)


def at_beginning(*step_funcs):
    closure = {'done': False}

    def _beg(sim, todo):
        if not closure['done']:
            for f in step_funcs:
                _eval_step_func(sim, f, todo)
            closure['done'] = True
    return _beg


def at_end(*step_funcs):
    def _end(sim, todo):
        if todo == 'finish':
            for func in step_funcs:
                _eval_step_func(sim, func, 'step')
            for func in step_funcs:
                _eval_step_func(sim, func, 'finish')
    return _end


def at_every(dt, *step_funcs):
    closure = {'tlast': 0.0}

    def _every(sim, todo):
        t = sim.round_time()
        if todo == 'finish' or t >= closure['tlast'] + dt + (-0.5 * sim.fields.dt):
            for func in step_funcs:
                _eval_step_func(sim, func, todo)
            closure['tlast'] = t
    return _every


def at_time(t, *step_funcs):
    closure = {'done': False}

    def _at_time(sim, todo):
        if not closure['done'] or todo == 'finish':
            for f in step_funcs:
                _eval_step_func(sim, f, todo)
        closure['done'] = closure['done'] or todo == 'step'
    return after_time(t, _at_time)


def before_time(t, *step_funcs):
    def _before_t(sim):
        return sim.round_time() < t
    return _when_true_funcs(_before_t, *step_funcs)


def during_sources(*step_funcs):
    closure = {'finished': False}

    def _during_sources(sim, todo):
        time = sim.fields.last_source_time()
        if sim.round_time() < time:
            for func in step_funcs:
                _eval_step_func(sim, func, 'step')
        elif closure['finished'] is False:
            for func in step_funcs:
                _eval_step_func(sim, func, 'finish')
            closure['finished'] = True
    return _during_sources


def in_volume(v, *step_funcs):
    closure = {'cur_eps': ''}

    def _in_volume(sim, todo):
        v_save = sim.output_volume
        eps_save = sim.last_eps_filename

        sim.output_volume = sim._fit_volume_to_simulation(v).swigobj

        if closure['cur_eps']:
            sim.last_eps_filename = closure['cur_eps']
        for func in step_funcs:
            _eval_step_func(sim, func, todo)

        closure['cur_eps'] = sim.last_eps_filename
        sim.output_volume = v_save
        if eps_save:
            sim.last_eps_filename = eps_save
    return _in_volume


def in_point(pt, *step_funcs):
    v = Volume(pt)
    return in_volume(v, *step_funcs)


def to_appended(fname, *step_funcs):
    closure = {'h5': None}

    def _to_appended(sim, todo):
        if closure['h5'] is None:
            closure['h5'] = sim.fields.open_h5file(fname, mp.h5file.WRITE, sim.get_filename_prefix())
        h5save = sim.output_append_h5
        sim.output_append_h5 = closure['h5']

        for func in step_funcs:
            _eval_step_func(sim, func, todo)

        if todo == 'finish':
            closure['h5'] = None
            sim.output_h5_hook(sim.fields.h5file_name(fname, sim.get_filename_prefix()))
        sim.output_append_h5 = h5save
    return _to_appended


def stop_when_fields_decayed(dt, c, pt, decay_by):

    closure = {
        'max_abs': 0,
        'cur_max': 0,
        't0': 0,
    }

    def _stop(sim):
        fabs = abs(sim.get_field_point(c, pt)) * abs(sim.get_field_point(c, pt))
        closure['cur_max'] = max(closure['cur_max'], fabs)

        if sim.round_time() <= dt + closure['t0']:
            return False
        else:
            old_cur = closure['cur_max']
            closure['cur_max'] = 0
            closure['t0'] = sim.round_time()
            closure['max_abs'] = max(closure['max_abs'], old_cur)
            if closure['max_abs'] != 0:
                fmt = "field decay(t = {}): {} / {} = {}"
                print(fmt.format(sim.meep_time(), old_cur, closure['max_abs'], old_cur / closure['max_abs']))
            return old_cur <= closure['max_abs'] * decay_by
    return _stop


def synchronized_magnetic(*step_funcs):
    def _sync(sim, todo):
        sim.fields.synchronize_magnetic_fields()
        for f in step_funcs:
            _eval_step_func(sim, f, todo)
        sim.fields.restore_magnetic_fields()
    return _sync


def when_true(cond, *step_funcs):
    return _when_true_funcs(cond, *step_funcs)


def when_false(cond, *step_funcs):
    return _when_true_funcs(lambda: not cond, *step_funcs)


def with_prefix(pre, *step_funcs):
    def _with_prefix(sim, todo):
        saved_pre = sim.filename_prefix
        sim.filename_prefix = pre + sim.get_filename_prefix()

        for f in step_funcs:
            _eval_step_func(sim, f, todo)
        sim.filename_prefix = saved_pre
    return _with_prefix


def display_csv(sim, name, data):
    for d in data:
        display_run_data(sim, name, d)


def display_progress(t0, t, dt):
    t_0 = mp.wall_time()
    closure = {'tlast': mp.wall_time()}

    def _disp(sim):
        t1 = mp.wall_time()
        if t1 - closure['tlast'] >= dt:
            msg_fmt = "Meep progress: {}/{} = {:.1f}% done in {:.1f}s, {:.1f}s to go"
            val1 = sim.meep_time() - t0
            val2 = val1 / (0.01 * t)
            val3 = t1 - t_0
            val4 = (val3 * (t / val1) - val3) if val1 != 0 else 0
            print(msg_fmt.format(val1, t, val2, val3, val4))
            closure['tlast'] = t1
    return _disp


def data_to_str(d):
    if type(d) is complex:
        sign = '+' if d.imag >= 0 else ''
        return "{}{}{}i".format(d.real, sign, d.imag)
    else:
        return str(d)


def display_run_data(sim, data_name, data):
    if isinstance(data, Sequence):
        data_str = [data_to_str(f) for f in data]
    else:
        data_str = [data_to_str(data)]
    print("{}{}:, {}".format(data_name, sim.run_index, ', '.join(data_str)))


def convert_h5(rm_h5, convert_cmd, *step_funcs):

    def convert(fname):
        if mp.my_rank() == 0:
            cmd = convert_cmd.split()
            cmd.append(fname)
            ret = subprocess.call(cmd)
            if ret == 0 and rm_h5:
                os.remove(fname)

    def _convert_h5(sim, todo):
        hooksave = sim.output_h5_hook
        sim.output_h5_hook = convert

        for f in step_funcs:
            _eval_step_func(sim, f, todo)

        sim.output_h5_hook = hooksave

    return _convert_h5


def output_png(compnt, options, rm_h5=True):
    closure = {'maxabs': 0.0}

    def _output_png(sim, todo):
        if todo == 'step':
            if sim.output_volume is None:
                ov = sim.fields.total_volume()
            else:
                ov = sim.output_volume

            closure['maxabs'] = max(closure['maxabs'],
                                    sim.fields.max_abs(compnt, ov))
            convert = sim.h5topng(rm_h5, "-M {} {}".format(closure['maxabs'], options),
                                  lambda sim: sim.output_component(compnt))
            convert(sim, todo)
    return _output_png


def output_epsilon(sim):
    sim.output_component(mp.Dielectric)


def output_mu(sim):
    sim.output_component(mp.Permeability)


def output_hpwr(sim):
    sim.output_component(mp.H_EnergyDensity)


def output_dpwr(sim):
    sim.output_component(mp.D_EnergyDensity)


def output_tot_pwr(sim):
    sim.output_component(mp.EnergyDensity)


def output_hfield(sim):
    sim.output_components('h', mp.Hx, mp.Hy, mp.Hz, mp.Hr, mp.Hp)


def output_hfield_x(sim):
    sim.output_component(mp.Hx)


def output_hfield_y(sim):
    sim.output_component(mp.Hy)


def output_hfield_z(sim):
    sim.output_component(mp.Hz)


def output_hfield_r(sim):
    sim.output_component(mp.Hr)


def output_hfield_p(sim):
    sim.output_component(mp.Hp)


def output_bfield(sim):
    sim.output_components('b', mp.Bx, mp.By, mp.Bz, mp.Br, mp.Bp)


def output_bfield_x(sim):
    sim.output_component(mp.Bx)


def output_bfield_y(sim):
    sim.output_component(mp.By)


def output_bfield_z(sim):
    sim.output_component(mp.Bz)


def output_bfield_r(sim):
    sim.output_component(mp.Br)


def output_bfield_p(sim):
    sim.output_component(mp.Bp)


def output_efield(sim):
    sim.output_components('e', mp.Ex, mp.Ey, mp.Ez, mp.Er, mp.Ep)


def output_efield_x(sim):
    sim.output_component(mp.Ex)


def output_efield_y(sim):
    sim.output_component(mp.Ey)


def output_efield_z(sim):
    sim.output_component(mp.Ez)


def output_efield_r(sim):
    sim.output_component(mp.Er)


def output_efield_p(sim):
    sim.output_component(mp.Ep)


def output_dfield(sim):
    sim.output_components('d', mp.Dx, mp.Dy, mp.Dz, mp.Dr, mp.Dp)


def output_dfield_x(sim):
    sim.output_component(mp.Dx)


def output_dfield_y(sim):
    sim.output_component(mp.Dy)


def output_dfield_z(sim):
    sim.output_component(mp.Dz)


def output_dfield_r(sim):
    sim.output_component(mp.Dr)


def output_dfield_p(sim):
    sim.output_component(mp.Dp)


# MPB compatibility
def output_poynting(sim):
    sim.output_components('s', mp.Sx, mp.Sy, mp.Sz, mp.Sr, mp.Sp)


def output_poynting_x(sim):
    sim.output_component(mp.Sx)


def output_poynting_y(sim):
    sim.output_component(mp.Sy)


def output_poynting_z(sim):
    sim.output_component(mp.Sz)


def output_poynting_r(sim):
    sim.output_component(mp.Sr)


def output_poynting_p(sim):
    sim.output_component(mp.Sp)


def output_sfield(sim):
    sim.output_components('s', mp.Sx, mp.Sy, mp.Sz, mp.Sr, mp.Sp)


def output_sfield_x(sim):
    sim.output_component(mp.Sx)


def output_sfield_y(sim):
    sim.output_component(mp.Sy)


def output_sfield_z(sim):
    sim.output_component(mp.Sz)


def output_sfield_r(sim):
    sim.output_component(mp.Sr)


def output_sfield_p(sim):
    sim.output_component(mp.Sp)


def get_ldos_freqs(f):
    start = f.omega_min / (2 * math.pi)
    stop = start + (f.domega / (2 * math.pi)) * f.Nomega
    return np.linspace(start, stop, num=f.Nomega, endpoint=False).tolist()


def dft_ldos(fcen, df, nfreq):
    ldos = mp._dft_ldos(fcen - df / 2, fcen + df / 2, nfreq)

    def _ldos(sim, todo):
        if todo == 'step':
            ldos.update(sim.fields)
        else:
            sim.ldos_data = mp._dft_ldos_ldos(ldos)
            sim.ldos_Fdata = mp._dft_ldos_F(ldos)
            sim.ldos_Jdata = mp._dft_ldos_J(ldos)
            display_csv(sim, 'ldos', zip(get_ldos_freqs(ldos), sim.ldos_data))
    return _ldos


def scale_flux_fields(s, flux):
    flux.scale_dfts(s)


def get_flux_freqs(f):
    return np.linspace(f.freq_min, f.freq_min + f.dfreq * f.Nfreq, num=f.Nfreq, endpoint=False).tolist()


def get_fluxes(f):
    return f.flux()


def scale_force_fields(s, force):
    force.scale_dfts(s)


def get_eigenmode_freqs(f):
    return np.linspace(f.freq_min, f.freq_min + f.dfreq * f.Nfreq, num=f.Nfreq, endpoint=False).tolist()


def get_force_freqs(f):
    return np.linspace(f.freq_min, f.freq_min + f.dfreq * f.Nfreq, num=f.Nfreq, endpoint=False).tolist()


def get_forces(f):
    return f.force()


def scale_near2far_fields(s, n2f):
    n2f.scale_dfts(s)


def get_near2far_freqs(f):
    return np.linspace(f.freq_min, f.freq_min + f.dfreq * f.Nfreq, num=f.Nfreq, endpoint=False).tolist()


def interpolate(n, nums):
    res = []
    if isinstance(nums[0], mp.Vector3):
        for low, high in zip(nums, nums[1:]):
            x = np.linspace(low.x, high.x, n + 1, endpoint=False).tolist()
            y = np.linspace(low.y, high.y, n + 1, endpoint=False).tolist()
            z = np.linspace(low.z, high.z, n + 1, endpoint=False).tolist()

            for i in range(len(x)):
                res.append(mp.Vector3(x[i], y[i], z[i]))
    else:
        for low, high in zip(nums, nums[1:]):
            res.extend(np.linspace(low, high, n + 1, endpoint=False).tolist())

    return res + [nums[-1]]


