from __future__ import division, print_function

import numbers
import os
import re
import sys

import meep as mp
from meep.geom import Vector3
from meep.source import EigenModeSource, check_positive


try:
    basestring
except NameError:
    basestring = str


CYLINDRICAL = -2


def get_num_args(func):
    return func.__code__.co_argcount


def is_bound(f):
    return hasattr(f, '__self__') and f.__self__ is not None


def py_v3_to_vec(dims, v3):
    if dims == 1:
        if v3.x == v3.y == v3.z:
            return mp.vec(dims, v3.x)
        else:
            return mp.vec(dims)
    elif dims == 2:
        return mp.vec(v3.x, v3.y)
    elif dims == 3:
        return mp.vec(v3.x, v3.y, v3.z)
    else:
        raise ValueError("Invalid dimensions in Volume: {}".format(dims))


class Pml(object):

    def __init__(self, thickness,
                 direction=-1,
                 side=-1,
                 r_asymptotic=1e-15,
                 mean_stretch=1.0,
                 pml_profile=lambda u: u * u):

        self.thickness = thickness
        self.direction = direction
        self.side = side
        self.r_asymptotic = r_asymptotic
        self.mean_stretch = mean_stretch
        self.pml_profile = pml_profile

        if direction == -1 and side == -1:
            self.swigobj = mp.pml(thickness, r_asymptotic, mean_stretch)
        elif direction == -1:
            self.swigobj = mp.pml(thickness, side, r_asymptotic, mean_stretch)
        else:
            self.swigobj = mp.pml(thickness, direction, side, r_asymptotic, mean_stretch)

    @property
    def r_asymptotic(self):
        return self._r_asymptotic

    @r_asymptotic.setter
    def r_asymptotic(self, val):
        self._r_asymptotic = check_positive('Pml.r_asymptotic', val)

    @property
    def mean_stretch(self):
        return self._mean_stretch

    @mean_stretch.setter
    def mean_stretch(self, val):
        if val >= 1:
            self._mean_stretch = val
        else:
            raise ValueError("Pml.mean_stretch must be >= 1. Got {}".format(val))


class Absorber(Pml):
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

    def __init__(self, center, size=Vector3(), dims=2,):
        self.center = center
        self.size = size
        self.dims = dims

        v1 = center - size.scale(0.5)
        v2 = center + size.scale(0.5)

        vec1 = py_v3_to_vec(self.dims, v1)
        vec2 = py_v3_to_vec(self.dims, v2)

        self.swigobj = mp.volume(vec1, vec2)


class Simulation(object):

    def __init__(self, cell_size, geometry, sources, resolution, eps_averaging=True,
                 dimensions=2, boundary_layers=[], symmetries=[], verbose=False):
        self.cell_size = cell_size
        self.geometry = geometry
        self.sources = sources
        self.resolution = resolution
        self.dimensions = dimensions
        self.boundary_layers = boundary_layers
        self.symmetries = symmetries
        self.geometry_center = Vector3()
        self.eps_averaging = eps_averaging
        self.subpixel_tol = 1e-4
        self.subpixel_maxeval = 100000
        self.ensure_periodicity = False
        self.extra_materials = []
        self.default_material = None
        self.epsion_input_file = ''
        self.num_chunks = 0
        self.courant = 0.5
        self.global_d_conductivity = 0
        self.global_b_conductivity = 0
        self.special_kz = False
        self.k_point = False
        self.fields = None
        self.structure = None
        self.accurate_fields_near_cylorigin = False
        self.m = 0
        self.force_complex_fields = False
        self.verbose = verbose
        self.progress_interval = 4
        self.init_fields_hooks = []
        self.progress_interval = 4
        self.run_index = 0
        self.filename_prefix = ''
        self.include_files = []
        self.output_append_h5 = None
        self.output_single_precision = False
        self.output_volume = []
        self.last_eps_filename = ''
        self.output_h5_hook = lambda fname: False
        self.harminv_data = []
        self.harminv_data_dt = 0
        self.harminv_results = []
        self.harminv_spectral_density = 1.1
        self.harminv_Q_thresh = 50.0
        self.harminv_rel_err_thresh = 1e20
        self.harminv_err_thresh = 0.01
        self.harminv_rel_amp_thresh = -1.0
        self.harminv_amp_thresh = -1.0
        self.interactive = False

        self._init_include_files()

    def _init_include_files(self):
        for arg in sys.argv:
            if isinstance(arg, basestring) and arg.endswith('.py'):
                self.include_files.append(arg)

    def _infer_dimensions(self, k):
        if k and self.dimensions == 3:

            def requires_2d(self, k):
                cond1 = False if k is None else not k
                cond2 = self.cell_size.z == 0
                cond3 = cond1 or self.special_kz or k.z == 0
                return cond2 and cond3

            if requires_2d(self, k):
                return 2
            else:
                return 3
        return self.dimensions

    def _init_structure(self, k=None):
        print('-' * 11)
        print('Initializing structure...')

        dims = self._infer_dimensions(k)

        if dims == 0 or dims == 1:
            gv = mp.vol1d(self.cell_size.z, self.resolution)
        elif dims == 2:
            gv = mp.vol2d(self.cell_size.x, self.cell_size.y, self.resolution)
        elif dims == 3:
            gv = mp.vol3d(self.cell_size.x, self.cell_size.y, self.cell_size.z, self.resolution)
        elif dims == CYLINDRICAL:
            gv = mp.volcyl(self.cell_size.x, self.cell_size.z, self.resolution)
        else:
            raise ValueError("Unsupported dimentionality: {}".format(dims))

        gv.center_origin()

        def dummy_eps(v):
            return 1

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

        br = _create_boundary_region_from_boundary_layers(self.boundary_layers, gv)

        self.structure = mp.structure(gv, dummy_eps, br, sym, self.num_chunks, self.courant,
                                      self.eps_averaging, self.subpixel_tol, self.subpixel_maxeval)

        mp.set_materials_from_geometry(self.structure, self.geometry, self.eps_averaging,
                                       self.subpixel_tol, self.subpixel_maxeval, self.ensure_periodicity)

    def _init_fields(self):
        is_cylindrical = self.dimensions == CYLINDRICAL

        if self.structure is None:
            self._init_structure(self.k_point)

        self.fields = mp.fields(
            self.structure,
            self.m if is_cylindrical else 0,
            self.k_point.z if self.special_kz and self.k_point else 0,
            not self.accurate_fields_near_cylorigin
        )

        if self.verbose:
            self.fields.verbose()

        def use_real(self):
            cond1 = is_cylindrical and self.m != 0
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
            self.fields.use_bloch(v)

        for s in self.sources:
            self.add_source(s)

        for hook in self.init_fields_hooks:
            hook()

    def meep_time(self):
        if self.fields is None:
            self._init_fields()
        return self.fields.time()

    def _round_time(self):
        if self.fields is None:
            self._init_fields()

        return self.fields.round_time()

    def _get_field_point(self, c, pt):
        v3 = py_v3_to_vec(self.dimensions, pt)
        return self.fields.get_field_from_comp(c, v3)

    def _get_filename_prefix(self):
        if self.include_files and self.filename_prefix == '':
            _, filename = os.path.split(self.include_files[0])
            return re.sub(r'\.py', '', filename)
        else:
            return self.filename_prefix

    def _eval_step_func(self, func, todo):
        num_args = get_num_args(func)
        if num_args == 0:
            raise ValueError("Step function '{}'' requires at least 1 argument".format(func.__name__))
        elif num_args == 1:
            if todo == 'step':
                if is_bound(func):
                    # Bound method already has self as first argument
                    func()
                else:
                    func(self)
        elif num_args == 2:
            if is_bound(func):
                # Bound method already has self as first argument
                func(todo)
            else:
                func(self, todo)
        else:
            raise ValueError("Step function '{}' can only take 1 or 2 arguments".format(func.__name__))

    def _combine_step_funcs(self, *step_funcs):
        closure = {'step_funcs': step_funcs}

        def _combine(sim, todo):
            for func in closure['step_funcs']:
                self._eval_step_func(func, todo)
        return _combine

    def _run_until(self, cond, step_funcs):
        self.interactive = False
        if self.fields is None:
            self._init_fields()

        if isinstance(cond, numbers.Number):
            closure = {
                'stop_time': cond,
                't0': self._round_time(),
            }

            def stop_cond(sim):
                return self._round_time() >= closure['t0'] + closure['stop_time']

            cond = stop_cond

            step_funcs = list(step_funcs)
            step_funcs.append(self.display_progress(
                closure['t0'], closure['t0'] + closure['stop_time'], self.progress_interval
            ))

        while not cond(self):
            for func in step_funcs:
                self._eval_step_func(func, 'step')
            self.fields.step()

        # Translating the recursive scheme version of run-until into an iterative version
        # (because python isn't tail-call-optimized) means we need one extra iteration to
        # be the same as scheme.
        for func in step_funcs:
            self._eval_step_func(func, 'step')

        for func in step_funcs:
            self._eval_step_func(func, 'finish')

        print("run {} finished at t = {} ({} timesteps)".format(self.run_index, self.meep_time(), self.fields.t))
        self.run_index += 1

    def _run_sources_until(self, cond, step_funcs):
        if self.fields is None:
            self._init_fields()

        ts = self.fields.last_source_time()

        if isinstance(cond, numbers.Number):
            new_cond = (ts - self._round_time()) + cond
        else:
            def f(sim):
                cond() and self._round_time() >= ts
            new_cond = f

        self._run_until(new_cond, step_funcs)

    def _run_sources(self, step_funcs):
        self._run_sources_until(self, 0, step_funcs)

    def _collect_harminv(self, data, data_dt):

        closure = {'data': data, 'data_dt': data_dt}

        def _collect1(c, pt):
            closure['data'] = []

            closure2 = {'t0': 0}

            def _collect2(sim):
                closure['data_dt'] = self.meep_time() - closure2['t0']
                closure2['t0'] = self.meep_time()
                self.harminv_data.append(self._get_field_point(c, pt))
            return _collect2
        return _collect1

    def _display_run_data(self, data_name, data):
        data_str = [str(f) for f in data]
        print("{}{}:, {}".format(data_name, self.run_index, ', '.join(data_str)))

    def _analyze_harminv(self, data, fcen, df, maxbands, dt=None):
        self._display_run_data('harminv', ['frequency', 'imag.', 'freq.', 'Q', '|amp|', 'amplitude', 'error'])

        bands = mp.py_do_harminv(self.harminv_data, dt if dt else self.fields.dt, fcen - df / 2, fcen + df / 2,
                                 maxbands, self.harminv_spectral_density, self.harminv_Q_thresh,
                                 self.harminv_rel_err_thresh, self.harminv_err_thresh,
                                 self.harminv_rel_amp_thresh, self.harminv_amp_thresh)

        for freq, amp, err in bands:
            Q = freq.real / (-2 * freq.imag)
            self._display_run_data('harminv', [freq.real, freq.imag, Q, abs(amp), amp, err])

        return bands

    def _harminv(self, data, dt, results, c, pt, fcen, df, maxbands):
        closure = {
            '_data': [],
            '_dt': 0,
            '_c': c,
            '_pt': pt,
            '_fcen': fcen,
            '_df': df,
            '_maxbands': maxbands,
            'data': data,
            'dt': dt,
            'results': results,
        }

        def _harm(sim):
            closure['data'] = closure['_data']
            closure['dt'] = closure['_dt']
            if closure['_maxbands'] is None or closure['_maxbands'] == 0:
                mb = 100
            else:
                mb = closure['_maxbands']

            closure['results'] = self._analyze_harminv(closure['data'], closure['_fcen'],
                                                       closure['_df'], mb, closure['_dt'])

        f1 = self._collect_harminv(closure['_data'], closure['_dt'])

        return self._combine_step_funcs(self.at_end(_harm), f1(closure['_c'], closure['_pt']))

    def harminv(self, c, pt, fcen, df, mxbands=None):
        return self._harminv(self.harminv_data, self.harminv_data_dt, self.harminv_results,
                             c, pt, fcen, df, mxbands)

    def add_source(self, src):
        if self.fields is None:
            self._init_fields()

        where = Volume(src.center, src.size).swigobj

        if isinstance(src, EigenModeSource):
            if src.direction < 0:
                direction = self.fields.normal_direction(where)
            else:
                direction = src.direction

            eig_vol = Volume(src.eig_lattice_center, src.eig_lattice_size).swigobj

            if src.amp_func is None:
                self.fields.add_eigenmode_src(
                    src.component,
                    src.src.swigobj,
                    direction,
                    where,
                    eig_vol,
                    src.eig_band,
                    src.eig_kpoint,
                    src.eig_match_freq,
                    src.eig_parity,
                    src.eig_resolution,
                    src.eig_tolerance,
                    src.amplitude
                )
            else:
                self.fields.add_eigenmode_src(
                    src.component,
                    src.src.swigobj,
                    direction,
                    where,
                    eig_vol,
                    src.eig_band,
                    src.eig_kpoint,
                    src.eig_parity,
                    src.eig_resolution,
                    src.eig_tolerance,
                    src.amplitude,
                    src.amp_func

                )
        else:
            if src.amp_func is None:
                self.fields.add_volume_source(
                    src.component,
                    src.src.swigobj,
                    where,
                    src.amplitude * 1.0,
                )
            else:
                self.fields.add_volume_source(
                    src.component,
                    src.src.swigobj,
                    where,
                    src.amp_func,
                    src.amplitude * 1.0
                )

    def add_flux(self, flux):
        pass

    def display_progress(self, t0, t, dt):
        closure = {
            't_0': mp.wall_time(),
            'tlast': mp.wall_time(),
        }

        def _disp(sim):
            t1 = mp.wall_time()
            if t1 - closure['tlast'] >= dt:
                msg_fmt = "Meep progress: {}/{} = {:.1f}% done in {:.1f}s, {:.1f}s to go"
                val1 = self.meep_time() - t0
                val2 = t
                val3 = (self.meep_time() - t0) / (0.01 * t)
                val4 = t1 - closure['t_0']
                val5 = ((t1 - closure['t_0']) * (t / (self.meep_time() - t0)) - (t1 - closure['t_0']))
                print(msg_fmt.format(val1, val2, val3, val4, val5))
                closure['tlast'] = t1
        return _disp

    def output_component(self, c, h5file=None):
        if self.fields is None:
            raise RuntimeError("Fields must be initialized before calling output_component")

        vol = self.fields.total_volume() if not self.output_volume else self.output_volume
        h5 = self.output_append_h5 if h5file is None else h5file
        append = h5file is None and self.output_append_h5 is not None

        self.fields.output_hdf5(c, vol, h5, append, self.output_single_precision, self._get_filename_prefix())

        if h5file is None:
            nm = self.fields.h5file_name(mp.component_name(c), self._get_filename_prefix(), True)
            if c == mp.Dielectric:
                self.last_eps_filename = nm
            self.output_h5_hook(nm)

    def output_epsilon(self):
        self.output_component(mp.Dielectric)

    # TODO(chogan): Write rest of convenience functions generated in define-output-field
    def output_efield_z(self):
        self.output_component(mp.Ez)

    def when_true_funcs(self, cond, *step_funcs):
        def _true(sim, todo):
            if todo == 'finish' or cond(sim):
                for f in step_funcs:
                    self._eval_step_func(f, todo)
        return _true

    def at_beginning(self, *step_funcs):
        # Work around python 2's lack of 'nonlocal' keyword
        closure = {
            'done': False,
            'step_funcs': step_funcs
        }

        def _beg(sim, todo):
            if not closure['done']:
                for f in closure['step_funcs']:
                    self._eval_step_func(f, todo)
                closure['done'] = True
        return _beg

    def at_every(self, dt, *step_funcs):
        if self.fields is None:
            self._init_fields()

        closure = {
            'tlast': self._round_time(),
            'step_funcs': step_funcs
        }

        def _every(sim, todo):
            t = self._round_time()
            if todo == 'finish' or t >= closure['tlast'] + dt + (-0.5 * self.fields.dt):
                for func in closure['step_funcs']:
                    self._eval_step_func(func, todo)
                closure['tlast'] = t
        return _every

    def at_end(self, *step_funcs):
        def _end(sim, todo):
            if todo == 'finish':
                for func in step_funcs:
                    self._eval_step_func(func, 'step')
                for func in step_funcs:
                    self._eval_step_func(func, 'finish')
        return _end

    def after_time(self, t, *step_funcs):
        if self.fields is None:
            self._init_fields()

        closure = {
            't0': self._round_time(),
            't': t
        }

        def _after_t(sim):
            return self._round_time() >= closure['t0'] + closure['t']

        return self.when_true_funcs(_after_t, *step_funcs)

    def after_sources(self, *step_funcs):
        if self.fields is None:
            self._init_fields()

        time = self.fields.last_source_time() - self._round_time()

        return self.after_time(time, *step_funcs)

    def run(self, *step_funcs, **kwargs):
        until = kwargs.pop('until', None)
        until_after_sources = kwargs.pop('sources', None)

        if kwargs:
            raise ValueError("Unrecognized keyword arguments: {}".format(kwargs.keys()))

        if until_after_sources:
            self._run_sources_until(until_after_sources, step_funcs)
        elif until:
            self._run_until(until, step_funcs)
        else:
            raise ValueError("Invalid run configuration")


def _create_boundary_region_from_boundary_layers(boundary_layers, gv):
    br = mp.boundary_region()

    for layer in boundary_layers:

        if isinstance(layer, Absorber):
            continue

        boundary_region_args = [
            mp.boundary_region.PML,
            layer.thickness,
            layer.r_asymptotic,
            layer.mean_stretch,
            mp.py_pml_profile,
            layer.pml_profile,
            1 / 3,  # TODO(chogan): Call adaptive_integration instead of hard-coding integral
            1 / 4,  # TODO(chogan): Call adaptive_integration instead of hard-coding integral
        ]

        if layer.direction == -1:
            d = mp.start_at_direction(gv.dim)
            loop_stop_directi = mp.stop_at_direction(gv.dim)

            while d < loop_stop_directi:
                if layer.side == -1:
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
            if layer.side == -1:
                b = mp.High
                loop_stop_bi = mp.Low

                while b != loop_stop_bi:
                    br += mp.boundary_region(*(boundary_region_args + [layer.direction]))
                    b = (b + 1) % 2
                    loop_stop_bi = mp.High
            else:
                br += mp.boundary_region(*(boundary_region_args + [layer.direction, layer.side]))
    return br
