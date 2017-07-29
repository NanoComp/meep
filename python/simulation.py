from __future__ import division, print_function

import numbers

import meep as mp
from meep.geom import Vector3
from meep.source import EigenModeSource, check_positive

# TODO(chogan): Hopefully we don't have to require funcsigs
try:
    from inspect import signature
except ImportError:
    from funcsigs import signature

CYLINDRICAL = -2
no_size = 1e-20


def get_num_args(func):
    sig = signature(func)
    return len(sig.parameters)


class Pml(object):

    def __init__(self, thickness,
                 direction=-1,
                 side=-1,
                 strength=1.0,
                 r_asymptotic=1e-15,
                 mean_stretch=1.0,
                 pml_profile=lambda u: u * u):

        self.direction = direction
        self.side = side
        self.strength = strength
        self.r_asymptotic = r_asymptotic
        self.mean_stretch = mean_stretch
        self.pml_profile = pml_profile

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


class Volume(object):

    def __init__(self, center, size=Vector3()):
        self.center = center
        self.size = size

        v1 = center - size.scale(0.5)
        v2 = center + size.scale(0.5)

        self.swigobj = mp.volume(v1, v2)


class Simulation(object):

    def __init__(self, cell, geometry, sources, resolution, dimensions=2, pml_layers=[], verbose=False):
        self.cell = cell
        self.geometry_center = None
        self.resolution = resolution
        self.eps_averaging = True
        self.subpixel_tol = 1e-4
        self.subpixel_maxeval = 100000
        self.ensure_periodicity = False
        self.geometry = geometry
        self.extra_materials = []
        self.default_material = None
        self.epsion_input_file = ''
        self.pml_layers = pml_layers
        self.symmetries = []
        self.num_chunks = 0
        self.courant = 0.5
        self.global_d_conductivity = 0
        self.global_b_conductivity = 0
        self.special_kz = False
        self.k_point = False
        self.fields = None
        self.structure = None
        self.accurate_fields_near_cylorigin = False
        self.sources = []
        self.m = 0
        self.force_complex_fields = False
        self.verbose = verbose
        self.progress_interval = 4
        self.dimensions = dimensions
        self.init_fields_hooks = []
        self.progress_interval = 4
        self.run_index = 0

    def _infer_dimensions(self, k):
        if k and self.dimensions == 3:

            def requires_2d(self, k):
                cond1 = False if k is None else not k
                cond2 = self.cell.size.z == self.no_size
                cond3 = cond1 or self.special_kz or k.z == 0
                return cond2 and cond3

            if requires_2d(self, k):
                return 2
            else:
                return 3
        return self.dimensions

    def _init_structure(self, k=None):
        # TODO(chogan): More descriptive name for not_not_k?
        # not_not_k = True if k is None else k

        # self.structure = mp.structure(
        #     self._infer_dimensions(k),
        #     self.cell.size,
        #     self.geometry_center,
        #     self.resolution,
        #     self.eps_averaging,
        #     self.subpixel_tol,
        #     self.subpixel_maxeval,
        #     self.ensure_periodicity and not_not_k,
        #     self.geometry,
        #     self.extra_materials,
        #     self.default_material,
        #     self.epsion_input_file,
        #     self.pml_layers,
        #     self.symmetries,
        #     self.num_chunks,
        #     self.courant,
        #     self.global_d_conductivity,
        #     self.global_b_conductivity
        # )
        dims = self._infer_dimensions(k)

        if dims == 0 or dims == 1:
            gv = mp.vol1d(self.cell.size.z, self.resolution)
        elif dims == 2:
            gv = mp.vol2d(self.cell.size.x, self.cell.size.y, self.resolution)
        elif dims == 3:
            gv = mp.vol3d(self.cell.size.x, self.cell.size.y, self.cell.size.z, self.resolution)
        elif dims == CYLINDRICAL:
            gv = mp.volcyl(self.cell.size.x, self.cell.size.z, self.resolution)
        else:
            raise ValueError("Unsupported dimentionality: {}".format(dims))

        gv.center_origin()

        def dummy_eps(v):
            return 1

        self.structure = mp.structure(gv, dummy_eps)
        mp.set_materials_from_geometry(self.structure, self.geometry)

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
            cond2 = not all([s.sphase.imag for s in self.symmetries])
            cond3 = not self.k_point
            cond4 = self.special_kz and self.k_point.x == 0 and self.k_point.y == 0
            cond5 = not (cond3 or cond4 or self.k_point == Vector3())
            return not (self.force_complex_fields or cond1 or cond2 or cond5)

        if use_real(self):
            self.fields.use_real_fields()
        else:
            print("Meep: using comlex fields.")

        if self.k_point:
            v = Vector3(self.k_point.x, self.k_point.y) if self.special_kz else self.k_point
            self.fields.use_bloch(v)

        for s in self.sources:
            self.add_source(s)

        for hook in self.init_fields_hooks:
            hook()

    def display_progress(self, t0, t, dt):
        t_0 = mp.wall_time()
        tlast = mp.wall_time()

        t1 = mp.wall_time()

        if t1 - tlast >= dt:
            msg_fmt = "Meep progress: {}/{} = {:.1g}%% done in {:.1g}s, {} s to go"
            val1 = mp.time() - t0
            val2 = t
            val3 = (mp.time() - t0) / (0.01 * t)
            val4 = t1 - t_0
            val5 = ((t1 - t_0) * (t / (mp.time() - t0)) - (t1 - t_0))
            print(msg_fmt.format(val1, val2, val3, val4, val5))
            tlast = t1

    def _round_time(self):
        if self.fields is None:
            self._init_fields()

        return self.fields.round_time()

    def _time(self):
        if self.fields is None:
            self._init_fields()

        return self.fields.time()

    def _eval_step_func(self, func, todo):
        # TODO(chogan): Should func have a 'self' param?
        num_args = get_num_args(func)
        if num_args == 0:
            if todo == 'step':
                func()
        else:
            func(todo)

    def _combine_step_funcs(self, step_funcs):
        closure = {'step_funcs': step_funcs}

        def f(todo):
            for f in closure['step_funcs']:
                self._eval_step_func(f, todo)
        return f

    def _call_step_func(self, name, funcs_dict):
        if name in funcs_dict:
            getattr(self, name)(funcs_dict[name])

    def _run_until(self, cond, step_funcs):
        # TODO(chogan): Interactive?
        # TODO(chogan): while loops make more sense than recursion here.
        if self.fields is None:
            self._init_fields()

        if isinstance(cond, numbers.Number):
            closure = {
                'stop_time': cond,
                't0': self._round_time()
            }

            def stop_cond():
                return self._round_time() >= closure['t0'] + closure['stop_time']

            # new_step_funcs = step_funcs.copy()
            # new_step_funcs.update({'display_progress': (self, t0, t0 + stop_time, self.progress_interval)})
            self._run_until(stop_cond, step_funcs)
        else:
            for func in step_funcs:
                self._eval_step_func(func, 'step')

            if cond():
                for func in step_funcs:
                    self._eval_step_func(func, 'finish')

                print("run {} finished at t = {} ({} timesteps)".format(self.run_index, self._time(), self.fields.t))
                self.run_index += 1
            else:
                self.fields.step()
                self._run_until(cond, step_funcs)

    def _run_sources_until(self, cond, step_funcs):
        if self.fields is None:
            self._init_fields()

        ts = self.fields.last_source_time()

        if isinstance(cond, numbers.Number):
            arg = (ts - self._round_time()) + cond
        else:
            def f():
                cond and self._round_time() >= ts
            arg = f

        self._run_until(step_funcs, arg)

    def _harminv(self, data, dt, results, c, pt, fcen, df, maxbands):
        pass

    def harminv(self, c, pt, fcen, df, *mxbands):
        return self._harminv([], 0, [], c, pt, fcen, df, mxbands)

    def add_source(self, src):
        if self.fields is None:
            self._init_fields()

        eig_vol = Volume(src.eigen_lattice_center, src.eigen_lattice_size).swigobj
        where = Volume(src.center, src.size).swigobj
        direction = self.fields.normal_direction(where)

        if isinstance(src, EigenModeSource):
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

    def output_component(self, c, h5file):
        if self.fields is None:
            raise RuntimeError("Fields must be initialized before calling output_component")

        vol = self.fields.total_volume() if self.output_volume is None else self.output_volume

        # self.fields.output_hdf5(c, vol)

    def output_epsilon(self):
        self.output_component(mp.Dielectric)

    def when_true_funcs(self, cond, *step_funcs):
        def f(todo):
            if todo == 'finish' or cond():
                for f in step_funcs:
                    self._eval_step_func(f, todo)
        return f

    def at_beginning(self, *step_funcs):
        # Work around python 2's lack of 'nonlocal' keyword
        closure = {
            'done': False,
            'step_funcs': step_funcs
        }

        def step_func(todo):
            if not closure['done']:
                for f in closure['step_funcs']:
                    self._eval_step_func(f, todo)
                closure['done'] = True
        return step_func

    def after_time(self, t, *step_funcs):
        if self.fields is None:
            self._init_fields()

        closure = {
            't0': self._round_time(),
            't': t
        }

        def f():
            return self._round_time() >= closure['t0'] + closure['t']

        self.when_true_funcs(f, *step_funcs)

    def after_sources(self, *step_funcs):
        if self.fields is None:
            self._init_fields()

        time = self.fields.last_source_time() - self._round_time()

        self.after_time(time, *step_funcs)

    def run(self, *step_funcs, **kwargs):
        until = kwargs.pop('until', None)
        sources = kwargs.pop('sources', None)

        if kwargs:
            raise ValueError("Unrecognized keyword arguments: {}".format(kwargs.keys()))

        if sources and until:
            self._run_sources_until(until, step_funcs)
        elif until:
            self._run_until(until, step_funcs)
        elif until is None and sources:
            self._run_sources(step_funcs)
        else:
            raise ValueError("Invalid run configuration")
