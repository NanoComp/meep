import meep as mp
from meep.geom import Vector3
from meep.source import EigenModeSource, check_positive

no_size = 1e-20


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

    # TODO(chogan): Look for scheme defaults for all properties
    def __init__(self, cell=None, geometry=None, source=None, verbose=False):
        self.cell = cell
        self.geometry_center = None
        self.resolution = 0
        self.eps_averaging = True
        self.subpixel_tol = 1e-4
        self.subpixel_maxeval = 100000
        self.ensure_periodicity = False
        self.geometry = geometry
        self.extra_materials = []
        self.default_material = None
        self.epsion_input_file = ''
        self.pml_layers = []
        self.symmetries = []
        self.num_chunks = 0
        self.courant = 0
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

    def _infer_dimensions(self, k):
        if k and self.dimensions == 3:
            # TODO(chogan): more descriptive name for not_k?
            not_k = False if k is None else not k
            if self.geometry_lattice.size.z == no_size and (not_k or self.special_kz or k.z == 0):
                return 2
            else:
                return 3
        return self.dimensions

    def _init_structure(self, k=None):
        # TODO(chogan): More descriptive name for not_not_k
        not_not_k = True if k is None else k
        self.structure = mp.structure(
            self._infer_dimensions(k),
            self.geometry_lattice.size,
            self.geometry_center,
            self.resolution,
            self.eps_averaging,
            self.subpixel_tol,
            self.subpixel_maxeval,
            self.ensure_periodicity and not_not_k,
            self.geometry,
            self.extra_materials,
            self.default_material,
            self.epsion_input_file,
            self.pml_layers,
            self.symmetries,
            self.num_chunks,
            self.courant,
            self.global_d_conductivity,
            self.global_b_conductivity
        )

    def _init_fields(self):
        is_cylindrical = self.dimensions == mp.CYLINDRICAL

        if self.structure is None:
            self._init_structure(self.k_point)

        self.fields = mp.fields(
            self.structure,
            m if is_cylindrical else 0,
            self.k_point.z if self.special_kz and k_point else 0,
            not self.accurate_fields_near_cylorigin
        )

        if self.verbose:
            self.fields.verbose()

        # TODO(chogan): Complete huge if condition
        if not (self.force_complex_fields or (is_cylindrical and m != 0)):
            self.fields.use_real_fields()
        else:
            print("Meep: using comlex fields.")

        if k_point:
            v = Vector3(k_point.x, k_point.y) if self.special_kz else k_point
            self.fields.use_bloch(v)

        for source in self.sources:
            self.add_source(s)

        # for hook in self.init_fields_hooks:
        #     hook()

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
                    src.eig_tolerance),
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

    def run(self):
        pass
