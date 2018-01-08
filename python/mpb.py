from __future__ import division

import numpy as np

import meep as mp
from meep.simulation import get_num_args
from meep.source import check_positive

U_MIN = 0
U_PROD = 1
U_SUM = 2


class MaterialGrid(object):

    def __init__(self,
                 epsilon_min,
                 epsilon_max,
                 size,
                 mu_min=1.0,
                 mu_max=1.0,
                 material_grid_kind=U_MIN,
                 matgrid_init=lambda x, y, z: 0.5):

        self.epsilon_min = epsilon_min
        self.epsilon_max = epsilon_max
        self.size = size
        self.mu_min = mu_min
        self.mu_max = mu_max
        self.material_grid_kind = material_grid_kind
        self.matgrid_init = matgrid_init

        # TODO: inexact->exact?
        g = np.array([0.333] + self.size)

    @property
    def material_grid_kind(self):
        return self._material_grid_kind

    @material_grid_kind.setter
    def material_grid_kind(self, val):
        if val < U_MIN or val > U_SUM:
            fmt = "material_grid_kind must be >= {} and <= {}: Got {}"
            raise ValueError(fmt.format(U_MIN, U_SUM, val))
        self._material_grid_kind = val

    @property
    def size(self):
        return self._size

    @size.setter
    def size(self, val):
        for v in [val.x, val.y, val.z]:
            check_positive("MaterialGrid.size", v)
            if v % 1 != 0:
                raise ValueError("MaterialGrid.size must be a Vector3 of integers")
        self._size = val

    @property
    def matgrid_init(self):
        return self._matgrid_init

    @matgrid_init.setter
    def matgrid_init(self, val):
        num_args = get_num_args(val)
        if num_args != 3:
            raise ValueError("Expected a function that takes 3 arguments. Got {}".format(num_args))
        self._matgrid_init = val


class Lattice(object):

    def __init__(self,
                 basis1=mp.Vector3(1, 0, 0),
                 basis2=mp.Vector3(0, 1, 0),
                 basis3=mp.Vector3(0, 0, 1),
                 size=mp.Vector3(1, 1, 1),
                 basis_size=mp.Vector3(1, 1, 1)):

        self.basis1 = basis1
        self.basis2 = basis2
        self.basis3 = basis3
        self.size = size
        self.basis_size = basis_size


# Class ellipsoid:
#     Class block:
#         Class geometric-object:
#             material-type material = ((material-type))
#             vector3 center
#         vector3 e1 = #(1.0 0.0 0.0)
#         vector3 e2 = #(0.0 1.0 0.0)
#         vector3 e3 = #(0.0 0.0 1.0)
#         vector3 size
# Class block:
#     Class geometric-object:
#         material-type material = ((material-type))
#         vector3 center
#     vector3 e1 = #(1.0 0.0 0.0)
#     vector3 e2 = #(0.0 1.0 0.0)
#     vector3 e3 = #(0.0 0.0 1.0)
#     vector3 size
# Class sphere:
#     Class geometric-object:
#         material-type material = ((material-type))
#         vector3 center
#     number radius
# Class wedge:
#     Class cylinder:
#         Class geometric-object:
#             material-type material = ((material-type))
#             vector3 center
#         vector3 axis = #(0.0 0.0 1.0)
#         number radius
#         number height
#     number wedge-angle = 6.283185307179586
#     vector3 wedge-start = #(1.0 0.0 0.0)
# Class cone:
#     Class cylinder:
#         Class geometric-object:
#             material-type material = ((material-type))
#             vector3 center
#         vector3 axis = #(0.0 0.0 1.0)
#         number radius
#         number height
#     number radius2 = 0
# Class cylinder:
#     Class geometric-object:
#         material-type material = ((material-type))
#         vector3 center
#     vector3 axis = #(0.0 0.0 1.0)
#     number radius
#     number height
# Class compound-geometric-object:
#     Class geometric-object:
#         material-type material = ((material-type))
#         vector3 center
#     geometric-object list component-objects = ()
# Class geometric-object:
#     material-type material = ((material-type))
#     vector3 center
# Class material-grid:
#     Class material-type:
#     integer material-grid-kind = 0
#     number epsilon-min
#     number epsilon-max
#     number mu-min = 1.0
#     number mu-max = 1.0
#     vector3 size
#     function matgrid-init = #<procedure 28b3c00 at ice-9/eval.scm:416:20 (a b c)>
# Class material-function:
#     Class material-type:
#     function material-func
# Class medium-anisotropic:
#     Class material-type:
#     vector3 epsilon-diag = #(1.0 1.0 1.0)
#     cvector3 epsilon-offdiag = #(0.0 0.0 0.0)
#     vector3 epsilon-offdiag-imag = #(0.0 0.0 0.0)
#     vector3 mu-diag = #(1.0 1.0 1.0)
#     cvector3 mu-offdiag = #(0.0 0.0 0.0)
#     vector3 mu-offdiag-imag = #(0.0 0.0 0.0)
# Class medium:
#     Class material-type:
#     number epsilon = 1.0
#     number mu = 1.0
# Class material-type:

class Matrix(object):

    def __init__(self):
        pass


class ModeSolver(object):

    def __init__(self,
                 resolution=10,
                 is_negative_epsilon_ok=False,
                 eigensolver_flops=0,
                 is_eigensolver_davidson=False,
                 eigensolver_nwork=3,
                 eigensolver_block_size=-11,
                 eigensolver_flags=68,
                 is_simple_preconditioner=False,
                 is_deterministic=False,
                 force_mu=False,
                 mu_input_file='',
                 epsilon_input_file='',
                 mesh_size=3,
                 target_freq=0.0,
                 tolerance=1.0e-7,
                 num_bands=1,
                 k_points=[],
                 ensure_periodicity=True,
                 geometry=[],
                 geometry_lattice=Lattice(),
                 geometry_center=mp.Vector3(0, 0, 0),
                 default_material=mp.Medium(epsilon=1),
                 dimensions=3,
                 randomize_fields=False):

        self.resolution = resolution
        self.is_negative_epsilon_ok = is_negative_epsilon_ok
        self.eigensolver_flops = eigensolver_flops
        self.is_eigensolver_davidson = is_eigensolver_davidson
        self.eigensolver_nwork = eigensolver_nwork
        self.eigensolver_block_size = eigensolver_block_size
        self.eigensolver_flags = eigensolver_flags
        self.is_simple_preconditioner = is_simple_preconditioner
        self.is_deterministic = is_deterministic
        self.force_mu = force_mu
        self.mu_input_file = mu_input_file
        self.epsilon_input_file = epsilon_input_file
        self.mesh_size = mesh_size
        self.target_freq = target_freq
        self.tolerance = tolerance
        self.num_bands = num_bands
        self.k_points = k_points
        self.ensure_periodicity = ensure_periodicity
        self.geometry = geometry
        self.geometry_lattice = geometry_lattice
        self.geometry_center = geometry_center
        self.default_material = default_material
        self.dimensions = dimensions
        self.randomize_fields = randomize_fields
        self.parity = None
        self.iterations = 0
        self.freqs = []
        self.eigensolver_flops = 0

    def run_parity(self, p, reset_fields, *band_functions):
        if self.randomize_fields and mpb.randomize_fields not in band_functions:
            band_functions.append(mpb.randomize_fields)

        #  (set! total-run-time (+ total-run-time
        #   (begin-time "total elapsed time for run: "
        #    (set! all-freqs '())
        #    (set! band-range-data '())
        #    (set! interactive? false)  ; don't be interactive if we call (run)
        #    (begin-time "elapsed time for initialization: "
        #            (init-params p (if reset-fields true false))
        #            (if (string? reset-fields) (load-eigenvectors reset-fields)))
        #    (let ((k-split (list-split k-points k-split-num k-split-index)))
        #      (set-kpoint-index (car k-split))
        #      (if (zero? (car k-split))
        #      (begin
        #            (output-epsilon) ; output epsilon immediately for 1st k block
        #            (if (using-mu?) (output-mu)))) ; and mu too, if we have it
        #      (if (> num-bands 0)
        #      (begin
        #        (map (lambda (k)
        #           (set! current-k k)
        #           (begin-time "elapsed time for k point: " (solve-kpoint k))
        #           (set! all-freqs (cons freqs all-freqs))
        #           (set! band-range-data
        #             (update-band-range-data band-range-data freqs k))
        #           (set! eigensolver-iters
        #             (append eigensolver-iters
        #                 (list (/ iterations num-bands))))
        #           (map (lambda (f)
        #              (if (zero? (procedure-num-args f))
        #                  (f) ; f is a thunk: evaluate once per k-point
        #                  (do ((band 1 (+ band 1))) ((> band num-bands))
        #                    (f band))))
        #                band-functions))
        #         (cdr k-split))
        #        (if (> (length (cdr k-split)) 1)
        #            (begin
        #          (output-band-range-data band-range-data)
        #          (set! gap-list (output-gaps band-range-data)))
        #            (set! gap-list '()))))))))
        #  (set! all-freqs (reverse all-freqs)) ; put them in the right order
        print("done")
