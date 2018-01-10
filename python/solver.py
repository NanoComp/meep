from __future__ import division

import time

# import numpy as np

import meep as mp
from meep import mpb
from meep.simulation import get_num_args
# from meep.source import check_positive

try:
    basestring
except NameError:
    basestring = str

U_MIN = 0
U_PROD = 1
U_SUM = 2


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


# TODO: Placeholders
def init_params(p, reset_fields):
    pass


def load_eigenvectors(fields):
    pass


def randomize_fields():
    pass


def list_split(l, num, index):
    pass


def set_kpoint_index(index):
    pass


def output_epsilon():
    pass


def using_mu():
    pass


def output_mu():
    pass


def solve_kpoint(k):
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
        self.parity = 0
        self.iterations = 0
        self.freqs = []
        self.all_freqs = []
        self.band_range_data = []
        self.eigensolver_flops = 0
        self.total_run_time = 0
        self.current_k = mp.Vector3()
        self.k_split_num = 1
        self.k_split_index = 0
        self.eigensolver_iters = []
        self.iterations = 0
        self.match_frequency = False

        self.mode_solver = mpb.mode_solver(self.num_bands, self.match_frequency, self.parity,
                                           self.resolution, self.geometry_lattice.size,
                                           self.tolerance)

    def update_band_range_data(self, brd, freqs, kpoint):
        pass

    def output_band_range_data(self):
        pass

    def output_gaps(self):
        pass

    def run_parity(self, p, reset_fields, *band_functions):
        if self.randomize_fields and randomize_fields not in band_functions:
            band_functions.append(randomize_fields)

        start = time.time()

        self.all_freqs = []
        self.band_range_data = []

        init_time = time.time()

        self.mode_solver.init(p, True if reset_fields else False)
        if isinstance(reset_fields, basestring):
            load_eigenvectors(reset_fields)

        print("elapsed time for initialization: {}".format(time.time() - init_time))
        ###### !!!!!!!!! #####
        return

        # TODO: Split over multiple processes
        k_split = list_split(self.k_points, self.k_split_num, self.k_split_index)
        set_kpoint_index(k_split[0])

        if k_split[0] == 0:
            output_epsilon()  # output epsilon immediately for 1st k block
            if using_mu():
                output_mu()  # and mu too, if we have it

        if self.num_bands > 0:
            for k in k_split[1]:
                self.current_k = k
                solve_kpoint_time = time.time()
                solve_kpoint(k)
                print("elapsed time for k point: {}".format(time.time() - solve_kpoint_time))

                # TODO: Make freqs an output var
                self.all_freqs.append(self.freqs)
                self.band_range_data = self.update_band_range_data(self.band_range_data,
                                                                   self.freqs, k)
                self.eigensolver_iters += [self.iterations / self.num_bands]

                for f in band_functions:
                    if get_num_args(f) == 0:
                        f()
                    else:
                        band = 1
                        while band < self.num_bands:
                            f(band)
                            band += 1

            if len(k_split[1]) > 1:
                self.output_band_range_data()
                self.gap_list = self.output_gaps()
            else:
                self.gap_list = []

        end = time.time() - start
        print("total elapsed time for run: {}".format(end))
        # TODO: Does this keep track of all runs? When should it be reset?
        self.total_run_time += end
        print("done")

    def run(self, *band_functions):
        self.run_parity(mp.NO_PARITY, True, band_functions)

    def run_zeven(self, *band_functions):
        self.run_parity(mp.EVEN_Z, True, band_functions)

    def run_zodd(self, *band_functions):
        self.run_parity(mp.ODD_Z, True, band_functions)

    def run_yeven(self, *band_functions):
        self.run_parity(mp.EVEN_Y, True, band_functions)

    def run_yodd(self, *band_functions):
        self.run_parity(mp.ODD_Y, True, band_functions)

    def run_yeven_zeven(self, *band_functions):
        self.run_parity(mp.EVEN_Y + mp.EVEN_Z, True, band_functions)

    def run_yeven_zodd(self, *band_functions):
        self.run_parity(mp.EVEN_Y + mp.ODD_Z, True, band_functions)

    def run_yodd_zeven(self, *band_functions):
        self.run_parity(mp.ODD_Y + mp.EVEN_Z, True, band_functions)

    def run_yodd_zodd(self, *band_functions):
        self.run_parity(mp.ODD_Y + mp.ODD_Z, True, band_functions)

    run_te = run_zeven
    run_tm = run_zodd
    run_te_yeven = run_yeven_zeven
    run_te_yodd = run_yodd_zeven
    run_tm_yeven = run_yeven_zodd
    run_tm_yodd = run_yodd_zodd
