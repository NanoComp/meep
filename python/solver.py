from __future__ import division

import functools
import operator
import time

import meep as mp
from meep import mpb
from meep.simulation import get_num_args

try:
    basestring
except NameError:
    basestring = str

U_MIN = 0
U_PROD = 1
U_SUM = 2


# TODO: Add tests
class Matrix(object):

    def __init__(self, c1=mp.Vector3(), c2=mp.Vector3(), c3=mp.Vector3()):
        self.c1 = c1
        self.c2 = c2
        self.c3 = c3

    def __getitem__(self, i):
        if i == 0:
            return self.c1
        elif i == 1:
            return self.c2
        elif i == 2:
            return self.c3
        else:
            raise IndexError("No value at index {}".format(i))

    def __mul__(self, m):
        # TODO: Matrix multiplication
        pass

    def __repr__(self):
        return "<{}\n {}\n {}>".format(self.row(0), self.row(1), self.row(2))

    def row(self, i):
        return mp.Vector3(self.c1[i], self.c2[i], self.c3[i])

    def determinant(self):
        sum1 = sum([
            functools.reduce(operator.mul, [self[x][x] for x in range(3)]),
            functools.reduce(operator.mul, [self[0][1], self[1][2], self[2][0]]),
            functools.reduce(operator.mul, [self[1][0], self[2][1], self[0][2]])
        ])
        sum2 = sum([
            functools.reduce(operator.mul, [self[0][2], self[1][1], self[2][0]]),
            functools.reduce(operator.mul, [self[0][1], self[1][0], self[2][2]]),
            functools.reduce(operator.mul, [self[1][2], self[2][1], self[0][0]])
        ])
        return sum1 - sum2

    def transpose(self):
        return Matrix(self.row(0), self.row(1), self.row(2))


# TODO: Add tests
class Lattice(object):

    def __init__(self,
                 size=mp.Vector3(1, 1, 1),
                 basis_size=mp.Vector3(1, 1, 1),
                 basis1=mp.Vector3(1, 0, 0),
                 basis2=mp.Vector3(0, 1, 0),
                 basis3=mp.Vector3(0, 0, 1)):

        self.size = size
        self.basis_size = basis_size
        self.basis1 = basis1
        self.basis2 = basis2
        self.basis3 = basis3

        @property
        def basis1(self):
            return self._basis1

        @basis1.setter
        def basis1(self, val):
            self._basis1 = val.unit()

        @property
        def basis2(self):
            return self._basis2

        @basis2.setter
        def basis2(self, val):
            self._basis2 = val.unit()

        @property
        def basis3(self):
            return self._basis3

        @basis3.setter
        def basis3(self, val):
            self._basis3 = val.unit()

        @property
        def b1(self):
            return self.basis1.scale(self.basis_size.x)

        @property
        def b2(self):
            return self.basis2.scale(self.basis_size.y)

        @property
        def b3(self):
            return self.basis3.scale(self.basis_size.z)

        @property
        def basis(self):
            B = Matrix(self.b1, self.b2, self.b3)

            if B.determinant() == 0:
                raise ValueError("Lattice basis vectors must be linearly independent.")

            return B

        @property
        def metric(self):
            B = self.basis
            return B.transpose() * B


# TODO: Placeholders
def load_eigenvectors(fields):
    pass


def output_mu():
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
                 randomize_fields=False,
                 filename_prefix=''):

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
        self.filename_prefix = filename_prefix
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
        self.mode_solver = None

    # The band-range-data is a list of tuples, each consisting of a (min, k-point)
    # tuple and a (max, k-point) tuple, with each min/max pair describing the
    # frequency range of a band and the k-points where it achieves its minimum/maximum.
    # Here, we update this data with a new list of band frequencies, and return the new
    # data.  If band-range-data is null or too short, the needed entries will be created.
    def update_band_range_data(self, brd, freqs, kpoint):

        def update_brd(brd, freqs, br_start):
            if not freqs:
                return br_start + brd
            else:
                br = [(mp.inf, -1), (-mp.inf, -1)] if not brd else brd[0]
                br_rest = [] if not brd else brd[1:]
                newmin = (freqs[0], kpoint) if freqs[0] < br[0][0] else br[0]
                newmax = (freqs[0], kpoint) if freqs[0] > br[0][1] else br[1]
                return update_brd(br_rest, freqs[1:], [(newmin, newmax)] + br_start)

        return update_brd(brd, freqs, [])

    def output_band_range_data(self, br_data):
        for tup, band in zip(br_data, range(len(br_data))):
            fmt = "Band {} range: {} at {} to {} at {}"
            min_freq, min_kpoint, max_freq, max_kpoint = tup
            print(fmt.format(band, min_freq, min_kpoint, max_freq, max_kpoint))

    # Output any gaps in the given band ranges, and return a list of the gaps as
    # a list of (percent, freq-min, freq-max) tuples.
    def output_gaps(self, br_data):

        def ogaps(br_cur, br_rest, i, gaps):
            if not br_rest:
                return gaps
            else:
                if br_cur[0][1] >= br_rest[0][0][0]:
                    return ogaps(br_rest[0], br_rest[1:], i + 1, gaps)
                else:
                    gap_size = 5
                    fmt = "Gap from band {} ({}) to band {} ({}), {}%"
                    print(fmt.format(i, br_cur[0][1], i + 1, br_rest[0][0][0], gap_size))
                    return ogaps(br_rest[0], br_rest[1:], i + 1, [gap_size, br_cur[0][1], br_rest[0][0][0]] + gaps)
        if not br_data:
            return []
        else:
            return ogaps(br_data[0], br_data[1:], 1, [])

    def output_epsilon(self):
        self.mode_solver.get_epsilon()
        # TODO
        # self.mode_solver.output_field_to_file(-1, self.get_filename_prefix)

    def output_mu(self):
        pass
        # TODO
        # self.mode_solver.get_mu()
        # TODO
        # self.mode_solver.output_field_to_file(-1, self.get_filename_prefix)

    def randomize_fields(self):
        self.mode_solver.randomize_fields()

    def run_parity(self, p, reset_fields, *band_functions):
        if self.randomize_fields and self.randomize_fields not in band_functions:
            band_functions.append(self.randomize_fields)

        start = time.time()

        self.all_freqs = []
        self.band_range_data = []

        init_time = time.time()

        print("Initializing eigensolver data")
        print("Computing {} bands with {} tolerance".format(self.num_bands, self.tolerance))

        # TODO: Can we keep the mode_solver around between runs, or does it need
        # to get created clean for each run?
        self.mode_solver = mpb.mode_solver(self.num_bands, p, self.resolution,
                                           self.geometry_lattice, self.tolerance,
                                           self.default_material, self.geometry,
                                           True if reset_fields else False)

        if isinstance(reset_fields, basestring):
            load_eigenvectors(reset_fields)

        print("{} k-points".format(len(self.k_points)))

        for kp in self.k_points:
            print("  {}".format(kp))

        print("elapsed time for initialization: {}".format(time.time() - init_time))

        # TODO: Split over multiple processes
        # k_split = list_split(self.k_points, self.k_split_num, self.k_split_index)
        k_split = (0, self.k_points)
        self.mode_solver.set_kpoint_index(k_split[0])

        if k_split[0] == 0:
            self.output_epsilon()  # output epsilon immediately for 1st k block
            if self.mode_solver.using_mu():
                self.output_mu()  # and mu too, if we have it

        if self.num_bands > 0:
            for k in k_split[1]:
                self.current_k = k
                solve_kpoint_time = time.time()
                self.mode_solver.solve_kpoint(k)
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
                        while band <= self.num_bands:
                            f(band)
                            band += 1

            if len(k_split[1]) > 1:
                self.output_band_range_data(self.band_range_data)
                self.gap_list = self.output_gaps(self.band_range_data)
            else:
                self.gap_list = []

        end = time.time() - start
        print("total elapsed time for run: {}".format(end))
        self.total_run_time += end
        print("done")

    def run(self, *band_functions):
        self.run_parity(mp.NO_PARITY, True, *band_functions)

    def run_zeven(self, *band_functions):
        self.run_parity(mp.EVEN_Z, True, *band_functions)

    def run_zodd(self, *band_functions):
        self.run_parity(mp.ODD_Z, True, *band_functions)

    def run_yeven(self, *band_functions):
        self.run_parity(mp.EVEN_Y, True, *band_functions)

    def run_yodd(self, *band_functions):
        self.run_parity(mp.ODD_Y, True, *band_functions)

    def run_yeven_zeven(self, *band_functions):
        self.run_parity(mp.EVEN_Y + mp.EVEN_Z, True, *band_functions)

    def run_yeven_zodd(self, *band_functions):
        self.run_parity(mp.EVEN_Y + mp.ODD_Z, True, *band_functions)

    def run_yodd_zeven(self, *band_functions):
        self.run_parity(mp.ODD_Y + mp.EVEN_Z, True, *band_functions)

    def run_yodd_zodd(self, *band_functions):
        self.run_parity(mp.ODD_Y + mp.ODD_Z, True, *band_functions)

    run_te = run_zeven
    run_tm = run_zodd
    run_te_yeven = run_yeven_zeven
    run_te_yodd = run_yodd_zeven
    run_tm_yeven = run_yeven_zodd
    run_tm_yodd = run_yodd_zodd
