from __future__ import division, print_function

import functools
import os
import numbers
import re
import sys
import time

import h5py
import numpy as np
import meep as mp
from . import mode_solver
from meep.simulation import get_num_args

try:
    basestring
except NameError:
    basestring = str

U_MIN = 0
U_PROD = 1
U_SUM = 2


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
                 geometry_lattice=mp.Lattice(),
                 geometry_center=mp.Vector3(0, 0, 0),
                 default_material=mp.Medium(epsilon=1),
                 dimensions=3,
                 random_fields=False,
                 filename_prefix='',
                 deterministic=False,
                 verbose=False):

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
        self.random_fields = random_fields
        self.filename_prefix = filename_prefix
        self.deterministic = deterministic
        self.verbose = verbose
        self.parity = ''
        self.iterations = 0
        self.all_freqs = []
        self.freqs = []
        self.band_range_data = []
        self.eigensolver_flops = 0
        self.total_run_time = 0
        self.current_k = mp.Vector3()
        self.k_split_num = 1
        self.k_split_index = 0
        self.eigensolver_iters = []
        self.mode_solver = None

    @property
    def resolution(self):
        return self._resolution

    @resolution.setter
    def resolution(self, val):
        if isinstance(val, numbers.Number):
            self._resolution = [val, val, val]
        elif isinstance(val, mp.Vector3):
            self._resolution = [val.x, val.y, val.z]
        else:
            t = type(val)
            raise TypeError("resolution must be a number or a Vector3: Got {}".format(t))

    def get_filename_prefix(self):
        if self.filename_prefix:
            return self.filename_prefix + '-'
        else:
            _, filename = os.path.split(sys.argv[0])

            if filename == 'ipykernel_launcher.py' or filename == '__main__.py':
                return ''
            else:
                return re.sub(r'\.py$', '', filename) + '-'

    def get_freqs(self):
        return self.mode_solver.get_freqs()

    def get_poynting(self, which_band):
        e = self.get_efield(which_band)
        h = self.get_hfield(which_band)
        # Reshape into rows of vector3s
        e = e.reshape((int(e.shape[0] / 3), 3))
        h = h.reshape((int(h.shape[0] / 3), 3))

        res = np.zeros(e.shape, dtype=np.complex128)

        def ExH(e, h):
            ev = mp.Vector3(e[0], e[1], e[2])
            hv = mp.Vector3(h[0], h[1], h[2])
            return ev.conj().cross(hv)

        for i in range(e.shape[0]):
            res[i] = np.array(ExH(e[i], h[i]))

        flat_res = res.ravel()
        # We only have to set curfield here to do `multiply_bloch_phase` in
        # _output_vector_field.
        self.mode_solver.set_curfield_cmplx(flat_res)
        self.mode_solver.set_curfield_type('v')

        return flat_res

    def get_epsilon(self):
        self.mode_solver.get_epsilon()

        dims = self.mode_solver.get_dims()
        eps = np.empty(np.prod(dims))
        self.mode_solver.get_curfield(eps)

        return np.reshape(eps, dims)

    def get_bfield(self, which_band):
        return self._get_field('b', which_band)

    def get_efield(self, which_band):
        return self._get_field('e', which_band)

    def get_dfield(self, which_band):
        return self._get_field('d', which_band)

    def get_hfield(self, which_band):
        return self._get_field('h', which_band)

    def get_charge_density(self, which_band):
        self.get_efield(which_band)
        self.mode_solver.compute_field_divergence()

    def _get_field(self, f, band):
        if self.mode_solver is None:
            raise ValueError("Must call a run function before attempting to get a field")

        size = self.mode_solver.get_field_size()
        field = np.zeros(size, dtype=np.complex128)

        if f == 'b':
            self.mode_solver.get_bfield(field, band)
        elif f == 'd':
            self.mode_solver.get_dfield(field, band)
        elif f == 'e':
            self.mode_solver.get_efield(field, band)
        elif f == 'h':
            self.mode_solver.get_hfield(field, band)

        return field

    def get_epsilon_point(self, p):
        return self.mode_solver.get_epsilon_point(p)

    def get_epsilon_inverse_tensor_point(self, p):
        return self.mode_solver.get_epsilon_inverse_tensor_point(p)

    def get_energy_point(self, p):
        return self.mode_solver.get_energy_point(p)

    def get_field_point(self, p):
        return self.mode_solver.get_field_point(p)

    def get_bloch_field_point(self, p):
        return self.mode_solver.get_bloch_field_point(p)

    def get_tot_pwr(self, which_band):
        self.get_dfield(which_band)
        self.compute_field_energy()

        dims = self.mode_solver.get_dims()
        epwr = np.zeros(np.prod(dims))
        self.mode_solver.get_curfield(epwr)
        epwr = np.reshape(epwr, dims)

        self.get_bfield(which_band)
        self.compute_field_energy()

        hpwr = np.zeros(np.prod(dims))
        self.mode_solver.get_curfield(hpwr)
        hpwr = np.reshape(hpwr, dims)

        tot_pwr = epwr + hpwr

        self.mode_solver.set_curfield(tot_pwr.ravel())
        self.mode_solver.set_curfield_type('R')

        return tot_pwr

    def get_eigenvectors(self, first_band, num_bands):
        dims = self.mode_solver.get_eigenvectors_slice_dims(num_bands)
        ev = np.zeros(np.prod(dims), dtype=np.complex128)
        self.mode_solver.get_eigenvectors(first_band - 1, num_bands, ev)
        return ev.reshape(dims)

    def set_eigenvectors(self, ev, first_band):
        self.mode_solver.set_eigenvectors(first_band - 1, ev.flatten())

    def save_eigenvectors(self, filename):
        with h5py.File(filename, 'w') as f:
            ev = self.get_eigenvectors(1, self.num_bands)
            f['rawdata'] = ev

    def load_eigenvectors(self, filename):
        with h5py.File(filename, 'r') as f:
            ev = f['rawdata'].value
            self.set_eigenvectors(ev, 1)
            self.mode_solver.curfield_reset()

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
                br = ((mp.inf, -1), (-mp.inf, -1)) if not brd else brd[0]
                br_rest = [] if not brd else brd[1:]
                newmin = (freqs[0], kpoint) if freqs[0] < br[0][0] else br[0]
                newmax = (freqs[0], kpoint) if freqs[0] > br[1][0] else br[1]
                new_start = br_start + [(newmin, newmax)]
                return update_brd(br_rest, freqs[1:], new_start)

        return update_brd(brd, freqs, [])

    def output_band_range_data(self, br_data):
        for tup, band in zip(br_data, range(1, len(br_data) + 1)):
            fmt = "Band {} range: {} at {} to {} at {}"
            min_band, max_band = tup
            min_freq, min_kpoint = min_band
            max_freq, max_kpoint = max_band
            print(fmt.format(band, min_freq, min_kpoint, max_freq, max_kpoint))

    # Output any gaps in the given band ranges, and return a list of the gaps as
    # a list of (percent, freq-min, freq-max) tuples.
    def output_gaps(self, br_data):

        def ogaps(br_cur, br_rest, i, gaps):
            if not br_rest:
                ordered_gaps = []
                gaps = list(reversed(gaps))
                for i in range(0, len(gaps), 3):
                    ordered_gaps.append((gaps[i + 2], gaps[i + 1], gaps[i]))
                return ordered_gaps
            else:
                br_rest_min_f = br_rest[0][0][0]
                br_cur_max_f = br_cur[1][0]
                if br_cur_max_f >= br_rest_min_f:
                    return ogaps(br_rest[0], br_rest[1:], i + 1, gaps)
                else:
                    gap_size = ((200 * (br_rest_min_f - br_cur_max_f)) /
                                (br_rest_min_f + br_cur_max_f))
                    fmt = "Gap from band {} ({}) to band {} ({}), {}%"
                    print(fmt.format(i, br_cur_max_f, i + 1, br_rest_min_f, gap_size))
                    return ogaps(br_rest[0], br_rest[1:], i + 1,
                                 [gap_size, br_cur_max_f, br_rest_min_f] + gaps)
        if not br_data:
            return []
        else:
            return ogaps(br_data[0], br_data[1:], 1, [])

    # Return the frequency gap from the band #lower-band to the band
    # #(lower-band+1), as a percentage of mid-gap frequency.  The "gap"
    # may be negative if the maximum of the lower band is higher than the
    # minimum of the upper band.  (The gap is computed from the
    # band-range-data of the previous run.)
    def retrieve_gap(self, lower_band):
        if lower_band + 1 > len(self.band_range_data):
            raise ValueError("retrieve-gap called for higher band than was calculated")

        f1 = self.band_range_data[lower_band - 1][1][0]
        f2 = self.band_range_data[lower_band][0][0]

        return (f2 - f1) / (0.005 * (f1 + f2))

    # Split a list L into num more-or-less equal pieces, returning the piece
    # given by index (in 0..num-1), along with the index in L of the first
    # element of the piece, as a list: [first-index, piece-of-L]
    def list_split(self, l, num, index):

        def list_sub(l, start, length, index, rest):
            if not l:
                return list(reversed(rest))
            if index >= start and index < (start + length):
                return list_sub(l[1:], start, length, index + 1, [l[0]] + rest)
            else:
                return list_sub(l[1:], start, length, index + 1, rest)

        if index >= num or index < 0:
            return (len(l), [])
        else:
            block_size = (len(l) + num - 1) // num
            start = index * block_size
            length = min(block_size, (len(l) - index * block_size))
            return (start, list_sub(l, start, length, 0, []))

    def output_field(self):
        self.output_field_to_file(mp.ALL, self.get_filename_prefix())

    def output_field_x(self):
        self.output_field_to_file(0, self.get_filename_prefix())

    def output_field_y(self):
        self.output_field_to_file(1, self.get_filename_prefix())

    def output_field_z(self):
        self.output_field_to_file(2, self.get_filename_prefix())

    def output_epsilon(self):
        self.mode_solver.get_epsilon()
        self.output_field_to_file(mp.ALL, self.get_filename_prefix())

    def output_mu(self):
        self.mode_solver.get_mu()
        self.output_field_to_file(mp.ALL, self.get_filename_prefix())

    def output_field_to_file(self, component, fname_prefix):
        curfield_type = self.mode_solver.get_curfield_type()
        output_k = self.mode_solver.get_output_k()

        if curfield_type in 'Rv':
            # Generic scalar/vector field. Don't know k
            output_k = [0, 0, 0]

        if curfield_type in 'dhbecv':
            self._output_vector_field(curfield_type, fname_prefix, output_k, component)
        elif curfield_type == 'C':
            self._output_complex_scalar_field(fname_prefix, output_k)
        elif curfield_type in 'DHBnmR':
            self._output_scalar_field(curfield_type, fname_prefix)
        else:
            raise ValueError("Unkown field type: {}".format(curfield_type))

        self.mode_solver.curfield_reset()

    def _output_complex_scalar_field(self, fname_prefix, output_k):
        curfield_type = 'C'
        kpoint_index = self.mode_solver.get_kpoint_index()
        curfield_band = self.mode_solver.curfield_band
        fname = "{}.k{:02d}.b{:02d}".format(curfield_type, kpoint_index, curfield_band)
        description = "{} field, kpoint {}, band {}, freq={:.6g}".format(
            curfield_type,
            kpoint_index,
            curfield_band,
            self.freqs[curfield_band - 1]
        )
        fname = self._create_fname(fname, fname_prefix, True)
        print("Outputting complex scalar field to {}...".format(fname))

        with h5py.File(fname, 'w') as f:
            f['description'] = description.encode()
            f['Bloch wavevector'] = np.array(output_k)
            self._write_lattice_vectors(f)

            dims = self.mode_solver.get_dims()
            field = np.empty(np.prod(dims), np.complex128)
            self.mode_solver.get_curfield_cmplx(field)

            reshaped_field = field.reshape(dims)

            f['c.r'] = np.real(reshaped_field)
            f['c.i'] = np.imag(reshaped_field)

    def _output_vector_field(self, curfield_type, fname_prefix, output_k, component):
        components = ['x', 'y', 'z']
        kpoint_index = self.mode_solver.get_kpoint_index()
        curfield_band = self.mode_solver.curfield_band
        fname = "{}.k{:02d}.b{:02d}".format(curfield_type, kpoint_index, curfield_band)

        if component >= 0:
            fname += ".{}".format(components[component])

        description = "{} field, kpoint {}, band {}, freq={:.6g}".format(
            curfield_type,
            kpoint_index,
            curfield_band,
            self.freqs[curfield_band - 1]
        )

        fname = self._create_fname(fname, fname_prefix, True)
        print("Outputting fields to {}...".format(fname))

        with h5py.File(fname, 'w') as f:
            f['description'] = description.encode()
            f['Bloch wavevector'] = np.array(output_k)
            self._write_lattice_vectors(f)

            self.mode_solver.multiply_bloch_phase()

            for c_idx, c in enumerate(components):
                if component >= 0 and c_idx != component:
                    continue

                dims = self.mode_solver.get_dims()
                field = np.empty(np.prod(dims) * 3, np.complex128)
                self.mode_solver.get_curfield_cmplx(field)
                component_field = field[c_idx::3].reshape(dims)

                name = "{}.r".format(c)
                f[name] = np.real(component_field)
                name = "{}.i".format(c)
                f[name] = np.imag(component_field)

    def _output_scalar_field(self, curfield_type, fname_prefix):
        components = ['x', 'y', 'z']

        if curfield_type == 'n':
            fname = 'epsilon'
            description = 'dielectric function, epsilon'
        elif curfield_type == 'm':
            fname = 'mu'
            description = 'permeability mu'
        else:
            kpoint_index = self.mode_solver.get_kpoint_index()
            curfield_band = self.mode_solver.curfield_band
            fname = "{}pwr.k{:02d}.b{:02d}".format(curfield_type.lower(),
                                                   kpoint_index, curfield_band)
            descr_fmt = "{} field energy density, kpoint {}, band {}, freq={:.6g}"
            description = descr_fmt.format(curfield_type, kpoint_index, curfield_band,
                                           self.freqs[curfield_band - 1])

        parity_suffix = False if curfield_type in 'mn' else True
        fname = self._create_fname(fname, fname_prefix, parity_suffix)
        print("Outputting {}...".format(fname))

        with h5py.File(fname, 'w') as f:
            f['description'] = description.encode()
            self._create_h5_dataset(f, 'data')
            self._write_lattice_vectors(f)

            if curfield_type == 'n':
                for inv in [False, True]:
                    inv_str = 'epsilon_inverse' if inv else 'epsilon'
                    for c1 in range(3):
                        for c2 in range(c1, 3):
                            self.mode_solver.get_epsilon_tensor(c1, c2, 0, inv)
                            dataname = "{}.{}{}".format(inv_str, components[c1],
                                                        components[c2])
                            self._create_h5_dataset(f, dataname)

                            if self.mode_solver.with_hermitian_epsilon() and c1 != c2:
                                self.mode_solver.get_epsilon_tensor(c1, c2, 1, inv)
                                dataname += '.i'
                                self._create_h5_dataset(f, dataname)

    def _write_lattice_vectors(self, h5file):
        lattice = np.zeros((3, 3))
        self.mode_solver.get_lattice(lattice)
        h5file['lattice vectors'] = lattice

    def _create_h5_dataset(self, h5file, key):
        dims = self.mode_solver.get_dims()
        arr = np.zeros(np.prod(dims))
        self.mode_solver.get_curfield(arr)
        h5file[key] = np.reshape(arr, dims)

    def _create_fname(self, fname, prefix, parity_suffix):
        parity_str = self.mode_solver.get_parity_string()
        if parity_suffix and parity_str:
            suffix = ".{}".format(parity_str)
        else:
            suffix = ''

        return prefix + fname + suffix + '.h5'

    def compute_field_energy(self):
        return self.mode_solver.compute_field_energy()

    def compute_field_divergence(self):
        mode_solver.compute_field_divergence()

    def compute_energy_in_objects(self, objs):
        return self.mode_solver.compute_energy_in_objects(objs)

    def compute_energy_in_dielectric(self, eps_low, eps_high):
        return self.mode_solver.compute_energy_in_dielectric(eps_low, eps_high)

    def compute_group_velocities(self):
        xarg = mp.cartesian_to_reciprocal(mp.Vector3(1), self.geometry_lattice)
        vx = self.mode_solver.compute_group_velocity_component(xarg)
        yarg = mp.cartesian_to_reciprocal(mp.Vector3(y=1), self.geometry_lattice)
        vy = self.mode_solver.compute_group_velocity_component(yarg)
        zarg = mp.cartesian_to_reciprocal(mp.Vector3(z=1), self.geometry_lattice)
        vz = self.mode_solver.compute_group_velocity_component(zarg)

        return [mp.Vector3(x, y, z) for x, y, z in zip(vx, vy, vz)]

    def compute_group_velocity_component(self, direction):
        return self.mode_solver.compute_group_velocity_component(direction)

    def compute_one_group_velocity(self, which_band):
        return self.mode_solver.compute_1_group_velocity(which_band)

    def compute_one_group_velocity_component(self, direction, which_band):
        return self.mode_solver.compute_1_group_velocity_component(direction,
                                                                   which_band)

    def randomize_fields(self):
        self.mode_solver.randomize_fields()

    def display_kpoint_data(self, name, data):
        k_index = self.mode_solver.get_kpoint_index()
        print("{}{}:, {}".format(self.parity, name, k_index), end='')

        for d in data:
            print(", {}".format(d), end='')
        print()

    def display_eigensolver_stats(self):
        num_runs = len(self.eigensolver_iters)

        if num_runs <= 0:
            return

        min_iters = min(self.eigensolver_iters)
        max_iters = max(self.eigensolver_iters)
        mean_iters = np.mean(self.eigensolver_iters)
        fmt = "eigensolver iterations for {} kpoints: {}-{}, mean = {}"
        print(fmt.format(num_runs, min_iters, max_iters, mean_iters), end='')

        sorted_iters = sorted(self.eigensolver_iters)
        idx1 = num_runs // 2
        idx2 = ((num_runs + 1) // 2) - 1
        median_iters = 0.5 * (sorted_iters[idx1] + sorted_iters[idx2])
        print(", median = {}".format(median_iters))

        mean_flops = self.eigensolver_flops / (num_runs * mean_iters)
        print("mean flops per iteration = {}".format(mean_flops))

        mean_time = self.total_run_time / (mean_iters * num_runs)
        print("mean time per iteration = {} s".format(mean_time))

    def run_parity(self, p, reset_fields, *band_functions):
        if self.random_fields and self.randomize_fields not in band_functions:
            band_functions.append(self.randomize_fields)

        start = time.time()

        self.all_freqs = []
        self.band_range_data = []

        init_time = time.time()

        print("Initializing eigensolver data")
        print("Computing {} bands with {} tolerance".format(self.num_bands, self.tolerance))

        if type(self.default_material) is not mp.Medium and callable(self.default_material):
            # TODO: Support epsilon_function user materials like meep?
            self.default_material.eps = False

        self.mode_solver = mode_solver(
            self.num_bands,
            p,
            self.resolution,
            self.geometry_lattice,
            self.tolerance,
            self.mesh_size,
            self.default_material,
            self.geometry,
            True if reset_fields else False,
            self.deterministic,
            self.target_freq,
            self.dimensions,
            self.verbose,
            self.ensure_periodicity,
            self.eigensolver_flops,
            self.is_negative_epsilon_ok,
            self.epsilon_input_file,
            self.mu_input_file,
            self.force_mu,
        )

        if isinstance(reset_fields, basestring):
            self.load_eigenvectors(reset_fields)

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
                self.iterations = self.mode_solver.get_iterations()
                print("elapsed time for k point: {}".format(time.time() - solve_kpoint_time))
                self.freqs = self.get_freqs()
                self.all_freqs.append(self.freqs)
                self.band_range_data = self.update_band_range_data(self.band_range_data,
                                                                   self.freqs, k)
                self.eigensolver_iters += [self.iterations / self.num_bands]

                for f in band_functions:
                    num_args = get_num_args(f)
                    if num_args == 1:
                        f(self)
                    elif num_args == 2:
                        band = 1
                        while band <= self.num_bands:
                            f(self, band)
                            band += 1
                    else:
                        raise ValueError("Band function should take 1 or 2 arguments. "
                                         "The first must be a ModeSolver instance")

            if len(k_split[1]) > 1:
                self.output_band_range_data(self.band_range_data)
                self.gap_list = self.output_gaps(self.band_range_data)
            else:
                self.gap_list = []

        end = time.time() - start
        print("total elapsed time for run: {}".format(end))
        self.total_run_time += end
        self.eigensolver_flops = self.mode_solver.get_eigensolver_flops()
        self.parity = self.mode_solver.get_parity_string()
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

    def find_k(self, p, omega, band_min, band_max, korig_and_kdir, tol,
               kmag_guess, kmag_min, kmag_max, *band_funcs):
        num_bands_save = self.num_bands
        kpoints_save = self.k_points
        nb = band_max - band_min + 1
        kdir = korig_and_kdir[1] if type(korig_and_kdir) is list else korig_and_kdir
        lat = self.geometry_lattice
        kdir1 = mp.cartesian_to_reciprocal(mp.reciprocal_to_cartesian(kdir, lat).unit(), lat)

        if type(korig_and_kdir) is list:
            korig = korig_and_kdir[0]
        else:
            korig = mp.Vector3()

        # k0s is a list caching the best k value found for each band:
        if type(kmag_guess) is list:
            k0s = kmag_guess
        else:
            k0s = [kmag_guess] * (band_max - band_min + 1)

        # dict to memoize all "band: k" results
        bktab = {}

        def rootfun(b):

            def _rootfun(k):
                # First, look in the cached table
                tab_val = bktab.get((b, k), None)
                if tab_val:
                    print("find-k {} at {}: {} (cached)".format(b, k, tab_val[0]))
                    return tab_val
                # Otherwise, compute bands and cache results
                else:
                    self.num_bands = b
                    self.k_points = [korig + kdir1.scale(k)]
                    self.run_parity(p, False)
                    v = self.mode_solver.compute_group_velocity_component(kdir1)

                    # Cache computed values
                    for _b, _f, _v in zip(range(band_min, b - band_min + 1, 1),
                                          self.freqs[band_min - 1:],
                                          v[band_min - 1:]):
                        tabval = bktab.get((_b, k0s[_b - band_min]), None)

                        if not tabval or abs(_f - omega) < abs(tabval[0]):
                            k0s[_b - band_min + 1] = k

                        bktab[(_b, k)] = (_f - omega, _v)

                    fun = self.freqs[-1] - omega
                    print("find-k {} at {}: {}".format(b, k, fun))
                    return (fun, v[-1])

            return _rootfun

        # Don't let previous computations interfere
        if self.mode_solver:
            self.randomize_fields()

        ks = []
        for b in range(band_max, band_max - nb, -1):
            ks.append(mp.find_root_deriv(rootfun(b), tol, kmag_min, kmag_max, k0s[b - band_min]))

        if band_funcs:
            for b, k in zip(range(band_max, band_max - nb, -1), reversed(ks)):
                self.num_bands = b
                self.k_points = [korig + kdir1.scale(k)]

                def bfunc(ms, b_prime):
                    if b_prime == b:
                        for f in band_funcs:
                            apply_band_func_thunk(ms, f, b, True)

                self.run_parity(p, False, bfunc)

        self.num_bands = num_bands_save
        self.k_points = kpoints_save
        ks = reversed(ks)
        print("{}kvals:, {}, {}, {}".format(self.parity, omega, band_min, band_max), end='')
        for k in korig:
            print(", {}".format(k), end='')
        for k in kdir1:
            print(", {}".format(k), end='')
        for k in ks:
            print(", {}".format(k), end='')
        print()

        return ks

    def first_brillouin_zone(self, k):
        """
        Function to convert a k-point k into an equivalent point in the
        first Brillouin zone (not necessarily the irreducible Brillouin zone)
        """
        def n(k):
            return mp.reciprocal_to_cartesian(k, self.geometry_lattice).norm()

        def try_plus(k, v):
            if n(k + v) < n(k):
                return try_plus(k + v, v)
            else:
                return k

        def _try(k, v):
            return try_plus(try_plus(k, v), mp.Vector3() - v)

        try_list = [
            mp.Vector3(1, 0, 0), mp.Vector3(0, 1, 0), mp.Vector3(0, 0, 1),
            mp.Vector3(0, 1, 1), mp.Vector3(1, 0, 1), mp.Vector3(1, 1, 0),
            mp.Vector3(0, 1, -1), mp.Vector3(1, 0, -1), mp.Vector3(1, -1, 0),
            mp.Vector3(1, 1, 1), mp.Vector3(-1, 1, 1), mp.Vector3(1, -1, 1),
            mp.Vector3(1, 1, -1),
        ]

        def try_all(k):
            return functools.reduce(_try, try_list, k)

        def try_all_and_repeat(k):
            knew = try_all(k)
            return try_all_and_repeat(knew) if n(knew) < n(k) else k

        k0 = k - mp.Vector3(*[round(x) for x in k])

        return try_all_and_repeat(k0) if n(k0) < n(k) else try_all_and_repeat(k)


# Predefined output functions (functions of the band index), for passing to `run`

def output_hfield(ms, which_band):
    ms.get_hfield(which_band)
    ms.output_field()


def output_hfield_x(ms, which_band):
    ms.get_hfield(which_band)
    ms.output_field_x()


def output_hfield_y(ms, which_band):
    ms.get_hfield(which_band)
    ms.output_field_y()


def output_hfield_z(ms, which_band):
    ms.get_hfield(which_band)
    ms.output_field_z()


def output_bfield(ms, which_band):
    ms.get_bfield(which_band)
    ms.output_field()


def output_bfield_x(ms, which_band):
    ms.get_bfield(which_band)
    ms.output_field_x()


def output_bfield_y(ms, which_band):
    ms.get_bfield(which_band)
    ms.output_field_y()


def output_bfield_z(ms, which_band):
    ms.get_bfield(which_band)
    ms.output_field_z()


def output_dfield(ms, which_band):
    ms.get_dfield(which_band)
    ms.output_field()


def output_dfield_x(ms, which_band):
    ms.get_dfield(which_band)
    ms.output_field_x()


def output_dfield_y(ms, which_band):
    ms.get_dfield(which_band)
    ms.output_field_y()


def output_dfield_z(ms, which_band):
    ms.get_dfield(which_band)
    ms.output_field_z()


def output_efield(ms, which_band):
    ms.get_efield(which_band)
    ms.output_field()


def output_efield_x(ms, which_band):
    ms.get_efield(which_band)
    ms.output_field_x()


def output_efield_y(ms, which_band):
    ms.get_efield(which_band)
    ms.output_field_y()


def output_efield_z(ms, which_band):
    ms.get_efield(which_band)
    ms.output_field_z()


def output_bpwr(ms, which_band):
    ms.get_bfield(which_band)
    ms.compute_field_energy()
    ms.output_field()


def output_dpwr(ms, which_band):
    ms.get_dfield(which_band)
    ms.compute_field_energy()
    ms.output_field()


def output_tot_pwr(ms, which_band):
    ms.get_tot_pwr(which_band)
    ms.output_field_to_file(-1, ms.get_filename_prefix() + 'tot.')


def output_dpwr_in_objects(output_func, min_energy, objects=[]):
    """
    The following function returns an output function that calls output_func for
    bands with D energy in objects > min-energy. For example,
    output_dpwr_in_objects(output_dfield, 0.20, some_object) would return an
    output function that would spit out the D field for bands with at least %20
    of their D energy in some-object.
    """

    def _output(ms, which_band):
        ms.get_dfield(which_band)
        ms.compute_field_energy()
        energy = ms.compute_energy_in_objects(objects)
        fmt = "dpwr:, {}, {}, {} "
        print(fmt.format(which_band, ms.freqs[which_band - 1], energy))
        if energy >= min_energy:
            apply_band_func(ms, output_func, which_band)

    return _output


def output_charge_density(ms, which_band):
    ms.get_charge_density(which_band)
    ms.output_field_to_file(-1, ms.get_filename_prefix())


def output_poynting(ms, which_band):
    ms.get_poynting(which_band)
    ms.output_field_to_file(-1, ms.get_filename_prefix() + 'flux.')


def output_poynting_x(ms, which_band):
    ms.get_poynting(which_band)
    ms.output_field_to_file(0, ms.get_filename_prefix() + 'flux.')


def output_poynting_y(ms, which_band):
    ms.get_poynting(which_band)
    ms.output_field_to_file(1, ms.get_filename_prefix() + 'flux.')


def output_poynting_z(ms, which_band):
    ms.get_poynting(which_band)
    ms.output_field_to_file(2, ms.get_filename_prefix() + 'flux.')


def display_yparities(ms):
    ms.display_kpoint_data('yparity', ms.mode_solver.compute_yparities())


def display_zparities(ms):
    ms.display_kpoint_data('zparity', ms.mode_solver.compute_zparities())


def display_group_velocities(ms):
    ms.display_kpoint_data('velocity', ms.compute_group_velocities())


# Band functions to pick a canonical phase for the eigenstate of the
# given band based upon the spatial representation of the given field

def fix_hfield_phase(ms, which_band):
    ms.get_hfield(which_band)
    ms.mode_solver.fix_field_phase()


def fix_bfield_phase(ms, which_band):
    ms.get_bfield(which_band)
    ms.mode_solver.fix_field_phase()


def fix_dfield_phase(ms, which_band):
    ms.get_dfield(which_band)
    ms.mode_solver.fix_field_phase()


def fix_efield_phase(ms, which_band):
    ms.get_efield(which_band)
    ms.mode_solver.fix_field_phase()


def apply_band_func_thunk(ms, band_func, which_band, eval_thunk):
    """
    We need a special function to evaluate band functions, since band functions
    can either be a function of the band number or a thunk (function of no arguments,
    evaluated once per k-point).
    """
    if get_num_args(band_func) == 1:
        if eval_thunk:
            band_func(ms)  # evaluate thunks once per k-point
    else:
        band_func(ms, which_band)


def apply_band_func(ms, band_func, which_band):
    apply_band_func_thunk(ms, band_func, which_band, which_band == 1)


def combine_band_functions(*band_funcs):
    """Combines zero or more band functions into one"""

    def _combine(ms, which_band):
        for f in band_funcs:
            apply_band_func(ms, f, which_band)
    return _combine


def output_at_kpoint(kpoint, *band_funcs):
    """Only invoke the given band functions for the specified k-point"""

    band_func = combine_band_functions(*band_funcs)

    def _output_at_kpoint(ms, which_band):
        if ms.current_k.close(kpoint, tol=1e-8 * kpoint.norm()):
            band_func(ms, which_band)

    return _output_at_kpoint
