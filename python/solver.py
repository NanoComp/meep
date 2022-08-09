import functools
import math
import numbers
import os
import re
import sys
import time

import h5py
import numpy as np
from meep.geom import init_do_averaging
from meep.simulation import get_num_args
from meep.verbosity_mgr import Verbosity

import meep as mp

from . import mode_solver, with_hermitian_epsilon

try:
    basestring
except NameError:
    basestring = str

U_MIN = 0
U_PROD = 1
U_MEAN = 2

verbosity = Verbosity(mp.cvar, "meep", 1)


class MPBArray(np.ndarray):
    def __new__(cls, input_array, lattice, kpoint=None, bloch_phase=False):
        # Input array is an already formed ndarray instance
        # We first cast to be our class type
        obj = np.asarray(input_array).view(cls)
        # add the new properties to the created instance
        obj.lattice = lattice
        obj.kpoint = kpoint
        obj.bloch_phase = bloch_phase
        # Finally, we must return the newly created object:
        return obj

    def __array_finalize__(self, obj):
        # ``self`` is a new object resulting from
        # ndarray.__new__(MPBArray, ...), therefore it only has
        # attributes that the ndarray.__new__ constructor gave it -
        # i.e. those of a standard ndarray.

        # We could have got to the ndarray.__new__ call in 3 ways:
        # From an explicit constructor - e.g. MPBArray(lattice):
        #    obj is None
        #    (we're in the middle of the MPBArray.__new__
        #    constructor, and self.lattice will be set when we return to
        #    MPBArray.__new__)
        if obj is None:
            return

        # From view casting - e.g arr.view(MPBArray):
        #    obj is arr
        #    (type(obj) can be MPBArray)
        # From new-from-template - e.g mpbarr[:3]
        #    type(obj) is MPBArray
        #
        # Note that it is here, rather than in the __new__ method,
        # that we set the default value for 'lattice', because this
        # method sees all creation of default objects - with the
        # MPBArray.__new__ constructor, but also with
        # arr.view(MPBArray).
        self.lattice = getattr(obj, "lattice", None)
        self.kpoint = getattr(obj, "kpoint", None)
        self.bloch_phase = getattr(obj, "bloch_phase", False)


class ModeSolver:
    def __init__(
        self,
        resolution=10,
        is_negative_epsilon_ok=False,
        eigensolver_flops=0,
        eigensolver_flags=68,
        use_simple_preconditioner=False,
        force_mu=False,
        mu_input_file="",
        epsilon_input_file="",
        mesh_size=3,
        target_freq=0.0,
        tolerance=1.0e-7,
        num_bands=1,
        k_points=None,
        ensure_periodicity=True,
        geometry=None,
        geometry_lattice=mp.Lattice(),
        geometry_center=mp.Vector3(0, 0, 0),
        default_material=mp.Medium(epsilon=1),
        dimensions=3,
        random_fields=False,
        filename_prefix="",
        deterministic=False,
        verbose=False,
        optimize_grid_size=True,
        eigensolver_nwork=3,
        eigensolver_block_size=-11,
    ):

        self.mode_solver = None
        self.resolution = resolution
        self.eigensolver_flags = eigensolver_flags
        self.k_points = k_points or []
        self.geometry = geometry or []
        self.geometry_lattice = geometry_lattice
        self.geometry_center = mp.Vector3(*geometry_center)
        self.default_material = default_material
        self.random_fields = random_fields
        self.filename_prefix = filename_prefix
        self.optimize_grid_size = optimize_grid_size
        self.parity = ""
        self.iterations = 0
        self.all_freqs = None
        self.freqs = []
        self.band_range_data = []
        self.total_run_time = 0
        self.current_k = mp.Vector3()
        self.k_split_num = 1
        self.k_split_index = 0
        self.eigensolver_iters = []

        grid_size = self._adjust_grid_size()

        if type(self.default_material) is not mp.Medium and callable(
            self.default_material
        ):
            init_do_averaging(self.default_material)
            self.default_material.eps = False

        self.mode_solver = mode_solver(
            num_bands,
            self.resolution,
            self.geometry_lattice,
            tolerance,
            mesh_size,
            self.default_material,
            deterministic,
            target_freq,
            dimensions,
            verbose,
            ensure_periodicity,
            eigensolver_flops,
            is_negative_epsilon_ok,
            epsilon_input_file,
            mu_input_file,
            force_mu,
            use_simple_preconditioner,
            grid_size,
            eigensolver_nwork,
            eigensolver_block_size,
        )

    @property
    def num_bands(self):
        return self.mode_solver.num_bands

    @num_bands.setter
    def num_bands(self, val):
        self.mode_solver.set_num_bands(val)

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
            raise TypeError(f"resolution must be a number or a Vector3: Got {t}")

        if self.mode_solver:
            self.mode_solver.resolution = self._resolution
            grid_size = self._adjust_grid_size()
            self.mode_solver.set_grid_size(grid_size)

    @property
    def geometry_lattice(self):
        return self._geometry_lattice

    @geometry_lattice.setter
    def geometry_lattice(self, val):
        self._geometry_lattice = val
        if self.mode_solver:
            self.mode_solver.set_libctl_geometry_lattice(val)
            grid_size = self._adjust_grid_size()
            self.mode_solver.set_grid_size(grid_size)

    @property
    def tolerance(self):
        return self.mode_solver.tolerance

    @tolerance.setter
    def tolerance(self, val):
        self.mode_solver.tolerance = val

    @property
    def mesh_size(self):
        return self.mode_solver.mesh_size

    @mesh_size.setter
    def mesh_size(self, val):
        self.mode_solver.mesh_size = val

    @property
    def deterministic(self):
        return self.mode_solver.deterministic

    @deterministic.setter
    def deterministic(self, val):
        self.mode_solver.deterministic = val

    @property
    def target_freq(self):
        return self.mode_solver.target_freq

    @target_freq.setter
    def target_freq(self, val):
        self.mode_solver.target_freq = val

    @property
    def dimensions(self):
        return self.mode_solver.get_libctl_dimensions()

    @dimensions.setter
    def dimensions(self, val):
        self.mode_solver.set_libctl_dimensions(val)

    @property
    def verbose(self):
        return self.mode_solver.verbose

    @verbose.setter
    def verbose(self, val):
        self.mode_solver.verbose = val

    @property
    def ensure_periodicity(self):
        return self.mode_solver.get_libctl_ensure_periodicity()

    @ensure_periodicity.setter
    def ensure_periodicity(self, val):
        self.mode_solver.set_libctl_ensure_periodicity(val)

    @property
    def eigensolver_flops(self):
        return self.mode_solver.eigensolver_flops

    @eigensolver_flops.setter
    def eigensolver_flops(self, val):
        self.mode_solver.eigensolver_flops = val

    @property
    def is_negative_epsilon_ok(self):
        return self.mode_solver.negative_epsilon_ok

    @is_negative_epsilon_ok.setter
    def is_negative_epsilon_ok(self, val):
        self.mode_solver.negative_epsilon_ok = val

    @property
    def epsilon_input_file(self):
        return self.mode_solver.epsilon_input_file

    @epsilon_input_file.setter
    def epsilon_input_file(self, val):
        self.mode_solver.epsilon_input_file = val

    @property
    def mu_input_file(self):
        return self.mode_solver.mu_input_file

    @mu_input_file.setter
    def mu_input_file(self, val):
        self.mode_solver.mu_input_file = val

    @property
    def force_mu(self):
        return self.mode_solver.force_mu

    @force_mu.setter
    def force_mu(self, val):
        self.mode_solver.force_mu = val

    @property
    def use_simple_preconditioner(self):
        return self.mode_solver.use_simple_preconditioner

    @use_simple_preconditioner.setter
    def use_simple_preconditioner(self, val):
        self.mode_solver.use_simple_preconditioner = val

    @property
    def eigensolver_nwork(self):
        return self.mode_solver.eigensolver_nwork

    @eigensolver_nwork.setter
    def eigensolver_nwork(self, val):
        self.mode_solver.eigensolver_nwork = val

    @property
    def eigensolver_block_size(self):
        return self.mode_solver.eigensolver_block_size

    @eigensolver_block_size.setter
    def eigensolver_block_size(self, val):
        self.mode_solver.eigensolver_block_size = val

    def _adjust_grid_size(self):
        grid_size = self._get_grid_size()

        if self.optimize_grid_size:
            grid_size = self._optimize_grid_size(grid_size)

        return grid_size

    def allow_negative_epsilon(self):
        self.is_negative_epsilon_ok = True
        self.target_freq = 1 / mp.inf

    def get_filename_prefix(self):
        if self.filename_prefix:
            return f"{self.filename_prefix}-"
        _, filename = os.path.split(sys.argv[0])

        return (
            ""
            if filename in ["ipykernel_launcher.py", "__main__.py"]
            else re.sub(r"\.py$", "", filename) + "-"
        )

    def get_freqs(self):
        return self.mode_solver.get_freqs()

    def multiply_bloch_phase(self, arr):
        dims = arr.shape
        arr = arr.ravel()
        self.mode_solver.multiply_bloch_phase(arr)
        return np.reshape(arr, dims)

    def get_poynting(self, which_band):
        e = self.get_efield(which_band, False)
        dims = e.shape
        e = e.ravel()
        h = self.get_hfield(which_band, False).ravel()
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
        self.mode_solver.set_curfield_cmplx(flat_res)
        self.mode_solver.set_curfield_type("v")

        arr = np.reshape(res, dims)
        return MPBArray(arr, self.get_lattice(), self.current_k)

    def get_epsilon(self):
        self.mode_solver.get_epsilon()
        return self.get_curfield_as_array(False)

    def get_mu(self):
        self.mode_solver.get_mu()
        return self.get_curfield_as_array(False)

    def get_bfield(self, which_band, bloch_phase=True):
        return self._get_field("b", which_band, bloch_phase)

    def get_efield(self, which_band, bloch_phase=True):
        return self._get_field("e", which_band, bloch_phase)

    def get_dfield(self, which_band, bloch_phase=True):
        return self._get_field("d", which_band, bloch_phase)

    def get_hfield(self, which_band, bloch_phase=True):
        return self._get_field("h", which_band, bloch_phase)

    def get_charge_density(self, which_band, bloch_phase=True):
        self.get_efield(which_band, bloch_phase)
        self.mode_solver.compute_field_divergence()

    def _get_field(self, f, band, bloch_phase):
        if self.mode_solver is None:
            raise ValueError(
                "Must call a run function before attempting to get a field"
            )

        if f == "b":
            self.mode_solver.get_bfield(band)
        elif f == "d":
            self.mode_solver.get_dfield(band)
        elif f == "e":
            self.mode_solver.get_efield(band)
        elif f == "h":
            self.mode_solver.get_hfield(band)

        dims = self.mode_solver.get_dims()

        while len(dims) < 3:
            dims += [1]

        dims += [3]
        arr = np.zeros(np.prod(dims), np.complex128)

        if bloch_phase:
            self.mode_solver.multiply_bloch_phase()

        self.mode_solver.get_curfield_cmplx(arr)

        arr = np.reshape(arr, dims)
        return MPBArray(
            arr, self.get_lattice(), self.current_k, bloch_phase=bloch_phase
        )

    def get_curfield_as_array(self, bloch_phase=True):
        dims = self.mode_solver.get_dims()
        arr = np.zeros(np.prod(dims))
        self.mode_solver.get_curfield(arr)
        arr = np.reshape(arr, dims)

        return MPBArray(
            arr, self.get_lattice(), self.current_k, bloch_phase=bloch_phase
        )

    def get_dpwr(self, band):
        self.get_dfield(band, False)
        self.compute_field_energy()
        return self.get_curfield_as_array(False)

    def get_bpwr(self, band):
        self.get_bfield(band, False)
        self.compute_field_energy()
        return self.get_curfield_as_array(False)

    def fix_field_phase(self):
        self.mode_solver.fix_field_phase()

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
        epwr = self.get_dpwr(which_band)
        hpwr = self.get_bpwr(which_band)

        tot_pwr = epwr + hpwr

        self.mode_solver.set_curfield(tot_pwr.ravel())
        self.mode_solver.set_curfield_type("R")

        return MPBArray(tot_pwr, self.get_lattice(), self.current_k, bloch_phase=False)

    def get_eigenvectors(self, first_band, num_bands):
        dims = self.mode_solver.get_eigenvectors_slice_dims(num_bands)
        ev = np.zeros(np.prod(dims), dtype=np.complex128)
        self.mode_solver.get_eigenvectors(first_band - 1, num_bands, ev)
        return MPBArray(ev.reshape(dims), self.get_lattice(), self.current_k)

    def set_eigenvectors(self, ev, first_band):
        self.mode_solver.set_eigenvectors(first_band - 1, ev.flatten())

    def save_eigenvectors(self, filename):
        with h5py.File(filename, "w") as f:
            ev = self.get_eigenvectors(1, self.num_bands)
            f["rawdata"] = ev

    def load_eigenvectors(self, filename):
        with h5py.File(filename, "r") as f:
            ev = f["rawdata"][()]
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
            br = brd[0] if brd else ((mp.inf, -1), (-mp.inf, -1))
            br_rest = brd[1:] if brd else []
            newmin = (freqs[0], kpoint) if freqs[0] < br[0][0] else br[0]
            newmax = (freqs[0], kpoint) if freqs[0] > br[1][0] else br[1]
            new_start = br_start + [(newmin, newmax)]
            return update_brd(br_rest, freqs[1:], new_start)

        return update_brd(brd, freqs, [])

    def output_band_range_data(self, br_data):
        if verbosity.mpb > 0:
            fmt = "Band {} range: {} at {} to {} at {}"
            for tup, band in zip(br_data, range(1, len(br_data) + 1)):
                min_band, max_band = tup
                min_freq, min_kpoint = min_band
                max_freq, max_kpoint = max_band
                print(fmt.format(band, min_freq, min_kpoint, max_freq, max_kpoint))

    # Output any gaps in the given band ranges, and return a list of the gaps as
    # a list of (percent, freq-min, freq-max) tuples.
    def output_gaps(self, br_data):
        def ogaps(br_cur, br_rest, i, gaps):
            if not br_rest:
                gaps = list(reversed(gaps))
                return [
                    (gaps[i + 2], gaps[i + 1], gaps[i]) for i in range(0, len(gaps), 3)
                ]
            else:
                br_rest_min_f = br_rest[0][0][0]
                br_cur_max_f = br_cur[1][0]
                if br_cur_max_f >= br_rest_min_f:
                    return ogaps(br_rest[0], br_rest[1:], i + 1, gaps)
                gap_size = (200 * (br_rest_min_f - br_cur_max_f)) / (
                    br_rest_min_f + br_cur_max_f
                )
                if verbosity.mpb > 0:
                    fmt = "Gap from band {} ({}) to band {} ({}), {}%"
                    print(fmt.format(i, br_cur_max_f, i + 1, br_rest_min_f, gap_size))
                return ogaps(
                    br_rest[0],
                    br_rest[1:],
                    i + 1,
                    [gap_size, br_cur_max_f, br_rest_min_f] + gaps,
                )

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

    def get_lattice(self):
        if self.mode_solver is None:
            raise RuntimeError("Must call ModeSolver.run before getting the lattice.")

        lattice = np.zeros((3, 3))
        self.mode_solver.get_lattice(lattice)

        return lattice

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

        if curfield_type in "Rv":
            # Generic scalar/vector field. Don't know k
            output_k = [0, 0, 0]

        if curfield_type in "dhbecv":
            self._output_vector_field(curfield_type, fname_prefix, output_k, component)
        elif curfield_type == "C":
            self._output_complex_scalar_field(fname_prefix, output_k)
        elif curfield_type in "DHBnmR":
            self._output_scalar_field(curfield_type, fname_prefix)
        else:
            raise ValueError(f"Unkown field type: {curfield_type}")

        self.mode_solver.curfield_reset()

    def _output_complex_scalar_field(self, fname_prefix, output_k):
        curfield_type = "C"
        kpoint_index = self.mode_solver.get_kpoint_index()
        curfield_band = self.mode_solver.curfield_band
        fname = f"{curfield_type}.k{kpoint_index:02d}.b{curfield_band:02d}"
        description = "{} field, kpoint {}, band {}, freq={:.6g}".format(
            curfield_type, kpoint_index, curfield_band, self.freqs[curfield_band - 1]
        )
        fname = self._create_fname(fname, fname_prefix, True)
        if verbosity.mpb > 0:
            print(f"Outputting complex scalar field to {fname}...")

        with h5py.File(fname, "w") as f:
            f["description"] = description.encode()
            f["Bloch wavevector"] = np.array(output_k)
            self._write_lattice_vectors(f)

            dims = self.mode_solver.get_dims()
            field = np.empty(np.prod(dims), np.complex128)
            self.mode_solver.get_curfield_cmplx(field)

            reshaped_field = field.reshape(dims)

            f["c.r"] = np.real(reshaped_field)
            f["c.i"] = np.imag(reshaped_field)

    def _output_vector_field(self, curfield_type, fname_prefix, output_k, component):
        components = ["x", "y", "z"]
        kpoint_index = self.mode_solver.get_kpoint_index()
        curfield_band = self.mode_solver.curfield_band
        fname = f"{curfield_type}.k{kpoint_index:02d}.b{curfield_band:02d}"

        if component >= 0:
            fname += f".{components[component]}"

        description = "{} field, kpoint {}, band {}, freq={:.6g}".format(
            curfield_type, kpoint_index, curfield_band, self.freqs[curfield_band - 1]
        )

        fname = self._create_fname(fname, fname_prefix, True)
        if verbosity.mpb > 0:
            print(f"Outputting fields to {fname}...")

        with h5py.File(fname, "w") as f:
            f["description"] = description.encode()
            f["Bloch wavevector"] = np.array(output_k)
            self._write_lattice_vectors(f)

            if curfield_type != "v":
                self.mode_solver.multiply_bloch_phase()

            for c_idx, c in enumerate(components):
                if component >= 0 and c_idx != component:
                    continue

                dims = self.mode_solver.get_dims()
                field = np.empty(np.prod(dims) * 3, np.complex128)
                self.mode_solver.get_curfield_cmplx(field)
                component_field = field[c_idx::3].reshape(dims)

                name = f"{c}.r"
                f[name] = np.real(component_field)
                name = f"{c}.i"
                f[name] = np.imag(component_field)

    def _output_scalar_field(self, curfield_type, fname_prefix):
        components = ["x", "y", "z"]

        if curfield_type == "n":
            fname = "epsilon"
            description = "dielectric function, epsilon"
        elif curfield_type == "m":
            fname = "mu"
            description = "permeability mu"
        else:
            kpoint_index = self.mode_solver.get_kpoint_index()
            curfield_band = self.mode_solver.curfield_band
            fname = "{}pwr.k{:02d}.b{:02d}".format(
                curfield_type.lower(), kpoint_index, curfield_band
            )
            descr_fmt = "{} field energy density, kpoint {}, band {}, freq={:.6g}"
            description = descr_fmt.format(
                curfield_type,
                kpoint_index,
                curfield_band,
                self.freqs[curfield_band - 1],
            )

        parity_suffix = curfield_type not in "mn"
        fname = self._create_fname(fname, fname_prefix, parity_suffix)
        if verbosity.mpb > 0:
            print(f"Outputting {fname}...")

        with h5py.File(fname, "w") as f:
            f["description"] = description.encode()
            self._create_h5_dataset(f, "data")
            self._write_lattice_vectors(f)

            if curfield_type == "n":
                for inv in [False, True]:
                    inv_str = "epsilon_inverse" if inv else "epsilon"
                    for c1 in range(3):
                        for c2 in range(c1, 3):
                            self.mode_solver.get_epsilon_tensor(c1, c2, 0, inv)
                            dataname = f"{inv_str}.{components[c1]}{components[c2]}"
                            self._create_h5_dataset(f, dataname)

                            if with_hermitian_epsilon() and c1 != c2:
                                self.mode_solver.get_epsilon_tensor(c1, c2, 1, inv)
                                dataname += ".i"
                                self._create_h5_dataset(f, dataname)

    def _write_lattice_vectors(self, h5file):
        lattice = np.zeros((3, 3))
        self.mode_solver.get_lattice(lattice)
        h5file["lattice vectors"] = lattice

    def _create_h5_dataset(self, h5file, key):
        h5file[key] = self.get_curfield_as_array(False)

    def _create_fname(self, fname, prefix, parity_suffix):
        parity_str = self.mode_solver.get_parity_string()
        suffix = f".{parity_str}" if parity_suffix and parity_str else ""
        return prefix + fname + suffix + ".h5"

    def compute_field_energy(self):
        return self.mode_solver.compute_field_energy()

    def compute_field_divergence(self):
        return self.mode_solver.compute_field_divergence()

    def compute_energy_in_objects(self, objs):
        return self.mode_solver.compute_energy_in_objects(objs)

    def compute_energy_in_dielectric(self, eps_low, eps_high):
        return self.mode_solver.compute_energy_in_dielectric(eps_low, eps_high)

    def compute_energy_integral(self, f):
        return self.mode_solver.compute_energy_integral(f)

    def compute_field_integral(self, f):
        return self.mode_solver.compute_field_integral(f)

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
        return self.mode_solver.compute_1_group_velocity_component(
            direction, which_band
        )

    def compute_zparities(self):
        return self.mode_solver.compute_zparities()

    def compute_yparities(self):
        return self.mode_solver.compute_yparities()

    def randomize_fields(self):
        self.mode_solver.randomize_fields()

    def display_kpoint_data(self, name, data):
        if verbosity.mpb > 0:
            k_index = self.mode_solver.get_kpoint_index()
            print(f"{self.parity}{name}:, {k_index}", end="")

            for d in data:
                print(f", {d}", end="")
            print()

    def display_eigensolver_stats(self):
        num_runs = len(self.eigensolver_iters)

        if num_runs <= 0:
            return

        min_iters = min(self.eigensolver_iters)
        max_iters = max(self.eigensolver_iters)
        mean_iters = np.mean(self.eigensolver_iters)
        if verbosity.mpb > 0:
            fmt = "eigensolver iterations for {} kpoints: {}-{}, mean = {}"
            print(fmt.format(num_runs, min_iters, max_iters, mean_iters), end="")

        sorted_iters = sorted(self.eigensolver_iters)
        idx1 = num_runs // 2
        idx2 = ((num_runs + 1) // 2) - 1
        median_iters = 0.5 * (sorted_iters[idx1] + sorted_iters[idx2])
        if verbosity.mpb > 0:
            print(f", median = {median_iters}")

        mean_flops = self.eigensolver_flops / (num_runs * mean_iters)
        if verbosity.mpb > 0:
            print(f"mean flops per iteration = {mean_flops}")

        mean_time = self.total_run_time / (mean_iters * num_runs)
        if verbosity.mpb > 0:
            print(f"mean time per iteration = {mean_time} s")

    def _get_grid_size(self):
        grid_size = mp.Vector3(
            self.resolution[0] * self.geometry_lattice.size.x,
            self.resolution[1] * self.geometry_lattice.size.y,
            self.resolution[2] * self.geometry_lattice.size.z,
        )

        grid_size.x = max(math.ceil(grid_size.x), 1)
        grid_size.y = max(math.ceil(grid_size.y), 1)
        grid_size.z = max(math.ceil(grid_size.z), 1)

        return grid_size

    def _optimize_grid_size(self, grid_size):
        grid_size.x = self.next_factor2357(grid_size.x)
        grid_size.y = self.next_factor2357(grid_size.y)
        grid_size.z = self.next_factor2357(grid_size.z)
        return grid_size

    def next_factor2357(self, n):
        def is_factor2357(n):
            def divby(n, p):
                return divby(n // p, p) if n % p == 0 else n

            return divby(divby(divby(divby(n, 2), 3), 5), 7) == 1

        if is_factor2357(n):
            return n
        return self.next_factor2357(n + 1)

    def init_params(self, p, reset_fields):
        self.mode_solver.init(p, reset_fields, self.geometry, self.default_material)

    def set_parity(self, p):
        self.mode_solver.set_parity(p)

    def solve_kpoint(self, k):
        self.mode_solver.solve_kpoint(k)
        self.freqs = self.get_freqs()

    def run_parity(self, p, reset_fields, *band_functions):
        if self.random_fields and self.randomize_fields not in band_functions:
            band_functions.append(self.randomize_fields)

        start = time.time()

        self.all_freqs = np.zeros((len(self.k_points), self.num_bands))
        self.band_range_data = []

        init_time = time.time()

        if verbosity.mpb > 0:
            print("Initializing eigensolver data")
            print(f"Computing {self.num_bands} bands with {self.tolerance} tolerance")

        self.init_params(p, reset_fields)

        if isinstance(reset_fields, basestring):
            self.load_eigenvectors(reset_fields)

        if verbosity.mpb > 0:
            print(f"{len(self.k_points)} k-points")

            for kp in self.k_points:
                print(f"  {kp}")

            print(f"elapsed time for initialization: {time.time() - init_time}")

        # TODO: Split over multiple processes
        # k_split = list_split(self.k_points, self.k_split_num, self.k_split_index)
        k_split = (0, self.k_points)
        self.mode_solver.set_kpoint_index(k_split[0])

        if self.num_bands > 0:
            for i, k in enumerate(k_split[1]):
                self.current_k = k
                solve_kpoint_time = time.time()
                self.solve_kpoint(k)
                self.iterations = self.mode_solver.get_iterations()
                if verbosity.mpb > 0:
                    print(
                        f"elapsed time for k point: {time.time() - solve_kpoint_time}"
                    )
                self.all_freqs[i, :] = np.array(self.freqs)
                self.band_range_data = self.update_band_range_data(
                    self.band_range_data, self.freqs, k
                )
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
                        raise ValueError(
                            "Band function should take 1 or 2 arguments. "
                            "The first must be a ModeSolver instance"
                        )

            if len(k_split[1]) > 1:
                self.output_band_range_data(self.band_range_data)
                self.gap_list = self.output_gaps(self.band_range_data)
            else:
                self.gap_list = []

        end = time.time() - start
        if verbosity.mpb > 0:
            print(f"total elapsed time for run: {end}")
        self.total_run_time += end
        self.eigensolver_flops = self.mode_solver.get_eigensolver_flops()
        self.parity = self.mode_solver.get_parity_string()
        if verbosity.mpb > 0:
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

    def find_k(
        self,
        p,
        omega,
        band_min,
        band_max,
        korig_and_kdir,
        tol,
        kmag_guess,
        kmag_min,
        kmag_max,
        *band_funcs,
    ):

        num_bands_save = self.num_bands
        kpoints_save = self.k_points
        nb = band_max - band_min + 1
        kdir = korig_and_kdir[1] if type(korig_and_kdir) is list else korig_and_kdir
        lat = self.geometry_lattice
        kdir1 = mp.cartesian_to_reciprocal(
            mp.reciprocal_to_cartesian(kdir, lat).unit(), lat
        )

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
                    if verbosity.mpb > 0:
                        print(f"find-k {b} at {k}: {tab_val[0]} (cached)")
                    return tab_val
                else:
                    self.num_bands = b
                    self.k_points = [korig + kdir1.scale(k)]
                    self.run_parity(p, False)
                    v = self.mode_solver.compute_group_velocity_component(kdir1)

                    # Cache computed values
                    for _b, _f, _v in zip(
                        range(band_min, b - band_min + 1),
                        self.freqs[band_min - 1 :],
                        v[band_min - 1 :],
                    ):
                        tabval = bktab.get((_b, k0s[_b - band_min]), None)

                        if not tabval or abs(_f - omega) < abs(tabval[0]):
                            k0s[_b - band_min + 1] = k

                        bktab[(_b, k)] = (_f - omega, _v)

                    fun = self.freqs[-1] - omega
                    if verbosity.mpb > 0:
                        print(f"find-k {b} at {k}: {fun}")
                    return (fun, v[-1])

            return _rootfun

        # Don't let previous computations interfere
        if self.mode_solver:
            self.randomize_fields()

        ks = []
        for b in range(band_max, band_max - nb, -1):
            ks.append(
                mp.find_root_deriv(
                    rootfun(b), tol, kmag_min, kmag_max, k0s[b - band_min]
                )
            )

        if band_funcs:
            for b, k in zip(range(1, band_max + 1), reversed(ks)):
                self.num_bands = b
                self.k_points = [korig + kdir1.scale(k)]

                def bfunc(ms, b_prime):
                    if b_prime == b:
                        for f in band_funcs:
                            apply_band_func_thunk(ms, f, b, True)

                self.run_parity(p, False, bfunc)

        self.num_bands = num_bands_save
        self.k_points = kpoints_save
        ks = list(reversed(ks))
        if verbosity.mpb > 0:
            print(
                f"{self.parity}kvals:, {omega}, {band_min}, {band_max}",
                end="",
            )
            for k in korig:
                print(f", {k}", end="")
            for k in kdir1:
                print(f", {k}", end="")
            for k in ks:
                print(f", {k}", end="")
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
            return try_plus(k + v, v) if n(k + v) < n(k) else k

        def _try(k, v):
            return try_plus(try_plus(k, v), mp.Vector3() - v)

        try_list = [
            mp.Vector3(1, 0, 0),
            mp.Vector3(0, 1, 0),
            mp.Vector3(0, 0, 1),
            mp.Vector3(0, 1, 1),
            mp.Vector3(1, 0, 1),
            mp.Vector3(1, 1, 0),
            mp.Vector3(0, 1, -1),
            mp.Vector3(1, 0, -1),
            mp.Vector3(1, -1, 0),
            mp.Vector3(1, 1, 1),
            mp.Vector3(-1, 1, 1),
            mp.Vector3(1, -1, 1),
            mp.Vector3(1, 1, -1),
        ]

        def try_all(k):
            return functools.reduce(_try, try_list, k)

        def try_all_and_repeat(k):
            knew = try_all(k)
            return try_all_and_repeat(knew) if n(knew) < n(k) else k

        k0 = k - mp.Vector3(*[round(x) for x in k])

        return try_all_and_repeat(k0) if n(k0) < n(k) else try_all_and_repeat(k)

    def get_dominant_planewave(self, band):
        return self.mode_solver.get_dominant_planewave(band)

    def transformed_overlap(self, W, w):
        return self.mode_solver.transformed_overlap(W, w)

    def compute_symmetry(self, band, W, w):
        return self.mode_solver.compute_symmetry(band, W, w)

    def compute_symmetries(self, W, w):
        return [
            self.mode_solver.compute_symmetry(band, W, w)
            for band in range(1, self.num_bands + 1)
        ]


# Predefined output functions (functions of the band index), for passing to `run`


def output_hfield(ms, which_band):
    ms.get_hfield(which_band, False)
    ms.output_field()


def output_hfield_x(ms, which_band):
    ms.get_hfield(which_band, False)
    ms.output_field_x()


def output_hfield_y(ms, which_band):
    ms.get_hfield(which_band, False)
    ms.output_field_y()


def output_hfield_z(ms, which_band):
    ms.get_hfield(which_band, False)
    ms.output_field_z()


def output_bfield(ms, which_band):
    ms.get_bfield(which_band, False)
    ms.output_field()


def output_bfield_x(ms, which_band):
    ms.get_bfield(which_band, False)
    ms.output_field_x()


def output_bfield_y(ms, which_band):
    ms.get_bfield(which_band, False)
    ms.output_field_y()


def output_bfield_z(ms, which_band):
    ms.get_bfield(which_band, False)
    ms.output_field_z()


def output_dfield(ms, which_band):
    ms.get_dfield(which_band, False)
    ms.output_field()


def output_dfield_x(ms, which_band):
    ms.get_dfield(which_band, False)
    ms.output_field_x()


def output_dfield_y(ms, which_band):
    ms.get_dfield(which_band, False)
    ms.output_field_y()


def output_dfield_z(ms, which_band):
    ms.get_dfield(which_band, False)
    ms.output_field_z()


def output_efield(ms, which_band):
    ms.get_efield(which_band, False)
    ms.output_field()


def output_efield_x(ms, which_band):
    ms.get_efield(which_band, False)
    ms.output_field_x()


def output_efield_y(ms, which_band):
    ms.get_efield(which_band, False)
    ms.output_field_y()


def output_efield_z(ms, which_band):
    ms.get_efield(which_band, False)
    ms.output_field_z()


def output_bpwr(ms, which_band):
    ms.get_bfield(which_band, False)
    ms.compute_field_energy()
    ms.output_field()


def output_dpwr(ms, which_band):
    ms.get_dfield(which_band, False)
    ms.compute_field_energy()
    ms.output_field()


def output_tot_pwr(ms, which_band):
    ms.get_tot_pwr(which_band)
    ms.output_field_to_file(-1, f"{ms.get_filename_prefix()}tot.")


def output_dpwr_in_objects(output_func, min_energy, objects=[]):
    """
    The following function returns an output function that calls output_func for
    bands with D energy in objects > min-energy. For example,
    output_dpwr_in_objects(output_dfield, 0.20, some_object) would return an
    output function that would spit out the D field for bands with at least %20
    of their D energy in some-object.
    """

    def _output(ms, which_band):
        ms.get_dfield(which_band, False)
        ms.compute_field_energy()
        energy = ms.compute_energy_in_objects(objects)
        if verbosity.mpb > 0:
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
    ms.output_field_to_file(-1, f"{ms.get_filename_prefix()}flux.")


def output_poynting_x(ms, which_band):
    ms.get_poynting(which_band)
    ms.output_field_to_file(0, f"{ms.get_filename_prefix()}flux.")


def output_poynting_y(ms, which_band):
    ms.get_poynting(which_band)
    ms.output_field_to_file(1, f"{ms.get_filename_prefix()}flux.")


def output_poynting_z(ms, which_band):
    ms.get_poynting(which_band)
    ms.output_field_to_file(2, f"{ms.get_filename_prefix()}flux.")


def display_yparities(ms):
    ms.display_kpoint_data("yparity", ms.mode_solver.compute_yparities())


def display_zparities(ms):
    ms.display_kpoint_data("zparity", ms.mode_solver.compute_zparities())


def display_group_velocities(ms):
    ms.display_kpoint_data("velocity", ms.compute_group_velocities())


# Band functions to pick a canonical phase for the eigenstate of the
# given band based upon the spatial representation of the given field


def fix_hfield_phase(ms, which_band):
    ms.get_hfield(which_band, False)
    ms.mode_solver.fix_field_phase()


def fix_bfield_phase(ms, which_band):
    ms.get_bfield(which_band, False)
    ms.mode_solver.fix_field_phase()


def fix_dfield_phase(ms, which_band):
    ms.get_dfield(which_band, False)
    ms.mode_solver.fix_field_phase()


def fix_efield_phase(ms, which_band):
    ms.get_efield(which_band, False)
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
