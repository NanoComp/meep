import math

import numpy as np

import meep as mp

from . import MPBArray, map_data


class MPBData:

    TWOPI = 6.2831853071795864769252867665590057683943388

    def __init__(
        self,
        lattice=None,
        kpoint=None,
        rectify=False,
        x=0,
        y=0,
        z=0,
        periods=0,
        resolution=0,
        phase_angle=0,
        pick_nearest=False,
        ve=None,
        verbose=False,
    ):

        self.lattice = lattice
        self.kpoint = kpoint
        self.rectify = rectify

        if periods:
            self.multiply_size = [periods, periods, periods]
        else:
            self.multiply_size = [x or 1, y or 1, z or 1]

        self.resolution = resolution
        self.phase_angle = phase_angle
        self.pick_nearest = pick_nearest
        self.ve = ve

        if self.ve:
            self.have_ve = True
            self.rectify = True
        else:
            self.have_ve = False
            self.ve = mp.Vector3(1, 0, 0)

        self.verbose = verbose
        self.scaleby = complex(1, 0)

        self.phase = complex(
            math.cos(self.TWOPI * self.phase_angle / 360.0),
            math.sin(self.TWOPI * self.phase_angle / 360.0),
        )
        self.scaleby *= self.phase

    def handle_dataset(self, in_arr):

        out_dims = [1, 1, 1]
        rank = len(in_arr.shape)
        num_ones = 3 - rank
        in_dims = list(in_arr.shape) + [1] * num_ones

        if np.iscomplexobj(in_arr):
            in_arr_re = np.real(in_arr)
            in_arr_im = np.imag(in_arr)
        else:
            in_arr_re = in_arr
            in_arr_im = None

        if self.verbose:
            fmt = "Input data is rank {}, size {}x{}x{}."
            print(fmt.format(rank, in_dims[0], in_dims[1], in_dims[2]))

        if self.resolution > 0:
            out_dims[0] = math.floor(self.Rout.c1.norm() * self.resolution + 0.5)
            out_dims[1] = math.floor(self.Rout.c2.norm() * self.resolution + 0.5)
            out_dims[2] = math.floor(self.Rout.c3.norm() * self.resolution + 0.5)
        else:
            for i in range(3):
                out_dims[i] = in_dims[i] * self.multiply_size[i]

        for i in range(rank, 3):
            out_dims[i] = 1

        N = 1
        for i in range(3):
            out_dims[i] = int(max(out_dims[i], 1))
            N *= out_dims[i]

        if self.verbose:
            print(f"Output data {out_dims[0]}x{out_dims[1]}x{out_dims[2]}")

        out_arr_re = np.zeros(int(N))

        if isinstance(in_arr_im, np.ndarray):
            out_arr_im = np.zeros(int(N))
        else:
            out_arr_im = np.array([])

        flat_in_arr_re = in_arr_re.ravel()
        flat_in_arr_im = (
            in_arr_im.ravel() if isinstance(in_arr_im, np.ndarray) else np.array([])
        )

        kvector = [self.kpoint.x, self.kpoint.y, self.kpoint.z] if self.kpoint else []
        map_data(
            flat_in_arr_re,
            flat_in_arr_im,
            np.array(in_dims, dtype=np.intc),
            out_arr_re,
            out_arr_im,
            np.array(out_dims, dtype=np.intc),
            self.coord_map,
            kvector,
            self.pick_nearest,
            self.verbose,
            False,
        )

        if np.iscomplexobj(in_arr):
            # multiply * scaleby for complex data
            complex_out = np.vectorize(complex)(out_arr_re, out_arr_im)
            complex_out *= self.scaleby

            return np.reshape(complex_out, out_dims[:rank])

        return np.reshape(out_arr_re, out_dims[:rank])

    def handle_cvector_dataset(self, in_arr, multiply_bloch_phase):
        in_x_re = np.real(in_arr[:, :, :, 0]).ravel()
        in_x_im = np.imag(in_arr[:, :, :, 0]).ravel()
        in_y_re = np.real(in_arr[:, :, :, 1]).ravel()
        in_y_im = np.imag(in_arr[:, :, :, 1]).ravel()
        in_z_re = np.real(in_arr[:, :, :, 2]).ravel()
        in_z_im = np.imag(in_arr[:, :, :, 2]).ravel()

        d_in = [[in_x_re, in_x_im], [in_y_re, in_y_im], [in_z_re, in_z_im]]
        in_dims = [in_arr.shape[0], in_arr.shape[1], 1]
        if self.verbose:
            print("Found complex vector dataset...")

        if self.verbose:
            fmt = "Input data is rank {}, size {}x{}x{}."
            rank = 2

            print(fmt.format(rank, in_dims[0], in_dims[1], in_dims[2]))

        # rotate vector field according to cart_map
        if self.verbose:
            fmt1 = "Rotating vectors by matrix [ {:.10g}, {:.10g}, {:.10g}"
            fmt2 = "                             {:.10g}, {:.10g}, {:.10g}"
            fmt3 = "                             {:.10g}, {:.10g}, {:.10g} ]"
            print(
                fmt1.format(self.cart_map.c1.x, self.cart_map.c2.x, self.cart_map.c3.x)
            )
            print(
                fmt2.format(self.cart_map.c1.y, self.cart_map.c2.y, self.cart_map.c3.y)
            )
            print(
                fmt3.format(self.cart_map.c1.z, self.cart_map.c2.z, self.cart_map.c3.z)
            )

        N = in_dims[0] * in_dims[1]
        for ri in range(2):
            for i in range(N):
                v = mp.Vector3(d_in[0][ri][i], d_in[1][ri][i], d_in[2][ri][i])
                v = self.cart_map * v
                d_in[0][ri][i] = v.x
                d_in[1][ri][i] = v.y
                d_in[2][ri][i] = v.z

        out_dims = [1, 1, 1]

        if self.resolution > 0:
            out_dims[0] = self.Rout.c1.norm() * self.resolution + 0.5
            out_dims[1] = self.Rout.c2.norm() * self.resolution + 0.5
            out_dims[2] = self.Rout.c3.norm() * self.resolution + 0.5
        else:
            for i in range(3):
                out_dims[i] = in_dims[i] * self.multiply_size[i]

        out_dims[2] = 1

        N = 1
        for i in range(3):
            out_dims[i] = int(max(out_dims[i], 1))
            N *= out_dims[i]

        if self.verbose:
            fmt = "Output data {}x{}x{}."
            print(fmt.format(out_dims[0], out_dims[1], out_dims[2]))

        kvector = [self.kpoint.x, self.kpoint.y, self.kpoint.z] if self.kpoint else []
        converted = []
        for dim in range(3):
            out_arr_re = np.zeros(int(N))
            out_arr_im = np.zeros(int(N))

            map_data(
                d_in[dim][0].ravel(),
                d_in[dim][1].ravel(),
                np.array(in_dims, dtype=np.intc),
                out_arr_re,
                out_arr_im,
                np.array(out_dims, dtype=np.intc),
                self.coord_map,
                kvector,
                self.pick_nearest,
                self.verbose,
                multiply_bloch_phase,
            )

            # multiply * scaleby
            complex_out = np.vectorize(complex)(out_arr_re, out_arr_im)
            complex_out *= self.scaleby
            converted.append(complex_out)

        result = np.zeros(np.prod(out_dims) * 3, np.complex128)
        result[0::3] = converted[0]
        result[1::3] = converted[1]
        result[2::3] = converted[2]

        return np.reshape(result, (out_dims[0], out_dims[1], 3))

    def init_output_lattice(self):

        cart_map = mp.Matrix(
            mp.Vector3(1, 0, 0), mp.Vector3(0, 1, 0), mp.Vector3(0, 0, 1)
        )

        Rin = mp.Matrix(
            mp.Vector3(*self.lattice[0]),
            mp.Vector3(*self.lattice[1]),
            mp.Vector3(*self.lattice[2]),
        )

        if self.verbose:
            print("Read lattice vectors")
            if self.kpoint:
                fmt = "Read Bloch wavevector ({:.6g}, {:.6g}, {:.6g})"
                print(fmt.format(self.kpoint.x, self.kpoint.y, self.kpoint.z))
            fmt = "Input lattice = ({:.6g}, {:.6g}, {:.6g}), ({:.6g}, {:.6g}, {:.6g}), ({:.6g}, {:.6g}, {:.6g})"
            print(
                fmt.format(
                    Rin.c1.x,
                    Rin.c1.y,
                    Rin.c1.z,
                    Rin.c2.x,
                    Rin.c2.y,
                    Rin.c2.z,
                    Rin.c3.x,
                    Rin.c3.y,
                    Rin.c3.z,
                )
            )

        Rout = mp.Matrix(Rin.c1, Rin.c2, Rin.c3)

        if self.rectify:
            # Orthogonalize the output lattice vectors.  If have_ve is true,
            # then the first new lattice vector should be in the direction
            # of the ve unit vector; otherwise, the first new lattice vector
            # is the first original lattice vector. Note that we do this in
            # such a way as to preserve the volume of the unit cell, and so
            # that our first vector (in the direction of ve) smoothly
            # interpolates between the original lattice vectors.

            ve = self.ve.unit() if self.have_ve else Rout.c1.unit()
            # First, compute c1 in the direction of ve by smoothly
            # interpolating the old c1/c2/c3 (formula is slightly tricky)
            V = Rout.c1.cross(Rout.c2).dot(Rout.c3)
            Rout.c2 = Rout.c2 - Rout.c1
            Rout.c3 = Rout.c3 - Rout.c1
            Rout.c1 = ve.scale(V / Rout.c2.cross(Rout.c3).dot(ve))

            # Now, orthogonalize c2 and c3
            Rout.c2 = Rout.c2 - ve.scale(ve.dot(Rout.c2))
            Rout.c3 = Rout.c3 - ve.scale(ve.dot(Rout.c3))
            Rout.c3 = Rout.c3 - Rout.c2.scale(
                Rout.c2.dot(Rout.c3) / Rout.c2.dot(Rout.c2)
            )

            cart_map.c1 = Rout.c1.unit()
            cart_map.c2 = Rout.c2.unit()
            cart_map.c3 = Rout.c3.unit()
            cart_map = cart_map.inverse()

        Rout.c1 = Rout.c1.scale(self.multiply_size[0])
        Rout.c2 = Rout.c2.scale(self.multiply_size[1])
        Rout.c3 = Rout.c3.scale(self.multiply_size[2])

        if self.verbose:
            fmt = "Output lattice = ({:.6g}, {:.6g}, {:.6g}), ({:.6g}, {:.6g}, {:.6g}), ({:.6g}, {:.6g}, {:.6g})"
            print(
                fmt.format(
                    Rout.c1.x,
                    Rout.c1.y,
                    Rout.c1.z,
                    Rout.c2.x,
                    Rout.c2.y,
                    Rout.c2.z,
                    Rout.c3.x,
                    Rout.c3.y,
                    Rout.c3.z,
                )
            )

        self.coord_map = Rin.inverse() * Rout
        self.Rout = Rout
        self.cart_map = cart_map

    def convert(self, arr, kpoint=None):
        if isinstance(arr, MPBArray):
            self.lattice = arr.lattice
            self.kpoint = arr.kpoint

        if self.lattice is None:
            err = (
                "Couldn't find 'lattice.' You must do one of the following:\n"
                + "  1. Pass the ModeSolver lattice to the MPBData constructor\n"
                + "     i.e., MPBData(lattice=ms.get_lattice())\n"
                + "  2. Create an MPBArray to pass to MPBData.convert()\n"
                + "     i.e., mpb_arr = MPBArray(arr, ms.get_lattice(), ... ); mpb_data.convert(mpb_arr))"
            )
            raise ValueError(err)

        if kpoint:
            self.kpoint = kpoint

        self.init_output_lattice()

        if len(arr.shape) == 4:
            return self.handle_cvector_dataset(arr, not arr.bloch_phase)
        else:
            return self.handle_dataset(arr)
