from __future__ import division

import math
# import re

import numpy as np
import meep as mp
from . import map_data


class MPBData(object):

    TWOPI = 6.2831853071795864769252867665590057683943388

    def __init__(self,
                 lattice,
                 rectify=False,
                 x=0,
                 y=0,
                 z=0,
                 periods=0,
                 resolution=0,
                 phase_angle=0,
                 pick_nearest=False,
                 transpose=False,
                 ve=None,
                 verbose=False):

        self.lattice = lattice
        self.rectify = rectify

        if periods:
            self.multiply_size = [periods, periods, periods]
        else:
            self.multiply_size = [
                x if x else 1,
                y if y else 1,
                z if z else 1
            ]

        self.resolution = resolution
        self.phase_angle = phase_angle
        self.pick_nearest = pick_nearest
        self.transpose = transpose
        self.ve = ve

        if self.ve:
            self.have_ve = True
            self.rectify = True
        else:
            self.have_ve = False
            self.ve = mp.Vector3(1, 0, 0)

        self.verbose = verbose
        self.scaleby = complex(1, 0)

        self.phase = complex(math.cos(self.TWOPI * self.phase_angle / 360.0),
                             math.sin(self.TWOPI * self.phase_angle / 360.0))
        self.scaleby *= self.phase

    def handle_dataset(self, in_arr, kvector):

        out_dims = [1, 1, 1]
        rank = len(in_arr.shape)
        num_ones = 3 - rank
        in_dims = [x for x in in_arr.shape] + [1] * num_ones

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
            out_dims[i] = max(out_dims[i], 1)
            N *= out_dims[i]

        out_dims2 = [1, 1, 1]

        if self.transpose:
            out_dims2[0] = int(out_dims[1])
            out_dims2[1] = int(out_dims[0])
            out_dims2[2] = int(out_dims[2])
        else:
            out_dims2[0] = int(out_dims[0])
            out_dims2[1] = int(out_dims[1])
            out_dims2[2] = int(out_dims[2])

        if self.verbose:
            print("Output data {}x{}x{}".format(out_dims2[0], out_dims2[1], out_dims2[2]))

        out_arr_re = np.zeros(int(N))

        if isinstance(in_arr_im, np.ndarray):
            out_arr_im = np.zeros(int(N))
        else:
            out_arr_im = np.array([])

        flat_in_arr_re = in_arr_re.ravel()
        flat_in_arr_im = in_arr_im.ravel() if isinstance(in_arr_im, np.ndarray) else np.array([])

        if kvector:
            kvector = [kvector.x, kvector.y, kvector.z]
        else:
            kvector = []

        map_data(flat_in_arr_re, flat_in_arr_im, np.array(in_dims, dtype=np.intc),
                 out_arr_re, out_arr_im, np.array(out_dims2, dtype=np.intc), self.coord_map,
                 kvector, self.pick_nearest, self.verbose)

        if np.iscomplexobj(in_arr):
            # multiply * scaleby for complex data
            complex_out = np.vectorize(complex)(out_arr_re, out_arr_im)
            complex_out *= self.scaleby

            return np.reshape(complex_out, out_dims2[:rank])

        return np.reshape(out_arr_re, out_dims2[:rank])

    # def handle_cvector_dataset(self, in_handle, out_handle, Rout, cart_map, kvector):
    #     d_in = [[0, 0], [0, 0], [0, 0]]
    #     in_dims = [1, 1, 1]
    #     out_dims = [1, 1, 1]
    #     rank = 3

    #     components = ['x', 'y', 'z']

    #     def try_individual_datasets(in_handle, out_handle, Rout, kvector):
    #         for dim in range(3):
    #             namr = "{}.r".format(components[dim])
    #             nami = "{}.i".format(components[dim])
    #             self.handle_dataset(in_handle, out_handle, namr, nami, Rout, kvector)

    #             namr = re.sub(r'\.r', '', namr)
    #             self.handle_dataset(in_handle, out_handle, namr, None, Rout, kvector)

    #     for dim in range(3):
    #         for ri in range(2):
    #             dims = [1, 1, 1]
    #             rnk = 3

    #             nam = components[dim]
    #             nam += '.i' if ri else '.r'
    #             d_in[dim][ri] = in_handle.get(nam, None)

    #             if d_in[dim][ri] is None:
    #                 try_individual_datasets(in_handle, out_handle, Rout, kvector)
    #                 return
    #             else:
    #                 d_in[dim][ri] = d_in[dim][ri].value

    #             if not dim and not ri:
    #                 rank = rnk
    #                 for i in range(3):
    #                     in_dims[i] = dims[i]
    #             else:
    #                 dims_not_equal = (in_dims[0] != dims[0] or
    #                                   in_dims[1] != dims[1] or
    #                                   in_dims[2] != dims[2])
    #                 if rank != rnk or dims_not_equal:
    #                     try_individual_datasets(in_handle, out_handle, Rout, kvector)
    #                     return

    #     if self.verbose:
    #         print("Found complex vector dataset...")

    #     if self.verbose:
    #         fmt = "Input data is rank {}, size {}x{}x{}."
    #         print(fmt.format(rank, in_dims[0], in_dims[1], in_dims[2]))

    #     # rotate vector field according to cart_map
    #     if self.verbose:
    #         fmt1 = "Rotating vectors by matrix [ {:.10g}, {:.10g}, {:.10g}"
    #         fmt2 = "                             {:.10g}, {:.10g}, {:.10g}"
    #         fmt3 = "                             {:.10g}, {:.10g}, {:.10g} ]"
    #         print(fmt1.format(cart_map.c1.x, cart_map.c2.x, cart_map.c3.x))
    #         print(fmt2.format(cart_map.c1.y, cart_map.c2.y, cart_map.c3.y))
    #         print(fmt3.format(cart_map.c1.z, cart_map.c2.z, cart_map.c3.z))

    #     N = in_dims[0] * in_dims[1] * in_dims[2]
    #     for ri in range(2):
    #         for i in range(N):
    #             v = mp.Vector3(d_in[0][ri][i], d_in[1][ri][i], d_in[2][ri][i])
    #             v = cart_map * v
    #             d_in[0][ri][i] = v.x
    #             d_in[1][ri][i] = v.y
    #             d_in[2][ri][i] = v.z

    #     if self.resolution > 0:
    #         out_dims[0] = Rout.c0.norm() * self.resolution + 0.5
    #         out_dims[1] = Rout.c1.norm() * self.resolution + 0.5
    #         out_dims[2] = Rout.c2.norm() * self.resolution + 0.5
    #     else:
    #         for i in range(3):
    #             out_dims[i] = in_dims[i] * self.multiply_size[i]

    #     for i in range(rank, 3):
    #         out_dims[i] = 1

    #     N = 1
    #     for i in range(3):
    #         out_dims[i] = max(out_dims[i], 1)
    #         N *= out_dims[i]

    #     out_dims2 = [0] * 3

    #     if self.transpose:
    #         out_dims2[0] = out_dims[1]
    #         out_dims2[1] = out_dims[0]
    #         out_dims2[2] = out_dims[2]
    #     else:
    #         out_dims2[0] = out_dims[0]
    #         out_dims2[1] = out_dims[1]
    #         out_dims2[2] = out_dims[2]

    #     if self.verbose:
    #         fmt = "Output data {}x{}x{}."
    #         print(fmt.format(out_dims2[0], out_dims2[1], out_dims2[2]))

    #     for dim in range(3):
    #         out_arr_re = np.zeros(int(N))
    #         d_out_im = np.zeros(int(N))

    #         self.map_data(d_in[dim][0], d_in[dim][1], out_arr_re, d_out_im, out_dims2, kvector)

    #         # multiply * scaleby
    #         complex_out = np.vectorize(complex)(out_arr_re, d_out_im)
    #         complex_out *= self.scaleby
    #         out_arr_re = complex_out[0:2]
    #         d_out_im = complex_out[1:2]

    #         nam = "{}.r-new".format(components[dim])
    #         if out_handle != in_handle:
    #             nam = re.sub('-new', '', nam)

    #         if self.verbose:
    #             print("Writing dataset to {}...".format(nam))

    #         out_handle[nam] = np.reshape(out_arr_re, out_dims)

    #         nam = re.sub('r', 'i', nam)

    #         if self.verbose:
    #             print("Writing dataset to {}...".format(nam))

    #         out_handle[nam] = np.reshape(d_out_im, out_dims)

    #         if self.verbose:
    #             print("Successfully wrote out data.")

    #     return

    def init_output_lattice(self, kvector):

        cart_map = mp.Matrix(
            mp.Vector3(1, 0, 0),
            mp.Vector3(0, 1, 0),
            mp.Vector3(0, 0, 1)
        )

        Rin = mp.Matrix(
            mp.Vector3(*self.lattice[0]),
            mp.Vector3(*self.lattice[1]),
            mp.Vector3(*self.lattice[2])
        )

        if self.verbose:
            print("Read lattice vectors")
            if kvector:
                print("Read Bloch wavevector ({:.6g}, {:.6g}, {:.6g})".format(kvector.x, kvector.y, kvector.z))
            fmt = "Input lattice = ({:.6g}, {:.6g}, {:.6g}), ({:.6g}, {:.6g}, {:.6g}), ({:.6g}, {:.6g}, {:.6g})"
            print(fmt.format(Rin.c1.x, Rin.c1.y, Rin.c1.z,
                             Rin.c2.x, Rin.c2.y, Rin.c2.z,
                             Rin.c3.x, Rin.c3.y, Rin.c3.z))

        Rout = mp.Matrix(Rin.c1, Rin.c2, Rin.c3)

        if self.rectify:
            # Orthogonalize the output lattice vectors.  If have_ve is true,
            # then the first new lattice vector should be in the direction
            # of the ve unit vector; otherwise, the first new lattice vector
            # is the first original lattice vector. Note that we do this in
            # such a way as to preserve the volume of the unit cell, and so
            # that our first vector (in the direction of ve) smoothly
            # interpolates between the original lattice vectors.

            if self.have_ve:
                ve = self.ve.unit()
            else:
                ve = Rout.c1.unit()

            # First, compute c1 in the direction of ve by smoothly
            # interpolating the old c1/c2/c3 (formula is slightly tricky)
            V = Rout.c1.cross(Rout.c2).dot(Rout.c3)
            Rout.c2 = Rout.c2 - Rout.c1
            Rout.c3 = Rout.c3 - Rout.c1
            Rout.c1 = ve.scale(V / Rout.c2.cross(Rout.c3).dot(ve))

            # Now, orthogonalize c2 and c3
            Rout.c2 = Rout.c2 - ve.scale(ve.dot(Rout.c2))
            Rout.c3 = Rout.c3 - ve.scale(ve.dot(Rout.c3))
            Rout.c3 = Rout.c3 - Rout.c2.scale(Rout.c2.dot(Rout.c3) /
                                              Rout.c2.dot(Rout.c2))

            cart_map.c1 = Rout.c1.unit()
            cart_map.c2 = Rout.c2.unit()
            cart_map.c3 = Rout.c3.unit()
            cart_map = cart_map.inverse()

        if self.transpose:
            cart_map = cart_map.transpose()
            v = cart_map.c1
            cart_map.c1 = cart_map.c2
            cart_map.c2 = v
            cart_map = cart_map.transpose()

        Rout.c1 = Rout.c1.scale(self.multiply_size[0])
        Rout.c2 = Rout.c2.scale(self.multiply_size[1])
        Rout.c3 = Rout.c3.scale(self.multiply_size[2])

        if self.verbose:
            fmt = "Output lattice = ({:.6g}, {:.6g}, {:.6g}), ({:.6g}, {:.6g}, {:.6g}), ({:.6g}, {:.6g}, {:.6g})"
            print(fmt.format(Rout.c1.x, Rout.c1.y, Rout.c1.z,
                             Rout.c2.x, Rout.c2.y, Rout.c2.z,
                             Rout.c3.x, Rout.c3.y, Rout.c3.z))

        self.coord_map = Rin.inverse() * Rout
        self.Rout = Rout
        self.cart_map = cart_map

    def convert(self, arr, kpoint=None):
        self.init_output_lattice(kpoint)

        return self.handle_dataset(arr, kpoint)

        # if cvector:
        #     self.handle_cvector_dataset(arr)
