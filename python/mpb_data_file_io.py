from __future__ import division

import math
import re
import sys

import h5py
import numpy as np
import meep as mp


def modf_positive(x):
    f, i = math.modf(x)
    if f < 0:
        f += 1.0
        if f >= 1.0:
            f = 0
        else:
            i -= 1.0
    return f, i


def adj_point(i1, nx, dx, xi):
    if dx >= 0:
        i2 = i1 + 1
        if i2 >= nx:
            i2 -= nx
            xi2 = xi + 1
        else:
            xi2 = xi
    else:
        i2 = i1 - 1
        if i2 < 0:
            i2 += nx
            xi2 = xi - 1
        else:
            xi2 = xi
        dx = -dx

    return i2, dx, xi2


def add_cmplx_times_phase(sum_re, sum_im, d_re, d_im, ix, iy, iz, s, scale_by):
    phase = 0.0
    p_re = 1.0
    p_im = 0.0

    new_phase = ix * s[0] + iy * s[1] + iz * s[2]

    if new_phase != phase:
        phase = new_phase
        p_re = math.cos(phase)
        p_im = math.sin(phase)

    sum_re += (d_re * p_re - d_im * p_im) * scale_by
    sum_im += (d_re * p_im + d_im * p_re) * scale_by

    return sum_re, sum_im


class MPBData(object):

    TWOPI = 6.2831853071795864769252867665590057683943388

    def __init__(self,
                 rectify=False,
                 x=0,
                 y=0,
                 z=0,
                 m=0,
                 resolution=0,
                 out_fname='',
                 data_name='',
                 phase_angle=0,
                 pick_nearest=False,
                 transpose=False,
                 ve=None,
                 verbose=False):

        self.rectify = rectify

        if m:
            self.multiply_size = [m, m, m]
        else:
            self.multiply_size = [
                x if x else 1,
                y if y else 1,
                z if z else 1
            ]

        self.resolution = resolution
        self.out_fname = out_fname
        self.data_name = data_name
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

    def map_data(self, d_in_re, d_in_im, d_out_re, d_out_im, n_out, kvector):
        rank = len(d_in_re.shape)
        num_ones = 3 - rank
        n_in = [x for x in d_in_re.shape] + [1] * num_ones

        flat_d_in_re = d_in_re.ravel()
        flat_d_in_im = d_in_im.ravel() if d_in_im else None

        min_out_re = 1e20
        max_out_re = -1e20
        min_out_im = 1e20
        max_out_im = -1e20

        self.coord_map.c1 = self.coord_map.c1.scale(1 / n_out[0])
        self.coord_map.c2 = self.coord_map.c2.scale(1 / n_out[1])
        self.coord_map.c3 = self.coord_map.c3.scale(1 / n_out[2])

        if kvector:
            s = [x * self.TWOPI for x in kvector]
        else:
            s = [0, 0, 0]

        # Compute shift so that the origin of the output cell is mapped to the origin
        # of the original primitive cell
        shiftx = 0.5 - (self.coord_map.c1.x * 0.5 * n_out[0] +
                        self.coord_map.c2.x * 0.5 * n_out[1] +
                        self.coord_map.c3.x * 0.5 * n_out[2])
        shifty = 0.5 - (self.coord_map.c1.y * 0.5 * n_out[0] +
                        self.coord_map.c2.y * 0.5 * n_out[1] +
                        self.coord_map.c3.y * 0.5 * n_out[2])
        shiftz = 0.5 - (self.coord_map.c1.z * 0.5 * n_out[0] +
                        self.coord_map.c2.z * 0.5 * n_out[1] +
                        self.coord_map.c3.z * 0.5 * n_out[2])

        for i in range(n_out[0]):
            for j in range(n_out[1]):
                for k in range(n_out[2]):
                    if self.transpose:
                        ijk = (j * n_out[0] + i) * n_out[2] + k
                    else:
                        ijk = (i * n_out[1] + j) * n_out[2] + k

                    # find the point corresponding to d_out[i,j,k] in the input array,
                    # and also find the next-nearest points.
                    x = (self.coord_map.c1.x * i + self.coord_map.c2.x * j +
                         self.coord_map.c3.x * k + shiftx)
                    y = (self.coord_map.c1.y * i + self.coord_map.c2.y * j +
                         self.coord_map.c3.y * k + shifty)
                    z = (self.coord_map.c1.z * i + self.coord_map.c2.z * j +
                         self.coord_map.c3.z * k + shiftz)

                    x, xi = modf_positive(x)
                    y, yi = modf_positive(y)
                    z, zi = modf_positive(z)

                    i1 = int(x * n_in[0])
                    j1 = int(y * n_in[1])
                    k1 = int(z * n_in[2])
                    dx = x * n_in[0] - i1
                    dy = y * n_in[1] - j1
                    dz = z * n_in[2] - k1

                    i2, dx, xi2 = adj_point(i1, n_in[0], dx, xi)
                    j2, dy, yi2 = adj_point(j1, n_in[1], dy, yi)
                    k2, dz, zi2 = adj_point(k1, n_in[2], dz, zi)

                    # dx, mdx, etcetera, are the weights for the various points in
                    # the input data, which we use for linearly interpolating to get
                    # the output point.
                    if self.pick_nearest:
                        # don't interpolate
                        dx = 0 if dx <= 0.5 else 1
                        dy = 0 if dy <= 0.5 else 1
                        dz = 0 if dz <= 0.5 else 1

                    mdx = 1.0 - dx
                    mdy = 1.0 - dy
                    mdz = 1.0 - dz

                    # Now, linearly interpolate the input to get the output. If the
                    # input/output are complex, we also need to multiply by the
                    # appropriate phase factor, depending upon which unit cell we
                    # are in.
                    def in_index(i, j, k):
                        return int((i * n_in[1] + j) * n_in[2] + k)

                    if isinstance(d_out_im, np.ndarray):
                        d_out_re[ijk] = 0
                        d_out_im[ijk] = 0

                        d_out_re[ijk], d_out_im[ijk] = add_cmplx_times_phase(
                            d_out_re[ijk], d_out_im[ijk],
                            flat_d_in_re[in_index(i1, j1, k1)],
                            flat_d_in_im[in_index(i1, j1, k1)],
                            xi, yi, zi, s, mdx * mdy * mdz
                        )

                        d_out_re[ijk], d_out_im[ijk] = add_cmplx_times_phase(
                            d_out_re[ijk], d_out_im[ijk],
                            flat_d_in_re[in_index(i1, j1, k2)],
                            flat_d_in_im[in_index(i1, j1, k2)],
                            xi, yi, zi2, s, mdx * mdy * dz
                        )
                        d_out_re[ijk], d_out_im[ijk] = add_cmplx_times_phase(
                            d_out_re[ijk], d_out_im[ijk],
                            flat_d_in_re[in_index(i1, j2, k1)],
                            flat_d_in_im[in_index(i1, j2, k1)],
                            xi, yi2, zi, s, mdx * dy * mdz
                        )

                        d_out_re[ijk], d_out_im[ijk] = add_cmplx_times_phase(
                            d_out_re[ijk], d_out_im[ijk],
                            flat_d_in_re[in_index(i1, j2, k2)],
                            flat_d_in_im[in_index(i1, j2, k2)],
                            xi, yi2, zi2, s, mdx * dy * dz
                        )

                        d_out_re[ijk], d_out_im[ijk] = add_cmplx_times_phase(
                            d_out_re[ijk], d_out_im[ijk],
                            flat_d_in_re[in_index(i2, j1, k1)],
                            flat_d_in_im[in_index(i2, j1, k1)],
                            xi2, yi, zi, s, dx * mdy * mdz
                        )

                        d_out_re[ijk], d_out_im[ijk] = add_cmplx_times_phase(
                            d_out_re[ijk], d_out_im[ijk],
                            flat_d_in_re[in_index(i2, j1, k2)],
                            flat_d_in_im[in_index(i2, j1, k2)],
                            xi2, yi, zi2, s, dx * mdy * dz
                        )

                        d_out_re[ijk], d_out_im[ijk] = add_cmplx_times_phase(
                            d_out_re[ijk], d_out_im[ijk],
                            flat_d_in_re[in_index(i2, j2, k1)],
                            flat_d_in_im[in_index(i2, j2, k1)],
                            xi2, yi2, zi, s, dx * dy * mdz
                        )

                        d_out_re[ijk], d_out_im[ijk] = add_cmplx_times_phase(
                            d_out_re[ijk], d_out_im[ijk],
                            flat_d_in_re[in_index(i2, j2, k2)],
                            flat_d_in_im[in_index(i2, j2, k2)],
                            xi2, yi2, zi2, s, dx * dy * dz
                        )

                        min_out_im = min(min_out_im, d_out_im[ijk])
                        max_out_im = max(max_out_im, d_out_im[ijk])
                    else:
                        d_out_re[ijk] = (flat_d_in_re[in_index(i1, j1, k1)] *
                                         mdx * mdy * mdz +
                                         flat_d_in_re[in_index(i1, j1, k2)] *
                                         mdx * mdy * dz +
                                         flat_d_in_re[in_index(i1, j2, k1)] *
                                         mdx * dy * mdz +
                                         flat_d_in_re[in_index(i1, j2, k2)] *
                                         mdx * dy * dz +
                                         flat_d_in_re[in_index(i2, j1, k1)] *
                                         dx * mdy * mdz +
                                         flat_d_in_re[in_index(i2, j1, k2)] *
                                         dx * mdy * dz +
                                         flat_d_in_re[in_index(i2, j2, k1)] *
                                         dx * dy * mdz +
                                         flat_d_in_re[in_index(i2, j2, k2)] *
                                         dx * dy * dz)

                    min_out_re = min(min_out_re, d_out_re[ijk])
                    max_out_re = max(max_out_re, d_out_re[ijk])

        if self.verbose:
            print("real part range: {:g} .. {:g}".format(min_out_re, max_out_re))
            if d_out_im:
                print("imag part range: {:g} .. {:g}".format(min_out_im, max_out_im))

    def handle_dataset(self, in_handle, out_handle, name_re, name_im, Rout, kvector):

        d_in_re = in_handle.get(name_re, None)
        if d_in_re is None:
            return
        else:
            d_in_re = d_in_re.value

        d_in_im = None
        out_dims = [1, 1, 1]

        rank = len(d_in_re.shape)
        num_ones = 3 - rank
        in_dims = [x for x in d_in_re.shape] + [1] * num_ones

        if self.verbose:
            print("Found dataset {}...".format(name_re))

        if name_im:
            d_in_im = in_handle[name_im].value

            if not d_in_im:
                print("Found {} dataset but not {}".format(name_re, name_im))
                return

            out_dims = [x for x in d_in_im.shape] + [1] * num_ones

            if out_dims != in_dims:
                sys.exit("re/im datasets must have same size.")

            if self.verbose:
                print("   and imaginary part dataset {}...".format(name_im))

        if self.verbose:
            fmt = "Input data is rank {}, size {}x{}x{}."
            print(fmt.format(rank, in_dims[0], in_dims[1], in_dims[2]))

        if self.resolution > 0:
            out_dims[0] = math.floor(Rout.c1.norm() * self.resolution + 0.5)
            out_dims[1] = math.floor(Rout.c2.norm() * self.resolution + 0.5)
            out_dims[2] = math.floor(Rout.c3.norm() * self.resolution + 0.5)
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
            out_dims2[0] = out_dims[1]
            out_dims2[1] = out_dims[0]
            out_dims2[2] = out_dims[2]
        else:
            out_dims2[0] = out_dims[0]
            out_dims2[1] = out_dims[1]
            out_dims2[2] = out_dims[2]

        if self.verbose:
            print("Output data {}x{}x{}".format(out_dims2[0], out_dims2[1], out_dims2[2]))

        d_out_re = np.zeros(int(N))
        if name_im:
            d_out_im = np.zeros(int(N))
        else:
            d_out_im = None

        self.map_data(d_in_re, d_in_im, d_out_re, d_out_im, out_dims2, kvector)

        if d_out_im:
            # multiply * scaleby for complex data
            complex_out = np.vectorize(complex)(d_out_re, d_out_im)
            complex_out *= self.scaleby
            d_out_re = complex_out[0:2]
            d_out_im = complex_out[1:2]

        if in_handle == out_handle:
            out_name = name_re + '-new'
        else:
            out_name = name_re

        if self.verbose:
            print("Writing dataset to {}...".format(out_name))

        d_out_re = np.reshape(d_out_re, out_dims2[:rank])

        out_handle[out_name] = d_out_re

        if d_out_im:
            if in_handle == out_handle:
                out_name_im = name_im + '-new'
            else:
                out_name_im = name_im

            if self.verbose:
                print("Writing dataset to {}...".format(out_name_im))

            d_out_im = np.reshape(d_out_im, out_dims2[:rank])
            out_handle[out_name_im] = d_out_im

        if self.verbose:
            print("Successfully wrote out data.")

    def handle_cvector_dataset(self, in_handle, out_handle, Rout, cart_map, kvector):
        d_in = [[0, 0], [0, 0], [0, 0]]
        in_dims = [1, 1, 1]
        out_dims = [1, 1, 1]
        rank = 3

        components = ['x', 'y', 'z']

        def try_individual_datasets(in_handle, out_handle, Rout, kvector):
            for dim in range(3):
                namr = "{}.r".format(components[dim])
                nami = "{}.i".format(components[dim])
                self.handle_dataset(in_handle, out_handle, namr, nami, Rout, kvector)

                namr = re.sub(r'\.r', '', namr)
                self.handle_dataset(in_handle, out_handle, namr, None, Rout, kvector)

        for dim in range(3):
            for ri in range(2):
                dims = [1, 1, 1]
                rnk = 3

                nam = components[dim]
                nam += '.i' if ri else '.r'
                d_in[dim][ri] = in_handle.get(nam, None)

                if d_in[dim][ri] is None:
                    try_individual_datasets(in_handle, out_handle, Rout, kvector)
                    return
                else:
                    d_in[dim][ri] = d_in[dim][ri].value

                if not dim and not ri:
                    rank = rnk
                    for i in range(3):
                        in_dims[i] = dims[i]
                else:
                    dims_not_equal = (in_dims[0] != dims[0] or
                                      in_dims[1] != dims[1] or
                                      in_dims[2] != dims[2])
                    if rank != rnk or dims_not_equal:
                        try_individual_datasets(in_handle, out_handle, Rout, kvector)
                        return

        if self.verbose:
            print("Found complex vector dataset...")

        if self.verbose:
            fmt = "Input data is rank {}, size {}x{}x{}."
            print(fmt.format(rank, in_dims[0], in_dims[1], in_dims[2]))

        # rotate vector field according to cart_map
        if self.verbose:
            fmt1 = "Rotating vectors by matrix [ {:.10g}, {:.10g}, {:.10g}"
            fmt2 = "                             {:.10g}, {:.10g}, {:.10g}"
            fmt3 = "                             {:.10g}, {:.10g}, {:.10g} ]"
            print(fmt1.format(cart_map.c1.x, cart_map.c2.x, cart_map.c3.x))
            print(fmt2.format(cart_map.c1.y, cart_map.c2.y, cart_map.c3.y))
            print(fmt3.format(cart_map.c1.z, cart_map.c2.z, cart_map.c3.z))

        N = in_dims[0] * in_dims[1] * in_dims[2]
        for ri in range(2):
            for i in range(N):
                v = mp.Vector3(d_in[0][ri][i], d_in[1][ri][i], d_in[2][ri][i])
                v = cart_map * v
                d_in[0][ri][i] = v.x
                d_in[1][ri][i] = v.y
                d_in[2][ri][i] = v.z

        if self.resolution > 0:
            out_dims[0] = Rout.c0.norm() * self.resolution + 0.5
            out_dims[1] = Rout.c1.norm() * self.resolution + 0.5
            out_dims[2] = Rout.c2.norm() * self.resolution + 0.5
        else:
            for i in range(3):
                out_dims[i] = in_dims[i] * self.multiply_size[i]

        for i in range(rank, 3):
            out_dims[i] = 1

        N = 1
        for i in range(3):
            out_dims[i] = max(out_dims[i], 1)
            N *= out_dims[i]

        out_dims2 = [0] * 3

        if self.transpose:
            out_dims2[0] = out_dims[1]
            out_dims2[1] = out_dims[0]
            out_dims2[2] = out_dims[2]
        else:
            out_dims2[0] = out_dims[0]
            out_dims2[1] = out_dims[1]
            out_dims2[2] = out_dims[2]

        if self.verbose:
            fmt = "Output data {}x{}x{}."
            print(fmt.format(out_dims2[0], out_dims2[1], out_dims2[2]))

        for dim in range(3):
            d_out_re = np.zeros(int(N))
            d_out_im = np.zeros(int(N))

            self.map_data(d_in[dim][0], d_in[dim][1], d_out_re, d_out_im, out_dims2, kvector)

            # multiply * scaleby
            complex_out = np.vectorize(complex)(d_out_re, d_out_im)
            complex_out *= self.scaleby
            d_out_re = complex_out[0:2]
            d_out_im = complex_out[1:2]

            nam = "{}.r-new".format(components[dim])
            if out_handle != in_handle:
                nam = re.sub('-new', '', nam)

            if self.verbose:
                print("Writing dataset to {}...".format(nam))

            out_handle[nam] = np.reshape(d_out_re, out_dims)

            nam = re.sub('r', 'i', nam)

            if self.verbose:
                print("Writing dataset to {}...".format(nam))

            out_handle[nam] = np.reshape(d_out_im, out_dims)

            if self.verbose:
                print("Successfully wrote out data.")

        return

    def handle_file(self, h5_fname, dname):

        Rin = mp.Matrix(
            mp.Vector3(1, 0, 0),
            mp.Vector3(0, 1, 0),
            mp.Vector3(0, 0, 1)
        )

        cart_map = mp.Matrix(
            mp.Vector3(1, 0, 0),
            mp.Vector3(0, 1, 0),
            mp.Vector3(0, 0, 1)
        )

        datanames = [
            "data",
            "epsilon.xx",
            "epsilon.xy",
            "epsilon.xz",
            "epsilon.yy",
            "epsilon.yz",
            "epsilon.zz",
            "epsilon_inverse.xx",
            "epsilon_inverse.xy",
            "epsilon_inverse.xz",
            "epsilon_inverse.yy",
            "epsilon_inverse.yz",
            "epsilon_inverse.zz"
        ]

        file_flag = 'r' if self.out_fname else 'r+'

        if self.verbose:
            print("Reading file {}...".format(h5_fname))

        in_handle = h5py.File(h5_fname, file_flag)

        if 'lattice vectors' in in_handle.keys():
            R = in_handle['lattice vectors'].value

            if R.shape == (3, 3):
                Rin.c1 = mp.Vector3(*R[0])
                Rin.c2 = mp.Vector3(*R[1])
                Rin.c3 = mp.Vector3(*R[2])

                if self.verbose:
                    print("Read lattice vectors")

        if 'Bloch wavevector' in in_handle.keys():
            kvector = in_handle['Bloch wavevector'].value
            if self.verbose:
                fmt = "Read Bloch wavevector ({}, {}, {})"
                print(fmt.format(kvector[0], kvector[1], kvector[2]))
        else:
            kvector = None

        if 'copies' in in_handle.keys():
            copies = in_handle['copies'].value
            Rin.c1 = Rin.c1.scale(copies[0])
            Rin.c2 = Rin.c2.scale(copies[1])
            Rin.c3 = Rin.c3.scale(copies[2])

            if kvector:
                kvector *= copies
            if self.verbose:
                fmt = "Read lattice copies ({}, {}, {})"
                print(fmt.format(copies[0], copies[1], copies[2]))
        else:
            copies = None

        if self.verbose:
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

        if self.out_fname:
            if self.verbose:
                print("Creating output file {}...".format(self.out_fname))
            out_handle = h5py.File(self.out_fname, 'x')
        else:
            if self.verbose:
                fmt = "Writing output datasets to input file {}..."
                print(fmt.format(h5_fname))
            out_handle = in_handle

        for dn in [dname] if dname else datanames:
            self.handle_dataset(in_handle, out_handle, dn, None, Rout, kvector)
            name_re = dn + '.r'
            name_im = dn + '.i'
            self.handle_dataset(in_handle, out_handle, name_re, name_im, Rout, kvector)

        self.handle_cvector_dataset(in_handle, out_handle, Rout, cart_map, kvector)

        out_handle.close()
        in_handle.close()

    def run(self, *fnames):

        for fn in fnames:
            if ':' in fn:
                h5_fname, dname = fn.split(':')
            else:
                h5_fname = fn
                dname = self.data_name

            self.handle_file(h5_fname, dname)
