from __future__ import division

import math
import os
import sys

import h5py

import meep as mp

examples_dir = os.path.realpath(os.path.dirname(__file__))
sys.path.insert(0, examples_dir)

# from mpb_tri_rods import main as run_tri_rods
# run_tri_rods()

TWOPI = 6.2831853071795864769252867665590057683943388

have_ve = 0
ve = mp.Vector3(1, 0, 0)
rectify = True
multiply_size = [3, 3, 3]
resolution = 32
transpose = False


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


def map_data(d_in_re, d_in_im, d_out_re, d_out_im, coord_map, kvector,
             pick_nearest, transpose):

    n_in = d_in_re.shape
    n_out = d_out_re.shape

    flat_d_out_re = d_out_re.ravel()
    flat_d_out_im = d_out_im.ravel()

    min_out_re = 1e20
    max_out_re = -1e20
    min_out_im = 1e20
    max_out_im = -1e20

    coord_map.c1 = coord_map.c1.scale(1 / n_out[0])
    coord_map.c2 = coord_map.c2.scale(1 / n_out[1])
    coord_map.c3 = coord_map.c3.scale(1 / n_out[2])

    if kvector:
        s = [x * TWOPI for x in kvector]
    else:
        s = [0, 0, 0]

    # Compute shift so that the origin of the output cell is mapped to the origin
    # of the original primitive cell
    shiftx = 0.5 - (coord_map.c1.x * 0.5 * n_out[0] +
                    coord_map.c2.x * 0.5 * n_out[1] +
                    coord_map.c3.x * 0.5 * n_out[2])
    shifty = 0.5 - (coord_map.c1.y * 0.5 * n_out[0] +
                    coord_map.c2.y * 0.5 * n_out[1] +
                    coord_map.c3.y * 0.5 * n_out[2])
    shiftz = 0.5 - (coord_map.c1.z * 0.5 * n_out[0] +
                    coord_map.c2.z * 0.5 * n_out[1] +
                    coord_map.c3.z * 0.5 * n_out[2])

    for i in range(n_out[0]):
        for j in range(n_out[1]):
            for k in range(n_out[2]):
                if transpose:
                    ijk = (j * n_out[0] + i) * n_out[2] + k
                else:
                    ijk = (i * n_out[1] + j) * n_out[2] + k

                # find the point corresponding to d_out[i,j,k] in the input array,
                # and also find the next-nearest points.
                x = coord_map.c1.x * i + coord_map.c2.x * j + coord_map.c3.x * k + shiftx
                y = coord_map.c1.y * i + coord_map.c2.y * j + coord_map.c3.y * k + shifty
                z = coord_map.c1.z * i + coord_map.c2.z * j + coord_map.c3.z * k + shiftz

                x, xi = modf_positive(x)
                y, yi = modf_positive(y)
                z, zi = modf_positive(z)

                i1 = x * n_in[0]
                j1 = y * n_in[1]
                k1 = z * n_in[2]
                dx = x * n_in[0] - i1
                dy = y * n_in[1] - j1
                dz = z * n_in[2] - k1

                i2, dx, xi2 = adj_point(i1, n_in[0], dx, xi)
                j2, dy, yi2 = adj_point(j1, n_in[1], dy, yi)
                k2, dz, zi2 = adj_point(k1, n_in[2], dz, zi)

                # dx, mdx, etcetera, are the weights for the various points in
                # the input data, which we use for linearly interpolating to get
                # the output point.
                if pick_nearest:
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
                    return (i * n_in[1] + j) * n_in[2] + k

                if d_out_im:
                    flat_d_out_re[ijk] = 0
                    flat_d_out_im[ijk] = 0

                    flat_d_out_re[ijk], flat_d_out_im[ijk] = add_cmplx_times_phase(
                        flat_d_out_re[ijk], flat_d_out_im[ijk],
                        d_in_re[in_index(i1, j1, k1)],
                        d_in_im[in_index(i1, j1, k1)],
                        xi, yi, zi, s, mdx * mdy * mdz
                    )

                    flat_d_out_re[ijk], flat_d_out_im[ijk] = add_cmplx_times_phase(
                        flat_d_out_re[ijk], flat_d_out_im[ijk],
                        d_in_re[in_index(i1, j1, k2)],
                        d_in_im[in_index(i1, j1, k2)],
                        xi, yi, zi2, s, mdx * mdy * dz
                    )
                    flat_d_out_re[ijk], flat_d_out_im[ijk] = add_cmplx_times_phase(
                        flat_d_out_re[ijk], flat_d_out_im[ijk],
                        d_in_re[in_index(i1, j2, k1)],
                        d_in_im[in_index(i1, j2, k1)],
                        xi, yi2, zi, s, mdx * dy * mdz
                    )

                    flat_d_out_re[ijk], flat_d_out_im[ijk] = add_cmplx_times_phase(
                        flat_d_out_re[ijk], flat_d_out_im[ijk],
                        d_in_re[in_index(i1, j2, k2)],
                        d_in_im[in_index(i1, j2, k2)],
                        xi, yi2, zi2, s, mdx * dy * dz
                    )

                    flat_d_out_re[ijk], flat_d_out_im[ijk] = add_cmplx_times_phase(
                        flat_d_out_re[ijk], flat_d_out_im[ijk],
                        d_in_re[in_index(i2, j1, k1)],
                        d_in_im[in_index(i2, j1, k1)],
                        xi2, yi, zi, s, dx * mdy * mdz
                    )

                    flat_d_out_re[ijk], flat_d_out_im[ijk] = add_cmplx_times_phase(
                        flat_d_out_re[ijk], flat_d_out_im[ijk],
                        d_in_re[in_index(i2, j1, k2)],
                        d_in_im[in_index(i2, j1, k2)],
                        xi2, yi, zi2, s, dx * mdy * dz
                    )

                    flat_d_out_re[ijk], flat_d_out_im[ijk] = add_cmplx_times_phase(
                        flat_d_out_re[ijk], flat_d_out_im[ijk],
                        d_in_re[in_index(i2, j2, k1)],
                        d_in_im[in_index(i2, j2, k1)],
                        xi2, yi2, zi, s, dx * dy * mdz
                    )

                    flat_d_out_re[ijk], flat_d_out_im[ijk] = add_cmplx_times_phase(
                        flat_d_out_re[ijk], flat_d_out_im[ijk],
                        d_in_re[in_index(i2, j2, k2)],
                        d_in_im[in_index(i2, j2, k2)],
                        xi2, yi2, zi2, s, dx * dy * dz
                    )

                    min_out_im = math.min(min_out_im, flat_d_out_im[ijk])
                    max_out_im = math.max(max_out_im, flat_d_out_im[ijk])
                else:
                    flat_d_out_re[ijk] = (d_in_re[in_index(i1, j1, k1)] * mdx * mdy * mdz +
                                          d_in_re[in_index(i1, j1, k2)] * mdx * mdy * dz +
                                          d_in_re[in_index(i1, j2, k1)] * mdx * dy * mdz +
                                          d_in_re[in_index(i1, j2, k2)] * mdx * dy * dz +
                                          d_in_re[in_index(i2, j1, k1)] * dx * mdy * mdz +
                                          d_in_re[in_index(i2, j1, k2)] * dx * mdy * dz +
                                          d_in_re[in_index(i2, j2, k1)] * dx * dy * mdz +
                                          d_in_re[in_index(i2, j2, k2)] * dx * dy * dz)

                min_out_re = math.min(min_out_re, flat_d_out_re[ijk])
                max_out_re = math.max(max_out_re, flat_d_out_re[ijk])

    # if verbose:
    #     print("real part range: {} .. {}".format(min_out_re, max_out_re))
    #     if d_out_im:
    #         print("imag part range: {} .. {}".format(min_out_im, max_out_im))


def handle_dataset(in_file, out_file, name_re, name_im, Rout, coord_map,
                   kvector, resolution, scaleby, multiply_size, pick_nearest,
                   transpose):

    out_dims = [1, 1, 1]
    out_dims2 = [1, 1, 1]

    with h5py.File(in_file, 'r') as fin, h5py.File(out_file, 'a') as fout:
        d_in_re = fin[name_re].value
        in_dims = d_in_re.shape

        if name_im:
            d_in_im = fin[name_im].value
            out_dims = d_in_im.shape

        if resolution > 0:
            out_dims[0] = Rout.c1.norm() * resolution + 0.5
            out_dims[1] = Rout.c2.norm() * resolution + 0.5
            out_dims[2] = Rout.c3.norm() * resolution + 0.5
        else:
            for i in range(3):
                out_dims[i] = in_dims[i] * multiply_size[i]

        out_dims[len(out_dims):] = 1
        # N (malloc size)

        if transpose:
            out_dims2[0] = out_dims[1]
            out_dims2[1] = out_dims[0]
            out_dims2[2] = out_dims[2]
        else:
            out_dims2[0] = out_dims[0]
            out_dims2[1] = out_dims[1]
            out_dims2[2] = out_dims[2]

        map_data(d_in_re, d_in_im, d_out_re, d_out_im, coord_map, kvector,
                 pick_nearest, transpose)

        if d_out_im:
            d_out_re *= scaleby
            d_out_re *= scaleby

        if in_file == out_file:
            out_name = name_re + '-new'

        fout[out_name] = d_out_re

        if d_out_im:
            if in_file == out_file:
                out_name_im = name_im + '-new'
            fout[out_name_im] = d_out_im


def handle_file(fname, out_fname, data_name, rectify, have_ve, ve, resolution,
                scaleby, multiply_size, pick_nearest, transpose):

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

    file_flag = 'r' if out_fname else 'r+'

    with h5py.File('mpb_data_analysis-epsilon.h5', file_flag) as f:

        if 'lattice vectors' in f.keys():
            R = f['lattice vectors'].value

            if R.shape == (3, 3):
                Rin.c1 = mp.Vector3(*R[0])
                Rin.c2 = mp.Vector3(*R[1])
                Rin.c3 = mp.Vector3(*R[2])

        if 'Bloch wavevector' in f.keys():
            kvector = f['Bloch wavevector'].value
        else:
            kvector = None

        if 'copies' in f.keys():
            copies = f['copies'].value
            Rin.c1 = Rin.c1.scale(copies[0])
            Rin.c2 = Rin.c2.scale(copies[1])
            Rin.c3 = Rin.c3.scale(copies[2])

            if kvector:
                kvector *= copies
        else:
            copies = None

        Rout = mp.Matrix(Rin.c1, Rin.c2, Rin.c3)

        if rectify:
            # Orthogonalize the output lattice vectors.  If have_ve is true, then the
            # first new lattice vector should be in the direction of the ve unit
            # vector; otherwise, the first new lattice vector is the first original
            # lattice vector. Note that we do this in such a way as to preserve the
            # volume of the unit cell, and so that our first vector (in the direction
            # of ve) smoothly interpolates between the original lattice vectors.

            if have_ve:
                ve = ve.unit()
            else:
                ve = Rout.c1.unit()

            # First, compute c1 in the direction of ve by smoothly interpolating the
            # old c1/c2/c3 (formula is slightly tricky)
            V = Rout.c1.cross(Rout.c2).dot(Rout.c3)
            Rout.c2 = Rout.c2 - Rout.c1
            Rout.c3 = Rout.c3 - Rout.c1
            Rout.c1 = ve.scale(V / Rout.c2.cross(Rout.c3).ve)

            # Now, orthogonalize c2 and c3
            Rout.c2 = Rout.c2 - ve.scale(ve.dot(Rout.c2))
            Rout.c3 = Rout.c3 - ve.scale(ve.dot(Rout.c3))
            Rout.c3 = Rout.c3 - Rout.c2.scale(Rout.c2.dot(Rout.c3) /
                                              Rout.c2.dot(Rout.c2))

            cart_map.c1 = Rout.c1.unit()
            cart_map.c2 = Rout.c2.unit()
            cart_map.c3 = Rout.c3.unit()
            cart_map = cart_map.inverse()

        if transpose:
            cart_map = cart_map.transpose()
            v = cart_map.c1
            cart_map.c1 = cart_map.c2
            cart_map.c2 = v
            cart_map = cart_map.transpose()

        Rout.c1 = Rout.c1.scale(multiply_size[0])
        Rout.c2 = Rout.c2.scale(multiply_size[1])
        Rout.c3 = Rout.c3.scale(multiply_size[2])

        coord_map = Rin.inverse() * Rout

        if out_fname:
            # if verbose:
            #     print("Creating output file {}...".format(out_fname))
            out_file = out_fname
        else:
            # if verbose:
            #     print("Writing output datasets to input file {}...".format(fname))
            out_file = fname

        for dname in [data_name] if data_name else datanames:

            handle_dataset(
                fname,
                out_fname,
                dname,
                None,
                Rout,
                coord_map,
                kvector,
                resolution,
                scaleby,
                multiply_size,
                pick_nearest,
                transpose
            )
