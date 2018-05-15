from __future__ import division

import functools
import math
import numbers
import operator
from collections import namedtuple
from copy import deepcopy
from numbers import Number

import numpy as np
import meep as mp


FreqRange = namedtuple('FreqRange', ['min', 'max'])


def check_nonnegative(prop, val):
    if val >= 0:
        return val
    else:
        raise ValueError("{} cannot be negative. Got {}".format(prop, val))


class Vector3(object):

    def __init__(self, x=0.0, y=0.0, z=0.0):
        self.x = float(x) if type(x) is int else x
        self.y = float(y) if type(y) is int else y
        self.z = float(z) if type(z) is int else z

    def __eq__(self, other):
        return self.x == other.x and self.y == other.y and self.z == other.z

    def __add__(self, other):
        x = self.x + other.x
        y = self.y + other.y
        z = self.z + other.z

        return Vector3(x, y, z)

    def __sub__(self, other):
        x = self.x - other.x
        y = self.y - other.y
        z = self.z - other.z

        return Vector3(x, y, z)

    def __mul__(self, other):
        if type(other) is Vector3:
            return self.dot(other)
        elif isinstance(other, Number):
            return self.scale(other)
        else:
            raise TypeError("No operation known for 'Vector3 * {}'".format(type(other)))

    def __rmul__(self, other):
        if isinstance(other, Number):
            return self.scale(other)
        else:
            raise TypeError("No operation known for '{} * Vector3'".format(type(other)))

    def __getitem__(self, i):
        if i == 0:
            return self.x
        elif i == 1:
            return self.y
        elif i == 2:
            return self.z
        else:
            raise IndexError("No value at index {}".format(i))

    def __repr__(self):
        return "Vector3<{}, {}, {}>".format(self.x, self.y, self.z)

    def __array__(self):
        return np.array([self.x, self.y, self.z])

    def conj(self):
        return Vector3(self.x.conjugate(), self.y.conjugate(), self.z.conjugate())

    def scale(self, s):
        x = self.x * s
        y = self.y * s
        z = self.z * s

        return Vector3(x, y, z)

    def dot(self, v):
        return self.x * v.x + self.y * v.y + self.z * v.z

    def cdot(self, v):
        return self.conj().dot(v)

    def cross(self, v):
        x = self.y * v.z - self.z * v.y
        y = self.z * v.x - self.x * v.z
        z = self.x * v.y - self.y * v.x

        return Vector3(x, y, z)

    def norm(self):
        return math.sqrt(abs(self.dot(self)))

    def unit(self):
        return self.scale(1 / self.norm())

    def close(self, v, tol=1.0e-7):
        return (abs(self.x - v.x) <= tol and
                abs(self.y - v.y) <= tol and
                abs(self.z - v.z) <= tol)

    def rotate(self, axis, theta):
        u = axis.unit()
        vpar = u.scale(u.dot(self))
        vcross = u.cross(self)
        vperp = self - vpar
        return vpar + (vperp.scale(math.cos(theta)) + vcross.scale(math.sin(theta)))

    # rotate vectors in lattice/reciprocal coords (note that the axis
    # is also given in the corresponding basis):

    def rotate_lattice(self, axis, theta, lat):
        a = lattice_to_cartesian(axis, lat)
        v = lattice_to_cartesian(self, lat)
        return cartesian_to_lattice(v.rotate(a, theta), lat)

    def rotate_reciprocal(self, axis, theta, lat):
        a = reciprocal_to_cartesian(axis, lat)
        v = reciprocal_to_cartesian(self, lat)
        return cartesian_to_reciprocal(v.rotate(a, theta), lat)


class Medium(object):

    def __init__(self, epsilon_diag=Vector3(1, 1, 1),
                 epsilon_offdiag=Vector3(0j, 0j, 0j),
                 mu_diag=Vector3(1, 1, 1),
                 mu_offdiag=Vector3(0j, 0j, 0j),
                 E_susceptibilities=[],
                 H_susceptibilities=[],
                 E_chi2_diag=Vector3(),
                 E_chi3_diag=Vector3(),
                 H_chi2_diag=Vector3(),
                 H_chi3_diag=Vector3(),
                 D_conductivity_diag=Vector3(),
                 B_conductivity_diag=Vector3(),
                 epsilon=None,
                 index=None,
                 mu=None,
                 chi2=None,
                 chi3=None,
                 D_conductivity=None,
                 B_conductivity=None,
                 E_chi2=None,
                 E_chi3=None,
                 H_chi2=None,
                 H_chi3=None,
                 valid_freq_range=None):

        if epsilon:
            epsilon_diag = Vector3(epsilon, epsilon, epsilon)
        elif index:
            i2 = index * index
            epsilon_diag = Vector3(i2, i2, i2)

        if mu:
            mu_diag = Vector3(mu, mu, mu)

        if D_conductivity:
            D_conductivity_diag = Vector3(D_conductivity, D_conductivity, D_conductivity)
        if B_conductivity:
            B_conductivity_diag = Vector3(B_conductivity, B_conductivity, B_conductivity)

        if E_chi2:
            E_chi2_diag = Vector3(E_chi2, E_chi2, E_chi2)
        if E_chi3:
            E_chi3_diag = Vector3(E_chi3, E_chi3, E_chi3)
        if H_chi2:
            H_chi2_diag = Vector3(H_chi2, H_chi2, H_chi2)
        if H_chi3:
            H_chi3_diag = Vector3(H_chi3, H_chi3, H_chi3)

        self.epsilon_diag = epsilon_diag
        self.epsilon_offdiag = epsilon_offdiag
        self.mu_diag = mu_diag
        self.mu_offdiag = mu_offdiag
        self.E_susceptibilities = E_susceptibilities
        self.H_susceptibilities = H_susceptibilities
        self.E_chi2_diag = Vector3(chi2, chi2, chi2) if chi2 else E_chi2_diag
        self.E_chi3_diag = Vector3(chi3, chi3, chi3) if chi3 else E_chi3_diag
        self.H_chi2_diag = H_chi2_diag
        self.H_chi3_diag = H_chi3_diag
        self.D_conductivity_diag = D_conductivity_diag
        self.B_conductivity_diag = B_conductivity_diag
        self.valid_freq_range = valid_freq_range


class Susceptibility(object):

    def __init__(self, sigma_diag=Vector3(), sigma_offdiag=Vector3(), sigma=None):
        self.sigma_diag = Vector3(sigma, sigma, sigma) if sigma else sigma_diag
        self.sigma_offdiag = sigma_offdiag


class LorentzianSusceptibility(Susceptibility):

    def __init__(self, frequency=0.0, gamma=0.0, **kwargs):
        super(LorentzianSusceptibility, self).__init__(**kwargs)
        self.frequency = frequency
        self.gamma = gamma


class DrudeSusceptibility(Susceptibility):

    def __init__(self, frequency=0.0, gamma=0.0, **kwargs):
        super(DrudeSusceptibility, self).__init__(**kwargs)
        self.frequency = frequency
        self.gamma = gamma


class NoisyLorentzianSusceptibility(LorentzianSusceptibility):

    def __init__(self, noise_amp=0.0, **kwargs):
        super(NoisyLorentzianSusceptibility, self).__init__(**kwargs)
        self.noise_amp = noise_amp


class NoisyDrudeSusceptibility(DrudeSusceptibility):

    def __init__(self, noise_amp=0.0, **kwargs):
        super(NoisyDrudeSusceptibility, self).__init__(**kwargs)
        self.noise_amp = noise_amp


class GeometricObject(object):

    def __init__(self, material=Medium(), center=Vector3(), epsilon_func=None):
        if type(material) is not Medium and callable(material):
            material.eps = False
        elif epsilon_func:
            epsilon_func.eps = True
            material = epsilon_func

        self.material = material
        self.center = center

    def __contains__(self, point):
        return mp.is_point_in_object(point, self)

    def __add__(self, vec):
        self.center += vec

    def shift(self, vec):
        c = deepcopy(self)
        c.center += vec
        return c

    def info(self, indent_by=0):
        mp.display_geometric_object_info(indent_by, self)


class Sphere(GeometricObject):

    def __init__(self, radius, **kwargs):
        self.radius = float(radius)
        super(Sphere, self).__init__(**kwargs)

    @property
    def radius(self):
        return self._radius

    @radius.setter
    def radius(self, val):
        self._radius = check_nonnegative("Sphere.radius", val)


class Cylinder(GeometricObject):

    def __init__(self, radius, axis=Vector3(0, 0, 1), height=1e20, **kwargs):
        self.axis = axis
        self.radius = float(radius)
        self.height = float(height)
        super(Cylinder, self).__init__(**kwargs)

    @property
    def radius(self):
        return self._radius

    @property
    def height(self):
        return self._height

    @radius.setter
    def radius(self, val):
        self._radius = check_nonnegative("Cylinder.radius", val)

    @height.setter
    def height(self, val):
        self._height = check_nonnegative("Cylinder.height", val)


class Wedge(Cylinder):

    def __init__(self, radius, wedge_angle=2 * math.pi, wedge_start=Vector3(1, 0, 0), **kwargs):
        self.wedge_angle = wedge_angle
        self.wedge_start = wedge_start
        super(Wedge, self).__init__(radius, **kwargs)


class Cone(Cylinder):

    def __init__(self, radius, radius2=0, **kwargs):
        self.radius2 = radius2
        super(Cone, self).__init__(radius, **kwargs)


class Block(GeometricObject):

    def __init__(self, size, e1=Vector3(1, 0, 0), e2=Vector3(0, 1, 0), e3=Vector3(0, 0, 1), **kwargs):
        self.size = size
        self.e1 = e1
        self.e2 = e2
        self.e3 = e3
        super(Block, self).__init__(**kwargs)


class Ellipsoid(Block):

    def __init__(self, **kwargs):
        super(Ellipsoid, self).__init__(**kwargs)


class Matrix(object):

    def __init__(self, c1=Vector3(), c2=Vector3(), c3=Vector3()):
        self.c1 = c1
        self.c2 = c2
        self.c3 = c3

    def __getitem__(self, i):
        return self.row(i)

    def __mul__(self, m):
        if type(m) is Matrix:
            return self.mm_mult(m)
        elif type(m) is Vector3:
            return self.mv_mult(m)
        elif isinstance(m, Number):
            return self.scale(m)
        else:
            raise TypeError("No operation known for 'Matrix * {}'".format(type(m)))

    def __repr__(self):
        return "<{}\n {}\n {}>".format(self.c1, self.c2, self.c3)

    def row(self, i):
        return Vector3(self.c1[i], self.c2[i], self.c3[i])

    def mm_mult(self, m):
        c1 = Vector3(self.row(0).dot(m.c1),
                     self.row(1).dot(m.c1),
                     self.row(2).dot(m.c1))
        c2 = Vector3(self.row(0).dot(m.c2),
                     self.row(1).dot(m.c2),
                     self.row(2).dot(m.c2))
        c3 = Vector3(self.row(0).dot(m.c3),
                     self.row(1).dot(m.c3),
                     self.row(2).dot(m.c3))

        return Matrix(c1, c2, c3)

    def mv_mult(self, v):
        return Vector3(*[self.row(i).dot(v) for i in range(3)])

    def scale(self, s):
        return Matrix(self.c1.scale(s), self.c2.scale(s), self.c3.scale(s))

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

    def inverse(self):
        v1x = self[1][1] * self[2][2] - self[1][2] * self[2][1]
        v1y = self[1][2] * self[2][0] - self[1][0] * self[2][2]
        v1z = self[1][0] * self[2][1] - self[1][1] * self[2][0]
        v1 = mp.Vector3(v1x, v1y, v1z)

        v2x = self[2][1] * self[0][2] - self[0][1] * self[2][2]
        v2y = self[0][0] * self[2][2] - self[0][2] * self[2][0]
        v2z = self[0][1] * self[2][0] - self[0][0] * self[2][1]
        v2 = mp.Vector3(v2x, v2y, v2z)

        v3x = self[0][1] * self[1][2] - self[1][1] * self[0][2]
        v3y = self[1][0] * self[0][2] - self[0][0] * self[1][2]
        v3z = self[1][1] * self[0][0] - self[1][0] * self[0][1]
        v3 = mp.Vector3(v3x, v3y, v3z)

        m = Matrix(v1, v2, v3)

        return m.scale(1 / self.determinant())


class Lattice(object):

    def __init__(self,
                 size=Vector3(1, 1, 1),
                 basis_size=Vector3(1, 1, 1),
                 basis1=Vector3(1, 0, 0),
                 basis2=Vector3(0, 1, 0),
                 basis3=Vector3(0, 0, 1)):

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


def lattice_to_cartesian(x, lat):
    if isinstance(x, Vector3):
        return lat.basis * x

    return (lat.basis * x) * lat.basis.inverse()


def cartesian_to_lattice(x, lat):
    if isinstance(x, Vector3):
        return lat.basis.inverse() * x

    return (lat.basis.inverse() * x) * lat.basis


def reciprocal_to_cartesian(x, lat):
    s = Vector3(*[1 if v == 0 else v for v in lat.size])

    m = Matrix(Vector3(s.x), Vector3(y=s.y), Vector3(z=s.z))
    Rst = (lat.basis * m).transpose()

    if isinstance(x, Vector3):
        return Rst.inverse() * x
    else:
        return (Rst.inverse() * x) * Rst


def cartesian_to_reciprocal(x, lat):
    s = Vector3(*[1 if v == 0 else v for v in lat.size])

    m = Matrix(Vector3(s.x), Vector3(y=s.y), Vector3(z=s.z))
    Rst = (lat.basis * m).transpose()

    if isinstance(x, Vector3):
        return Rst * x
    else:
        return (Rst * x) * Rst.inverse()


def lattice_to_reciprocal(x, lat):
    return cartesian_to_reciprocal(lattice_to_cartesian(x, lat), lat)


def reciprocal_to_lattice(x, lat):
    return cartesian_to_lattice(reciprocal_to_cartesian(x, lat), lat)


def geometric_object_duplicates(shift_vector, min_multiple, max_multiple, go):

    def _dup(min_multiple, lst):
        if min_multiple <= max_multiple:
            shifted = go.shift(shift_vector.scale(min_multiple))
            return _dup(min_multiple + 1, [shifted] + lst)
        else:
            return lst

    return _dup(min_multiple, [])


def geometric_objects_duplicates(shift_vector, min_multiple, max_multiple, go_list):
    dups = []
    for go in go_list:
        dups += geometric_object_duplicates(shift_vector, min_multiple, max_multiple, go)
    return dups


def geometric_objects_lattice_duplicates(lat, go_list, *usize):

    def lat_to_lattice(v):
        return cartesian_to_lattice(lat.basis * v, lat)

    u1 = usize[0] if usize else 1
    u2 = usize[1] if len(usize) >= 2 else 1
    u3 = usize[2] if len(usize) >= 3 else 1
    s = lat.size

    b1 = lat_to_lattice(mp.Vector3(u1))
    b2 = lat_to_lattice(mp.Vector3(0, u2, 0))
    b3 = lat_to_lattice(mp.Vector3(0, 0, u3))

    n1 = math.ceil((s.x if s.x else 1e-20) / u1)
    n2 = math.ceil((s.y if s.y else 1e-20) / u2)
    n3 = math.ceil((s.z if s.z else 1e-20) / u3)

    min3 = -math.floor((n3 - 1) / 2)
    max3 = math.ceil((n3 - 1) / 2)
    d3 = geometric_objects_duplicates(b3, int(min3), int(max3), go_list)

    min2 = -math.floor((n2 - 1) / 2)
    max2 = math.ceil((n2 - 1) / 2)
    d2 = geometric_objects_duplicates(b2, int(min2), int(max2), d3)

    min1 = -math.floor((n1 - 1) / 2)
    max1 = math.ceil((n1 - 1) / 2)

    return geometric_objects_duplicates(b1, int(min1), int(max1), d2)


# Return a 'memoized' version of the function f, which caches its
# arguments and return values so as never to compute the same thing twice.
def memoize(f):
    f_memo_tab = {}

    def _mem(y=None):
        tab_val = f_memo_tab.get(y, None)
        if tab_val:
            return tab_val

        fy = f(y)
        f_memo_tab[y] = fy
        return fy
    return _mem


# Find a root by Newton's method with bounds and bisection,
# given a function f that returns a pair of (value . derivative)
def find_root_deriv(f, tol, x_min, x_max, x_guess=None):
    # Some trickiness: we only need to evaluate the function at x_min and
    # x_max if a Newton step fails, and even then only if we haven't already
    # bracketed the root, so do this via lazy evaluation.
    f_memo = memoize(f)

    def lazy(x):
        if isinstance(x, numbers.Number):
            return x
        return x()

    def pick_bound(which):
        def _pb():
            fmin_tup = f_memo(x_min)
            fmax_tup = f_memo(x_max)
            fmin = fmin_tup[0]
            fmax = fmax_tup[0]

            if which(fmin):
                return x_min
            elif which(fmax):
                return x_max
            else:
                raise ValueError("failed to bracket the root in find_root_deriv")
        return _pb

    def in_bounds(x, f, df, a, b):
        return (f - (df * (x - a))) * (f - (df * (x - b))) < 0

    def newton(x, a, b, dx):
        if abs(dx) < abs(tol * x):
            return x

        fx_tup = f_memo(x)
        f = fx_tup[0]
        df = fx_tup[1]

        if f == 0:
            return x

        a_prime = x if f < 0 else a
        b_prime = x if f > 0 else b

        cond = f_memo(lazy(a_prime))[0] * f_memo(lazy(b_prime))[0] > 0
        if dx != x_max - x_min and dx * (f / df) < 0 and cond:
            raise ValueError("failed to bracket the root in find_root_deriv")

        if isinstance(a, numbers.Number) and isinstance(b, numbers.Number):
            is_in_bounds = in_bounds(x, f, df, a, b)
        else:
            is_in_bounds = in_bounds(x, f, df, x_min, x_max)

        if is_in_bounds:
            return newton(x - (f / df), a_prime, b_prime, f / df)

        av = lazy(a)
        bv = lazy(b)
        dx_prime = 0.5 * (bv - av)
        a_pp = av if a == a_prime else a_prime
        b_pp = bv if b == b_prime else b_prime

        return newton((av + bv) * 0.5, a_pp, b_pp, dx_prime)

    if x_guess is None:
        x_guess = (x_min + x_max) * 0.5

    return newton(x_guess, pick_bound(lambda aa: aa < 0),
                  pick_bound(lambda aa: aa > 0), x_max - x_min)
