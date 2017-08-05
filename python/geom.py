import functools
import math
import operator

import meep as mp


def check_nonnegative(prop, val):
    if val >= 0:
        return val
    else:
        raise ValueError("{} cannot be negative. Got {}".format(prop, val))


class Vector3(object):

    def __init__(self, x=0.0, y=0.0, z=0.0):
        self.x = float(x)
        self.y = float(y)
        self.z = float(z)

    def __eq__(self, other):
        return self.x == other.x and self.y == other.y and self.z == other.z

    def __add__(self, other):
        self.x += other.x
        self.y += other.y
        self.z += other.z

        return self

    def __sub__(self, other):
        self.x -= other.x
        self.y -= other.y
        self.z -= other.z

        return self

    def __getitem__(self, i):
        if i == 0:
            return self.x
        elif i == 1:
            return self.y
        elif i == 2:
            return self.z
        else:
            raise IndexError("No value at index {}".format(i))

    def scale(self, s):
        self.x *= s
        self.y *= s
        self.z *= s

        return self


# TODO(chogan): Write tests
class Matrix3x3(object):

    def __init__(self, c1=Vector3(), c2=Vector3(), c3=Vector3()):
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

    def row(self, i):
        return Vector3(self.c1[i], self.c2[i], self.c3[i])

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
        return Matrix3x3(self.row(0), self.row(1), self.row(2))


class Medium(object):

    def __init__(self, epsilon_diag=Vector3(1, 1, 1),
                 epsilon_offdiag=Vector3(),
                 mu_diag=Vector3(1, 1, 1),
                 mu_offdiag=Vector3(),
                 E_susceptibilities=[],
                 H_susceptibilities=[],
                 E_chi2_diag=Vector3(),
                 E_chi3_diag=Vector3(),
                 H_chi2_diag=Vector3(),
                 H_chi3_diag=Vector3(),
                 D_conductivity_diag=Vector3(),
                 B_conductivity_diag=Vector3()):

        self.epsilon_diag = epsilon_diag
        self.epsilon_offdiag = epsilon_offdiag
        self.mu_diag = mu_diag
        self.mu_offdiag = mu_offdiag
        self.E_susceptibilities = E_susceptibilities
        self.H_susceptibilities = H_susceptibilities
        self.E_chi2_diag = E_chi2_diag
        self.E_chi3_diag = E_chi3_diag
        self.H_chi2_diag = H_chi2_diag
        self.H_chi3_diag = H_chi3_diag
        self.D_conductivity_diag = D_conductivity_diag
        self.B_conductivity_diag = B_conductivity_diag


class Susceptibility:

    def __init__(self, sigma_diag=Vector3(), sigma_offdiag=Vector3()):
        self.sigma_diag = sigma_diag
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

    def __init__(self, material=Medium(), center=Vector3()):
        self.material = material
        self.center = center

    def __contains__(self, point):
        return mp.point_in_objectp(point, self)


class Sphere(GeometricObject):

    def __init__(self, radius, **kwargs):
        self.radius = radius
        super(Sphere, self).__init__(**kwargs)

    @property
    def radius(self):
        return self._radius

    @radius.setter
    def radius(self, val):
        self._radius = check_nonnegative("Sphere.radius", val)


class Cylinder(GeometricObject):

    def __init__(self, radius, axis=Vector3(0, 0, 1), height=float('inf'), **kwargs):
        self.axis = axis
        self.radius = radius
        self.height = height
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


# TODO(chogan): Write tests
class Lattice(object):

    def __init__(self, basis1=Vector3(1, 0, 0),
                 basis2=Vector3(0, 1, 0),
                 basis3=Vector3(0, 0, 1),
                 size=Vector3(1, 1, 1),
                 basis_size=Vector3(1, 1, 1)):

        self.basis1 = basis1
        self.basis2 = basis2
        self.basis3 = basis3
        self.size = size
        self.basis_size = basis_size

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
        B = Matrix3x3(self.b1, self.b2, self.b3)
        if B.determinant() == 0:
            raise ValueError("Lattice basis vectors must be linearly independent!")
        return B

    @property
    def metric(self, B):
        return B.transpose() * self.basis
