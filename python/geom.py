import math
import meep as mp


def non_negative(x):
    return x >= 0


class Vector3(object):

    def __init__(self, x=0.0, y=0.0, z=0.0):
        self.x = float(x)
        self.y = float(y)
        self.z = float(z)


class MaterialType(object):

    def __init__(self, data=None):
        self.data = data


class GeometricObject(object):

    def __init__(self, material=MaterialType(), center=None):
        self.material = material
        self.center = center

    def __contains__(self, point):
        return mp.point_in_objectp(point, self)


class Sphere(GeometricObject):

    def __init__(self, radius=None, **kwargs):
        self.radius = radius
        super(Sphere, self).__init__(**kwargs)

    @property
    def radius(self):
        return self._radius

    @radius.setter
    def radius(self, val):
        if not val or non_negative(val):
            self._radius = val
        else:
            raise ValueError("Radius cannot be negative. Got {}".format(val))


class Cylinder(GeometricObject):

    def __init__(self, axis=Vector3(0, 0, 1), radius=None, height=None, **kwargs):
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
        if non_negative(val) or not val:
            self._radius = val
        else:
            raise ValueError("Cylinder.radius cannot be negative: {}".format(val))

    @height.setter
    def height(self, val):
        if non_negative(val) or not val:
            self._height = val
        else:
            raise ValueError("Cylinder.height cannot be negative: {}".format(val))


class Wedge(Cylinder):

    def __init__(self, wedge_angle=8 * math.atan(1), wedge_start=Vector3(1, 0, 0), **kwargs):
        self.wedge_angle = wedge_angle
        self.wedge_start = wedge_start
        super(Wedge, self).__init__(**kwargs)


class Cone(Cylinder):

    def __init__(self, radius2=0, **kwargs):
        self.radius2 = radius2
        super(Cone, self).__init__(**kwargs)


# TODO(chogan): Write tests
class Block(GeometricObject):

    def __init__(self, size, e1=Vector3(1, 0, 0), e2=Vector3(0, 1, 0),
                 e3=Vector3(0, 0, 1), **kwargs):
        self.size = size
        self.e1 = e1
        self.e2 = e2
        self.e3 = e3
        super(Block, self).__init__(**kwargs)

    # @property
    # def e1(self):
    #     return self._e1

    # @property
    # def e2(self):
    #     return self._e2

    # @property
    # def e3(self):
    #     return self._e3

    # # TODO(chogan): Cache the result?
    # @property
    # def projection_matrix(self):
    #     return Matrix3x3(self.e1, self.e2, self.e3).inverse()

    # @e1.setter
    # def e1(self, val):
    #     self._e1 = self.process_func1(val)

    # @e2.setter
    # def e2(self, val):
    #     self._e2 = self.process_func2(val)

    # @e3.setter
    # def e3(self, val):
    #     self._e3 = self.process_func3(val)


class Ellipsoid(Block):

    def __init__(self, **kwargs):
        self._inverse_semi_axes = None
        super(Ellipsoid, self).__init__(**kwargs)

    @property
    def inverse_semi_axes(self):
        return self._inverse_semi_axes  # map(lambda x: 2.0 / x, self.size)


class CompoundGeometricObject(GeometricObject):

    def __init__(self, component_objects=[], **kwargs):
        super(CompoundGeometricObject, self).__init__(**kwargs)
        self.component_objects = component_objects
