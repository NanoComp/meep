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


# TODO(chogan): Add tests
class Cylinder(GeometricObject):

    def __init__(self, axis=Vector3(0, 0, 1), radius=None, height=None, **kwargs):
        # self.process_func = unit_vector3
        self.axis = axis
        self.radius = radius
        self.height = height
        super(Cylinder, self).__init__(**kwargs)

    # @property
    # def axis(self):
    #     return self._axis

    @property
    def radius(self):
        return self._radius

    @property
    def height(self):
        return self._height

    # @axis.setter
    # def axis(self, val):
    #     self._axis = self.process_func(val)

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
