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
