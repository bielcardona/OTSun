from FreeCAD import Base
from .math import one_orthogonal_vector


def axial_rotation_from_axis_and_angle(axis_origin, axis_dir, angle):
    z_axis = Base.Vector(0, 0, 1)
    local_cs = Base.Placement(axis_origin, Base.Rotation(z_axis, axis_dir))
    return local_cs.multiply(Base.Placement(Base.Vector(), Base.Rotation(angle, 0, 0)).multiply(local_cs.inverse()))


def axial_rotation_from_vector_and_image(origin, vector0, vector1):
    """
    Returns a rotation that transforms the ray with origin and vector equals to vector0 to the ray with same
    origin and vector vector1.

    Parameters
    ----------
    origin : Base.Vector
    vector0 : Base.Vector
    vector1 : Base.Vector

    Returns
    -------

    """
    try:
        normal = vector0.cross(vector1)
        normal.normalize()
    except Base.FreeCADError:
        # we assume now that the two vectors are collinear
        normal = one_orthogonal_vector(vector0)
    angle = vector0.getAngle(vector1)
    return axial_rotation_from_axis_and_angle(origin, normal, angle)


class Joint:
    pass


class AxialJoint(Joint):
    def __init__(self, axis_origin, axis_vector):
        self.axis_origin = axis_origin
        self.axis_vector = axis_vector


class CentralJoint(Joint):
    def __init__(self, center):
        self.center = center


class SolarTracking:
    """
    Class used to model the tracking of the light source by mobile elements in a scene.

    Parameters
    ----------
    target : point/line where light rays should be directed
    source_direction : direction of the sun rays
    scene : scene with mobile elements
    element_joint_map : dict with keys the mobile elements and values the corresponding joint
    """
    def __init__(self, target, source_direction, scene, element_joint_map):
        self.target = target
        self.source_direction = source_direction
        self.scene = scene
        self.joints = element_joint_map
        self.movements = []


class CentralTracking(SolarTracking):
    def __init__(self, *args):
        super(CentralTracking, self).__init__(*args)


class AxialTracking(SolarTracking):
    def __init__(self, *args):
        super(AxialTracking, self).__init__(*args)