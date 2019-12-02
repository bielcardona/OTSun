from FreeCAD import Base
from .math import one_orthogonal_vector, projection_on_vector, projection_on_orthogonal_of_vector
from numpy import pi


def axial_rotation_from_axis_and_angle(axis_origin, axis_dir, angle):
    """
    Returns a rotation with given axis and angle
    """
    z_axis = Base.Vector(0, 0, 1)
    local_cs = Base.Placement(axis_origin, Base.Rotation(z_axis, axis_dir))
    return local_cs.multiply(
        Base.Placement(Base.Vector(), Base.Rotation(angle * 180.0 / pi, 0, 0)).multiply(local_cs.inverse()))


def axial_rotation_from_vector_and_image(origin, vector0, vector1):
    """
    Returns a rotation that transforms the ray with origin and vector equals to vector0 to the ray with same
    origin and vector vector1.
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
    """
    Base class to represent joints
    """
    pass


class AxialJoint(Joint):
    """
    Class used to represent joints that rotate around an axis
    """
    def __init__(self, axis_origin, axis_vector, normal_vector):
        self.axis_origin = axis_origin
        self.axis_vector = axis_vector
        self.normal_vector = normal_vector

    def compute_rotation(self, target):
        pointing_vector = projection_on_orthogonal_of_vector(target-self.axis_origin, self.axis_vector)
        angle = self.normal_vector.getAngle(pointing_vector)
        return axial_rotation_from_axis_and_angle(self.axis_origin, self.axis_vector, angle)


class CentralJoint(Joint):
    """
     Class used to represent joints that rotate around a point
     """

    def __init__(self, center, normal_vector):
        self.center = center
        self.normal_vector = normal_vector

    def compute_rotation(self, target):
        pointing_vector = target - self.center
        return axial_rotation_from_vector_and_image(self.center, self.normal_vector, pointing_vector)


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
        self.element_joint_map = element_joint_map
        self.movements = None

    def compute_movements(self):
        self.movements = []
        for element, joint in self.element_joint_map.items():
            movement = joint.compute_rotation(self.target)
            self.movements.append((element, movement))

    def make_movements(self):
        if self.movements is None:
            self.compute_movements()
        for (element, movement) in self.movements:
            element.Placement = movement.multiply(element.Placement)

