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
    def __init__(self, axis_origin, axis_vector):
        self.axis_origin = axis_origin
        self.axis_vector = axis_vector

    def compute_rotation(self, target, normal, light_direction):
        pointing_vector = projection_on_orthogonal_of_vector(target-self.axis_origin, self.axis_vector)
        angle = normal.getAngle(pointing_vector)
        return axial_rotation_from_axis_and_angle(self.axis_origin, self.axis_vector, angle/2)


class CentralJoint(Joint):
    """
    Class used to represent joints that rotate around a point
    """

    def __init__(self, center):
        self.center = center

    def compute_rotation(self, target, normal, light_direction):
        pointing = target - self.center
        pointing.normalize()
        light_unit = light_direction * 1.0
        light_unit.normalize()
        desired_normal = pointing - light_unit
        return axial_rotation_from_vector_and_image(self.center, normal, desired_normal)


class SolarTracking:
    pass


class SourceTracking(SolarTracking):
    pass


class TargetTracking(SolarTracking):
    """
    Class used to model the tracking of the light source by mobile elements in a scene.

    Parameters
    ----------
    target : point/line where light rays should be directed
    source_direction : direction of the sun rays
    scene : scene with mobile elements
    object_joint_map : dict with keys the mobile elements and values the corresponding joint
    object_normal_map : dict with keys the mobile elements and values the corresponding main direction
    """
    def __init__(self, target, source_direction, scene, object_joint_map, object_normal_map):
        self.target = target
        self.source_direction = source_direction
        self.scene = scene
        self.object_joint_map = object_joint_map
        self.object_normal_map = object_normal_map
        self.movements = None

    def compute_movements(self):
        self.movements = []
        for fc_object, joint in self.object_joint_map.items():
            movement = joint.compute_rotation(self.target, self.object_normal_map[fc_object], self.source_direction)
            self.movements.append((fc_object, movement))

    def make_movements(self):
        if self.movements is None:
            self.compute_movements()
        for (fc_object, movement) in self.movements:
            fc_object.Placement = movement.multiply(fc_object.Placement)
            for element, obj in self.scene.element_object_dict.items():
                if obj == fc_object:
                    element.Placement = movement.multiply(element.Placement)
            # self.scene.element_object_dict[element].Placement = movement.multiply(self.scene.element_object_dict[element])
            # modify self.element_normal_map[element]
        self.scene.recompute_boundbox()

    def undo_movements(self):
        if self.movements is None:
            return
        for (fc_object, movement) in self.movements:
            inv_movement = movement.inverse()
            fc_object.Placement = inv_movement.multiply(fc_object.Placement)
            for obj, elem in self.scene.element_object_dict.items():
                if elem == fc_object:
                    obj.Placement = inv_movement.multiply(obj.Placement)
