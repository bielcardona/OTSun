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


def get_labels(obj):
    label = obj.Label
    start = label.find("(")
    end = label.find(")")
    if start < 0 or end < 0:
        return []
    string = label[start + 1:end]
    return string.split(',')


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

    def compute_rotation_to_point(self, target, normal, light_direction):
        pointing = projection_on_orthogonal_of_vector(target-self.axis_origin, self.axis_vector)
        pointing.normalize()
        light_unit = projection_on_orthogonal_of_vector(light_direction, self.axis_vector)
        light_unit.normalize()
        desired_normal = pointing - light_unit
        angle = normal.getAngle(desired_normal)
        # return axial_rotation_from_vector_and_image(self.axis_origin, normal, desired_normal)
        # angle = normal.getAngle(pointing_vector)
        return axial_rotation_from_axis_and_angle(self.axis_origin, self.axis_vector, angle)

    def compute_rotation_to_direction(self, normal, light_direction):
        light_unit = projection_on_orthogonal_of_vector(light_direction, self.axis_vector)
        light_unit.normalize()
        desired_normal = light_unit * (-1.0)
        angle = normal.getAngle(desired_normal)
        return axial_rotation_from_axis_and_angle(self.axis_origin, self.axis_vector, angle)


class CentralJoint(Joint):
    """
    Class used to represent joints that rotate around a point
    """

    def __init__(self, center):
        self.center = center

    def compute_rotation_to_point(self, target, normal, light_direction):
        pointing = target - self.center
        pointing.normalize()
        light_unit = light_direction * 1.0
        light_unit.normalize()
        desired_normal = pointing - light_unit
        return axial_rotation_from_vector_and_image(self.center, normal, desired_normal)

    def compute_rotation_to_direction(self, normal, light_direction):
        light_unit = light_direction * 1.0
        light_unit.normalize()
        desired_normal = light_unit * (-1.0)
        return axial_rotation_from_vector_and_image(self.center, normal, desired_normal)


class SolarTracking:
    def __init__(self):
        self.scene = None
        self.object_joint_map = None
        self.object_normal_map = None
        self.movements = None

    def get_joints_from_objects(self):
        for obj in self.scene.objects:
            try:
                joint_label = get_labels(obj)[1]
                joint_obj = [obj2 for obj2 in self.scene.objects if obj2.Label == joint_label][0]
                if joint_obj.Shape.ShapeType == "Vertex":
                    center = joint_obj.Shape.Vertexes[0].Point
                    self.object_joint_map[obj] = CentralJoint(center)
                elif joint_obj.Shape.ShapeType in ["Wire", "Edge"]:
                    start = joint_obj.Shape.Vertexes[0].Point
                    end = joint_obj.Shape.Vertexes[1].Point
                    self.object_joint_map[obj] = AxialJoint(start, end-start)
                pass
            except IndexError:
                pass

    def get_normals_from_objects(self):
        for obj in self.scene.objects:
            try:
                normal_label = get_labels(obj)[2]
                normal_obj = [obj2 for obj2 in self.scene.objects if obj2.Label == normal_label][0]
                if normal_obj.Shape.ShapeType in ["Wire", "Edge"]:
                    start = normal_obj.Shape.Vertexes[0].Point
                    end = normal_obj.Shape.Vertexes[1].Point
                    self.object_normal_map[obj] = end-start
            except IndexError:
                pass

    def compute_movements(self):
        # Must be subclassed
        return None

    def make_movements(self):
        if self.movements is None:
            self.compute_movements()
        for (fc_object, movement) in self.movements:
            fc_object.Placement = movement.multiply(fc_object.Placement)
            for element, obj in self.scene.element_object_dict.items():
                if obj == fc_object:
                    element.Placement = movement.multiply(element.Placement)
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


class SourceTracking(SolarTracking):
    def __init__(self, source_direction, scene, object_joint_map=None, object_normal_map=None):
        super().__init__()
        # self.target = target
        self.source_direction = source_direction
        self.scene = scene
        if object_joint_map is None:
            self.object_joint_map = {}
            self.get_joints_from_objects()
        else:
            self.object_joint_map = object_joint_map
        if object_normal_map is None:
            self.object_normal_map = {}
            self.get_normals_from_objects()
        else:
            self.object_normal_map = object_normal_map
        self.movements = None
    pass

    def compute_movements(self):
        self.movements = []
        for fc_object, joint in self.object_joint_map.items():
            movement = joint.compute_rotation_to_direction(self.object_normal_map[fc_object], self.source_direction)
            self.movements.append((fc_object, movement))


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
    def __init__(self, target, source_direction, scene, object_joint_map=None, object_normal_map=None):
        super().__init__()
        self.target = target
        self.source_direction = source_direction
        self.scene = scene
        if object_joint_map is None:
            self.object_joint_map = {}
            self.get_joints_from_objects()
        else:
            self.object_joint_map = object_joint_map
        if object_normal_map is None:
            self.object_normal_map = {}
            self.get_normals_from_objects()
        else:
            self.object_normal_map = object_normal_map
        self.movements = None

    def compute_movements(self):
        self.movements = []
        for fc_object, joint in self.object_joint_map.items():
            movement = joint.compute_rotation_to_point(self.target, self.object_normal_map[fc_object],
                                                       self.source_direction)
            self.movements.append((fc_object, movement))
