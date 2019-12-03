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
        pointing = projection_on_orthogonal_of_vector(target-self.axis_origin, self.axis_vector)
        pointing.normalize()
        light_unit = projection_on_orthogonal_of_vector(light_direction, self.axis_vector)
        light_unit.normalize()
        desired_normal = pointing - light_unit
        angle = normal.getAngle(desired_normal)
        # return axial_rotation_from_vector_and_image(self.axis_origin, normal, desired_normal)
        # angle = normal.getAngle(pointing_vector)
        return axial_rotation_from_axis_and_angle(self.axis_origin, self.axis_vector, angle)


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
    def __init__(self, target, source_direction, scene, object_joint_map=None, object_normal_map=None):
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

    def get_joints_from_objects(self):
        for obj in self.scene.objects:
            try:
                label = obj.Label
                start = label.find("(")
                end = label.find(")")
                string = label[start + 1:end]
                joint_label = string.split(',')[1]
                joint_obj = [obj2 for obj2 in self.scene.objects if obj2.Label == joint_label][0]
                if joint_obj.Shape.ShapeType == "Vertex":
                    center = joint_obj.Shape.Vertexes[0].Point
                    self.object_joint_map[obj] = CentralJoint(center)
                elif joint_obj.Shape.ShapeType == "Wire":
                    start = joint_obj.Shape.Vertexes[0].Point
                    end = joint_obj.Shape.Vertexes[1].Point
                    self.object_joint_map[obj] = AxialJoint(start, end-start)
                pass
            except:
                pass

    def get_normals_from_objects(self):
        for obj in self.scene.objects:
            try:
                label = obj.Label
                start = label.find("(")
                end = label.find(")")
                string = label[start + 1:end]
                normal_label = string.split(',')[2]
                normal_obj = [obj2 for obj2 in self.scene.objects if obj2.Label == normal_label][0]
                if normal_obj.Shape.ShapeType == "Wire":
                    start = normal_obj.Shape.Vertexes[0].Point
                    end = normal_obj.Shape.Vertexes[1].Point
                    self.object_normal_map[obj] = end-start
            except:
                pass

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
