"""Module otsun.movements for computing movements

The module defines the classes `Joint` (with subclasses `AxialJoint` and `CentralJoint`) and `MultiTracking`
"""
from FreeCAD import Base
from .math import one_orthogonal_vector, projection_on_orthogonal_of_vector, EPSILON
from numpy import pi

from autologging import traced
from .logging_unit import logger


def orientation(u, v, w):
    """Computes the orientation (+/-1) of a basis"""
    det = u.dot(v.cross(w))
    if det >= 0:
        return +1
    else:
        return -1


def signed_angle(axis, v1, v2):
    """
    Computes the (signed) angle from v1 to v2, wrt the axis vector
    """
    angle = v1.getAngle(v2)
    angle *= orientation(axis, v1, v2)
    return angle


def apply_movement_to_vector(movement, vector):
    pA = Base.Vector(0, 0, 0)
    pB = pA + vector
    pAprime = movement.multVec(pA)
    pBprime = movement.multVec(pB)
    return pBprime - pAprime


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


def axial_rotation_from_axis_vector_and_image(axis_origin, axis_dir, vector0, vector1):
    angle = vector0.getAngle(vector1)
    return axial_rotation_from_axis_and_angle(axis_origin, axis_dir, angle)


def biaxial_rotation_from_vector_and_image(origin, vector_v, vector_h, current_normal, desired_normal):
    # First rotate around vector_h, so that the image of vector_v is orthogonal
    # to desired_normal
    image_vector_horizontal = desired_normal.cross(vector_h)
    g1 = axial_rotation_from_vector_and_image(origin, vector_v, image_vector_horizontal)
    new_normal = apply_movement_to_vector(g1, current_normal)
    g2 = axial_rotation_from_vector_and_image(origin, new_normal, desired_normal)
    return g2.multiply(g1)


def get_labels(obj):
    """Gets the different labels of a given FreeCAD object"""
    label = obj.Label
    start = label.find("(")
    end = label.find(")")
    if start < 0 or end < 0:
        return []
    string = label[start + 1:end]
    return string.split(',')


@traced(logger)
class Joint:
    """
    Base class to represent joints
    """
    pass


@traced(logger)
class AxialJoint(Joint):
    """
    Class used to represent joints that rotate around an axis
    """

    def __init__(self, axis_origin, axis_vector):
        self.axis_origin = axis_origin
        self.axis_vector = axis_vector

    def compute_rotation_to_point(self, target, normal, light_direction):
        """
        Computes a rotation so that ray in the direction of light_direction is reflected and hits the target.
        """
        pointing = projection_on_orthogonal_of_vector(target - self.axis_origin, self.axis_vector)
        pointing.normalize()
        light_unit = projection_on_orthogonal_of_vector(light_direction, self.axis_vector)
        if light_unit.Length > EPSILON:
            light_unit.normalize()
            desired_normal = pointing - light_unit
            angle = signed_angle(self.axis_vector, normal, desired_normal)  # normal.getAngle(desired_normal)
        else:
            angle = 0
        return axial_rotation_from_axis_and_angle(self.axis_origin, self.axis_vector, angle)

    def compute_rotation_to_direction(self, normal, light_direction):
        """
        Computes a rotation so that ray in the direction of light_direction is returned to its origin.
        """
        light_unit = projection_on_orthogonal_of_vector(light_direction, self.axis_vector)
        light_unit.normalize()
        desired_normal = light_unit * (-1.0)
        angle = signed_angle(self.axis_vector, normal, desired_normal)  # normal.getAngle(desired_normal)
        return axial_rotation_from_axis_and_angle(self.axis_origin, self.axis_vector, angle)


@traced(logger)
class TwoOrthogonalAxialJoint(Joint):
    """
    Class used to represent joints with two orthogonal axis
    """

    def __init__(self, axis_origin, axis_vector_vertical, axis_vector_horizontal):
        self.axis_origin = axis_origin
        self.axis_vector_vertical = axis_vector_vertical
        self.axis_vector_horizontal = axis_vector_horizontal

    def compute_rotation_to_point(self, target, normal, light_direction):
        """
        Computes a rotation so that ray in the direction of light_direction is reflected and hits the target.
        """
        pointing = target - self.axis_origin
        pointing.normalize()
        light_unit = light_direction * 1.0
        light_unit.normalize()
        desired_normal = pointing - light_unit
        return biaxial_rotation_from_vector_and_image(self.axis_origin, self.axis_vector_vertical,
                                                      self.axis_vector_horizontal, normal, desired_normal)

    def compute_rotation_to_direction(self, normal, light_direction):
        """
        Computes a rotation so that ray in the direction of light_direction is returned to its origin.
        """
        light_unit = light_direction * 1.0
        light_unit.normalize()
        desired_normal = light_unit * (-1.0)
        return biaxial_rotation_from_vector_and_image(self.axis_origin, self.axis_vector_vertical,
                                                      self.axis_vector_horizontal, normal, desired_normal)


@traced(logger)
class CentralJoint(Joint):
    """
    Class used to represent joints that rotate around a point
    """

    def __init__(self, center):
        self.center = center

    def compute_rotation_to_point(self, target, normal, light_direction):
        """
        Computes a rotation so that ray in the direction of light_direction is reflected and hits the target.
        """
        pointing = target - self.center
        pointing.normalize()
        light_unit = light_direction * 1.0
        light_unit.normalize()
        desired_normal = pointing - light_unit
        return axial_rotation_from_vector_and_image(self.center, normal, desired_normal)

    def compute_rotation_to_direction(self, normal, light_direction):
        """
        Computes a rotation so that ray in the direction of light_direction is returned to its origin.
        """
        light_unit = light_direction * 1.0
        light_unit.normalize()
        desired_normal = light_unit * (-1.0)
        return axial_rotation_from_vector_and_image(self.center, normal, desired_normal)


@traced(logger)
class MultiTracking:
    """
    Class that represents how objects are linked to joints and implements their movement.
    """

    def __init__(self, source_direction, scene, analyze_labels=True):
        self.scene = scene
        self.source_direction = source_direction
        self.object_joint_map = {}
        self.object_normal_map = {}
        self.object_target_map = {}
        self.object_movements_map = {}
        if analyze_labels:
            self.get_data_from_objects()
            self.compute_movements()

    def get_object_from_label(self, label):
        return [obj2 for obj2 in self.scene.objects if obj2.Label == label][0]

    @staticmethod
    def get_dimension_of_object(obj):
        obj_type = obj.Shape.ShapeType
        if obj_type == 'Vertex':
            return 0
        if obj_type in ['Wire', 'Edge']:
            return 1
        return None

    @staticmethod
    def get_point_from_vertex(vertex):
        return vertex.Shape.Vertexes[0].Point

    @staticmethod
    def get_point_and_vector_from_edge(edge):
        start = edge.Shape.Vertexes[0].Point
        end = edge.Shape.Vertexes[1].Point
        return start, end - start

    @staticmethod
    def get_center_two_edges(ob1, ob2):
        sh1 = ob1.Shape
        sh2 = ob2.Shape
        d, pairs_data, _ = sh1.distToShape(sh2)
        pairs = pairs_data[0]
        return (pairs[0]+pairs[1])/2

    # def get_joint_from_label(self, label):
    #     """
    #     Builds a joint according to the label of the object
    #     """
    #     joint_obj = [obj2 for obj2 in self.scene.objects if obj2.Label == label][0]
    #     if joint_obj.Shape.ShapeType == "Vertex":
    #         center = joint_obj.Shape.Vertexes[0].Point
    #         return CentralJoint(center)
    #     elif joint_obj.Shape.ShapeType in ["Wire", "Edge"]:
    #         start = joint_obj.Shape.Vertexes[0].Point
    #         end = joint_obj.Shape.Vertexes[1].Point
    #         return AxialJoint(start, end - start)
    #
    # def get_normal_from_label(self, label):
    #     """
    #     Finds the principal vector of an object according to its label
    #     """
    #     normal_obj = [obj2 for obj2 in self.scene.objects if obj2.Label == label][0]
    #     if normal_obj.Shape.ShapeType in ["Wire", "Edge"]:
    #         start = normal_obj.Shape.Vertexes[0].Point
    #         end = normal_obj.Shape.Vertexes[1].Point
    #         return end - start

    def get_target_from_label(self, label):
        """
        Finds the target of an object according to its label
        """
        target_obj = [obj2 for obj2 in self.scene.objects if obj2.Label == label][0]
        if target_obj.Shape.ShapeType == "Vertex":
            target = target_obj.Shape.Vertexes[0].Point
            return target

    def get_data_from_objects(self):
        """
        Gets all data of joints from the objects in the scene
        """
        for obj in self.scene.objects:
            labels = get_labels(obj)
            if len(labels) < 2:
                # no movement
                continue
            objects = [self.get_object_from_label(label) for label in labels[1:]]
            dimensions = [self.get_dimension_of_object(ob) for ob in objects]
            if dimensions in [[1, 1, 1, 0], [1, 1, 1]]:  # Biaxial joint
                _, axis1_vector = self.get_point_and_vector_from_edge(objects[0])
                _, axis2_vector = self.get_point_and_vector_from_edge(objects[1])
                center = self.get_center_two_edges(objects[0], objects[1])
                _, normal = self.get_point_and_vector_from_edge(objects[2])
                self.object_joint_map[obj] = TwoOrthogonalAxialJoint(center, axis2_vector, axis1_vector)
                self.object_normal_map[obj] = normal
                pass
            elif dimensions in [[1, 1, 0], [1, 1]]:  # Axial joint
                axis_point, axis_vector = self.get_point_and_vector_from_edge(objects[0])
                _, normal = self.get_point_and_vector_from_edge(objects[1])
                self.object_joint_map[obj] = AxialJoint(axis_point, axis_vector)
                self.object_normal_map[obj] = normal
                pass
            elif dimensions in [[0, 1, 0], [0, 1]]:  # Central joint
                center = self.get_point_from_vertex(objects[0])
                _, normal = self.get_point_and_vector_from_edge(objects[1])
                self.object_joint_map[obj] = CentralJoint(center)
                self.object_normal_map[obj] = normal
            else:
                raise Exception("Malformed joint characterization")
            if dimensions[-1] == 0:
                self.object_target_map[obj] = self.get_point_from_vertex(objects[-1])
            else:
                self.object_target_map[obj] = None
            #
            # if len(labels) == 2:
            #     # gives joint but not normal
            #     raise Exception("Need to give joint and normal")
            # self.object_joint_map[obj] = self.get_joint_from_label(labels[1])
            # self.object_normal_map[obj] = self.get_normal_from_label(labels[2])
            # if len(labels) == 3:
            #     # no target given
            #     self.object_target_map[obj] = None
            # else:
            #     self.object_target_map[obj] = self.get_target_from_label(labels[3])

    def compute_movements(self):
        """
        Computes the movements so that the objects point in the right direction according to their target
        """
        for obj in self.object_joint_map:
            joint = self.object_joint_map[obj]
            normal = self.object_normal_map[obj]
            target = self.object_target_map[obj]
            if target is None:
                # source tracking
                movement = joint.compute_rotation_to_direction(normal, self.source_direction)
            else:
                movement = joint.compute_rotation_to_point(target, normal, self.source_direction)
            self.object_movements_map[obj] = movement

    def make_movements(self):
        """
        Makes the computed movements
        """
        for obj in self.object_movements_map:
            movement = self.object_movements_map[obj]
            obj.Placement = movement.multiply(obj.Placement)
            for element, el_obj in self.scene.element_object_dict.items():
                if obj == el_obj:
                    element.Placement = movement.multiply(element.Placement)
        self.scene.recompute_boundbox()

    def undo_movements(self):
        """
        Undo the movements so that they are in their original position
        """
        for obj in self.object_movements_map:
            movement = self.object_movements_map[obj]
            movement = movement.inverse()
            obj.Placement = movement.multiply(obj.Placement)
            for element, el_obj in self.scene.element_object_dict.items():
                if obj == el_obj:
                    element.Placement = movement.multiply(element.Placement)
        self.scene.recompute_boundbox()
