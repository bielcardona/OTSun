"""Module otsun.movements for computing movements

The module defines the classes `Joint` (with subclasses `AxialJoint` and `CentralJoint`) and `MultiTracking`
"""
from typing import Any

from FreeCAD import Base

from .scene import Scene
from .math import one_orthogonal_vector, projection_on_orthogonal_of_vector, EPSILON
from numpy import pi
from numpy.random import normal as random_normal

from autologging import traced
import logging
logger = logging.getLogger(__name__)


def orientation(u: Base.Vector, v: Base.Vector, w: Base.Vector) -> int:
    """Computes the orientation (+/-1) of a basis"""
    det = u.dot(v.cross(w))
    if det >= 0:
        return +1
    else:
        return -1


def signed_angle(axis: Base.Vector, v1: Base.Vector, v2: Base.Vector) -> float:
    """
    Computes the (signed) angle from v1 to v2, wrt the axis vector
    """
    angle = v1.getAngle(v2)
    angle *= orientation(axis, v1, v2)
    return angle


def apply_movement_to_vector(movement: Base.Placement, vector: Base.Vector) -> Base.Vector:
    pA = Base.Vector(0, 0, 0)
    pB = pA + vector
    pAprime = movement.multVec(pA)
    pBprime = movement.multVec(pB)
    return pBprime - pAprime


def axial_rotation_from_axis_and_angle(axis_origin: Base.Vector, axis_dir: Base.Vector,
                                       angle:float) -> Base.Placement:
    """
    Returns a rotation with given axis and angle
    """
    z_axis = Base.Vector(0, 0, 1)
    local_cs = Base.Placement(axis_origin, Base.Rotation(z_axis, axis_dir))
    return local_cs.multiply(
        Base.Placement(Base.Vector(), Base.Rotation(angle * 180.0 / pi, 0, 0)).multiply(local_cs.inverse()))


def axial_rotation_from_vector_and_image(origin: Base.Vector, vector0: Base.Vector,
                                         vector1: Base.Vector, error: float=0) -> Base.Placement:
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
    angle += error
    return axial_rotation_from_axis_and_angle(origin, normal, angle)


def biaxial_rotation_from_vector_and_image(origin: Base.Vector, vector_v: Base.Vector, vector_h: Base.Vector,
                                           current_normal: Base.Vector, desired_normal: Base.Vector,
                                           error_v: float=0, error_h: float=0) -> Base.Placement:
    # First rotate around vector_h, so that the image of vector_v is orthogonal
    # to desired_normal
    image_vector_horizontal = desired_normal.cross(vector_h)
    g1 = axial_rotation_from_vector_and_image(origin, vector_v, image_vector_horizontal, error_v)
    new_normal = apply_movement_to_vector(g1, current_normal)
    g2 = axial_rotation_from_vector_and_image(origin, new_normal, desired_normal, error_h)
    return g2.multiply(g1)


def get_labels(obj: Any) -> list[str]:
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
    def compute_rotation_to_point(self,
                                  target: Base.Vector,
                                  normal: Base.Vector,
                                  light_direction: Base.Vector) -> Base.Placement:
        pass

    def compute_rotation_to_direction(self,
                                      normal: Base.Vector,
                                      light_direction: Base.Vector) -> Base.Placement:
        pass

@traced(logger)
class AxialJoint(Joint):
    """
    Class used to represent joints that rotate around an axis
    """

    def __init__(self,
                 axis_origin: Base.Vector,
                 axis_vector: Base.Vector,
                 sigma: float | None = None):
        self.axis_origin: Base.Vector = axis_origin
        self.axis_vector: Base.Vector = axis_vector
        self.sigma: float | None = sigma
        self.error: float = 0
        if self.sigma is not None:
            self.error = random_normal(scale=self.sigma/1000)


    def compute_rotation_to_point(self,
                                  target: Base.Vector,
                                  normal: Base.Vector,
                                  light_direction: Base.Vector) -> Base.Placement:
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
        angle += self.error
        return axial_rotation_from_axis_and_angle(self.axis_origin, self.axis_vector, angle)

    def compute_rotation_to_direction(self,
                                      normal: Base.Vector,
                                      light_direction: Base.Vector) -> Base.Placement:
        """
        Computes a rotation so that ray in the direction of light_direction is returned to its origin.
        """
        light_unit = projection_on_orthogonal_of_vector(light_direction, self.axis_vector)
        light_unit.normalize()
        desired_normal = light_unit * (-1.0)
        angle = signed_angle(self.axis_vector, normal, desired_normal)  # normal.getAngle(desired_normal)
        angle += self.error
        return axial_rotation_from_axis_and_angle(self.axis_origin, self.axis_vector, angle)


@traced(logger)
class TwoOrthogonalAxialJoint(Joint):
    """
    Class used to represent joints with two orthogonal axis
    """

    def __init__(self,
                 axis_origin: Base.Vector,
                 axis_vector_vertical: Base.Vector,
                 axis_vector_horizontal: Base.Vector,
                 sigma_v: float|None = None,
                 sigma_h: float|None = None):
        self.axis_origin: Base.Vector = axis_origin
        self.axis_vector_vertical: Base.Vector = axis_vector_vertical
        self.axis_vector_horizontal: Base.Vector = axis_vector_horizontal
        self.sigma_v: float = sigma_v
        self.sigma_h: float = sigma_h
        self.error_v: float = 0
        if self.sigma_v is not None:
            self.error_v = random_normal(scale=self.sigma_v/1000)
        self.error_h: float = 0
        if self.sigma_h is not None:
            self.error_h = random_normal(scale=self.sigma_h/1000)


    def compute_rotation_to_point(self, target: Base.Vector, normal: Base.Vector,
                                  light_direction: Base.Vector) -> Base.Placement:
        """
        Computes a rotation so that ray in the direction of light_direction is reflected and hits the target.
        """
        pointing = target - self.axis_origin
        pointing.normalize()
        light_unit = light_direction * 1.0
        light_unit.normalize()
        desired_normal = pointing - light_unit
        return biaxial_rotation_from_vector_and_image(self.axis_origin, self.axis_vector_vertical,
                                                      self.axis_vector_horizontal, normal, desired_normal,
                                                      self.error_v, self.error_h)

    def compute_rotation_to_direction(self,
                                      normal: Base.Vector,
                                      light_direction: Base.Vector) -> Base.Placement:
        """
        Computes a rotation so that ray in the direction of light_direction is returned to its origin.
        """
        light_unit = light_direction * 1.0
        light_unit.normalize()
        desired_normal = light_unit * (-1.0)
        return biaxial_rotation_from_vector_and_image(self.axis_origin, self.axis_vector_vertical,
                                                      self.axis_vector_horizontal, normal, desired_normal,
                                                      self.error_v, self.error_h)


@traced(logger)
class CentralJoint(Joint):
    """
    Class used to represent joints that rotate around a point
    """

    def __init__(self, center: Base.Vector):
        self.center: Base.Vector = center

    def compute_rotation_to_point(self, target: Base.Vector, normal: Base.Vector,
                                  light_direction: Base.Vector) -> Base.Placement:
        """
        Computes a rotation so that ray in the direction of light_direction is reflected and hits the target.
        """
        pointing = target - self.center
        pointing.normalize()
        light_unit = light_direction * 1.0
        light_unit.normalize()
        desired_normal = pointing - light_unit
        return axial_rotation_from_vector_and_image(self.center, normal, desired_normal)

    def compute_rotation_to_direction(self, normal: Base.Vector, light_direction: Base.Vector) -> Base.Placement:
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

    def __init__(self,
                 source_direction: Base.Vector,
                 scene: Scene,
                 analyze_labels: bool=True):
        self.scene: Scene = scene
        self.source_direction: Base.Vector = source_direction
        self.object_joint_map: dict[Any, Joint] = {}
        self.object_normal_map: dict[Any, Base.Vector] = {}
        self.object_target_map: dict[Any, Base.Vector|None] = {}
        self.object_movements_map: dict[Any, Base.Placement] = {}
        if analyze_labels:
            self.get_data_from_objects()
            self.compute_movements()

    def get_object_from_label(self, label: str) -> Any:
        return [obj2 for obj2 in self.scene.objects if obj2.Label == label][0]

    @staticmethod
    def get_dimension_of_object(obj: Any) -> int | None:
        obj_type = obj.Shape.ShapeType
        if obj_type == 'Vertex':
            return 0
        if obj_type in ['Wire', 'Edge']:
            return 1
        return None

    @staticmethod
    def get_point_from_vertex(vertex: Any) -> Base.Vector:
        return vertex.Shape.Vertexes[0].Point

    @staticmethod
    def get_point_and_vector_from_edge(edge: Any) -> tuple[Base.Vector, Base.Vector]:
        start = edge.Shape.Vertexes[0].Point
        end = edge.Shape.Vertexes[1].Point
        return start, end - start

    @staticmethod
    def get_center_two_edges(ob1: Any, ob2: Any) -> Base.Vector:
        sh1 = ob1.Shape
        sh2 = ob2.Shape
        d, pairs_data, _ = sh1.distToShape(sh2)
        pairs = pairs_data[0]
        return (pairs[0]+pairs[1])/2

    def get_target_from_label(self, label: str) -> Base.Vector:
        """
        Finds the target of an object according to its label
        """
        target_obj = [obj2 for obj2 in self.scene.objects if obj2.Label == label][0]
        if target_obj.Shape.ShapeType == "Vertex":
            target = target_obj.Shape.Vertexes[0].Point
            return target

    def get_deviation(self, axis_label: str) -> float | None:
        try:
            return self.scene.extra_data['axis_deviation'][axis_label]
        except KeyError:
            return None

    def get_data_from_objects(self) -> None:
        """
        Gets all data of joints from the objects in the scene
        """
        cached_joints = dict()
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
                sigma1 = self.get_deviation(objects[0].Label)
                sigma2 = self.get_deviation(objects[1].Label)
                if (labels[1], labels[2]) not in cached_joints:
                    cached_joints[(labels[1], labels[2])] = TwoOrthogonalAxialJoint(
                        center, axis2_vector, axis1_vector, sigma2, sigma1)
                self.object_joint_map[obj] = cached_joints[(labels[1], labels[2])]
                self.object_normal_map[obj] = normal
                pass
            elif dimensions in [[1, 1, 0], [1, 1]]:  # Axial joint
                axis_point, axis_vector = self.get_point_and_vector_from_edge(objects[0])
                _, normal = self.get_point_and_vector_from_edge(objects[1])
                sigma = self.get_deviation(objects[0].Label)
                if labels[1] not in cached_joints:
                    cached_joints[labels[1]] = AxialJoint(axis_point, axis_vector, sigma)
                self.object_joint_map[obj] = cached_joints[labels[1]]
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

    def compute_movements(self) -> None:
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

    def make_movements(self) -> None:
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

    def undo_movements(self) -> None:
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
