"""Module otsun.scene for the modelization optical systems

The module defines the class `Scene` that models the elements in an optical system
"""

from .materials import Material, VolumeMaterial, SurfaceMaterial, TwoLayerMaterial
from .logging_unit import logger
from .math import correct_normal

EPSILON = 1E-6


class Scene:
    """
    Class used to define the Scene. It encodes all the objects
    that interfere with light rays.
    """

    def __init__(self, objects):
        self.objects = objects
        self.faces = []  # All the faces in the Scene
        self.solids = []  # All the solids in the Scene
        self.name_of_solid = {}
        self.materials = {}  # Assign the materials to objects
        # self.boundaries = {}  # Assign boundary faces to solids
        # self.sum_energy = 0.0 # Energy of the Scene
        self.epsilon = EPSILON # Tolerance for solid containment # 2 nm.
        self.boundbox = None
        self.element_object_dict = {}

        for obj in objects:
            # noinspection PyNoneFunctionAssignment
            logger.debug(f"loading object {obj.Label}")
            material = Material.get_from_label(obj.Label)
            if not material:
                logger.warning("Material not found for object %s", obj.Label)
                continue
            shape = obj.Shape
            faces = shape.Faces
            solids = shape.Solids
            if isinstance(material, VolumeMaterial):
                if not solids:
                    logger.warning(f"Volume material {material.name} associated to object {obj.Label} without solids")
                for solid in solids:
                    logger.debug(f"Assigning material {material.name} to solid {solid} in {obj.Label}")
                    self.name_of_solid[solid] = obj.Label
                    self.materials[solid] = material
                    self.element_object_dict[solid] = obj
                self.solids.extend(solids)
                self.faces.extend(faces)
            elif isinstance(material, SurfaceMaterial) or isinstance(material, TwoLayerMaterial):
                if solids:
                    logger.warning(f"Surface material {material.name} associated to object {obj.Label} with solids")
                if not faces:
                    logger.warning(f"Surface material {material.name} associated to object {obj.Label} without faces")
                for face in faces:
                    logger.debug(f"Assigning material {material.name} to face {face} in {obj.Label}")
                    self.materials[face] = material
                    self.element_object_dict[face] = obj
                self.faces.extend(faces)
            else:
                logger.warning(f"Material {material.name} associated to {obj.Label} is not Surface or Volume material")
                continue
            # solids = obj.Shape.Solids
            # faces = obj.Shape.Faces
            # if solids and isinstance(material, VolumeMaterial):  # Object is a solid
            #     for solid in solids:
            #         logger.debug(f"Assigning material {material.name} to solid {solid} in {obj.Label}")
            #         self.name_of_solid[solid] = obj.Label
            #         self.materials[solid] = material
            #         self.element_object_dict[solid] = obj
            # else:  # Object is made of faces
            #     for face in faces:
            #         logger.debug(f"Assigning material {material.name} to face {face} in {obj.Label}")
            #         self.materials[face] = material
            #         self.element_object_dict[face] = obj
            # self.solids.extend(solids)
            # self.faces.extend(faces)
            # # self.materials[obj] = material ## Cal?
            if not self.boundbox:
                self.boundbox = obj.Shape.BoundBox
            else:
                self.boundbox.add(obj.Shape.BoundBox)

        self.remove_duplicate_faces()

        self.diameter = self.boundbox.DiagonalLength

    def recompute_boundbox(self):
        """
        Recomputes the boundbox, so that all objects are contained in it
        """
        boundbox = None
        for elem in self.solids:
            if boundbox is None:
                boundbox = elem.BoundBox
            else:
                boundbox.add(elem.BoundBox)
        for elem in self.faces:
            if boundbox is None:
                boundbox = elem.BoundBox
            else:
                boundbox.add(elem.BoundBox)
        self.boundbox = boundbox

    def remove_duplicate_faces(self):
        """
        Removes redundant faces, so that the computation of next_intersection does not find duplicated points.
        """
        logger.debug("Removing duplicate faces")
        faces_with_material = [face for face in self.faces if
                               face in self.materials]
        faces_no_material = [face for face in self.faces if
                             face not in self.materials]
        if not faces_no_material:
            logger.debug("Done Removing duplicate faces (nothing to do)")
            return
        complex_no_material = faces_no_material[0]
        for (i,face) in enumerate(faces_no_material[1:]):
            logger.debug(f"Fusing face {i} of {len(faces_no_material)-1}")
            complex_no_material = complex_no_material.fuse(face)
        for (i,face) in enumerate(faces_with_material):
            logger.debug(f"Cutting face {i} of {len(faces_with_material) - 1}")
            complex_no_material = complex_no_material.cut(face)
        self.faces = faces_with_material
        self.faces.extend(complex_no_material.Faces)
        logger.debug("Done Removing duplicate faces")

    def solid_at_point(self, point):
        """
        Returns the solid that a point is inside.
        """
        for solid in self.solids:
            if solid.isInside(point, self.epsilon, False):
                return solid
        return None

    def face_at_point(self, point):
        """
        Returns the face that a point is inside.
        """
        for face in self.faces:
            if face.isInside(point, self.epsilon, True):
                return face
        return None

    def next_solid_at_point_in_direction(self, point, normal, direction):
        """
        Returns the next solid found in the given direction from a given point
        """
        external_normal = correct_normal(normal, direction)
        point_plus_epsilon = point + external_normal * (-2) * self.epsilon
        return self.solid_at_point(point_plus_epsilon)
