"""Module otsun.scene for the modelization optical systems

The module defines the class `Scene` that models the elements in an optical system
"""
from __future__ import annotations
from typing import Any

import FreeCAD
from Part import Face, Solid
from FreeCAD import Base

from .helpers import parameters_for_freecad_doc, get_document
from .materials import Material, VolumeMaterial, SurfaceMaterial, TwoLayerMaterial, MaterialNotFoundError
from .math import correct_normal
from autologging import traced
from pathlib import Path
import logging
logger = logging.getLogger(__name__)

EPSILON = 1E-6

@traced(logger)
class Scene:
    """
    Class used to define the Scene. It encodes all the objects
    that interfere with light rays.
    """

    def __init__(
            self, 
            objects: list[Any], 
            extra_data: dict[str, Any] | None = None
            ):
        self.objects: list[Any] = objects
        self.faces: list[Face] = []  # All the faces in the Scene
        self.solids: list[Solid] = []  # All the solids in the Scene
        self.name_of_solid: dict[Solid, str] = {}
        self.bb_of_solid: dict[Solid, Base.BoundBox] = {}  # Bounding box of solid
        self.materials: dict[int, Material] = {}  # Assign the materials to objects
        # self.boundaries = {}  # Assign boundary faces to solids
        # self.sum_energy = 0.0 # Energy of the Scene
        self.epsilon: float = EPSILON # Tolerance for solid containment # 2 nm.
        self.boundbox: Base.BoundBox | None = None
        self.element_object_dict: dict[ Face | Solid, Any] = {}
        self.aperture_pv: float | None = None
        self.aperture_th: float | None = None

        self.extra_data : dict[str, Any]
        if extra_data is None:
            self.extra_data = {}
        else:
            self.extra_data = extra_data
        # Format for extra_data:
        # extra_data['vacuum_material']: VolumeMaterial to be considered after 1st intersection
        if 'vacuum_material' in self.extra_data:
            vacuum_material_name = self.extra_data['vacuum_material'][0]
            material = Material.by_name.get(vacuum_material_name, None)
            if not material:
                logger.warning("Material %s not found for vacuum", vacuum_material_name)
            else:
                self.extra_data['vacuum_material'] = material
        if 'axis_deviation' in self.extra_data:
            pairs_list = self.extra_data['axis_deviation'][:]
            self.extra_data['axis_deviation'] = {}
            iterator = iter(pairs_list)
            pairs_iterator = zip(iterator, iterator)
            for axis, sigma in pairs_iterator:
                self.extra_data['axis_deviation'][axis] = sigma
        # extra_data['axis_deviation']: Dictionary where each key is the axis in the scene and the value its std dev
        if 'aperture_pv' in self.extra_data:
            self.aperture_pv = self.extra_data['aperture_pv']
        if 'aperture_th' in self.extra_data:
            self.aperture_th = self.extra_data['aperture_th']

        for obj in objects:
            # noinspection PyNoneFunctionAssignment
            logger.debug(f"loading object {obj.Label}")
            try:
                material = Material.get_from_label(obj.Label)
            except MaterialNotFoundError as e:
                logger.warning(f"{e}")
                continue
            if not material:
                logger.debug(f"Object {obj.Label} does not have a material")
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
                    self.materials[id(solid)] = material
                    self.element_object_dict[solid] = obj
                    self.bb_of_solid[solid] = solid.BoundBox
                self.solids.extend(solids)
                self.faces.extend(faces)
            elif isinstance(material, SurfaceMaterial) or isinstance(material, TwoLayerMaterial):
                if solids:
                    logger.warning(f"Surface material {material.name} associated to object {obj.Label} with solids")
                if not faces:
                    logger.warning(f"Surface material {material.name} associated to object {obj.Label} without faces")
                for face in faces:
                    logger.debug(f"Assigning material {material.name} to face {face} in {obj.Label}")
                    self.materials[id(face)] = material
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

        # self.remove_duplicate_faces()
        if not self.boundbox:
            raise ValueError("No objects in the scene")

        self.diameter : float = self.boundbox.DiagonalLength

    @classmethod
    def from_freecad_document(
            cls,
            doc: str | Path | FreeCAD.Document,
            extra_data: dict[str, Any] | None = None
    ) -> Scene:
        """
        Creates a Scene from a FreeCAD document
        """
        doc = get_document(doc)
        objects = doc.Objects
        extra_data_from_doc = parameters_for_freecad_doc(doc)
        if extra_data:
            extra_data_from_doc.update(extra_data)
            # Overwrite extra_data_from_doc with the one passed as argument if it exists
            # If aperture_pv is given in extra_data, it will overwrite the one from the document
        return cls(objects, extra_data_from_doc)

    def recompute_boundbox(self) -> None:
        """
        Recomputes the boundbox, so that all objects are contained in it
        """
        boundbox = None
        for solid in self.solids:
            if boundbox is None:
                boundbox = solid.BoundBox
            else:
                boundbox.add(solid.BoundBox)
        for face in self.faces:
            if boundbox is None:
                boundbox = face.BoundBox
            else:
                boundbox.add(face.BoundBox)
        self.boundbox = boundbox

    def remove_duplicate_faces(self) -> None:
        """
        Removes redundant faces, so that the computation of next_intersection does not find duplicated points.
        """
        logger.debug("Removing duplicate faces")
        faces_with_material = [face for face in self.faces if
                               id(face) in self.materials]
        faces_no_material = [face for face in self.faces if
                             id(face) not in self.materials]
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

    def solid_at_point(self, point: Base.Vector) -> Solid | None:
        """
        Returns the solid that a point is inside.
        """
        for solid in self.solids:
            bb = self.bb_of_solid[solid]
            if bb.isInside(point):
                if solid.isInside(point, self.epsilon, False):
                    return solid
        return None

    def face_at_point(self, point: Base.Vector) -> Face | None:
        """
        Returns the face that a point is inside.
        """
        for face in self.faces:
            if face.isInside(point, self.epsilon, True):
                return face
        return None

    def next_solid_at_point_in_direction(
            self, 
            point: Base.Vector, 
            normal: Base.Vector, 
            direction: Base.Vector
        ) -> Solid|None:
        """
        Returns the next solid found in the given direction from a given point
        """
        external_normal = correct_normal(normal, direction)
        point_plus_epsilon = point + external_normal * (-2) * self.epsilon
        return self.solid_at_point(point_plus_epsilon)
