from .materials import Material
from .logging_unit import logger

class Scene:
    """
    Class used to define the Scene. It encodes all the objects
    that interfere with light rays.
    """

    def __init__(self, objects):
        self.faces = []  # All the faces in the Scene
        self.solids = []  # All the solids in the Scene
        self.materials = {}  # Assign the materials to objects
        # self.sum_energy = 0.0 # Energy of the Scene
        self.epsilon = 40 * 10.0 ** (-6.0)  # Tolerance for solid containment # 10 nm.
        self.boundbox = None

        for obj in objects:
            # noinspection PyNoneFunctionAssignment
            material = Material.get_from_label(obj.Label)
            if not material:
                logger.info("Material not found for object %s", obj.Label)
                continue
            solids = obj.Shape.Solids
            faces = obj.Shape.Faces
            if solids:  # Object is a solid
                for solid in solids:
                    self.materials[solid] = material
            else:  # Object is made of faces
                for face in faces:
                    self.materials[face] = material
            self.solids.extend(solids)
            self.faces.extend(faces)
            # self.materials[obj] = material ## Cal?
            if not self.boundbox:
                self.boundbox = obj.Shape.BoundBox
            else:
                self.boundbox.add(obj.Shape.BoundBox)

        self.diameter = self.boundbox.DiagonalLength
        self.remove_duplicate_faces()

    def remove_duplicate_faces(self):
        faces_with_material = [face for face in self.faces if
                               face in self.materials]
        faces_no_material = [face for face in self.faces if
                             face not in self.materials]
        if not faces_no_material:
            return
        complex_no_material = faces_no_material[0]
        for face in faces_no_material[1:]:
            complex_no_material = complex_no_material.fuse(face)
        for face in faces_with_material:
            complex_no_material = complex_no_material.cut(face)
        self.faces = faces_with_material
        self.faces.extend(complex_no_material.Faces)

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
