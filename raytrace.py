# -*- coding: utf-8 -*-
"""
Created on Fri Nov 25 15:10:31 2016

"""

from FreeCAD import Base
import Part
import numpy as np
import itertools
import random
import math


# ---
# Helper functions for reflexion and refraction
# ---
def reflexion(incident, normal):
    """
    Implementation of law of reflexion
    """
    # We assume all vectors are normalized
    return incident - normal * 2.0 * normal.dot(incident), "Reflexion"


# noinspection PyUnresolvedReferences,PyUnresolvedReferences
def refraction(incident, normal, n1, n2):
    """
    Implementation of Snell's law of refraction
    """
    # https://en.wikipedia.org/wiki/Snell's_law#Vector_form
    mynormal = normal * 1.0
    if mynormal.dot(incident) > 0:
        # noinspection PyAugmentAssignment
        mynormal = mynormal * (-1.0)
    r = n1 / n2
    c1 = - mynormal.dot(incident)
    c2sq = 1.0 - r * r * (1.0 - c1 * c1)
    if c2sq < 0:
        return reflexion(incident, normal)
    c2 = c2sq ** 0.5
    return incident * r + mynormal * (r * c1 - c2), "Refraction"


def polar_to_cartesian(phi, theta):
    """
    Convert polar coordinates (given in degrees) to cartesian
    :param phi:
    :param theta:
    :return:
    """
    rad = math.acos(-1.0) / 180.0
    x = -math.sin(theta * rad) * math.cos(phi * rad)
    y = -math.sin(theta * rad) * math.sin(phi * rad)
    z = -math.cos(theta * rad)
    return Base.Vector(x, y, z)


# ---
# Helper functions for input of functions
# ---
def constant_function(c):
    return lambda x: c


def tabulated_function(xvalues, yvalues):
    return lambda x: np.interp(x, xvalues, yvalues)


# ---
# Classes for materials
# ---

# region Materials

class Material(object):
    """
    Class used to represent materials and their physical properties
    """
    by_name = {}

    def __init__(self, name, properties):
        self.by_name[name] = self
        self.name = name
        self.properties = properties

    @classmethod
    def get_from_label(cls, label):
        if ("(" not in label) or (")" not in label):
            return None
        start = label.find("(")
        end = label.find(")")
        name = label[start + 1:end]
        return cls.by_name.get(name, None)

    @classmethod
    def create(cls, name, properties):
        _ = cls(name, properties)

    def change_of_direction(self, ray, normal_vector):
        pass
        # Compute new direction (ray.current_material)


class VolumeMaterial(Material, object):
    def __init__(self, name, properties):
        """
        Initializes a Volume Material. The properties parameter must be a dict with the physical properties
        describing the material. At least, the following must be provided:
        'index_of_refraction': index of refraction of the material, as a function of its wavelength.
        """
        super(VolumeMaterial, self).__init__(name, properties)
        self.kind = 'Volume'

    def change_of_direction(self, ray, normal_vector):
        wavelength = ray.wavelength
        n1 = ray.current_medium.properties['index_of_refraction'](wavelength)
        n2 = self.properties['index_of_refraction'](wavelength)
        direction, phenomenon = refraction(ray.directions[-1], normal_vector, n1, n2)
        if phenomenon == "Refraction":
            ray.current_medium = self
        return direction, phenomenon
        pass


def create_simple_volume_material(name, index_of_refraction):
    VolumeMaterial.create(name, {'index_of_refraction': constant_function(index_of_refraction)})


create_simple_volume_material("Vacuum", 1.0)
# noinspection PyNoneFunctionAssignment
vacuum_medium = VolumeMaterial.by_name["Vacuum"]


class SurfaceMaterial(Material, object):
    def __init__(self, name, properties):
        """
        Initializes a Surface Material. The properties parameter must be a dict with the physical properties
        describing the material. At least, the following must be provided:
        'probability_of_reflexion': probability that a photon gets reflected, as a function of its wavelength.
        'probability_of_absortion': probability that a photon gets absorbed, as a function of its wavelength.
        'transmitance': probability that a photon passes through the material, as a function of its wavelength.
        """
        super(SurfaceMaterial, self).__init__(name, properties)
        self.kind = 'Surface'

    def change_of_direction(self, ray, normal_vector):
        # TODO: Implementar bé. El que hi ha és provisional
        phenomena = ["Reflexion", "Absortion"]
        probabilities = [self.properties['probability_of_reflexion'](ray.properties['wavelength']),
                         self.properties['probability_of_absortion'](ray.properties['wavelength'])]
        phenomenon = np.random.choice(phenomena, 1, p = probabilities)[0]
        if phenomenon == 'Reflexion':
            return reflexion(ray.directions[-1], normal_vector)
        if phenomenon == "Absortion":
            return Base.Vector(0.0, 0.0, 0.0), "Absortion"
        pass

    def scatter_direction(self, ray, direction):
        # TODO: pensar
        pass


def create_opaque_simple_material(name, por):
    SurfaceMaterial.create(name, {'probability_of_reflexion': constant_function(por),
                                  'probability_of_absortion': constant_function(1 - por),
                                  'transmitance': constant_function(0)})


def create_transparent_simple_material(name, por):
    SurfaceMaterial.create(name, {'probability_of_reflexion': constant_function(por),
                                  'probability_of_absortion': constant_function(0),
                                  'transmitance': constant_function(1 - por)})

def create_mirror_simple_material(name, properties):
    props = properties
    SurfaceMaterial.create(name, {'probability_of_reflexion': constant_function(props['por']),
                                  'probability_of_absortion': constant_function(1 - props['por']),
                                  'transmitance': constant_function(0),
                                  'energy_collector': False})
 
# endregion

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
        self.epsilon = 0.00001  # Tolerance for solid containment
        self.boundbox = None

        for obj in objects:
            # noinspection PyNoneFunctionAssignment
            material = Material.get_from_label(obj.Label)
            if not material:
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
            if solid.isInside(point, self.epsilon, True):
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


class SunWindow:
    """
    Class used to define a sun window, defined as a rectangle in space such
    that the projection of the scene on the plane is enclosed in the rectangle
    """

    @staticmethod
    def find_min_rectangle(points, normal):
        min_area = None
        for (p, q) in itertools.combinations(points, 2):
            v1 = q - p
            if v1.Length < 0.001:  # TODO: customize
                continue
            v1.normalize()
            v2 = v1.cross(normal)
            v2.normalize()
            xs = [v1.dot(r - p) for r in points]
            ys = [v2.dot(r - p) for r in points]
            minx = min(xs)
            maxx = max(xs)
            miny = min(ys)
            maxy = max(ys)
            length1 = maxx - minx
            length2 = maxy - miny
            area = length1 * length2
            if not min_area or area < min_area:
                # noinspection PyUnusedLocal
                min_area = area
                best_origin = p + v1 * minx + v2 * miny
                best_v1 = v1
                best_v2 = v2
                best_length1 = length1
                best_length2 = length2
            # noinspection PyUnboundLocalVariable
            return best_origin, best_v1, best_v2, best_length1, best_length2

    def __init__(self, scene, direction):
        xs = [scene.boundbox.XMin, scene.boundbox.XMax]
        ys = [scene.boundbox.YMin, scene.boundbox.YMax]
        zs = [scene.boundbox.ZMin, scene.boundbox.ZMax]
        coords = itertools.product(xs, ys, zs)
        points = [Base.Vector(c) for c in coords]

        point_of_plane = (scene.boundbox.Center -
                          direction * 0.5 * scene.boundbox.DiagonalLength)
        projected_points = [p.projectToPlane(point_of_plane, direction)
                            for p in points]
        (self.origin, self.v1, self.v2, self.length1, self.length2) = (
            SunWindow.find_min_rectangle(projected_points, direction))
        self.aperture = self.length1 * self.length2
        self.main_direction = direction

    def random_point(self):
        return (self.origin + self.v1 * self.length1 * random.random() +
                self.v2 * self.length2 * random.random())

    def random_direction(self):
        return self.main_direction

    def add_to_document(self, doc):
        sw = Part.makePolygon([self.origin,
                               self.origin + self.v1 * self.length1,
                               self.origin + self.v1 * self.length1 +
                               self.v2 * self.length2,
                               self.origin + self.v2 * self.length2
                               ], True)
        doc.addObject("Part::Feature", "SunWindow").Shape = sw


class SunWindowBuie(SunWindow):
    def random_direction(self):
        return self.main_direction + random.random()

    def __init__(self, scene, direction, circum_solar_ratio):
        super(SunWindowBuie,self).__init__(scene, direction)
        self.length1 *= 1.2
        self.length2 *= 1.2
        self.origin = self.origin - 0.1 * self.v1 * self.length1 - 0.1 * self.v2 * self.length2
        		
		
class Ray:
    """
    Class used to model a sun ray. It keeps information of the path it 
    describes as it interacts with the scene and its energy.
    """

    def __init__(self, scene, origin, direction, properties):
        self.scene = scene
        self.points = [origin]
        self.directions = [direction]
        self.materials = [vacuum_medium]
        self.properties = properties
        self.wavelength = properties['wavelength']
        self.energy = properties['energy']
        self.finished = False
        self.got_absorbed = False
        self.current_medium = vacuum_medium

    def next_intersection(self):
        """
        Computes the next intersection of the ray with the scene. 
        Returns a pair (point,face), where point is the point of intersection
        and face is the face that intersects.
        If no intersection is found, it returns a "point at infinity" and None.
        """
        max_distance = 5 * self.scene.diameter
        p0 = self.points[-1]
        direction = self.directions[-1]
        p0 = p0 + direction * self.scene.epsilon
        intersections = []
        p1 = p0 + direction * max_distance
        segment = Part.Line(p0, p1)
        shape_segment = segment.toShape()
        for face in self.scene.faces:
            punts = face.section(shape_segment).Vertexes
            for punt in punts:
                intersections.append([punt.Point, face])
        if not intersections:
            return p1, None
        closestpair = min(intersections,
                          key=lambda (pair): p0.distanceToPoint(pair[0]))
        return tuple(closestpair)

    def next_direction(self, face):
        """
        Computes the next direction of the ray as it interacts with the scene,
        the material where the ray will be travelling next and
        a string representing the physical phenomenon that took place
        """
        current_direction = self.directions[-1]
        current_point = self.points[-1]
        current_material = self.materials[-1]
        uv = face.Surface.parameter(current_point)
        normal = face.normalAt(uv[0], uv[1])
        normal.normalize()
        if face in self.scene.materials:
            # face is active
            material = self.scene.materials[face]
            direction, phenomenon = material.change_of_direction(self, normal)
        else:
            # face is not active
            point_plus_delta = current_point + current_direction * self.scene.epsilon
            next_solid = self.scene.solid_at_point(point_plus_delta)
            nearby_material = self.scene.materials.get(next_solid, vacuum_medium)
            direction, phenomenon = nearby_material.change_of_direction(self, normal)
        next_material = None
        if phenomenon == 'Refraction':
            # noinspection PyUnboundLocalVariable
            next_material = nearby_material
        elif phenomenon == 'Reflexion':
            next_material = current_material
        elif phenomenon == 'Absortion':
            next_material = None
        return direction, next_material, phenomenon

    def update_energy(self):
        # TODO: @Ramon
        pass

    def run(self, max_hops=20):
        """
        Makes a sun ray propagate until it gets absorbed, it exits the scene,
        or gets caught in multiple (> max_hops) reflections.
        """
        count = 0
        while (not self.finished) and (count < max_hops):
            count += 1
            point, face = self.next_intersection()
            self.points.append(point)
            if not face:
                self.finished = True
                self.got_absorbed = False
                break
            vector, material, phenomenon = self.next_direction(face)
            self.update_energy()
            self.directions.append(vector)
            self.materials.append(material)
            if phenomenon == 'Absortion':
                self.finished = True
                self.got_absorbed = True

    def add_to_document(self, doc):
        lshape_wire = Part.makePolygon(self.points)
        my_shape_ray = doc.addObject("Part::Feature", "Ray")
        my_shape_ray.Shape = lshape_wire


class Experiment:
    """
    Sets up and runs and experiment in a given scene with a sun window in a 
    given direction.
    """

    def __init__(self, scene, direction, number_of_rays, wavelength, initial_energy,
                 show_in_doc=None):
        self.scene = scene
        self.direction = direction
        self.sunwindow = SunWindow(scene, direction)
        if show_in_doc:
            self.sunwindow.add_to_document(show_in_doc)
        self.number_of_rays = number_of_rays
        self.wavelength = wavelength
        self.initial_energy = initial_energy
        self.captured_energy = 0

    def run(self, show_in_doc=None):
        for _ in xrange(self.number_of_rays):
            ray = Ray(self.scene, self.sunwindow.random_point(), self.direction,
                      {'wavelength': self.wavelength, 'energy': self.initial_energy})
            ray.run()
            if show_in_doc:
                ray.add_to_document(show_in_doc)
            if ray.got_absorbed:
                self.captured_energy += ray.energy
