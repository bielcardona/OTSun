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

try:
    import xml.etree.cElementTree as ET
except ImportError:
    import xml.etree.ElementTree as ET


# ---
# Helper functions
#
def reflexion(incident, normal):
    """
    Implementation of law of reflexion
    """
    # We assume all vectors are normalized
    return incident - normal * 2.0 * normal.dot(incident)

def refraction(incident, normal, n1, n2):
    """
    Implementation of Snell's law of refraction
    """
    # https://en.wikipedia.org/wiki/Snell's_law#Vector_form
    mynormal = normal * 1.0
    if mynormal.dot(incident) > 0:
        mynormal = mynormal * (-1.0)
    r = n1 / n2
    c1 = - mynormal.dot(incident)
    c2sq = 1.0 - r * r * (1.0 - c1 * c1)
    if c2sq < 0:
        return reflexion(incident, normal)
    c2 = c2sq ** 0.5
    return incident * r + mynormal * (r * c1 - c2)



class Material:
    """
    Class used to represent materials and their physical properties
    """
    by_name = {}

    def __init__(self, name, parameters):
        self.by_name[name] = self
        self.name = name
        self.parameters = parameters

    @classmethod
    def get_from_label(cls, label):
        if ("(" not in label) or (")" not in label):
            return None
        start = label.find("(")
        end = label.find(")")
        name = label[start + 1:end]
        return cls.by_name.get(name, None)

    def change_of_direction(self, ray, normal_vector):
        pass
        # Compute new direction (ray.current_material)


class VolumeMaterial(Material):
    def __init__(self, name, parameters):
        super(VolumeMaterial,self).__init__(name, parameters)
        pass

    def change_of_direction(self, ray, normal_vector):
        wavelength = ray.wavelength

class SimpleVolumeMaterial(VolumeMaterial):
    def __init__(self, name, parameters):
        super(SimpleVolumeMaterial,self).__init__(name, parameters)
        self.index_of_refraction = index_of_refraction

    def change_of_direction(self, ray, normal_vector):
        pass

glass1 = SimpleVolumeMaterial(1.7)


class SurfaceMaterial(Material):
    def __init__(self, name, parameters):
        super(SurfaceMaterial,self).__init__(name, parameters)
        pass

    def scatter_direction(self, ray, direction):
        # pensar!
        pass


class SimpleSurfaceMaterial(SurfaceMaterial):
    def __init__(self, name, parameters):
        super(SimpleSurfaceMaterial,self).__init__(name, parameters)
        self.probability_of_reflexion = probability_of_reflexion
        self.probability_of_absortion = probability_of_absortion
        self.transmitance = 1 - probability_of_absortion - probability_of_reflexion

    def change_of_direction(self, ray, normal_vector):
        pass

def opaque_simple_material(por):
    return SimpleSurfaceMaterial(por,1-por)

def transparent_simple_material(por):
    return SimpleSurfaceMaterial(por,0)

perfect_mirror = opaque_simple_material(1)
perfect_absorber = opaque_simple_material(0)

class TabulatedSurfaceMaterial(SurfaceMaterial):
    def __init__(self, reflactances_file):
        super(TabulatedSurfaceMaterial,self).__init__()
        self.reflactances = None # carregar dades del fitxer

    def get_reflactance(self,wavelength):
        pass

class ClosedFormSurfaceMaterial(SurfaceMaterial):
    def __init__(self, function_of_reflactance):
        super(ClosedFormSurfaceMaterial,self).__init__()
        self.function_of_reflactance = function_of_reflactance

    def get_reflactance(self,wavelength):
        return self.function_of_reflactance(wavelength)

def function1(wavelength):
    return math.exp(wavelength/1000)

strangemat1 = ClosedFormSurfaceMaterial(function1)



vacuum_medium = Material("Vacuum", "Solid",
                         {'index_of_refraction': 1.0,
                          'probability_of_reflexion': 0.0})


class Scene:
    """
    Class used to define the Scene. It encodes all the objects 
    that interfere with light rays. 
    """

    def __init__(self, objects):
        self.faces = []  # All the faces in the Scene
        self.solids = []  # All the solids in the Scene
        self.materials = {}  # Asign the materials to objects
        # self.sum_energy = 0.0 # Energy of the Scene
        self.epsilon = 0.00001  # Tolerance for solid containment
        self.boundbox = None

        for obj in objects:
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
                min_area = area
                best_origin = p + v1 * minx + v2 * miny
                best_v1 = v1
                best_v2 = v2
                best_length1 = length1
                best_length2 = length2
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
        return self.main_direction + random()

class Ray:
    """
    Class used to model a sun ray. It keeps information of the path it 
    describes as it interacts with the scene and its energy.
    """

    def __init__(self, wavelength, origin, direction, energy, scene):
        self.scene = scene
        self.wavelength = wavelength
        self.points = [origin]
        self.directions = [direction]
        self.energy = energy
        self.finished = False
        self.got_absorbed = False

    def next_intersection(self):
        """
        Computes the next intersection of the ray with the scene. 
        Returns a pair (point,face), where point is the point of intersection
        and face is the face that intersects.
        If no intersection is found, it returns None.
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
            self.points.append(p1)
            self.directions.append(direction)
            self.finished = True
            self.got_absorbed = False
            return None
        closestpair = min(intersections,
                          key=lambda (pair): p0.distanceToPoint(pair[0]))
        closest = closestpair[0]
        self.points.append(closest)

        middle_point = (p0 + closest) * 0.5
        dist = p0.distanceToPoint(closest)
        solid = self.scene.solid_at_point(middle_point)
        if solid:
            material = self.scene.materials[solid]
            exco = material.parameters['extinction_coefficient']
            self.energy = self.energy * np.exp(exco * dist)
            print self.energy, exco, dist
        return closestpair


    def next_direction(self, pair):
        """
        Computes the next direction of the ray as it interacts with the scene. 
        The direction is stored in the data of the ray.
        """
        p1 = self.points[-1]
        v0 = self.directions[-1]
        if not pair:
            self.directions.append(v0)
            return None
        face = pair[1]
        uv = face.Surface.parameter(p1)
        normal = face.normalAt(uv[0], uv[1])
        normal.normalize()
        if face in self.scene.materials:
            material = self.scene.materials[face]
            if material.kind == "Antireflectant":
                pass
            if material.kind == "Mirror":
                print "Mirror"
                por = material.parameters["probability_of_reflection"]
                if random.random() < por:
                    self.directions.append(reflexion(v0, normal))
                else:
                    self.energy = 0
                    self.finished = True
                    self.got_absorbed = False
            if material.kind == "Absorber":
                poa = material.parameters.get("probability_of_absortion", 1)
                if random.random() < poa:
                    self.finished = True
                    self.got_absorbed = True
                else:
                    self.energy = 0
                    self.finished = True
                    self.got_absorbed = False
        else:
            print "Refraction"
            sol1 = self.scene.solid_at_point(p1 - v0 * self.scene.epsilon)
            sol2 = self.scene.solid_at_point(p1 + v0 * self.scene.epsilon)
            mat1 = self.scene.materials.get(sol1, vacuum_medium)
            mat2 = self.scene.materials.get(sol2, vacuum_medium)
            por = mat2.parameters.get('probability_of_reflection', 0)
            rnd = random.random()
            if rnd < por:
                self.directions.append(reflexion(v0, normal))
            else:
                n1 = mat1.parameters.get('index_of_refraction', 1)
                n2 = mat2.parameters.get('index_of_refraction', 1)
                self.directions.append(refraction(v0, normal, n1, n2))

    def run(self, max_hops=20):
        """
        Makes a sun ray propagate until it gets absorbed, it exits the scene,
        or gets caught in multiple (> max_hops) reflections.
        """
        count = 0
        while not self.finished:
            pair = self.next_intersection()
            if not pair: break
            self.next_direction(pair)
            count += 1
            if count > max_hops: break

    def add_to_document(self, doc):
        lshape_wire = Part.makePolygon(self.points)
        MyShape_ray = doc.addObject("Part::Feature", "Ray")
        MyShape_ray.Shape = lshape_wire


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


class Experiment:
    """
    Sets up and runs and experiment in a given scene with a sun window in a 
    given direction.
    """

    def __init__(self, scene, direction, number_of_rays, wavelength, initial_energy,
                 show_in_doc=None):
        self.scene = scene
        ## self.direction = direction
        self.sunwindow = SunWindow(scene, direction)
        if show_in_doc:
            self.sunwindow.add_to_document(show_in_doc)
        self.number_of_rays = number_of_rays
        self.wavelength = wavelength
        self.initial_energy = initial_energy
        self.captured_energy = 0

    def run(self, show_in_doc=None):
        for _ in xrange(self.number_of_rays):
            ray = Ray(self.sunwindow.random_point(), self.direction,
                      self.initial_energy, self.scene)
            ray.run()
            if show_in_doc:
                ray.add_to_document(show_in_doc)
            if ray.got_absorbed:
                self.captured_energy += ray.energy


if __name__ == "__main__":
    import FreeCAD
    from os.path import expanduser, join

    home = expanduser("~")
    filename = join(home, "materials.xml")
    doc = FreeCAD.activeDocument()
    sel = doc.Objects
    Material.load_from_file(filename)
    scene = Scene(sel)
    exp = Experiment(scene, Base.Vector(1, 1, 1), 100, 1.0 ,1.0, doc)
    exp.run(doc)
