from autologging import traced
from .logging_unit import logger
from .materials import vacuum_medium
from .optics import Phenomenon
import numpy as np
import Part

@traced(logger)
class Ray(object):
    """
    Class used to model a sun ray

    It keeps information of the path it describes and its polarization as it interacts
    with the scene and its energy.

    Parameters
    ----------
    scene : otsun.Scene
        Scene where the ray moves
    origin : Base.Vector
        Point of emission of the ray
    direction : Base.Vector
        Direction of the ray when it was emitted
    properties : dict
        Dictionary with properties of the ray. Mandatory:
        'wavelength' : float
            Wavelength of the ray
        'energy' : float
            Initial enerty of the ray
        'polarization_vector' : Base.Vector
            Initial polarization vector of the ray
    #TODO: Make this properties parameters of the class (not a dict)

    Attributes
    ----------
    points : list of Base.Vector
        Points where the ray passes at each iteration
    directions : list of Base.Vector
        Directions of the ray at each iteration
    normals : list of Base.Vector
        Vectors normal to the surfaces where the ray hits (except the first)
        # TODO: Review... is it needed? strange!
    materials : list of otsun.Material
        List of materials where the ray is located at each iteration
    wavelength : float
        Wavelength of ray
    energy : float
        Current energy of ray
    polarization_vectors : list of Base.Vector
        Polarization of the ray at each iteration
    finished : bool
    got_absorbed : bool
    current_medium : otsun.Material
        Current material where the ray is
        # TODO: Find relation between current_medium and materials[-1] (!?)
    PV_energy : float
    PV_values : list of ?
    in_PV : bool
    PV_absorbed : list of ?
    """

    def __init__(self, scene, origin, direction, properties):
        self.scene = scene
        self.points = [origin]
        self.directions = [direction]
        self.normals = []
        self.materials = [vacuum_medium]
        self.properties = properties
        self.wavelength = properties['wavelength']
        self.energy = properties['energy']
        self.polarization_vectors = [properties['polarization_vector']]
        self.finished = False
        self.got_absorbed = False
        self.current_medium = vacuum_medium
        self.PV_energy = 0.0
        self.PV_values = []
        self.in_PV = False
        self.PV_absorbed = []

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

    def next_state_and_material(self, face):
        """
        Computes the next direction of the ray as it interacts with the scene,
        the material where the ray will be travelling next and
        a string representing the physical phenomenon that took place
        """
        current_direction = self.directions[-1]
        current_point = self.points[-1]
        current_material = self.materials[-1]
        polarization_vector = self.polarization_vectors[-1]
        uv = face.Surface.parameter(current_point)
        normal = face.normalAt(uv[0], uv[1])
        normal.normalize()
        if face in self.scene.materials:
            # face is active
            material = self.scene.materials[face]
            point_plus_delta = current_point + current_direction * self.scene.epsilon
            next_solid = self.scene.solid_at_point(point_plus_delta)
            nearby_material = self.scene.materials.get(next_solid, vacuum_medium)
            state = material.change_of_direction(self, normal, nearby_material)
            ### polarization_vector, direction, phenomenon = results[0], results[1], results[2]
            # TODO polarization_vector
        else:
            # face is not active
            point_plus_delta = current_point + current_direction * self.scene.epsilon
            next_solid = self.scene.solid_at_point(point_plus_delta)
            nearby_material = self.scene.materials.get(next_solid, vacuum_medium)
            if nearby_material.kind == 'Surface': # solid with surface material on all faces
                state = nearby_material.change_of_direction(self, normal, nearby_material)
            else:
                state = nearby_material.change_of_direction(self, normal)
            ### polarization_vector, direction, phenomenon = nearby_material.change_of_direction(self,
            ###                                                                                 normal)
            # TODO polarization_vector
        next_material = None
        if state.phenomenon == Phenomenon.REFRACTION:
            next_material = nearby_material
        elif state.phenomenon == Phenomenon.REFLEXION:
            next_material = current_material
        elif state.phenomenon == Phenomenon.ABSORPTION:
            next_material = None
        elif state.phenomenon == Phenomenon.GOT_ABSORBED:  # it is needed? Review
            next_material = None
        return state, next_material, normal
        ### polarization_vector, direction, next_material, phenomenon

    def update_energy(self, material):
        # TODO: @Ramon
        point_1 = self.points[-1]
        point_2 = self.points[-2]
        if 'extinction_coefficient' in material.properties:
            alpha = material.properties['extinction_coefficient'](self.wavelength) * 4 * np.pi / (self.wavelength / 1E6) # mm-1
            d = point_1.distanceToPoint(point_2)
            self.energy = self.energy * np.exp(- alpha * d)
        if 'attenuation_coefficient' in material.properties:
            alpha = material.properties['attenuation_coefficient'](self.wavelength) # mm-1
            d = point_1.distanceToPoint(point_2)
            self.energy = self.energy * np.exp(- alpha * d)

    def run(self, max_hops=100):
        """
        Makes a sun ray propagate until it gets absorbed, it exits the scene,
        or gets caught in multiple (> max_hops) reflections.
        """
        count = 0
        while (not self.finished) and (count < max_hops):
            assert self.current_medium == self.materials[-1]
            count += 1
            point, face = self.next_intersection()
            self.points.append(point)
            if not face:
                self.finished = True
                self.got_absorbed = False
                break
            energy_before = self.energy
            point_1 = self.points[-1]
            point_2 = self.points[-2]
            middle_point = point_1.add(point_2) * 0.5
            actual_solid = self.scene.solid_at_point(middle_point)
            if actual_solid:
                current_material = self.scene.materials[actual_solid]
                self.update_energy(current_material)
                if 'PV_material' in self.scene.materials[actual_solid].properties:
                    self.in_PV = True
                    self.PV_energy = energy_before
                    alpha = self.scene.materials[actual_solid].properties['extinction_coefficient'](
                    self.wavelength) * 4 * np.pi / (self.wavelength / 1E6) # mm-1
                    angle_incident = np.arccos (- self.normals[-1].dot(self.directions[-1])) * 180.0 / np.pi
                    self.PV_values.append((point_2.x,point_2.y,point_2.z,point_1.x,point_1.y,point_1.z,energy_before,self.energy,self.wavelength,alpha,angle_incident))
                    self.PV_absorbed.append(energy_before - self.energy)
            state, material, normal = self.next_state_and_material(face)  # TODO polarization_vector
            self.directions.append(state.direction)
            self.materials.append(material)
            self.normals.append(normal)
            self.polarization_vectors.append(state.polarization)
            if self.energy < 1E-4: # needed for PV calculations
                self.finished = True
            if state.phenomenon == Phenomenon.ABSORPTION:
                self.finished = True
            if state.phenomenon == Phenomenon.GOT_ABSORBED:
                self.got_absorbed = True
                self.finished = True

    def add_to_document(self, doc):
        lshape_wire = Part.makePolygon(self.points)
        my_shape_ray = doc.addObject("Part::Feature", "Ray")
        my_shape_ray.Shape = lshape_wire
#        lshape_wire = Part.makePolygon(self.points)
#        doc.addObject("Part::Feature", "Ray").Shape = lshape_wire


