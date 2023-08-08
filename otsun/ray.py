"""Module otsun.ray for the modelization of light rays

The module defines the class `Ray`
"""

from autologging import traced
from .logging_unit import logger
from .materials import vacuum_medium, PVMaterial, SurfaceMaterial, TwoLayerMaterial, PolarizedThinFilm
from .optics import Phenomenon, OpticalState
import numpy as np
import Part
from FreeCAD import Base

try:
    part_line = Part.LineSegment
except AttributeError:
    part_line = Part.Line

# Zero energy level
LOW_ENERGY = 1E-6


def _center(bb):
    return Base.Vector((bb.XMin+bb.XMax)/2, (bb.YMin+bb.YMax)/2, (bb.ZMin+bb.ZMax)/2)


def _interval_intersects(x1, x2, y1, y2):
    return x1 <= y2 and y1 <= x2


def _bb_intersects(bb1, bb2):
    return (_interval_intersects(bb1.XMin, bb1.XMax, bb2.XMin, bb2.XMax) and
            _interval_intersects(bb1.YMin, bb1.YMax, bb2.YMin, bb2.YMax) and
            _interval_intersects(bb1.ZMin, bb1.ZMax, bb2.ZMin, bb2.ZMax))


def _distance_point_to_line(point, line_point, line_vector):
    v = point-line_point
    cross = v.cross(line_vector)
    return cross.Length


def _distance_point_to_ray(point, ray_origin, ray_vector):
    v = point - ray_origin
    if v.dot(ray_vector) >= 0:
        return _distance_point_to_line(point, ray_origin, ray_vector)
    return v.Length


def _line_may_intersect_bb(bb, line_point, line_vector):
    center = _center(bb)
    d = _distance_point_to_line(center, line_point, line_vector)
    return d <= bb.DiagonalLength/2


def _ray_may_intersect_bb(bb, ray_origin, ray_vector):
    center = _center(bb)
    d = _distance_point_to_ray(center, ray_origin, ray_vector)
    return d <= bb.DiagonalLength/2


@traced(logger)
class Ray(object):
    """
    Class used to model a sun ray

    It keeps information of the path it describes and its polarization as it interacts
    with the scene and its energy.

    Parameters
    ----------
    scene : Scene
        Scene where the ray moves
    origin : Base.Vector
        Point of emission of the ray
    direction : Base.Vector
        Direction of the ray when it was emitted
    wavelength : float
        Wavelength of the ray
    energy : float
        Initial energy of the ray
    polarization_vector : Base.Vector
        Initial polarization vector of the ray

    Attributes
    ----------
    points : list of Base.Vector
        Points where the ray passes at each iteration
    optical_states : list of OpticalState
        OpticalStates of the rays at each iteration
    last_normal : Base.Vector
        Last vector normal to the surface where the ray hits (used for PV)
    wavelength : float
        Wavelength of ray
    energy : float
        Current energy of ray
    polarization_vectors : list of Base.Vector
        Polarization of the ray at each iteration
    finished : bool
    Th_absorbed : bool
    PV_values : list of tuple of floats
        List where each entry is a PV_value (a tuple of floats)
    PV_absorbed : list of float
        List of values of absorbed PV energy
    """

    def __init__(self, scene, origin, direction,
                 wavelength, energy, polarization_vector):
        self.scene = scene
        self.points = [origin]
        state = OpticalState(
            polarization_vector,
            direction,
            Phenomenon.JUST_STARTED,
            vacuum_medium
        )
        self.optical_states = [state]
        self.current_solid = None
        self.last_normal = None
        self.last_touched_face = None
        self.wavelength = wavelength
        self.energy = energy
        self.polarization_vectors = [polarization_vector]
        self.finished = False
        self.Th_absorbed = False
        self.PV_values = []
        self.PV_absorbed = []

    def __str__(self):
        return "Pos.: %s, OS: %s, Energy: %s" % (
            self.points[-1], self.optical_states[-1], self.energy
        )

    def current_medium(self):
        """
        Get medium where ray is currently traveling

        Returns
        -------
            otsun.VolumeMaterial
        """
        return self.optical_states[-1].material

    def current_direction(self):
        """
        Get current direction

        Returns
        -------
            Base.Vector
        """
        return self.optical_states[-1].direction

    def current_polarization(self):
        """
        Get current polarization

        Returns
        -------
            Base.Vector
        """
        return self.optical_states[-1].polarization

    def weighted_distance(self, p0, pair):
        if pair[1] in self.scene.materials:
            return p0.distanceToPoint(pair[0])
        else:
            return p0.distanceToPoint(pair[0]) + self.scene.epsilon

    def next_intersection(self):
        """
        Finds next intersection of the ray with the scene

        Returns a pair (point,face), where point is the point of intersection
        and face is the face that intersects.
        If no intersection is found, it returns a "point at infinity" and None.

        Returns
        -------
        point : Base.Vector
            Point of intersection
        face : Part.Face or None
            Face where it intersects
        """
        max_distance = 5 * self.scene.diameter
        p0 = self.points[-1]
        direction = self.current_direction()
        # p0 = p0 + direction * self.scene.epsilon
        intersections = []
        p1 = p0 + direction * max_distance
        segment = part_line(p0, p1)
        shape_segment = segment.toShape()
        bb1 = shape_segment.BoundBox
        candidates = self.scene.faces
        feasible_faces = 0
        filtered_faces_1 = 0
        filtered_faces_2 = 0
        feasible_but_empty = 0
        for face in candidates:
            bb2 = face.BoundBox
            if not _ray_may_intersect_bb(bb2, p0, direction):
                filtered_faces_1 += 1
                continue
            if not _bb_intersects(bb1, bb2):
                filtered_faces_2 += 1
                continue
            feasible_faces += 1
            punts = face.section(shape_segment).Vertexes
            # if face == self.last_touched_face:
            #     punts = [punt for punt in punts if p0.distanceToPoint(punt.Point) > self.scene.epsilon]
            if not punts:
                feasible_but_empty += 1
                continue
                # logger.debug("feasible face but empty intersection")
            # material = self.scene.materials[face]
            if (face in self.scene.materials) & (face != self.last_touched_face):
                for punt in punts:
                    intersections.append([punt.Point, face])
            else:
                for punt in punts:
                    if p0.distanceToPoint(punt.Point) > self.scene.epsilon:
                        intersections.append([punt.Point, face])

        logger.debug(f"Found {len(intersections)} points in {feasible_faces} faces. Filtered {filtered_faces_1}+{filtered_faces_2}. Feasible but empty {feasible_but_empty}")
#        intersections = [punt_cara for punt_cara in intersections if
#                         p0.distanceToPoint(punt_cara[0]) > self.scene.epsilon]
        if not intersections:
            logger.debug("No true intersection found with scene")
            return p1, None
        closest_pair = min(intersections,
                           key=lambda pair: self.weighted_distance(p0, pair))
        return tuple(closest_pair)

    def next_state_solid_and_normal(self, face):
        """
        Computes the next optical state after the ray hits the face, and the normal vector

        Parameters
        ----------
        face

        Returns
        -------
        state : OpticalState
            Optical State that the ray will have after hitting the face
        next_solid : Part.Solid
            Solid where the ray will travel after hitting the face
        normal : Base.Vector
            Normal vector to the face where the ray hits
        """
        current_direction = self.current_direction()
        current_point = self.points[-1]
        uv = face.Surface.parameter(current_point)
        normal = face.normalAt(uv[0], uv[1])
        normal.normalize()
        nearby_solid = self.scene.next_solid_at_point_in_direction(current_point, normal, current_direction)
        nearby_material = self.scene.materials.get(nearby_solid, vacuum_medium)
        if face in self.scene.materials:
            # face is active
            material = self.scene.materials[face]
            # point_plus_delta = current_point + current_direction * 2 * self.scene.epsilon
            state = material.change_of_optical_state(self, normal, nearby_material)
            # TODO polarization_vector
        else:
            # face is not active
            # point_plus_delta = current_point + current_direction * 2 * self.scene.epsilon
            # next_solid = self.scene.solid_at_point(point_plus_delta)
            # nearby_material = self.scene.materials.get(next_solid, vacuum_medium)
            if isinstance(nearby_material, (SurfaceMaterial, TwoLayerMaterial)):
                # solid with surface material on all faces
                state = nearby_material.change_of_optical_state(self, normal, nearby_material)
            else:
                state = nearby_material.change_of_optical_state(self, normal)
            # TODO polarization_vector
        logger.debug(state)
        next_solid = None
        if state.phenomenon in [Phenomenon.ENERGY_ABSORBED, Phenomenon.ABSORPTION]:
            pass
        elif state.phenomenon == Phenomenon.REFLEXION:
            next_solid = self.current_solid
        elif state.phenomenon in [Phenomenon.REFRACTION, Phenomenon.TRANSMITTANCE]:
            next_solid = nearby_solid
        return state, next_solid, normal

    def update_energy(self):
        material = self.current_medium()
        if (material.name == "Vacuum" and
                'vacuum_material' in self.scene.extra_data and
                len(self.points) > 2
        ):
            # implementation of dispersive vacuum after first interaction with scene
            material = self.scene.extra_data['vacuum_material']
        # TODO: @Ramon
        point_1 = self.points[-1]
        point_2 = self.points[-2]
        if material.properties.get('extinction_coefficient', None):
            alpha = material.properties['extinction_coefficient'](self.wavelength) * \
                    4 * np.pi / (self.wavelength / 1E6)  # mm-1
            if alpha:
                d = point_1.distanceToPoint(point_2)
                self.energy = self.energy * np.exp(- alpha * d)
        if material.properties.get('attenuation_coefficient', None):
            alpha = material.properties['attenuation_coefficient'](self.wavelength)  # mm-1
            if alpha:
                d = point_1.distanceToPoint(point_2)
                self.energy = self.energy * np.exp(- alpha * d)

    def run(self, max_hops=200):
        """
        Makes the ray propagate

        Makes the sun ray propagate until it gets absorbed, it exits the scene,
        or gets caught in multiple (> max_hops) reflections.

        Parameters
        ----------
        max_hops : int
            Maximum number of iterations
        """
        count = 0
        while (not self.finished) and (count < max_hops):
            logger.debug("Ray running. Hop %s, %s, Solid %s", count, self,
                        self.scene.name_of_solid.get(self.current_solid, "Void"))
            count += 1

            # Find next intersection
            # Update points
            point, face = self.next_intersection()
            self.last_touched_face = face
            self.points.append(point)
            if not face:
                self.finished = True
                self.Th_absorbed = False
                break

            # Update energy
            energy_before = self.energy
            self.update_energy()

            # Treat PV material
            current_material = self.current_medium()
            if isinstance(current_material, PVMaterial):
                PV_absorbed_energy, PV_value = current_material.get_PV_data(self, energy_before)
                self.PV_values.append(PV_value)
                self.PV_absorbed.append(PV_absorbed_energy)

            # Update optical state
            state, next_solid, normal = self.next_state_solid_and_normal(face)

            # TODO: The following line is not elegant.
            #  Provisional solution for updating energy when passed a thin film
            if 'factor_energy_absorbed' in state.extra_data:
                factor = state.extra_data['factor_energy_absorbed']
                self.energy = self.energy * (1 - factor)

            # Update optical_states
            self.optical_states.append(state)
            self.last_normal = normal
            self.current_solid = next_solid

            # Test for finish
            if state.phenomenon == Phenomenon.ABSORPTION:
                self.finished = True
            if state.phenomenon == Phenomenon.ENERGY_ABSORBED:
                self.Th_absorbed = True
                self.finished = True
            if self.energy < LOW_ENERGY:
                self.finished = True
        logger.debug("Ray stopped. Hop %s, %s, Solid %s", count, self,
                    self.scene.name_of_solid.get(self.current_solid, "Void"))

    def add_to_document(self, doc):
        """
        Draws the ray in a FreeCAD document

        Parameters
        ----------
        doc : App.Document
        """
        lshape_wire = Part.makePolygon(self.points)
        my_shape_ray = doc.addObject("Part::Feature", "Ray")
        my_shape_ray.Shape = lshape_wire
