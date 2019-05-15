from autologging import traced
from .logging_unit import logger
from .materials import vacuum_medium, PVMaterial, PolarizedThinFilm
from .optics import Phenomenon, OpticalState
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
    got_absorbed : bool
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
        self.last_normal = None
        self.wavelength = wavelength
        self.energy = energy
        self.polarization_vectors = [polarization_vector]
        self.finished = False
        self.got_absorbed = False
        self.PV_values = []
        self.PV_absorbed = []

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
        closest_pair = min(intersections,
                           key=lambda pair: p0.distanceToPoint(pair[0]))
        return tuple(closest_pair)

    def next_state_and_normal(self, face):
        """
        Computes the next optical state after the ray hits the face, and the normal vector

        Parameters
        ----------
        face

        Returns
        -------
        state : self.OpticalState
            Optical State that the ray will have after hitting the face
        normal : Base.Vector
            Normal vector to the face where the ray hits
        """
        current_direction = self.current_direction()
        current_point = self.points[-1]
        uv = face.Surface.parameter(current_point)
        normal = face.normalAt(uv[0], uv[1])
        normal.normalize()
        if face in self.scene.materials:
            # face is active
            material = self.scene.materials[face]
            point_plus_delta = current_point + current_direction * self.scene.epsilon
            next_solid = self.scene.solid_at_point(point_plus_delta)
            nearby_material = self.scene.materials.get(next_solid, vacuum_medium)
            state = material.change_of_optical_state(self, normal, nearby_material)
            # TODO polarization_vector
        else:
            # face is not active
            point_plus_delta = current_point + current_direction * self.scene.epsilon
            next_solid = self.scene.solid_at_point(point_plus_delta)
            nearby_material = self.scene.materials.get(next_solid, vacuum_medium)
            if nearby_material.kind == 'Surface': # solid with surface material on all faces
                state = nearby_material.change_of_optical_state(self, normal, nearby_material)
            else:
                state = nearby_material.change_of_optical_state(self, normal)
            # TODO polarization_vector
        return state, normal

    def update_energy(self):
        material = self.current_medium()
        # TODO: @Ramon
        point_1 = self.points[-1]
        point_2 = self.points[-2]
        if 'extinction_coefficient' in material.properties:
            alpha = material.properties['extinction_coefficient'](self.wavelength) * \
                    4 * np.pi / (self.wavelength / 1E6) # mm-1
            d = point_1.distanceToPoint(point_2)
            self.energy = self.energy * np.exp(- alpha * d)
        if 'attenuation_coefficient' in material.properties:
            alpha = material.properties['attenuation_coefficient'](self.wavelength) # mm-1
            d = point_1.distanceToPoint(point_2)
            self.energy = self.energy * np.exp(- alpha * d)

    def run(self, max_hops=100):
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
            count += 1

            # Find next intersection
            # Update points
            point, face = self.next_intersection()
            self.points.append(point)
            if not face:
                self.finished = True
                self.got_absorbed = False
                break

            # Update optical state
            state, normal = self.next_state_and_normal(face)

            # Update energy
            energy_before = self.energy
            self.update_energy()
            # TODO: The following line is not elegant.
            #  Provisional solution for updating energy when passed a thin film
            if 'factor_energy_absorbed_thin_film' in state.extra_data:
                factor = state.extra_data['factor_energy_absorbed_thin_film']
                self.energy = self.energy * (1 - factor)

            # Treat PV material
            current_material = self.current_medium()
            if isinstance(current_material, PVMaterial):
                PV_absorbed_energy, PV_value = current_material.get_PV_data(self, energy_before)
                self.PV_values.append(PV_value)
                self.PV_absorbed.append(PV_absorbed_energy)

            # Update optical_states
            self.optical_states.append(state)
            self.last_normal = normal

            # Test for finish
            if self.energy < 1E-4:
                # needed for PV calculations TODO @Ramon: Why?
                self.finished = True
            if state.phenomenon == Phenomenon.ABSORPTION:
                self.finished = True
            if state.phenomenon == Phenomenon.GOT_ABSORBED:
                self.got_absorbed = True
                self.finished = True

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

