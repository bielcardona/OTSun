import itertools
import Part
from .materials import vacuum_medium
from .optics import *
from .math import pick_random_from_CDF
from .logging_unit import logger
from autologging import traced


class SunWindow(object):
    """
    Class used to define a sun window, defined as a rectangle in space such
    that the projection of the scene on the plane is enclosed in the rectangle,
    the rectangle has an increase of 10% in each length
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
            Length1 = best_length1 * 1.04
            Length2 = best_length2 * 1.04
            origin = best_origin - best_v1 * Length1 * 0.02 - best_v2 * Length2 * 0.02
            return origin, best_v1, best_v2, Length1, Length2

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
        return (self.origin + self.v1 * self.length1 * myrandom() +
                self.v2 * self.length2 * myrandom())

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


class BuieDistribution:
    """
    A distribution for the direct light from the sun according to the Buie model:
    Buie, D., 2005. Corrigendum to "The effective size of the solar cone for
    solar concentrating systems" [Solar Energy 74 (2003) 417-427].
    Solar Energy 79, 568-570. 5.2
    """

    def __init__(self, CircumSolarRatio):
        self.CSR = CircumSolarRatio
        self.SD = 4.65  # Solar Disk in mrad
        self.SS = 43.6  # Solar Size in mrad
        self.aa1 = BuieDistribution.a1(self.CSR, self.SD)  # normalization constant for the disk region
        self.aa2 = BuieDistribution.a2(self.CSR, self.SD, self.SS)  # normalization constant for the circumsolar region
        self.CDF_Disk_Region = BuieDistribution.CDF_disk_region(self.aa1,
                                                                self.SD)  # Buie distribution for the disk region

    def angle_distribution(self, random_number):
        u = random_number
        if u < 1.0 - self.CSR:
            th_u = BuieDistribution.th_solar_disk_region(u, self.aa1, self.SD) / 1000.0 * 180.0 / np.pi
        else:
            th_u = BuieDistribution.th_circumsolar_region(u, self.CSR, self.SD, self.SS,
                                                          self.aa2) / 1000.0 * 180.0 / np.pi
        return th_u

    @staticmethod
    def solar_disk__density_distribution(angle):
        return 2.0 * np.pi * np.cos(0.326 * angle) / np.cos(0.305 * angle) * angle

    @staticmethod
    def circumsolar__density_distribution(angle, CSR):
        gamma = 2.2 * np.log(0.52 * CSR) * CSR ** (0.43) - 0.1
        kappa = 0.9 * np.log(13.5 * CSR) * CSR ** (-0.3)
        return 2.0 * np.pi * np.exp(kappa) * angle ** (gamma + 1.0)

    @staticmethod
    def a1(CSR, SD):
        """ Parameter a1 needed for the normalization of the probability distribution in the disk region"""
        f = BuieDistribution.solar_disk__density_distribution
        th = np.arange(0.0, SD, 0.001)
        f_th = np.vectorize(f)
        aa1 = (1.0 - CSR) / np.trapz(f_th(th), dx=0.001)
        return aa1

    @staticmethod
    def a2(CSR, SD, SS):
        """ Parameter a2 needed for the normalization of the probability distribution in the circumsolar region"""
        f = BuieDistribution.circumsolar__density_distribution
        th = np.arange(SD, SS, 0.001)
        f_th = np.vectorize(f)
        aa2 = CSR / np.trapz(f_th(th, CSR), dx=0.001)
        return aa2

    @staticmethod
    def CDF_disk_region(a1, SD):
        """ Cumulative Distribution Function in the solar disk region"""
        f = BuieDistribution.solar_disk__density_distribution
        th = np.arange(0.0, SD, 0.001)
        f_th = np.vectorize(f)
        cumulative = np.cumsum(f_th(th))
        CDF = (th, a1 * cumulative / 1000.0)
        return CDF

    @staticmethod
    def th_solar_disk_region(u, a1, SD):
        """ Random angle based on the probability distribution in the circumsolar region"""
        CDF = BuieDistribution.CDF_disk_region(a1, SD)
        idx = (np.abs(CDF[1] - u)).argmin()
        return CDF[0][idx]

    @staticmethod
    def th_circumsolar_region(u, CSR, SD, SS, aa2):
        """ Random angle based on the CDF in the circumsolar region"""
        gamma = 2.2 * np.log(0.52 * CSR) * CSR ** (0.43) - 0.1
        kappa = 0.9 * np.log(13.5 * CSR) * CSR ** (-0.3)
        u_csr = (u - 1.0) + CSR  # Since CSR-CDF starts at zero
        f1 = u_csr * (gamma + 2.0) / (aa2 * 2 * np.pi * np.exp(kappa))
        f2 = SD ** (gamma + 2.0)
        exp = (1.0 / (gamma + 2.0))
        th_u = np.power(f1 + f2, exp)
        return th_u


@traced(logger)
class Ray:
    """
    Class used to model a sun ray and its polarization vector. It keeps information of the path it
    describes as it interacts with the scene and its energy.
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
            # noinspection PyUnboundLocalVariable
            next_material = nearby_material
        elif state.phenomenon == Phenomenon.REFLEXION:
            next_material = current_material
        elif state.phenomenon == Phenomenon.ABSORTION:
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
            if state.phenomenon == Phenomenon.ABSORTION:
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


class LightSource:
    """
    Sets up a light source with a given scene, a given emitting region and a given light spectrum.
    The emitting region must provide the main direction.
    Light spectrum could be: a constant value (for a single wavelength in nanometers), or a spectrum distribution.
    The distribution (dispersion) for the main direction is provided in "direction_distribution".
    The polarization_vector is a Base.Vector for polarized light. If is not given unpolarized light is generated.
    """

    def __init__(self, scene, emitting_region, light_spectrum, initial_energy, direction_distribution=None,
                 polarization_vector=None):
        self.scene = scene
        self.emitting_region = emitting_region
        self.light_spectrum = light_spectrum
        self.initial_energy = initial_energy
        self.direction_distribution = direction_distribution
        self.polarization_vector = polarization_vector
        self.wavelengths = []

    def emit_ray(self):
        point = self.emitting_region.random_point()
        main_direction = self.emitting_region.main_direction  # emitting main direction
        direction = main_direction
        if self.direction_distribution is not None:  # main direction has a distribution
            theta = self.direction_distribution.angle_distribution(myrandom())
            phi = 360.0 * myrandom()
            direction = dispersion_from_main_direction(main_direction, theta, phi)
            if self.polarization_vector:  # single polarization vector is active
                polarization_vector = dispersion_polarization(main_direction, self.polarization_vector, theta, phi)
        if self.polarization_vector is None:  # unpolarization is active
            polarization_vector = random_polarization(direction)  # random polarization from light direction
        if np.isscalar(self.light_spectrum):
            wavelength = self.light_spectrum  # experiment with a single wavelength (nanometers)
        else:
            wavelength = pick_random_from_CDF(self.light_spectrum)  # light spectrum is active (nanometers)
        ray = Ray(self.scene, point, direction, {'wavelength': wavelength,
                                                 'energy': self.initial_energy,
                                                 'polarization_vector': polarization_vector})
        return ray

