"""
Module that implements rays and its sources
"""

import itertools
import Part
from .materials import vacuum_medium
from .ray import Ray
from .optics import *
from .math import pick_random_from_cdf
from .logging_unit import logger
from autologging import traced


class SunWindow(object):
    """
    Class that implements a Sun window (rectangle that emits rays)

    Source of light defined by a rectangle "at infinity" in space that emits rays perpendicular to it.

    Parameters
    ----------
    scene : Scene
        Scene that contains the sun window
    main_direction : Base.Vector
        Vector orthogonal to the emitting region

    Attributes
    ----------
    origin : Base.Vector
        Center of the rectangle
    v1, v2 : Base.Vector
        Unit vectors parallel to the sides of the rectangle
    length1, length2 : float
        Lengths of the sides of the rectangle
    aperture : float
        Area of the rectangle
    """

    def __init__(self, scene, main_direction):
        xs = [scene.boundbox.XMin, scene.boundbox.XMax]
        ys = [scene.boundbox.YMin, scene.boundbox.YMax]
        zs = [scene.boundbox.ZMin, scene.boundbox.ZMax]
        coords = itertools.product(xs, ys, zs)
        points = [Base.Vector(c) for c in coords]

        point_of_plane = (scene.boundbox.Center -
                          main_direction * 0.5 * scene.boundbox.DiagonalLength)
        projected_points = [p.projectToPlane(point_of_plane, main_direction)
                            for p in points]
        (self.origin, self.v1, self.v2, self.length1, self.length2) = (
            SunWindow.find_min_rectangle(projected_points, main_direction))
        self.aperture = self.length1 * self.length2
        self.main_direction = main_direction

    @staticmethod
    def find_min_rectangle(points, normal):
        """
        Computes the minimum rectangle covering points in a direction

        Given a list of `points`, take its projection in a `normal` direction,
        and the rectangle with minimum area that encloses this projections

        Parameters
        ----------
        points : list of Base.Vector
            List of points to enclose
        normal : Base.Vector
            Vector orthogonal to the rectangle to be found
        Returns
        -------
        origin : Base.Vector
            Center of the rectangle
        best_v1, best_v2 : Base.Vector
            Unit vector parallel to the sides of the rectangle
        best_v2 : Base.Vector
            Other unit vector parallel to the sides of the rectangle
        length1 : float
            Length of side of the rectangle
        length2 : float
            Length of other side of the rectangle
        """
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
        length1 = best_length1 * 1.04
        length2 = best_length2 * 1.04
        origin = best_origin - best_v1 * length1 * 0.02 - best_v2 * length2 * 0.02
        return origin, best_v1, best_v2, length1, length2

    def random_point(self):
        """
        Returns a random point on the rectangle

        Returns
        -------
        Base.Vector
        """
        return (self.origin + self.v1 * self.length1 * myrandom() +
                self.v2 * self.length2 * myrandom())

    def random_direction(self):
        """
        Returns the main direction

        Maybe in the future will return some random vector

        Returns
        -------
        Base.Vector
        """
        return self.main_direction

    def add_to_document(self, doc):
        """
        Adds the rectangle to the FreeCAD document

        Parameters
        ----------
        doc : App.Document
        """
        sw = Part.makePolygon([self.origin,
                               self.origin + self.v1 * self.length1,
                               self.origin + self.v1 * self.length1 +
                               self.v2 * self.length2,
                               self.origin + self.v2 * self.length2
                               ], True)
        doc.addObject("Part::Feature", "SunWindow").Shape = sw


class BuieDistribution(object):
    """
    Implementation of the Buie Distribution for Sun emission

    A distribution for the direct light from the sun according to the Buie model:
    Buie, D., 2005. Corrigendum to "The effective size of the solar cone for
    solar concentrating systems" [Solar Energy 74 (2003) 417-427].
    Solar Energy 79, 568-570. 5.2

    Parameters
    ----------
    CircumSolarRatio : float

    Attributes
    ----------
    CSR : float
    SD : float
    SS : float
    aa1 : float
    aa2 : float
        Parameters describing the distribution
    CDF_Disk_Region : tuple of np.ndarray
        Cumulative Distribution Function in the solar disk region
    """

    def __init__(self, CircumSolarRatio):
        self.CSR = CircumSolarRatio
        self.SD = 4.65  # Solar Disk in mrad
        self.SS = 43.6  # Solar Size in mrad
        # normalization constant for the disk region
        self.aa1 = BuieDistribution.a1(self.CSR, self.SD)
        # normalization constant for the circumsolar region
        self.aa2 = BuieDistribution.a2(self.CSR, self.SD, self.SS)
        self.CDF_Disk_Region = BuieDistribution.CDF_disk_region(self.aa1,
                                                                self.SD)  # Buie distribution for the disk region

    def angle_distribution(self, random_number):
        """
        #TODO: Say one sentence @Ramon

        Parameters
        ----------
        random_number : float

        Returns
        -------
        float
        """
        u = random_number
        if u < 1.0 - self.CSR:
            th_u = BuieDistribution.th_solar_disk_region(u, self.aa1, self.SD) / 1000.0 * 180.0 / np.pi
        else:
            th_u = BuieDistribution.th_circumsolar_region(u, self.CSR, self.SD, self.SS,
                                                          self.aa2) / 1000.0 * 180.0 / np.pi
        return th_u

    @staticmethod
    def solar_disk__density_distribution(angle):
        """
        #TODO: Say one sentence @Ramon

        Parameters
        ----------
        angle : float

        Returns
        -------
        float
        """
        return 2.0 * np.pi * np.cos(0.326 * angle) / np.cos(0.305 * angle) * angle

    @staticmethod
    def circumsolar__density_distribution(angle, CSR):
        """
        #TODO: Say one sentence @Ramon

        Parameters
        ----------
        angle : float
        CSR : float

        Returns
        -------
        float
        """
        gamma = 2.2 * np.log(0.52 * CSR) * CSR ** (0.43) - 0.1
        kappa = 0.9 * np.log(13.5 * CSR) * CSR ** (-0.3)
        return 2.0 * np.pi * np.exp(kappa) * angle ** (gamma + 1.0)

    @staticmethod
    def a1(CSR, SD):
        """ Parameter a1 needed for the normalization of the probability distribution in the disk region
        #TODO: Say one sentence @Ramon

        Parameters
        ----------
        CSR
        SD

        Returns
        -------

        """
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
        """
        Cumulative Distribution Function in the solar disk region

        Parameters
        ----------
        a1 : float
        SD : float
            Parameters of the model

        Returns
        -------
        tuple of np.ndarray
        """
        f = BuieDistribution.solar_disk__density_distribution
        th = np.arange(0.0, SD, 0.001)
        f_th = np.vectorize(f)
        cumulative = np.cumsum(f_th(th))
        CDF = (th, a1 * cumulative / 1000.0)
        return CDF

    @staticmethod
    def th_solar_disk_region(u, a1, SD):
        """ Random angle based on the probability distribution in the circumsolar region"""
        # TODO: It is recomputed every time? @Ramon
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


class LightSource(object):
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
            wavelength = pick_random_from_cdf(self.light_spectrum)  # light spectrum is active (nanometers)
        ray = Ray(self.scene, point, direction, wavelength, self.initial_energy, polarization_vector)
        return ray

