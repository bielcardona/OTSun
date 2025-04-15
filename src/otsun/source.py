"""
Module otsun.source that implements rays and its sources
"""

import itertools
from typing import Callable

import Part
import numpy as np
from FreeCAD import Base, Document
from numpy import ndarray

from abc import ABC, abstractmethod

from .math import pick_random_from_cdf, myrandom, tabulated_function, two_orthogonal_vectors, area_of_triangle, random_point_of_triangle
from .optics import dispersion_from_main_direction, random_polarization, dispersion_polarization
from .ray import Ray
from .scene import Scene
from random import choices

from scipy.spatial import ConvexHull

EPSILON = 1E-6
# Tolerance for considering equal to zero

class AbstractEmittingObject(ABC):
    @abstractmethod
    def random_point(self) -> Base.Vector:
        pass

    @abstractmethod
    def normal_direction(self, point: Base.Vector) -> Base.Vector:
        pass

    # def random_direction(self) -> Base.Vector:
    #     ...


class GeneralizedSunWindow(AbstractEmittingObject):
    def __init__(
            self, 
            scene: Scene, 
            main_direction: Base.Vector):
        bbs = []
        for shape in itertools.chain(scene.solids, scene.faces):
            bbs.append(shape.BoundBox)
        origin = (scene.boundbox.Center -
                          main_direction * 0.5 * scene.boundbox.DiagonalLength)
        u, v = two_orthogonal_vectors(main_direction)
        points = []
        plane_points = []
        for bb in bbs:
            xs = [bb.XMin, bb.XMax]
            ys = [bb.YMin, bb.YMax]
            zs = [bb.ZMin, bb.ZMax]
            coords = itertools.product(xs, ys, zs)
            for c in coords:
                point = Base.Vector(c)
                points.append(point)
                pu = u.dot(point - origin) * 1.02
                pv = v.dot(point - origin) * 1.02
                plane_points.append((pu,pv))
        hull = ConvexHull(plane_points)
        plane_vertices = [plane_points[i] for i in hull.vertices]
        self.vertices: list[Base.Vector] = [origin + u*pu + v*pv for (pu,pv) in plane_vertices]
        self.triangles: list[tuple[Base.Vector, Base.Vector, Base.Vector]] = [
            (self.vertices[0], self.vertices[i], self.vertices[i+1])
            for i in range(1,len(self.vertices)-1)]
        self.triangle_areas: list[float] = list(map(area_of_triangle, self.triangles))
        self.aperture: float = sum(self.triangle_areas)
        self.main_direction: Base.Vector = main_direction

    def add_to_document(self, doc: Document) -> None:
        sw = Part.makePolygon(self.vertices, True)
        doc.addObject("Part::Feature", "SunWindow").Shape = sw

    def random_point(self) -> Base.Vector:
        random_triangle = choices(self.triangles, self.triangle_areas)[0]
        return random_point_of_triangle(random_triangle)

    def normal_direction(self, point: Base.Vector) -> Base.Vector:
        return self.main_direction

    # def random_direction(self) -> Base.Vector:
    #     """
    #     Returns the main direction
    #
    #     Maybe in the future will return some random vector
    #
    #     Returns
    #     -------
    #     Base.Vector
    #     """
    #     return self.main_direction


class SunWindow(AbstractEmittingObject):
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

    def __init__(self, scene: Scene, main_direction: Base.Vector):
        bbs = []
        for shape in itertools.chain(scene.solids, scene.faces):
            bbs.append(shape.BoundBox)
        projected_points = []
        for bb in bbs:
            xs = [bb.XMin, bb.XMax]
            ys = [bb.YMin, bb.YMax]
            zs = [bb.ZMin, bb.ZMax]
            coords = itertools.product(xs, ys, zs)
            points = [Base.Vector(c) for c in coords]

            point_of_plane = (scene.boundbox.Center -
                              main_direction * 0.5 * scene.boundbox.DiagonalLength)
            projected_points.extend([p.projectToPlane(point_of_plane, main_direction)
                                for p in points])
        (self.origin, self.v1, self.v2, self.length1, self.length2) = (
            SunWindow.find_min_rectangle(projected_points, main_direction))
        self.aperture = self.length1 * self.length2
        self.main_direction = main_direction

    @staticmethod
    def find_min_rectangle(points: list[Base.Vector], normal: Base.Vector
        ) -> tuple[Base.Vector, Base.Vector, Base.Vector, float, float]:
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
            if v1.Length < EPSILON:  # TODO: customize
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

    def random_point(self) -> Base.Vector :
        """
        Returns a random point on the rectangle

        Returns
        -------
        Base.Vector
        """
        return (self.origin + self.v1 * self.length1 * myrandom() +
                self.v2 * self.length2 * myrandom())

    def normal_direction(self, point: Base.Vector) -> Base.Vector:
        return self.main_direction


    # def random_direction(self) -> Base.Vector :
    #     """
    #     Returns the main direction
    #
    #     Maybe in the future will return some random vector
    #
    #     Returns
    #     -------
    #     Base.Vector
    #     """
    #     return self.main_direction

    def add_to_document(self, doc: Document) -> None:
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


class EmittingSphere(AbstractEmittingObject):
    def __init__(self, sphere: Part.Feature):
        self.sphere = sphere

    def random_point(self) -> Base.Vector:
        u = 2.0 * np.pi * myrandom()
        v = np.arccos(2.0 * myrandom() - 1.0)
        r = self.sphere.Radius * (1.0 + EPSILON)
        x = r * np.sin(v) * np.cos(u)
        y = r * np.sin(v) * np.sin(u)
        z = r * np.cos(v)
        return Base.Vector(x, y, z) + self.sphere.Placement.Base

    def normal_direction(self, point: Base.Vector) -> Base.Vector:
        return point - self.sphere.Placement.Base

    # def random_direction(self) -> Base.Vector:
    #     ...


class LightSource(object):
    """
    Sets up a light source with a given scene, a given emitting region and a given light spectrum.
    The emitting region must provide the main direction.
    Light spectrum could be: a constant value (for a single wavelength in nanometers), or a spectrum distribution.
    The distribution (dispersion) for the main direction is provided in "direction_distribution".
    The polarization_vector is a Base.Vector for polarized light. If is not given unpolarized light is generated.
    """

    def __init__(
            self,
            scene: Scene,
            emitting_region: AbstractEmittingObject,
            light_spectrum: float | tuple[ndarray, ndarray],
            initial_energy: float,
            direction_distribution: Callable[[float],float]|None = None,
            polarization_vector: Base.Vector|None = None
    ):
        self.scene = scene
        self.emitting_region = emitting_region
        self.light_spectrum = light_spectrum
        self.initial_energy = initial_energy
        self.direction_distribution = direction_distribution
        self.polarization_vector = polarization_vector
        self.wavelengths = []

    def emit_ray(self) -> Ray:
        """
        Simulates the emission of a ray
        """
        point = self.emitting_region.random_point()
        ## main_direction = self.emitting_region.main_direction  # emitting main direction
        main_direction = self.emitting_region.normal_direction(point)
        direction = main_direction
        if self.direction_distribution is not None:  # main direction has a distribution
            theta = self.direction_distribution(myrandom())
            phi = 360.0 * myrandom()
            direction = dispersion_from_main_direction(main_direction, theta, phi)
            if self.polarization_vector:  # single polarization vector is active
                polarization_vector = dispersion_polarization(main_direction, self.polarization_vector, theta, phi)
        if self.polarization_vector is None:  # unpolarization is active
            polarization_vector = random_polarization(direction)  # random polarization from light direction
        else:
            polarization_vector = self.polarization_vector
            polarization_vector.normalize()
        if np.isscalar(self.light_spectrum):
            wavelength = self.light_spectrum  # experiment with a single wavelength (nanometers)
        else:
            wavelength = pick_random_from_cdf(self.light_spectrum)  # light spectrum is active (nanometers)
        ray = Ray(self.scene, point, direction, wavelength, self.initial_energy, polarization_vector)
        return ray


# Auxiliary functions for buie_distribution

def _calculate_a1(csr: float, sd: float) -> float:
    """ Parameter a1 needed for the normalization of the probability distribution in the disk region
    """

    th = np.arange(0.0, sd, 0.001)
    values = [2.0 * np.pi * np.cos(0.326 * angle) / np.cos(0.305 * angle) * angle
              for angle in th]
    a1 = (1.0 - csr) / np.trapz(values, dx=0.001)
    return a1


def _circumsolar__density_distribution(angle: float, csr: float) -> float:
    gamma = 2.2 * np.log(0.52 * csr) * csr ** 0.43 - 0.1
    kappa = 0.9 * np.log(13.5 * csr) * csr ** (-0.3)
    return 2.0 * np.pi * np.exp(kappa) * angle ** (gamma + 1.0)


def _calculate_a2(csr: float, sd: float, ss: float) -> float:
    """ Parameter a2 needed for the normalization of the probability distribution in the circumsolar region"""
    f = _circumsolar__density_distribution
    th = np.arange(sd, ss, 0.001)
    f_th = np.vectorize(f)
    a2 = csr / np.trapz(f_th(th, csr), dx=0.001)
    return a2


def _calculate_cdf_disk_region(a1: float, sd: float) -> tuple[np.ndarray, np.ndarray]:
    """
    Cumulative Distribution Function in the solar disk region
    """
    th = np.arange(0.0, sd, 0.001)
    values = [2.0 * np.pi * np.cos(0.326 * angle) / np.cos(0.305 * angle) * angle
              for angle in th]
    cumulative = np.cumsum(values)
    cdf = (th, a1 * cumulative / 1000.0)
    return cdf


def _th_solar_disk_region(u: float, cdf: tuple[ndarray, ndarray]) -> float:
    """ Random angle based on the probability distribution in the circumsolar region"""
    # TODO: It is recomputed every time? @Ramon
    # CDF = CDF_Disk_Region # calculate_CDF_disk_region(a1, SD)
    idx = (np.abs(cdf[1] - u)).argmin()
    return cdf[0][idx]


def _th_circumsolar_region(u: float, csr: float, sd: float, a2: float) -> float:
    """ Random angle based on the CDF in the circumsolar region"""
    gamma = 2.2 * np.log(0.52 * csr) * csr ** 0.43 - 0.1
    kappa = 0.9 * np.log(13.5 * csr) * csr ** (-0.3)
    u_csr = (u - 1.0) + csr  # Since CSR-CDF starts at zero
    f1 = u_csr * (gamma + 2.0) / (a2 * 2 * np.pi * np.exp(kappa))
    f2 = sd ** (gamma + 2.0)
    exp = (1.0 / (gamma + 2.0))
    th_u = np.power(f1 + f2, exp)
    return th_u


def buie_distribution(circum_solar_ratio: float) -> Callable[[float], float]:
    """
    Implementation of the Buie Distribution for Sun emission

    A distribution for the direct light from the sun according to the Buie model:
    Buie, D., 2005. Corrigendum to "The effective size of the solar cone for
    solar concentrating systems" [Solar Energy 74 (2003) 417-427].
    Solar Energy 79, 568-570. 5.2

    Parameters
    ----------
    circum_solar_ratio : float

    Returns
    -------
    angle distribution for random input: function
        Function that interpolates by straight line segments the input data
    """
    csr = circum_solar_ratio
    sd = 4.65
    # Solar Disk in mrad
    ss = 43.6
    # Solar Size in mrad
    a2 = None # _calculate_a2(CSR, SD, SS)
    cdf_disk_region = None # _calculate_CDF_disk_region(a1, SD)
    #    Buie distribution for the disk region
    u_values = np.arange(0.0, 1.001, 0.001)
    dist_values = []
    for u in u_values:
        if u <= 1.0 - csr:
            if cdf_disk_region is None:
                a1 = _calculate_a1(csr, sd)
                #     normalization constant for the disk region
                cdf_disk_region = _calculate_cdf_disk_region(a1, sd)
            dist_values.append(_th_solar_disk_region(u, cdf_disk_region) / 1000.0 * 180.0 / np.pi)
        else:
            if a2 is None:
                a2 = _calculate_a2(csr, sd, ss)
                #    normalization constant for the circumsolar region
            dist_values.append(_th_circumsolar_region(u, csr, sd, a2) / 1000.0 * 180.0 / np.pi)
    f = tabulated_function(u_values, dist_values)
    return f
