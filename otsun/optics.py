# -*- coding: utf-8 -*-

from FreeCAD import Base
import numpy as np
from .math import myrandom
from enum import Enum
from numpy.lib.scimath import sqrt
from autologging import traced
from .logging_unit import logger


class Phenomenon(Enum):
    """Enum for optical phenomena
    """
    REFLEXION = 1
    REFRACTION = 2
    ABSORPTION = 3
    TRANSMITTANCE = 4
    GOT_ABSORBED = 5


class OpticalState(object):
    """Optical state of a ray.

    The OpticalState class gathers together information about the optical
    state of a ray.
    """
    
    def __init__(self, polarization, direction, phenomenon):
        """Optical state of a ray.

        The OpticalState class gathers together information about the optical
        state of a ray.

        Parameters
        ----------
        polarization : Base.Vector
            polarization vector of the ray
        direction : Base.Vector
            direction vector of the ray
        phenomenon : Phenomenon
            last phenomenon that the ray experimented
        """
        self.polarization = polarization
        self.direction = direction
        self.phenomenon = phenomenon


# ---
# Helper functions for reflexion and refraction
# ---
@traced(logger)
def reflexion(incident, normal, polarization_vector, polarization_vector_calculated_before=False):
    """
    Implementation of law of reflexion for incident and polarization vector.

    Parameters
    ----------
    normal : Base.Vector
        normal vector
    polarization_vector :
    polarization_vector_calculated_before :
    incident : Base.Vector
        incident vector

    Returns
    -------
    Base.Vector
        reflected vector from incident and normal surface
    """
    mynormal = normal * 1.0
    if mynormal.dot(incident) > 0:  # Ray intercepted on the backside of the surface
        mynormal = mynormal * (-1.0)
    reflected = incident - mynormal * 2.0 * mynormal.dot(incident)  # refelcted vector
    c1 = - mynormal.dot(incident)  # cos (incidence_angle)
    angle = 2.0 * (np.pi - np.arccos(c1)) * 180.0 / np.pi  # angle for make a rotation to the polarizacion vector
    if not polarization_vector_calculated_before:
        normal_parallel_plane = incident.cross(
            mynormal)  # axis for the rotation, only parallel component must be rotated
        if normal_parallel_plane == Base.Vector(0, 0, 0):  # to avoid null vector in a common case
            normal_parallel_plane = Base.Vector(1, 0, 0)
        rotation = Base.Rotation(normal_parallel_plane, angle)
        polarization_vector = rotation.multVec(polarization_vector)
    else:
        polarization_vector = polarization_vector
    return OpticalState(polarization_vector, reflected, Phenomenon.REFLEXION)


@traced(logger)
def lambertian_reflexion(incident, normal):
    """
    Implementation of lambertian reflection for diffusely reflecting surface

    Parameters
    ----------
    incident : Base.Vector
        direction vector of the incident ray
    normal : Base.Vector
        normal vector of the surface at the point of incidence

    Returns
    -------
    OpticalState
    """
    my_normal = normal * 1.0
    if my_normal.dot(incident) > 0:  # Ray intercepted on the backside of the surface
        my_normal = my_normal * (-1.0)
    dot = - 1.0
    while dot <= 0.01:
        random_vector = Base.Vector(myrandom() - 0.5, myrandom() - 0.5, myrandom() - 0.5)
        random_vector.normalize()
        dot = my_normal.dot(random_vector)
    new_direction = random_vector
    random_polarization_vector = random_polarization(new_direction)  ## To be tested!!!!!
    return OpticalState(random_polarization_vector, new_direction, Phenomenon.REFLEXION)


@traced(logger)
def single_gaussian_dispersion(normal, state, sigma_1):
    """
    Single gaussian dispersion based on the ideal reflected direction

    Computes

    Parameters
    ----------
    normal : Base.Vector
        normal vector of the surface at the point of incidence
    state : OpticalState
        optical state of the ray
    sigma_1 : float
        dispersion coefficient

    Returns
    -------
    OpticalState

    """
    # TODO: @Ramon: Review
    rad = np.pi / 180.0
    u = myrandom()
    theta = (-2. * sigma_1 ** 2. * np.log(u)) ** 0.5 / 1000.0 / rad
    v = state.direction
    axis_1 = normal.cross(state.direction)
    rotation_1 = Base.Rotation(axis_1, theta)
    new_v1 = rotation_1.multVec(v)
    u = myrandom()
    phi = 360. * u
    axis_2 = v
    rotation_2 = Base.Rotation(axis_2, phi)
    new_v2 = rotation_2.multVec(new_v1)
    polarization_vector = state.polarization
    new_pol_1 = rotation_1.multVec(polarization_vector)
    new_polarization_vector = rotation_2.multVec(new_pol_1)
    return OpticalState(new_polarization_vector, new_v2, Phenomenon.REFLEXION)


@traced(logger)
def double_gaussian_dispersion(normal, state, sigma_1, sigma_2, k):
    """
    Double gaussian dispersion based on the ideal reflected direction

    Computes

    Parameters
    ----------
    normal : Base.Vector
        normal vector of the surface at the point of incidence
    state : OpticalState
        optical state of the ray
    sigma_1 : float
        dispersion coefficient (first case)
    sigma_2 : float
        dispersion coefficient (second case)
    k : float
        threshold for randomly applying first or second case

    Returns
    -------
    OpticalState

    """
    rad = np.pi / 180.0
    k_ran = myrandom()
    u = myrandom()
    if k_ran < k:
        theta = (-2. * sigma_1 ** 2. * np.log(u)) ** 0.5 / 1000.0 / rad
    else:
        theta = (-2. * sigma_2 ** 2. * np.log(u)) ** 0.5 / 1000.0 / rad
    v = state.direction
    axis_1 = normal.cross(state.direction)
    rotation_1 = Base.Rotation(axis_1, theta)
    new_v1 = rotation_1.multVec(v)
    u = myrandom()
    phi = 360. * u
    axis_2 = v
    rotation_2 = Base.Rotation(axis_2, phi)
    new_v2 = rotation_2.multVec(new_v1)
    polarization_vector = state.polarization
    new_pol_1 = rotation_1.multVec(polarization_vector)
    new_polarization_vector = rotation_2.multVec(new_pol_1)
    return OpticalState(new_polarization_vector, new_v2, Phenomenon.REFLEXION)


@traced(logger)
def refraction(incident, normal, n1, n2, polarization_vector):
    """Implementation of Snell's law of refraction

    Parameters
    ----------
    incident
    normal
    n1
    n2
    polarization_vector

    Returns
    -------

    """
    # TODO: document
    # https://en.wikipedia.org/wiki/Snell's_law#Vector_form
    my_normal = normal * 1.0
    if my_normal.dot(incident) > 0:  # Ray intercepted on the backside of the surface
        my_normal = my_normal * (-1.0)
    r = n1 / n2
    c1 = - my_normal.dot(incident)  # cos (incidence_angle)
    c2sq = 1.0 - r * r * (1.0 - c1 * c1)  # cos (refracted_angle) ** 2
    if c2sq.real < 0:  # total internal reflection
        return reflexion(incident, normal, polarization_vector)
    c2 = sqrt(c2sq)  # cos (refracted_angle)
    normal_parallel_plane = incident.cross(my_normal)  # normal vector of the parallel plane
    if normal_parallel_plane == Base.Vector(0, 0, 0):  # to avoid null vector at my_normal and incident parallel vectors
        normal_parallel_plane = Base.Vector(1, 0, 0)
    normal_parallel_plane.normalize()
    normal_perpendicular_plane = normal_parallel_plane.cross(incident)  # normal vector of the perpendicular plane
    # http://www.maplesoft.com/support/help/Maple/view.aspx?path=MathApps/ProjectionOfVectorOntoPlane
    parallel_v = polarization_vector - normal_parallel_plane * polarization_vector.dot(normal_parallel_plane)
    parallel_component = parallel_v.Length
    perpendicular_v = polarization_vector - normal_perpendicular_plane * polarization_vector.dot(
        normal_perpendicular_plane)
    perpendicular_component = perpendicular_v.Length
    ref_per = perpendicular_component / (perpendicular_component + parallel_component)
    perpendicular_polarized = False
    # https://en.wikipedia.org/wiki/Fresnel_equations # Fresnel equations
    if myrandom() < ref_per:
        a = (n1 * c1 - n2 * c2) / (n1 * c1 + n2 * c2)
        r = a * a.conjugate()  # reflectance for s-polarized (perpendicular) light
        perpendicular_polarized = True
        polarization_vector = perpendicular_v.normalize()
    else:
        a = (n1 * c2 - n2 * c1) / (n1 * c2 + n2 * c1)
        r = a * a.conjugate()  # reflectance for p-polarized (parallel) light
        polarization_vector = parallel_v.normalize()
    if myrandom() < r.real:  # ray reflected
        return reflexion(incident, normal, polarization_vector, True)
    else:  # ray refracted
        refracted_direction = incident * r.real + my_normal * (r.real * c1.real - c2.real)
        perp_v = refracted_direction.cross(my_normal)
        if perp_v == Base.Vector(0, 0, 0):  # to avoid null vector at my_normal and incident parallel vectors
            perp_v = Base.Vector(1, 0, 0)
        para_v = refracted_direction.cross(perp_v)
        if perpendicular_polarized:
            return OpticalState(perp_v, refracted_direction, Phenomenon.REFRACTION)
        else:
            return OpticalState(para_v, refracted_direction, Phenomenon.REFRACTION)


@traced(logger)
def shure_refraction(incident, normal, n1, n2, polarization_vector, perpendicular_polarized):
    """Implementation of Snell's law of refraction

    Parameters
    ----------
    incident
    normal
    n1
    n2
    polarization_vector
    perpendicular_polarized

    Returns
    -------

    """
    # TODO: document
    # https://en.wikipedia.org/wiki/Snell's_law#Vector_form
    mynormal = normal * 1.0
    if mynormal.dot(incident) > 0:  # Ray intercepted on the backside of the surface
        mynormal = mynormal * (-1.0)
    r = n1.real / n2.real
    c1 = - mynormal.dot(incident)  # cos (incidence_angle)
    c2sq = 1.0 - r * r * (1.0 - c1 * c1)  # cos (refracted_angle) ** 2
    c2 = c2sq ** 0.5  # cos (refracted_angle)
    refracted_direction = incident * r + mynormal * (r * c1 - c2)
    if perpendicular_polarized == 'thin film':  # refraction in thin film
        return OpticalState(polarization_vector, refracted_direction, Phenomenon.REFRACTION)
    perp_v = refracted_direction.cross(mynormal)
    if perp_v == Base.Vector(0, 0, 0):  # to avoid null vector at mynormal and incident parallel vectors
        perp_v = Base.Vector(1, 0, 0)
    para_v = refracted_direction.cross(perp_v)
    if perpendicular_polarized:
        return OpticalState(perp_v, refracted_direction, Phenomenon.REFRACTION)
    else:
        return OpticalState(para_v, refracted_direction, Phenomenon.REFRACTION)


# ---
# Helper function for dispersions and polarization vector
# ---
@traced(logger)
def dispersion_from_main_direction(main_direction, theta, phi):
    """

    Parameters
    ----------
    main_direction
    theta
    phi

    Returns
    -------

    """
    # TODO: document
    v = main_direction
    v_p = Base.Vector(v[1], -v[0], 0)
    if v_p == Base.Vector(0, 0, 0):  # to avoid null vector at mynormal and incident parallel vectors
        v_p = Base.Vector(1, 0, 0)
    v_p.normalize()
    rotation_1 = Base.Rotation(v_p, theta)
    new_v1 = rotation_1.multVec(v)
    axis_2 = main_direction
    rotation_2 = Base.Rotation(axis_2, phi)
    new_v2 = rotation_2.multVec(new_v1)
    return new_v2


@traced(logger)
def dispersion_polarization(main_direction, polarization_vector, theta, phi):
    """

    Parameters
    ----------
    main_direction
    polarization_vector
    theta
    phi

    Returns
    -------

    """
    # TODO: document
    v = main_direction
    v_p = Base.Vector(v[1], -v[0], 0)
    if v_p == Base.Vector(0, 0, 0):  # to avoid null vector at mynormal and incident parallel vectors
        v_p = Base.Vector(1, 0, 0)
    v_p.normalize()
    rotation_1 = Base.Rotation(v_p, theta)
    new_v1 = rotation_1.multVec(polarization_vector)
    axis_2 = main_direction
    rotation_2 = Base.Rotation(axis_2, phi)
    new_v2 = rotation_2.multVec(new_v1)
    return new_v2


@traced(logger)
def random_polarization(direction):
    """

    Parameters
    ----------
    direction

    Returns
    -------

    """
    # TODO: document
    v_p = Base.Vector(direction[1], -direction[0], 0)
    if v_p == Base.Vector(0, 0, 0):
        v_p = Base.Vector(1, 0, 0)
    v_p.normalize()
    phi = 360.0 * myrandom()
    rotation = Base.Rotation(direction, phi)
    new_v_p = rotation.multVec(v_p)
    return new_v_p

