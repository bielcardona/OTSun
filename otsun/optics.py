# -*- coding: utf-8 -*-

from FreeCAD import Base
import numpy as np
from .math import myrandom
from enum import Enum
from numpy.lib.scimath import sqrt



class Phenomenon(Enum):
    REFLEXION = 1
    REFRACTION = 2
    ABSORTION = 3
    TRANSMITANCE = 4
    GOT_ABSORBED = 5


class OpticalState(object):
    def __init__(self, polarization, direction, phenomenon):
        self.polarization = polarization
        self.direction = direction
        self.phenomenon = phenomenon


# ---
# Helper functions for reflexion and refraction
# ---
def reflexion(incident, normal, polarization_vector, polarization_vector_calculated_before=False):
    """
    Implementation of law of reflexion for incident and polarization vector.

    Parameters
    ----------
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


def lambertian_reflexion(incident, normal):
    """
    Implementation of lambertian reflection for diffusely reflecting surface
    """
    mynormal = normal * 1.0
    if mynormal.dot(incident) > 0:  # Ray intercepted on the backside of the surface
        mynormal = mynormal * (-1.0)
    dot = - 1.0
    while (dot <= 0.01):
        random_vector = Base.Vector(myrandom() - 0.5, myrandom() - 0.5, myrandom() - 0.5)
        random_vector.normalize()
        dot = mynormal.dot(random_vector)
    new_direction = random_vector
    random_polarization_vector = random_polarization(new_direction)  ## To be tested!!!!!
    return OpticalState(random_polarization_vector, new_direction, Phenomenon.REFLEXION)


def single_gaussian_dispersion(normal, state, sigma_1):
    """
    Implementation of single gaussian dispersion based on the ideal reflected direction
    """
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


def double_gaussian_dispersion(normal, state, sigma_1, sigma_2, k):
    """
    Implementation of double gaussian dispersion based on the ideal reflected direction
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


def TW_absorptance_ratio(normal, b_constant, c_constant, incident):
    """
    Implementation of the angular Solar Absorptance model for selective absorber material.
    1 - b * (1/cos - 1) ** c
    Model based on the study:
    Tesfamichael, T., and Wackelgard, E., 2000, "Angular Solar Absorptance and
    Incident Angle Modifier of Selective Absorbers for Solar Thermal Collectors,"
    Sol. Energy, 68, pp. 335–341.
    """
    # We assume the normal is normalized.
    mynormal = normal * 1.0
    if mynormal.dot(incident) > 0:  # Ray intercepted on the backside of the surface
        mynormal = mynormal * (-1.0)
    incidence_angle = np.arccos(mynormal.dot(incident) * (-1.0))
    incidence_angle_deg = incidence_angle * 180.0 / np.pi
    if incidence_angle_deg < 80.0:
        absortion_ratio = 1.0 - b_constant * (1.0 / np.cos(incidence_angle) - 1.0) ** c_constant
    else:
        y0 = 1.0 - b_constant * (1.0 / np.cos(80.0 * np.pi / 180.0) - 1.0) ** c_constant
        m = y0 / 10.0
        absortion_ratio = y0 - m * (incidence_angle_deg - 80.0)
    return absortion_ratio


def calculate_reflexion_metallic(incident, normal, n1, n2, polarization_vector):
    """
    Implementation of Fresnel equations for metallic materials
    """
    mynormal = normal * 1.0
    if mynormal.dot(incident) > 0:  # Ray intercepted on the backside of the surface
        # noinspection PyAugmentAssignment
        mynormal = mynormal * (-1.0)
    r = n1 / n2
    c1 = - mynormal.dot(incident)  # cos (incidence_angle)
    c2 = sqrt(1.0 - r * r * (1.0 - c1 * c1))  # cos (refracted_angle)

    normal_parallel_plane = incident.cross(mynormal)  # normal vector of the parallel plane
    if normal_parallel_plane == Base.Vector(0, 0, 0):  # to avoid null vector at mynormal and incident parallel vectors
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

    if myrandom() < ref_per:
        a = (n1 * c1 - n2 * c2) / (n1 * c1 + n2 * c2)
        R = a * a.conjugate()  # reflectance for s-polarized (perpendicular) light
        perpendicular_polarized = True
        polarization_vector = perpendicular_v.normalize()
    else:
        a = (n1 * c2 - n2 * c1) / (n1 * c2 + n2 * c1)
        R = a * a.conjugate()  # reflectance for p-polarized (parallel) light
        polarization_vector = parallel_v.normalize()
    if myrandom() < R.real:  # ray reflected
        return 1, 0, 0, polarization_vector, perpendicular_polarized, True
    else:  # ray refracted
        return 0, 1, 0, polarization_vector, perpendicular_polarized, True


# noinspection PyUnresolvedReferences,PyUnresolvedReferences
def refraction(incident, normal, n1, n2, polarization_vector):
    """
    Implementation of Snell's law of refraction
    """
    # https://en.wikipedia.org/wiki/Snell's_law#Vector_form
    mynormal = normal * 1.0
    if mynormal.dot(incident) > 0:  # Ray intercepted on the backside of the surface
        mynormal = mynormal * (-1.0)
    r = n1 / n2
    c1 = - mynormal.dot(incident)  # cos (incidence_angle)
    c2sq = 1.0 - r * r * (1.0 - c1 * c1)  # cos (refracted_angle) ** 2
    if c2sq.real < 0:  # total internal reflection
        return reflexion(incident, normal, polarization_vector)
    c2 = sqrt(c2sq)  # cos (refracted_angle)
    normal_parallel_plane = incident.cross(mynormal)  # normal vector of the parallel plane
    if normal_parallel_plane == Base.Vector(0, 0, 0):  # to avoid null vector at mynormal and incident parallel vectors
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
        R = a * a.conjugate()  # reflectance for s-polarized (perpendicular) light
        perpendicular_polarized = True
        polarization_vector = perpendicular_v.normalize()
    else:
        a = (n1 * c2 - n2 * c1) / (n1 * c2 + n2 * c1)
        R = a * a.conjugate()  # reflectance for p-polarized (parallel) light
        polarization_vector = parallel_v.normalize()
    if myrandom() < R.real:  # ray reflected
        return reflexion(incident, normal, polarization_vector, True)
    else:  # ray refracted
        refracted_direction = incident * r.real + mynormal * (r.real * c1.real - c2.real)
        perp_v = refracted_direction.cross(mynormal)
        if perp_v == Base.Vector(0, 0, 0):  # to avoid null vector at mynormal and incident parallel vectors
            perp_v = Base.Vector(1, 0, 0)
        para_v = refracted_direction.cross(perp_v)
        if perpendicular_polarized:
            return OpticalState(perp_v, refracted_direction, Phenomenon.REFRACTION)
        else:
            return OpticalState(para_v, refracted_direction, Phenomenon.REFRACTION)


def calculate_probabilities_polarizaton_coating(incident, normal, n1, n2, polarization_vector, properties, wavelength):
    # returns probability of Reflexion, probability of Absortion, probability of Transmitance, polarization_vector
    mynormal = normal * 1.0
    backside = False
    if mynormal.dot(incident) > 0:  # Ray intercepted on the backside of the surface
        # noinspection PyAugmentAssignment
        mynormal = mynormal * (-1.0)
        backside = True
    r = n1 / n2
    c1 = - mynormal.dot(incident)  # cos (incidence_angle)
    c2sq = 1.0 - r * r * (1.0 - c1 * c1)  # cos (refracted_angle) ** 2
    if properties['transparent_material']:  # transparent coating
        if c2sq.real < 0:  # total internal reflection
            return reflexion(incident, normal, polarization_vector)
    c2 = sqrt(c2sq)  # cos (refracted_angle)
    normal_parallel_plane = incident.cross(mynormal)  # normal vector of the parallel plane
    if normal_parallel_plane == Base.Vector(0, 0, 0):  # to avoid null vector at mynormal and incident parallel vectors
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

    if backside == True and properties['transparent_material']:  # Ray intercepted on the backside of the surface
        angle = np.arccos(c2.real) * 180.0 / np.pi
    else:
        angle = np.arccos(c1) * 180.0 / np.pi
    Matrix_Reflectance = properties['Matrix_polarized_reflectance_coating']
    Matrix_R = Matrix_Reflectance(angle, wavelength)
    if myrandom() < ref_per:
        R = calculate_reflectance(Matrix_R, angle, wavelength)[1]  # reflectance for s-polarized (perpendicular) light
        perpendicular_polarized = True
        polarization_vector = perpendicular_v.normalize()
    else:
        angle = np.arccos(c1) * 180.0 / np.pi
        R = calculate_reflectance(Matrix_R, angle, wavelength)[3]  # reflectance for p-polarized (parallel) light
        polarization_vector = parallel_v.normalize()
    if myrandom() < R:  # ray reflected
        return 1, 0, 0, polarization_vector, perpendicular_polarized
    else:  # ray refracted or absorbed
        if properties['energy_collector']:  # absorber coating
            return 0, 1, 0, polarization_vector, perpendicular_polarized
        if properties['specular_material']:  # reflector coating
            return 0, 1, 0, polarization_vector, perpendicular_polarized
        if properties['transparent_material']:  # transparent coating
            return 0, 0, 1, polarization_vector, perpendicular_polarized


def shure_refraction(incident, normal, n1, n2, polarization_vector, perpendicular_polarized):
    """
    Implementation of Snell's law of refraction
    """
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


def calculate_state_thin_film(incident, normal, n1, n2, polarization_vector, properties, wavelength):
    # returns optical state of the ray in thin film material
    mynormal = normal * 1.0
    backside = False
    if mynormal.dot(incident) > 0:  # Ray intercepted on the backside of the surface
        mynormal = mynormal * (-1.0)
        backside = True
    r = n1.real / n2.real
    c1 = - mynormal.dot(incident)  # cos (incidence_angle)
    c2sq = 1.0 - r * r * (1.0 - c1 * c1)  # cos (refracted_angle) ** 2
    if c2sq.real < 0:  # total internal reflection
        return 0.0, reflexion(incident, normal, polarization_vector)
    c2 = sqrt(c2sq)  # cos (refracted_angle)

    normal_parallel_plane = incident.cross(mynormal)  # normal vector of the parallel plane
    if normal_parallel_plane == Base.Vector(0, 0, 0):  # to avoid null vector at mynormal and incident parallel vectors
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

    if backside == True:  # Ray intercepted on the backside of the transparent surface
        angle = np.arccos(c2.real) * 180.0 / np.pi
    else:
        angle = np.arccos(c1) * 180.0 / np.pi
    Matrix_Reflectance = properties['Matrix_reflectance_thin_film']
    Matrix_R = Matrix_Reflectance(angle, wavelength)
    if myrandom() < ref_per:
        R = calculate_reflectance(Matrix_R, angle, wavelength)[1]  # reflectance for s-polarized (perpendicular) light
        perpendicular_polarized = True
        polarization_vector = perpendicular_v.normalize()
    else:
        angle = np.arccos(c1) * 180.0 / np.pi
        R = calculate_reflectance(Matrix_R, angle, wavelength)[3]  # reflectance for p-polarized (parallel) light
        polarization_vector = parallel_v.normalize()
    if myrandom() < R:  # ray reflected
        return 0.0, reflexion(incident, normal, polarization_vector)
    else:
        Matrix_Transmittance = properties['Matrix_transmittance_thin_film']
        Matrix_T = Matrix_Transmittance(angle, wavelength)
        if perpendicular_polarized:
            T = calculate_reflectance(Matrix_T, angle, wavelength)[1]
        else:
            T = calculate_reflectance(Matrix_T, angle, wavelength)[3]
        energy_absorbed_thin_film = (1 - R - T) / (1 - R)
        refracted_direction = incident * r.real + mynormal * (r.real * c1.real - c2.real)
        return energy_absorbed_thin_film, OpticalState(polarization_vector, refracted_direction, Phenomenon.REFRACTION)


# ---
# Helper function for dispersions and polarization vector
# ---
def dispersion_from_main_direction(main_direction, theta, phi):
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


def dispersion_polarization(main_direction, polarization_vector, theta, phi):
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


def random_polarization(direction):
    v_p = Base.Vector(direction[1], -direction[0], 0)
    if v_p == Base.Vector(0, 0, 0):
        v_p = Base.Vector(1, 0, 0)
    v_p.normalize()
    phi = 360.0 * myrandom()
    rotation = Base.Rotation(direction, phi)
    new_v_p = rotation.multVec(v_p)
    return new_v_p


# ---
# Helper function for reflectance depending on the wavelength for coatings layers.
# ---

def matrix_reflectance(data_material):
    """
    data_material: wavelenth in nm, angle in deg., reflectance s-polarized (perpendicular), reflectance p-polarized (parallel)
    the values should be in the corresponding order columns with constants steps
	We generate a matrix for the interpolation depending on the x=angle and y=wavelength values (lambda values)
    """
    steps = np.diff(data_material[:, 1])
    try:
        a = min(steps[steps > 0.0])
    except:
        a = 1E-4
    steps = np.diff(data_material[:, 0])
    try:
        w = min(steps[steps > 0.0])
    except:
        w = 1E-4
    return lambda x, y: [row for row in data_material if
                         (x + a > row[1] and x - a < row[1]) and (y + w > row[0] and y - w < row[0])]


def calculate_reflectance(matrix_reflectance, angle, wavelength):
    """
    data_array: values needed for the interpolation
	We generate an interpolation to find the perpendicular and parallel reflectance depending on the angle and wavelength
    """
    if len(matrix_reflectance) == 0:  # wavelength experiments out of range of the data_material
        return 'R_per', 0.0, 'R_par', 0.0
    R_Matrix = np.mat(matrix_reflectance)
    if len(R_Matrix) == 1:  # interpolation is not needed
        return 'R_per', R_Matrix[0, 2], 'R_par', R_Matrix[0, 3]

    if len(R_Matrix) == 2:  # simple interpolation
        if R_Matrix[0, 0] == R_Matrix[1, 0]:  # identical wavelengths; hence angle interpolation (column 1)
            xvalues = [R_Matrix[0, 1], R_Matrix[1, 1]]
            yvalues = [R_Matrix[0, 2], R_Matrix[1, 2]]
            R_per = np.interp(angle, xvalues, yvalues)

            xvalues = [R_Matrix[0, 1], R_Matrix[1, 1]]
            yvalues = [R_Matrix[0, 3], R_Matrix[1, 3]]
            R_par = np.interp(angle, xvalues, yvalues)
            return 'R_per', R_per, 'R_par', R_par

        if R_Matrix[0, 1] == R_Matrix[1, 1]:  # identical angles; hence wavelength interpolation (column 0)
            xvalues = [R_Matrix[0, 0], R_Matrix[1, 0]]
            yvalues = [R_Matrix[0, 2], R_Matrix[1, 2]]
            R_per = np.interp(wavelength, xvalues, yvalues)

            xvalues = [R_Matrix[0, 0], R_Matrix[1, 0]]
            yvalues = [R_Matrix[0, 3], R_Matrix[1, 3]]
            R_par = np.interp(wavelength, xvalues, yvalues)
            return 'R_per', R_per, 'R_par', R_par

    if len(R_Matrix) == 4:  # double interpolation
        xvalues = [R_Matrix[0, 1], R_Matrix[1, 1]]
        yvalues = [R_Matrix[0, 2], R_Matrix[1, 2]]
        R1 = np.interp(angle, xvalues, yvalues)
        xvalues = [R_Matrix[2, 1], R_Matrix[3, 1]]
        yvalues = [R_Matrix[2, 2], R_Matrix[3, 2]]
        R2 = np.interp(angle, xvalues, yvalues)
        xvalues = [R_Matrix[0, 0], R_Matrix[3, 0]]
        yvalues = [R1, R2]
        R_per = np.interp(wavelength, xvalues, yvalues)

        xvalues = [R_Matrix[0, 1], R_Matrix[1, 1]]
        yvalues = [R_Matrix[0, 3], R_Matrix[1, 3]]
        R1 = np.interp(angle, xvalues, yvalues)
        xvalues = [R_Matrix[2, 1], R_Matrix[3, 1]]
        yvalues = [R_Matrix[2, 3], R_Matrix[3, 3]]
        R2 = np.interp(angle, xvalues, yvalues)
        xvalues = [R_Matrix[0, 0], R_Matrix[3, 0]]
        yvalues = [R1, R2]
        R_par = np.interp(wavelength, xvalues, yvalues)
        return 'R_per', R_per, 'R_par', R_par

# endregion

