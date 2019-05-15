"""
Implementation of optical effects on rays
"""

from FreeCAD import Base
import numpy as np
from .math import myrandom
from enum import Enum
from numpy.lib.scimath import sqrt
from autologging import traced
from .logging_unit import logger
# from .materials import Material, vacuum_medium


class Phenomenon(Enum):
    """
    Enum for optical phenomena.

    This enum is used to describe which optical phenomenon affected a ray
    """
    JUST_STARTED = 0
    REFLEXION = 1
    REFRACTION = 2
    ABSORPTION = 3
    TRANSMITTANCE = 4
    GOT_ABSORBED = 5


class OpticalState(object):
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
    material : Material
        Material where the ray is located
    extra_data : dict
        Dictionary where materials can put extra data
    """
    
    def __init__(self, polarization, direction, phenomenon, material = None, extra_data = None):
        self.polarization = polarization
        self.direction = direction
        self.phenomenon = phenomenon
        self.material = material
        if extra_data is None:
            extra_data = {}
        self.extra_data = extra_data

    def apply_single_gaussian_dispersion(self, normal, sigma_1):
        """
        Apply a single gaussian dispersion to the optical state

        Parameters
        ----------
        normal : Base.Vector
            normal vector of the surface at the point of incidence
        sigma_1 : float
            dispersion coefficient
         """
        # TODO: @Ramon: Review
        rad = np.pi / 180.0
        u = myrandom()
        theta = (-2. * sigma_1 ** 2. * np.log(u)) ** 0.5 / 1000.0 / rad
        v = self.direction
        axis_1 = normal.cross(self.direction)
        rotation_1 = Base.Rotation(axis_1, theta)
        new_v1 = rotation_1.multVec(v)
        u = myrandom()
        phi = 360. * u
        axis_2 = v
        rotation_2 = Base.Rotation(axis_2, phi)
        new_v2 = rotation_2.multVec(new_v1)
        polarization_vector = self.polarization
        new_pol_1 = rotation_1.multVec(polarization_vector)
        new_polarization_vector = rotation_2.multVec(new_pol_1)
        self.polarization = new_polarization_vector
        self.direction = new_v2

    def apply_double_gaussian_dispersion(self, normal, sigma_1, sigma_2, k):
        """
        Apply a double gaussian dispersion to the state


        Parameters
        ----------
        normal : Base.Vector
            normal vector of the surface at the point of incidence
        sigma_1 : float
            dispersion coefficient (first case)
        sigma_2 : float
            dispersion coefficient (second case)
        k : float
            threshold for randomly applying first or second case
        """
        rad = np.pi / 180.0
        k_ran = myrandom()
        u = myrandom()
        if k_ran < k:
            theta = (-2. * sigma_1 ** 2. * np.log(u)) ** 0.5 / 1000.0 / rad
        else:
            theta = (-2. * sigma_2 ** 2. * np.log(u)) ** 0.5 / 1000.0 / rad
        v = self.direction
        axis_1 = normal.cross(self.direction)
        rotation_1 = Base.Rotation(axis_1, theta)
        new_v1 = rotation_1.multVec(v)
        u = myrandom()
        phi = 360. * u
        axis_2 = v
        rotation_2 = Base.Rotation(axis_2, phi)
        new_v2 = rotation_2.multVec(new_v1)
        polarization_vector = self.polarization
        new_pol_1 = rotation_1.multVec(polarization_vector)
        new_polarization_vector = rotation_2.multVec(new_pol_1)
        self.direction = new_v2
        self.polarization = new_polarization_vector


# ---
# Helper functions for reflexion and refraction
# ---
@traced(logger)
def reflexion(incident, normal, polarization_vector,
              polarization_vector_calculated_before=False):
    """
    Implementation of law of reflexion for incident and polarization vector.

    Parameters
    ----------
    incident : Base.Vector
        vector incident
    normal : Base.Vector
        vector normal to the surface
    polarization_vector :
    polarization_vector_calculated_before : bool
        True if the polarization vector was computed before
    incident : Base.Vector
        incident vector

    Returns
    -------
    OpticalState
        optical state of the reflected ray
    """
    my_normal = normal * 1.0
    if my_normal.dot(incident) > 0:
        # Ray intercepted on the backside of the surface
        my_normal = my_normal * (-1.0)
    reflected = incident - my_normal * 2.0 * my_normal.dot(incident)
    c1 = - my_normal.dot(incident)  # cos (incidence_angle)
    # angle for making a rotation to the polarizacion vector
    # TODO: @Ramon Review if formula for angle is OK
    angle = 2.0 * (np.pi - np.arccos(c1)) * 180.0 / np.pi
    if not polarization_vector_calculated_before:
        # axis for the rotation, only parallel component must be rotated
        normal_parallel_plane = incident.cross(my_normal)
        if normal_parallel_plane == Base.Vector(0, 0, 0):
            # to avoid null vector in a common case
            normal_parallel_plane = Base.Vector(1, 0, 0)
        rotation = Base.Rotation(normal_parallel_plane, angle)
        polarization_vector = rotation.multVec(polarization_vector)
    # else:
    #     polarization_vector = polarization_vector
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
        # TODO: this is not uniform!
        random_vector = Base.Vector(myrandom() - 0.5, myrandom() - 0.5, myrandom() - 0.5)
        random_vector.normalize()
        dot = my_normal.dot(random_vector)
    new_direction = random_vector
    random_polarization_vector = random_polarization(new_direction)  ## To be tested!!!!!
    return OpticalState(random_polarization_vector, new_direction, Phenomenon.REFLEXION)


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
    if my_normal.dot(incident) > 0:
        # Ray intercepted on the backside of the surface
        my_normal = my_normal * (-1.0)
    r = n1 / n2
    c1 = - my_normal.dot(incident)  # cos (incidence_angle)
    c2sq = 1.0 - r * r * (1.0 - c1 * c1)  # cos (refracted_angle) ** 2
    if c2sq.real < 0:
        # total internal reflection
        return reflexion(incident, normal, polarization_vector)
    c2 = sqrt(c2sq)  # cos (refracted_angle)
    normal_parallel_plane = incident.cross(my_normal)
    # vector perpendicular to incident and parallel to plane
    if normal_parallel_plane == Base.Vector(0, 0, 0):
        # to avoid null vector at my_normal and incident parallel vectors
        normal_parallel_plane = Base.Vector(1, 0, 0)
    normal_parallel_plane.normalize()
    parallel_v = polarization_vector - \
                 normal_parallel_plane * polarization_vector.dot(normal_parallel_plane)
    # TODO: Review
    # parallel_v is the projection of polarization_vector on the plane
    # containing the normal to the surface and the incident vector
    perpendicular_v = polarization_vector - my_normal * polarization_vector.dot(my_normal)
    # TODO: Review
    # perpendicular_v is the projection of polarization_vector on the plane
    # parallel to the surface
    parallel_component = parallel_v.Length
    perpendicular_component = perpendicular_v.Length
    ref_per = perpendicular_component / (perpendicular_component + parallel_component)
    # https://en.wikipedia.org/wiki/Fresnel_equations # Fresnel equations
    if myrandom() < ref_per:
        # Perpendicular polarization
        a = (n1 * c1 - n2 * c2) / (n1 * c1 + n2 * c2)
        r = a * a.conjugate()  # reflectance for s-polarized (perpendicular) light
        perpendicular_polarized = True
        polarization_vector = perpendicular_v.normalize()
    else:
        # Parallel polarization
        a = (n1 * c2 - n2 * c1) / (n1 * c2 + n2 * c1)
        r = a * a.conjugate()  # reflectance for p-polarized (parallel) light
        perpendicular_polarized = False
        polarization_vector = parallel_v.normalize()
    if myrandom() < r.real:
        # ray reflected
        return reflexion(incident, normal, polarization_vector, True)
    else:
        # ray refracted
        refracted_direction = incident * r.real + \
                              my_normal * (r.real * c1.real - c2.real)
        perp_v = refracted_direction.cross(my_normal)
        if perp_v == Base.Vector(0, 0, 0):
            # to avoid null vector at my_normal and incident parallel vectors
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
    if mynormal.dot(incident) > 0:
        # Ray intercepted on the backside of the surface
        mynormal = mynormal * (-1.0)
    r = n1.real / n2.real
    c1 = - mynormal.dot(incident)  # cos (incidence_angle)
    c2sq = 1.0 - r * r * (1.0 - c1 * c1)  # cos (refracted_angle) ** 2
    c2 = c2sq ** 0.5  # cos (refracted_angle)
    refracted_direction = incident * r + mynormal * (r * c1 - c2)
    if perpendicular_polarized == 'thin film':
        # refraction in thin film
        return OpticalState(polarization_vector, refracted_direction, Phenomenon.REFRACTION)
    perp_v = refracted_direction.cross(mynormal)
    if perp_v == Base.Vector(0, 0, 0):
        # to avoid null vector at mynormal and incident parallel vectors
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
    if v_p == Base.Vector(0, 0, 0):
        # to avoid null vector at mynormal and incident parallel vectors
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
    if v_p == Base.Vector(0, 0, 0):
        # to avoid null vector at mynormal and incident parallel vectors
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


# ---
# Helper function for reflectance depending on the wavelength for coatings layers.
# ---

@traced(logger)
def matrix_reflectance(data_material):
    """
    Computes the matrix reflectance

    Parameters
    ----------
    data_material : array_like
        columns are: wavelenth in nm, angle in deg., reflectance s-polarized (perpendicular),
        reflectance p-polarized (parallel) the values should be in the corresponding order columns
        with constants steps

    Returns
    -------
        function
    """
    # TODO: document
    steps = np.diff(data_material[:, 1])
    try:
        a = min(steps[steps > 0.0])
    except ValueError:
        a = 1E-4
    steps = np.diff(data_material[:, 0])
    try:
        w = min(steps[steps > 0.0])
    except ValueError:
        w = 1E-4
    return lambda x, y: [row for row in data_material if
                         (x + a > row[1] > x - a) and (y + w > row[0] > y - w)]


@traced(logger)
def calculate_reflectance(m_reflectance, angle, wavelength):
    """Compute perpendicular and parallel components of reflectance

    Interpolates the value of the perperdicular and parallel reflectance from
    the matrix of values, depending on the angle and the wavelength

    Parameters
    ----------
    m_reflectance : float or tuple of list of floats
    angle : float
    wavelength : float

    Returns
    -------
    tuple of float or tuple of np.complex
    """
    if len(m_reflectance) == 0:  # wavelength experiments out of range of the data_material
        return  0.0, 0.0
    r_matrix = np.mat(m_reflectance)
    if len(r_matrix) == 1:
        # interpolation is not needed
        return r_matrix[0, 2], r_matrix[0, 3]

    if len(r_matrix) == 2:  # simple interpolation
        if r_matrix[0, 0] == r_matrix[1, 0]:
            # identical wavelengths; hence angle interpolation (column 1)
            x_values = [r_matrix[0, 1], r_matrix[1, 1]]
            y_values = [r_matrix[0, 2], r_matrix[1, 2]]
            r_per = np.interp(angle, x_values, y_values)

            x_values = [r_matrix[0, 1], r_matrix[1, 1]]
            y_values = [r_matrix[0, 3], r_matrix[1, 3]]
            r_par = np.interp(angle, x_values, y_values)
            return r_per, r_par

        if r_matrix[0, 1] == r_matrix[1, 1]:
            # identical angles; hence wavelength interpolation (column 0)
            x_values = [r_matrix[0, 0], r_matrix[1, 0]]
            y_values = [r_matrix[0, 2], r_matrix[1, 2]]
            r_per = np.interp(wavelength, x_values, y_values)

            x_values = [r_matrix[0, 0], r_matrix[1, 0]]
            y_values = [r_matrix[0, 3], r_matrix[1, 3]]
            r_par = np.interp(wavelength, x_values, y_values)
            return r_per, r_par

    if len(r_matrix) == 4:
        # double interpolation
        x_values = [r_matrix[0, 1], r_matrix[1, 1]]
        y_values = [r_matrix[0, 2], r_matrix[1, 2]]
        r1 = np.interp(angle, x_values, y_values)
        x_values = [r_matrix[2, 1], r_matrix[3, 1]]
        y_values = [r_matrix[2, 2], r_matrix[3, 2]]
        r2 = np.interp(angle, x_values, y_values)
        x_values = [r_matrix[0, 0], r_matrix[3, 0]]
        y_values = [r1, r2]
        r_per = np.interp(wavelength, x_values, y_values)

        x_values = [r_matrix[0, 1], r_matrix[1, 1]]
        y_values = [r_matrix[0, 3], r_matrix[1, 3]]
        r1 = np.interp(angle, x_values, y_values)
        x_values = [r_matrix[2, 1], r_matrix[3, 1]]
        y_values = [r_matrix[2, 3], r_matrix[3, 3]]
        r2 = np.interp(angle, x_values, y_values)
        x_values = [r_matrix[0, 0], r_matrix[3, 0]]
        y_values = [r1, r2]
        r_par = np.interp(wavelength, x_values, y_values)
        return r_per, r_par


