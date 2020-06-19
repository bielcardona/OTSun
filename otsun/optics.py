"""Module otsun.optics

Implementation of optical effects on rays
"""

from FreeCAD import Base
import numpy as np
from .math import arccos, myrandom, one_orthogonal_vector, correct_normal, \
    parallel_orthogonal_components, rad_to_deg, normalize
from enum import Enum
from numpy.lib.scimath import sqrt
from autologging import traced
from .logging_unit import logger

# from .materials import Material, vacuum_medium

EPSILON = 1E-6
# Tolerance for considering equal to zero
INF = 1E20



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
    ENERGY_ABSORBED = 5


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
    solid : Shape
        solid where the ray is located
    material : Material
        Material where the ray is located
    extra_data : dict
        Dictionary where materials can put extra data
    """

    def __init__(self, polarization, direction, phenomenon, material=None, extra_data=None):
        self.polarization = polarization
        self.direction = direction
        self.phenomenon = phenomenon
        self.material = material  # TODO: Set solid
        if extra_data is None:
            extra_data = {}
        self.extra_data = extra_data

    def __str__(self):
        return "Dir.: %s, Pol.: %s, Phen.: %s" % (
            self.direction, self.polarization, self.phenomenon
        )

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

    def apply_dispersion(self, properties, normal_vector):
        if properties.get('sigma_1', None):
            sigma_1 = properties['sigma_1']
            if properties.get('sigma_2', None):
                sigma_2 = properties['sigma_2']
                k = properties.get('k', None) or 0.5
                self.apply_double_gaussian_dispersion(
                    normal_vector, sigma_1, sigma_2, k)
            else:
                self.apply_single_gaussian_dispersion(normal_vector, sigma_1)


# ---
# Helper functions for reflection and refraction
# ---


@traced(logger)
def simple_reflection(incident, normal):
    normal = correct_normal(normal, incident)
    c1 = - normal.dot(incident)
    reflected = incident + normal * 2.0 * c1
    return reflected.normalize()


@traced(logger)
def simple_polarization_reflection(incident, normal, normal_parallel_plane, polarization):
    c1 = - normal.dot(incident)
    angle = rad_to_deg(np.pi - 2.0 * arccos(c1))
    rotation = Base.Rotation(normal_parallel_plane, angle)
    return rotation.multVec(polarization)


@traced(logger)
def simple_polarization_refraction(incident, normal, normal_parallel_plane, c2, polarization_vector):
    c1 = - normal.dot(incident)
    angle = rad_to_deg(arccos(c2.real) - arccos(c1))
    rotation = Base.Rotation(normal_parallel_plane, angle)
    return rotation.multVec(polarization_vector)


@traced(logger)
def reflection(incident, normal_vector, polarization_vector,
               polarization_vector_calculated_before=False):
    """
    Implementation of law of reflection for incident and polarization vector.

    Parameters
    ----------
    incident : Base.Vector
        vector incident
    normal_vector : Base.Vector
        vector normal to the surface
    polarization_vector: Base.Vector
        Polarization vector of the ray
    polarization_vector_calculated_before : bool
        True if the polarization vector was computed before

    Returns
    -------
    OpticalState
        optical state of the reflected ray
    """
    normal = correct_normal(normal_vector, incident)
    reflected = simple_reflection(incident, normal)
    # reflection changes the polarization vector
    if not polarization_vector_calculated_before:
        # we calculate the new polarization vector
        parallel_v, perpendicular_v, normal_parallel_plane = \
            parallel_orthogonal_components(polarization_vector, incident, normal)
        polarization_vector = simple_polarization_reflection(incident, normal, normal_parallel_plane,
                                                             polarization_vector)
    return OpticalState(polarization_vector, reflected, Phenomenon.REFLEXION)  # TODO: Set solid


@traced(logger)
def lambertian_reflection(incident, normal_vector):
    """
    Implementation of lambertian reflection for diffusely reflecting surface

    Parameters
    ----------
    incident : Base.Vector
        vector incident
    normal_vector : Base.Vector
        vector normal to the surface

    Returns
    -------
    OpticalState
        optical state of the reflected ray
    """
    normal = correct_normal(normal_vector, incident)
    dot = - 1.0
    while dot <= 0.01:
        random_vector = Base.Vector(myrandom() - 0.5, myrandom() - 0.5, myrandom() - 0.5)
        if random_vector.Length < 1:
            random_vector.normalize()
            dot = normal.dot(random_vector)
    new_direction = random_vector
    random_polarization_vector = random_polarization(new_direction)
    return OpticalState(random_polarization_vector,
                        new_direction, Phenomenon.REFLEXION)  # TODO: Set solid


@traced(logger)
def refraction(incident, normal_vector, n1, n2, polarization_vector, lambertian_surface=False):
    """Implementation of Fresnel equations of refraction

    Parameters
    ----------
    incident : Base.Vector
        direction vector of the incident ray	
    normal_vector: Base.Vector
        normal vector of the surface at the point of incidence
    n1: complex
        complex refractive index where ray is currently traveling
    n2: complex
        complex refractive index of nearby material
    polarization_vector: Base.Vector
        Polarization vector of the ray
    lambertian_surface: Bool
        Indicates if the surface has lambertian reflection

    Returns
    -------
    OpticalState
    """
    # Fresnel equations. Oblique incident in absorbing media.
    # See Chapter 2 of the book Thin-films optical filters
    # See Fresnel equations in WikipediA
    normal = correct_normal(normal_vector, incident)
    r = n1 / n2
    c1 = - normal.dot(incident)
    # cos (incident_angle)
    c2sq = 1.0 - r * r * (1.0 - c1 * c1)
    # cos (refracted_angle) ** 2
    if c2sq.real < 0:
        # total internal reflection
        if not lambertian_surface:
            return reflection(incident, normal, polarization_vector)
        else:
            return lambertian_reflection(incident, normal)
    c2 = sqrt(c2sq)
    # cos (refracted_angle)
    parallel_v, perpendicular_v, normal_parallel_plane = parallel_orthogonal_components(polarization_vector, incident,
                                                                                        normal)
    # parallel and perpendicular components of polarization vector and orthogonal vector of the parallel plane
    ref_per = perpendicular_v.Length ** 2.0 / polarization_vector.Length ** 2.0
    # weight of perpendicular component: 0 < ref_per < 1
    # We decide the polarization projection onto the parallel / perpendicular plane
    if myrandom() < ref_per:
        # reflectance for s-polarized (perpendicular) light
        a = (n1 * c1 - n2 * c2) / (n1 * c1 + n2 * c2)
        reflectance = a * a.conjugate()
        perpendicular_polarized = True
        polarization_vector = normalize(perpendicular_v)
    else:
        # reflectance for p-polarized (parallel) light
        a = (n1 * c2 - n2 * c1) / (n1 * c2 + n2 * c1)
        reflectance = a * a.conjugate()
        perpendicular_polarized = False
        polarization_vector = normalize(parallel_v)
    # The ray can be reflected or refracted 
    if myrandom() < reflectance.real:
        # ray reflected
        if not lambertian_surface:
            reflected_direction = simple_reflection(incident, normal)
            if not perpendicular_polarized:
                # reflection changes the parallel component of incident polarization
                polarization_vector = simple_polarization_reflection(incident, normal, normal_parallel_plane,
                                                                     polarization_vector)
            return OpticalState(polarization_vector, reflected_direction,
                                Phenomenon.REFLEXION)  # TODO: Set solid
        else:
            return lambertian_reflection(incident, normal)
    else:
        # ray refracted: computing the refracted direction
        if c2.real > 1:
            # avoiding invalid solutions for metallic materials
            c2 = 1
        refracted_direction = incident * r.real + \
                              normal * (r.real * c1 - c2.real)
        refracted_direction.normalize()
        if not perpendicular_polarized:
            # refraction changes the parallel component of incident polarization
            polarization_vector = simple_polarization_refraction(incident, normal, normal_parallel_plane, c2,
                                                                 polarization_vector)
        return OpticalState(polarization_vector, refracted_direction,
                            Phenomenon.REFRACTION)  # TODO: Set solid


@traced(logger)
def shure_refraction(incident, normal_vector, n1, n2, polarization_vector, lambertian_surface=False):
    """Implementation of Snell's law of refraction

    Parameters
    ----------
    incident : Base.Vector
        direction vector of the incident ray
    normal_vector: Base.Vector
        normal vector of the surface at the point of incidence
    n1: complex
        complex refractive index where ray is currently traveling
    n2: complex
        complex refractive index of nearby material
    polarization_vector: Base.Vector
        Polarization vector of the ray
    lambertian_surface: Bool
        Indicates if the surface has lambertian reflection

    Returns
    -------
    OpticalState
    """
    # TODO: document

    normal = correct_normal(normal_vector, incident)
    r = n1 / n2
    c1 = - normal.dot(incident)
    # cos (incident_angle)
    c2sq = 1.0 - r * r * (1.0 - c1 * c1)
    # cos (refracted_angle) ** 2
    if c2sq.real < 0:
        # total internal reflection
        if not lambertian_surface:
            return reflection(incident, normal, polarization_vector)
        else:
            return lambertian_reflection(incident, normal)
    c2 = sqrt(c2sq)
    # cos (refracted_angle)
    if c2.real > 1:
        # avoiding invalid solutions for metallic materials
        c2 = 1
    parallel_v, perpendicular_v, normal_parallel_plane = parallel_orthogonal_components(polarization_vector, incident,
                                                                                        normal)
    # parallel and perpendicular components of polarization vector and orthogonal vector of the parallel plane
    # ray refracted: computing the refracted direction
    refracted_direction = incident * r.real + normal * (r.real * c1 - c2.real)
    refracted_direction.normalize()
    polarization_vector = simple_polarization_refraction(incident, normal, normal_parallel_plane, c2,
                                                         polarization_vector)
    return OpticalState(polarization_vector, refracted_direction,
                        Phenomenon.REFRACTION)  # TODO: Set solid


# ---
# Helper function for dispersions and polarization vector
# ---
@traced(logger)
def dispersion_from_main_direction(main_direction, theta, phi):
    """
    Computes dispersion from the main direction in terms of angles theta and phi
    """
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
    Computes dispersion of polarization vector in terms of angles theta and phi
    """
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
    Returns a random polarization orthogonal to the given direction
    """
    orthogonal_vector = one_orthogonal_vector(direction)
    phi = 360.0 * myrandom()
    rotation = Base.Rotation(direction, phi)
    random_polarization_v = rotation.multVec(orthogonal_vector)
    return random_polarization_v


# ---
# Helper function for reflectance depending on the wavelength for coatings layers.
# ---

def _round_or_floor_ceil(x):
    xround = int(round(x))
    if abs(x - xround) < 1E-4:
        return [xround]
    xfloor = int(np.floor(x))
    return [xfloor, xfloor + 1]


@traced(logger)
def matrix_reflectance(data_material):
    """
    Computes the matrix of reflectances of a material

    data_material: wavelenth in nm, angle in deg., reflectance s-polarized (perpendicular),
    reflectance p-polarized (parallel)
    the values should be in the corresponding order columns with constants steps
    We generate a matrix for the interpolation depending on the x=angle and
    y=wavelength values (lambda values)

    Parameters
    ----------
    data_material
    """
    steps = np.diff(data_material[:, 1])
    try:
        delta_a = min(steps[steps > 0.0])
    except ValueError:
        delta_a = 1E-4
    min_a = min(data_material[:, 1])
    steps = np.diff(data_material[:, 0])
    try:
        delta_w = min(steps[steps > 0.0])
    except ValueError:
        delta_w = 1E-4
    min_w = min(data_material[:, 0])
    data_dict = {}
    for row in data_material:
        data_dict[(int(round((row[0] - min_w) / delta_w)),
                   int(round((row[1] - min_a) / delta_a)))] = row

    #    @memoize
    # @lru_cache(maxsize=None)
    def internal_matrix_reflectance(x, y):
        index_a = (x - min_a) / delta_a
        index_w = (y - min_w) / delta_w
        index_a_ints = _round_or_floor_ceil(index_a)
        index_w_ints = _round_or_floor_ceil(index_w)
        rows = []
        for i_a in index_a_ints:
            for i_w in index_w_ints:
                try:
                    rows.append(data_dict[(i_w, i_a)])
                except KeyError:
                    pass
        return rows

    return internal_matrix_reflectance


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
        return 0.0, 0.0
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
