"""
Mathematical functions that are used as helper functions
"""

import numpy as np
from FreeCAD import Base
import random
import time
import collections

EPSILON = 1E-6 
# Tolerance for considering equal to zero
INF = 1E20
# Infinite

def polar_to_cartesian(phi, theta):
    """Convert polar coordinates of unit vector to cartesian

    Parameters
    ----------
    phi : float
        phi angle (ISO 31-11) in degrees
    theta : float
        theta angle (ISO 31-11) in degrees

    Returns
    -------
    Base.Vector

    """
    rad = np.pi / 180.0
    x = np.sin(theta * rad) * np.cos(phi * rad)
    y = np.sin(theta * rad) * np.sin(phi * rad)
    z = np.cos(theta * rad)
    return Base.Vector(x, y, z)


def rad_to_deg(angle):
    return angle * 180.0 / np.pi

# ---
# Helper functions for input of functions
# ---
def constant_function(c):
    """Create a constant function

    Parameters
    ----------
    c : float
        constant to return
    Returns
    -------
    function
        Constant function equal to `c`
    """
    return lambda x: c


def tabulated_function(xvalues, yvalues):
    """Create a linear interpolating function from tabulated values

    Parameters
    ----------
    xvalues : list of float
        x coordinates of the tabulated values
    yvalues : list of float
        y coordinates of the tabulated values

    Returns
    -------
    function
        Function that interpolates by straight line segments the input data
    """
    cached = {}
    def function_to_cache(x):
        if x in cached:
            return cached[x]
        value = np.interp(x, xvalues, yvalues)
        cached[x] = value
        return value
    return function_to_cache


# ---
# Helper function for random Linear congruential generator
# ---
_previous = None
def random_congruential(seed=None):
    """Random Linear congruential generator based on  MTH$RANDOM

    Parameters
    ----------
    seed : float
        seed to use in the generation of random numbers

    Returns
    -------
    float

    """
    # http://daviddeley.com/random/random4.htm
    a = 69069.0
    c = 1.0
    m = 2.0 ** 32.0
    rm = 1.0 / m
    global _previous
    if not seed:
        if not _previous:
            _previous = time.time()
    else:
        _previous = seed
    _previous = np.remainder(_previous * a + c, m)
    return _previous * rm


# ---
# Define the random algorithm
# ---
# myrandom = random_congruential
myrandom = random.random

# ---
# Helper function for Cumulative Function Distribution and Randomly generatting distribution
# ---

def cdf_from_pdf_file(data_file):
    """
    Computes CDF from PDF values stored in a file

    Creates a Cumulative Distribution Function from Probability Density
    Function data file. Each line must be a pair of numbers x y=pdf(x).
    It returns the CDF as two lists; first on is the list of x-values,
    second one is the list of corresponding CDF values.

    Parameters
    ----------

    data_file: file or str
        file or filename where PDF values are stored

    Returns
    -------
    list of float, list of float
        x-values and y-values of CDF
    """
    data_array = np.loadtxt(data_file, usecols=(0, 1))
    x = data_array[:, 0]
    y = data_array[:, 1]
    x_cdf = x
    n = np.size(y)
    y_ii = []
    for i in np.arange(n - 1):
        y_i = (y[i + 1] + y[i]) / 2.0 * (x[i + 1] - x[i])
        y_ii = np.append(y_ii, y_i)
    y_ii = np.append(y_ii, y_ii[-1])
    k_integration = np.trapz(y_ii, x_cdf)
    y_cdf = np.cumsum(y_ii) / k_integration
    return x_cdf, y_cdf / y_cdf[-1]


def pick_random_from_cdf(cdf):
    """
    Pick a random value according to a given CDF.

    We apply the Inverse transform sampling: https://en.wikipedia.org/wiki/Inverse_transform_sampling

    Parameters
    ----------
    cdf : tuple of list of float
        First list is list of x-values; second one is list of values of CDF

    Returns
    -------
    float
    """
    return np.interp(random.random(), cdf[1], cdf[0])

def parallel_orthogonal_components(vector, incident, normal):
    """Decomposition of vector in components

    Given `vector` (a polarization),
    `incident` (direction of a ray) and
    `normal` (vector orthogonal to a plane),
    decompose `vector` it in
    a component contained in the reflexion (parallel) plane (determined by normal and incident): p-polarized (parallel) light
    a component contained in the orthogonal plane to the reflexion plane: s-polarized (perpendicular) light
    also returns the normal vector to the reflexion plane

    Parameters
    ----------
    vector : Base.Vector
    incident : Base.Vector
    normal : Base.Vector

    Returns
    -------
    parallel : Base.Vector
    orthogonal : Base.Vector
    normal_of_parallel_plane: Base.Vector
    """
    polarization_vector = vector
    normal_parallel_plane = incident.cross(normal)
    # orthogonal vector to reflexion plane (parallel_plane)
    if normal_parallel_plane.Length < EPSILON:
        normal_parallel_plane = one_orthogonal_vector(normal)
    normal_parallel_plane.normalize()
    normal_perpendicular_plane = incident.cross(normal_parallel_plane)
    # orthogonal vector to perpendicular_plane
    parallel_v = polarization_vector - \
                 normal_parallel_plane * polarization_vector.dot(normal_parallel_plane)
    # parallel_v is the projection of polarization_vector onto parallel_plane
    perpendicular_v = polarization_vector - \
                 normal_perpendicular_plane * polarization_vector.dot(normal_perpendicular_plane)
    # perpendicular_v is the projection of polarization_vector onto the perpendicular_plane
    return parallel_v, perpendicular_v, normal_parallel_plane
	
def two_orthogonal_vectors(vector):
    """Gives two orthogonal vectors of a vector

    Given `vector` find two orthogonal vectors

    Parameters
    ----------
    vector : Base.Vector

    Returns
    -------
    orthogonal_1 : Base.Vector
    orthogonal_2 : Base.Vector
    """
    orthogonal_1 = one_orthogonal_vector(vector)
    orthogonal_2 = vector.cross(orthogonal_1)
    return orthogonal_1.normalize(), orthogonal_2.normalize()

def one_orthogonal_vector(vector):
    """Gives one orthogonal vector of a vector

    Given `vector` find one orthogonal vector

    Parameters
    ----------
    vector : Base.Vector

    Returns
    -------
    orthogonal : Base.Vector
    """
    min_pos = np.argmin([abs(vector[0]), abs(vector[1]), abs(vector[2])])
    if min_pos == 0:
        orthogonal = Base.Vector(0, vector[2], -vector[1])
    elif min_pos == 1:
        orthogonal = Base.Vector(vector[2], 0, -vector[0])
    else:
        orthogonal = Base.Vector(vector[1], -vector[0], 0)
    return orthogonal.normalize()

def correct_normal(normal, incident):
    if normal.dot(incident) > 0:
        return normal * (-1)
    else:
        return normal

def normalize(vector):
    if vector.Length < EPSILON:
        vector = vector * INF
    return vector.normalize()
	
def arccos(x):
    assert abs(x) < 1 + EPSILON
    if abs(x) < 1 - EPSILON:
        return np.arccos(x)
    if x > 0:
        return 0
    return np.pi

def buie_distribution(CircumSolarRatio):
    """
    Implementation of the Buie Distribution for Sun emission

    A distribution for the direct light from the sun according to the Buie model:
    Buie, D., 2005. Corrigendum to "The effective size of the solar cone for
    solar concentrating systems" [Solar Energy 74 (2003) 417-427].
    Solar Energy 79, 568-570. 5.2

    Parameters
    ----------
    CircumSolarRatio : float

    Returns
    -------
    angle distribution for random input: function
        Function that interpolates by straight line segments the input data
    """
    def calculate_a1(CSR, SD):
        """ Parameter a1 needed for the normalization of the probability distribution in the disk region
        #TODO: Say one sentence @Ramon

        Parameters
        ----------
        CSR
        SD

        Returns
        -------

        """

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
        f = solar_disk__density_distribution
        th = np.arange(0.0, SD, 0.001)
        f_th = np.vectorize(f)
        a1 = (1.0 - CSR) / np.trapz(f_th(th), dx=0.001)
        return a1

    def calculate_a2(CSR, SD, SS):
        """ Parameter a2 needed for the normalization of the probability distribution in the circumsolar region"""
        def circumsolar__density_distribution(angle, CSR):
            gamma = 2.2 * np.log(0.52 * CSR) * CSR ** (0.43) - 0.1
            kappa = 0.9 * np.log(13.5 * CSR) * CSR ** (-0.3)
            return 2.0 * np.pi * np.exp(kappa) * angle ** (gamma + 1.0)
        f = circumsolar__density_distribution
        th = np.arange(SD, SS, 0.001)
        f_th = np.vectorize(f)
        a2 = CSR / np.trapz(f_th(th, CSR), dx=0.001)
        return a2

    def calculate_CDF_disk_region(a1, SD):
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
        f = solar_disk__density_distribution
        th = np.arange(0.0, SD, 0.001)
        f_th = np.vectorize(f)
        cumulative = np.cumsum(f_th(th))
        CDF = (th, a1 * cumulative / 1000.0)
        return CDF

    CSR = CircumSolarRatio
    SD = 4.65
    # Solar Disk in mrad
    SS = 43.6
    # Solar Size in mrad
    a1 = calculate_a1(CSR, SD)
#     normalization constant for the disk region
    a2 = calculate_a2(CSR, SD, SS)
#    normalization constant for the circumsolar region
    CDF_Disk_Region = calculate_CDF_disk_region(a1,SD)
#    Buie distribution for the disk region

    def th_solar_disk_region(u, a1, SD):
        """ Random angle based on the probability distribution in the circumsolar region"""
        # TODO: It is recomputed every time? @Ramon
        CDF = calculate_CDF_disk_region(a1, SD)
        idx = (np.abs(CDF[1] - u)).argmin()
        return CDF[0][idx]

    def th_circumsolar_region(u, CSR, SD, a2):
        """ Random angle based on the CDF in the circumsolar region"""
        gamma = 2.2 * np.log(0.52 * CSR) * CSR ** (0.43) - 0.1
        kappa = 0.9 * np.log(13.5 * CSR) * CSR ** (-0.3)
        u_csr = (u - 1.0) + CSR  # Since CSR-CDF starts at zero
        f1 = u_csr * (gamma + 2.0) / (a2 * 2 * np.pi * np.exp(kappa))
        f2 = SD ** (gamma + 2.0)
        exp = (1.0 / (gamma + 2.0))
        th_u = np.power(f1 + f2, exp)
        return th_u
    random_array = np.arange(0.0, 1.001, 0.001)
    ang_distr = {}
    for u in random_array:
        if u < 1.0 - CSR:
            ang_distr[u] = th_solar_disk_region(u, a1, SD) / 1000.0 * 180.0 / np.pi
        else:
            ang_distr[u] = th_circumsolar_region(u, CSR, SD, a2) / 1000.0 * 180.0 / np.pi
    ang_distribution = collections.OrderedDict(sorted(ang_distr.items()))
    f = tabulated_function(ang_distribution.keys(), ang_distribution.values())
    return f
