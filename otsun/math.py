"""
Mathematical functions that are used as helper functions
"""

import numpy as np
from FreeCAD import Base
import random
import time

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
    return lambda x: np.interp(x, xvalues, yvalues)


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
    output: list of float, list of float
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
    a component parallel to the plane (and orthogonal to incident)
    a component contained in the plane determined by normal and incident

    Parameters
    ----------
    vector : Base.Vector
    incident : Base.Vector
    normal : Base.Vector

    Returns
    -------
    parallel : Base.Vector
    orthogonal : Base.Vector
    """
    normal_parallel_plane = incident.cross(normal)
    # normal vector of the parallel plane
    if normal_parallel_plane == Base.Vector(0, 0, 0):
        # to avoid null vector at mynormal and incident parallel vectors
        normal_parallel_plane = Base.Vector(1, 0, 0)
    normal_parallel_plane.normalize()
    normal_perpendicular_plane = normal_parallel_plane.cross(incident)
    # normal vector of the perpendicular plane
    parallel_v = vector - normal_parallel_plane * \
                 vector.dot(normal_parallel_plane)
    perpendicular_v = vector - \
                      normal_perpendicular_plane * \
                      vector.dot(normal_perpendicular_plane)
    return parallel_v, perpendicular_v
