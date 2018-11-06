import numpy as np
from FreeCAD import Base
import random

def polar_to_cartesian(phi, theta):
    """
    Convert polar coordinates (given in degrees) to cartesian
    :param phi:
    :param theta:
    :return:
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
    return lambda x: c


def tabulated_function(xvalues, yvalues):
    return lambda x: np.interp(x, xvalues, yvalues)


# ---
# Helper function for random Linear congruential generator
# ---
def random_congruential(seed=None):
    """
    Implementation of a random Linear congruential generator based on the MTH$RANDOM
    """
    # http://daviddeley.com/random/random4.htm
    a = 69069.0
    c = 1.0
    m = 2.0 ** 32.0
    rm = 1.0 / m
    if seed is not None:
        random_congruential.previous = seed
    seed = np.remainder(random_congruential.previous * a + c, m)
    random_congruential.previous = seed
    return seed * rm


# ---
# Define the random algorithm
# ---
myrandom = random_congruential


# myrandom = random.random

# ---
# Helper function for Cumulative Function Distribution and Randomly generatting distribution
# ---

def create_CDF_from_PDF(data_file):
    """
    Creates a Cumulative Distribution Function from Probability Density Function data file (x,y)
    data_file: x,y values
    output: x_CDF, y_CDF
    """
    data_array = np.loadtxt(data_file, usecols=(0, 1))
    x = data_array[:, 0]
    y = data_array[:, 1]
    x_CDF = x
    n = np.size(y)
    y_ii = []
    for i in np.arange(n - 1):
        y_i = (y[i + 1] + y[i]) / 2.0 * (x[i + 1] - x[i])
        y_ii = np.append(y_ii, y_i)
    y_ii = np.append(y_ii, y_ii[-1])
    k_integration = np.trapz(y_ii, x_CDF)
    y_CDF = np.cumsum(y_ii) / k_integration
    return x_CDF, y_CDF / y_CDF[-1]


def pick_random_from_CDF(cumulative_distribution_function):
    """
    Randomly generatting distribution acording to a Cumulative Distribution Function
    We apply the Inverse transform sampling: https://en.wikipedia.org/wiki/Inverse_transform_sampling
    """
    CDF = cumulative_distribution_function
    return np.interp(random.random(), CDF[1], CDF[0])
