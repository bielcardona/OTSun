# -*- coding: utf-8 -*-
"""
Created on Fri Nov 25 15:10:31 2016

"""

from FreeCAD import Base
import Part
import numpy as np
import itertools
import random
import time
import dill
import zipfile
from enum import Enum
from numpy.lib.scimath import sqrt

# ---
# Helper functions for outputs in Spectral Analysis
# ---

def spectrum_to_constant_step(file_in, wavelength_step, wavelength_min, wavelength_max):
    data_array = np.loadtxt(file_in, usecols=(0,1))
    wl_spectrum = data_array[:,0]
    I_spectrum = data_array[:,1]
    array_inter = [[x, np.interp(x,wl_spectrum,I_spectrum)] for x in np.arange(wavelength_min, wavelength_max + wavelength_step / 2.0, wavelength_step)]
    return np.array(array_inter)	

def make_histogram_from_experiment_results(reults_wavelength, resutls_energy, step_wavelength, aperture_collector, aperture_source):
    y_energy = np.array(np.concatenate(resutls_energy))
    y_energy = (y_energy / aperture_collector) / (1.0 / aperture_source)
    x_wavelength = np.array(np.concatenate(reults_wavelength))
    data_= np.array([x_wavelength, y_energy])
    data_ = data_.T
    bins_ = np.arange(int(np.amin(x_wavelength)), np.amax(x_wavelength)+step_wavelength*1.1, step_wavelength)
    table_ = np.histogram(data_[:,0], bins = bins_, weights = data_[:,1])
    norm = np.histogram(data_[:,0], bins = bins_)
    u = np.divide(table_[0],norm[0])
    bins = np.arange(int(np.amin(x_wavelength)), np.amax(x_wavelength)+step_wavelength, step_wavelength)
    table_ = np.column_stack((bins, u))
    return table_
	
def twoD_array_to_constant_step(twoD_array, step, wavelength_min, wavelength_max):
    wl_spectrum = twoD_array[:,0]
    I_spectrum = twoD_array[:,1]
    array_inter = [[x, np.interp(x,wl_spectrum,I_spectrum)] for x in np.arange(wavelength_min,wavelength_max+step/2.0,step)]
    return np.array(array_inter)
	
def spectral_response(optical_absortion_wavelength, iqe):
    q_e = 1.60217662E-19
    h = 6.62607E-34
    c = 299792458.0
    hc = h*c
    SR = []
    opt_wavelength = optical_absortion_wavelength	
    if np.isreal(iqe):
        SR = [[opt[0], iqe * opt[0] * opt[1] * q_e * 1E-9/ hc, ] for opt in opt_wavelength]
    else:
        data_array = np.loadtxt(iqe, usecols=(0,1))
        wl_ = data_array[:,0]
        iqe_ = data_array[:,1]        
        SR = [[opt[0], np.interp(opt[0],wl_,iqe_) * opt[0] * opt[1] * q_e * 1E-9/ hc, ] for opt in opt_wavelength]		
    return np.array(SR)

def photo_current(spectral_response, source_spectrum):
    wavelengths = source_spectrum[:,0]
    SR_by_spectrum = spectral_response[:,1] * source_spectrum[:,1]
    photo_current = np.trapz(SR_by_spectrum, x = wavelengths)
    return photo_current
	
# end region

# ---
# Helper functions for outputs in Total Analysis
# ---

def integral_from_data_file(file_in):
    source_spectrum = np.loadtxt(file_in, usecols=(0,1))
    integral = np.trapz(source_spectrum[:,1], x = source_spectrum[:,0])
    return integral	

# end region

# Base class for optical phonomena

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
def reflexion(incident, normal, polarization_vector, polarization_vector_calculated_before = False):
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
        normal_parallel_plane = incident.cross(mynormal)  # axis for the rotation, only parallel component must be rotated
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
    Tesfamichael, T., and Wäckelgard, E., 2000, “Angular Solar Absorptance and
    Incident Angle Modifier of Selective Absorbers for Solar Thermal Collectors,”
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
    c2 = sqrt( 1.0 - r * r * (1.0 - c1 * c1) )  # cos (refracted_angle)

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
        R = a * a.conjugate() # reflectance for s-polarized (perpendicular) light
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
        R = a * a.conjugate() # reflectance for s-polarized (perpendicular) light
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
    if mynormal.dot(incident) > 0: # Ray intercepted on the backside of the surface
        # noinspection PyAugmentAssignment
        mynormal = mynormal * (-1.0)
        backside = True
    r = n1 / n2
    c1 = - mynormal.dot(incident) # cos (incidence_angle) 
    c2sq = 1.0 - r * r * (1.0 - c1 * c1) # cos (refracted_angle) ** 2
    if properties['transparent_material']: # transparent coating
        if c2sq.real < 0: # total internal reflection
            return reflexion(incident, normal, polarization_vector)
    c2 = sqrt(c2sq) # cos (refracted_angle)
    normal_parallel_plane = incident.cross(mynormal) # normal vector of the parallel plane
    if normal_parallel_plane == Base.Vector(0,0,0): # to avoid null vector at mynormal and incident parallel vectors
        normal_parallel_plane = Base.Vector(1,0,0)
    normal_parallel_plane.normalize()
    normal_perpendicular_plane = normal_parallel_plane.cross(incident) # normal vector of the perpendicular plane 
    # http://www.maplesoft.com/support/help/Maple/view.aspx?path=MathApps/ProjectionOfVectorOntoPlane
    parallel_v = polarization_vector - normal_parallel_plane * polarization_vector.dot(normal_parallel_plane)
    parallel_component = parallel_v.Length
    perpendicular_v = polarization_vector - normal_perpendicular_plane * polarization_vector.dot(normal_perpendicular_plane)
    perpendicular_component = perpendicular_v.Length
    ref_per = perpendicular_component /(perpendicular_component + parallel_component)
    perpendicular_polarized = False
    # https://en.wikipedia.org/wiki/Fresnel_equations # Fresnel equations

    if backside == True and properties['transparent_material']:  # Ray intercepted on the backside of the surface
        angle = np.arccos(c2.real) * 180.0 / np.pi	
    else:
        angle = np.arccos(c1) * 180.0 / np.pi	    
    Matrix_Reflectance = properties['Matrix_polarized_reflectance_coating']
    Matrix_R = Matrix_Reflectance(angle,wavelength)
    if myrandom() < ref_per:
        R = calculate_reflectance(Matrix_R, angle, wavelength)[1] # reflectance for s-polarized (perpendicular) light
        perpendicular_polarized = True
        polarization_vector = perpendicular_v.normalize()
    else:
        angle = np.arccos(c1) * 180.0 / np.pi
        R = calculate_reflectance(Matrix_R, angle, wavelength)[3]  # reflectance for p-polarized (parallel) light
        polarization_vector = parallel_v.normalize()
    if myrandom() < R: # ray reflected    
        return 1, 0, 0, polarization_vector, perpendicular_polarized
    else: # ray refracted or absorbed		
        if properties['energy_collector']: # absorber coating
            return 0, 1, 0, polarization_vector, perpendicular_polarized
        if properties['specular_material']: # reflector coating
                return 0, 1, 0, polarization_vector, perpendicular_polarized
        if properties['transparent_material']: # transparent coating
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
    if perpendicular_polarized == 'thin film': # refraction in thin film
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
    if mynormal.dot(incident) > 0: # Ray intercepted on the backside of the surface
        mynormal = mynormal * (-1.0)
        backside = True
    r = n1.real / n2.real
    c1 = - mynormal.dot(incident) # cos (incidence_angle) 
    c2sq = 1.0 - r * r * (1.0 - c1 * c1) # cos (refracted_angle) ** 2
    if c2sq.real < 0: # total internal reflection
        return 0.0, reflexion(incident, normal, polarization_vector)
    c2 = sqrt(c2sq) # cos (refracted_angle)

    normal_parallel_plane = incident.cross(mynormal) # normal vector of the parallel plane
    if normal_parallel_plane == Base.Vector(0,0,0): # to avoid null vector at mynormal and incident parallel vectors
        normal_parallel_plane = Base.Vector(1,0,0)
    normal_parallel_plane.normalize()
    normal_perpendicular_plane = normal_parallel_plane.cross(incident) # normal vector of the perpendicular plane 
    # http://www.maplesoft.com/support/help/Maple/view.aspx?path=MathApps/ProjectionOfVectorOntoPlane
    parallel_v = polarization_vector - normal_parallel_plane * polarization_vector.dot(normal_parallel_plane)
    parallel_component = parallel_v.Length
    perpendicular_v = polarization_vector - normal_perpendicular_plane * polarization_vector.dot(normal_perpendicular_plane)
    perpendicular_component = perpendicular_v.Length
    ref_per = perpendicular_component /(perpendicular_component + parallel_component)
    perpendicular_polarized = False
    # https://en.wikipedia.org/wiki/Fresnel_equations # Fresnel equations

    if backside == True:  # Ray intercepted on the backside of the transparent surface
        angle = np.arccos(c2.real) * 180.0 / np.pi	
    else:
        angle = np.arccos(c1) * 180.0 / np.pi	    
    Matrix_Reflectance = properties['Matrix_reflectance_thin_film']
    Matrix_R = Matrix_Reflectance(angle,wavelength)
    if myrandom() < ref_per:
        R = calculate_reflectance(Matrix_R, angle, wavelength)[1] # reflectance for s-polarized (perpendicular) light
        perpendicular_polarized = True
        polarization_vector = perpendicular_v.normalize()
    else:
        angle = np.arccos(c1) * 180.0 / np.pi
        R = calculate_reflectance(Matrix_R, angle, wavelength)[3]  # reflectance for p-polarized (parallel) light
        polarization_vector = parallel_v.normalize()
    if myrandom() < R: # ray reflected
        return 0.0, reflexion(incident, normal, polarization_vector)
    else:
        Matrix_Transmittance = properties['Matrix_transmittance_thin_film']
        Matrix_T = Matrix_Transmittance(angle,wavelength)
        if perpendicular_polarized:
            T = calculate_reflectance(Matrix_T, angle, wavelength)[1]
        else:
            T = calculate_reflectance(Matrix_T, angle, wavelength)[3]
        energy_absorbed_thin_film = (1 - R - T) / (1 - R)
        refracted_direction = incident * r.real + mynormal * (r.real * c1.real - c2.real)
        return energy_absorbed_thin_film, OpticalState(polarization_vector, refracted_direction, Phenomenon.REFRACTION)

		
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


# ---
# Helper function for reflectance depending on the wavelength for coatings layers.
# ---	

def matrix_reflectance(data_material):
    """
    data_material: wavelenth in nm, angle in deg., reflectance s-polarized (perpendicular), reflectance p-polarized (parallel)
    the values should be in the corresponding order columns with constants steps
	We generate a matrix for the interpolation depending on the x=angle and y=wavelength values (lambda values) 
    """
    steps = np.diff(data_material[:,1])
    try:
        a = min(steps[steps>0.0])
    except:
        a = 1E-4
    steps = np.diff(data_material[:,0])
    try:
        w = min(steps[steps>0.0])
    except:
        w = 1E-4
    return lambda x,y: [row for row in data_material if (x+a > row[1] and x-a < row[1]) and (y+w > row[0] and y-w < row[0])]


def calculate_reflectance(matrix_reflectance, angle, wavelength):
    """
    data_array: values needed for the interpolation 
	We generate an interpolation to find the perpendicular and parallel reflectance depending on the angle and wavelength
    """
    if len(matrix_reflectance) == 0: # wavelength experiments out of range of the data_material
        return 'R_per',0.0,'R_par',0.0
    R_Matrix = np.mat(matrix_reflectance)
    if len(R_Matrix) == 1: # interpolation is not needed
        return 'R_per',R_Matrix[0,2],'R_par',R_Matrix[0,3]


    if len(R_Matrix) == 2: # simple interpolation
        if R_Matrix[0,0] == R_Matrix[1,0]: # identical wavelengths; hence angle interpolation (column 1)
            xvalues = [R_Matrix[0,1],R_Matrix[1,1]]
            yvalues = [R_Matrix[0,2],R_Matrix[1,2]]
            R_per = np.interp(angle, xvalues, yvalues)

            xvalues = [R_Matrix[0,1],R_Matrix[1,1]]
            yvalues = [R_Matrix[0,3],R_Matrix[1,3]]
            R_par = np.interp(angle, xvalues, yvalues)
            return 'R_per',R_per,'R_par',R_par
			
        if R_Matrix[0,1] == R_Matrix[1,1]: # identical angles; hence wavelength interpolation (column 0)
            xvalues = [R_Matrix[0,0],R_Matrix[1,0]]
            yvalues = [R_Matrix[0,2],R_Matrix[1,2]]
            R_per = np.interp(wavelength, xvalues, yvalues)

            xvalues = [R_Matrix[0,0],R_Matrix[1,0]]
            yvalues = [R_Matrix[0,3],R_Matrix[1,3]]
            R_par = np.interp(wavelength, xvalues, yvalues)
            return 'R_per',R_per,'R_par',R_par

			
    if len(R_Matrix) == 4: # double interpolation
        xvalues = [R_Matrix[0,1],R_Matrix[1,1]]
        yvalues = [R_Matrix[0,2],R_Matrix[1,2]]
        R1 = np.interp(angle, xvalues, yvalues)
        xvalues = [R_Matrix[2,1],R_Matrix[3,1]]
        yvalues = [R_Matrix[2,2],R_Matrix[3,2]]
        R2 = np.interp(angle, xvalues, yvalues)
        xvalues = [R_Matrix[0,0],R_Matrix[3,0]]
        yvalues = [R1,R2]
        R_per = np.interp(wavelength, xvalues, yvalues)
        
        xvalues = [R_Matrix[0,1],R_Matrix[1,1]]
        yvalues = [R_Matrix[0,3],R_Matrix[1,3]]
        R1 = np.interp(angle, xvalues, yvalues)
        xvalues = [R_Matrix[2,1],R_Matrix[3,1]]
        yvalues = [R_Matrix[2,3],R_Matrix[3,3]]
        R2 = np.interp(angle, xvalues, yvalues)
        xvalues = [R_Matrix[0,0],R_Matrix[3,0]]
        yvalues = [R1,R2]
        R_par = np.interp(wavelength, xvalues, yvalues)
        return 'R_per',R_per,'R_par',R_par
	
	
# ---
# Classes for materials
# ---

# region Materials

class Material(object):
    """
    Class used to represent materials and their physical properties

    Attributes
    ----------
    name : str
        Holds the name of the material
    properties : dict
        Dictionary with physical properties of the material
    """
    by_name = {}

    def __init__(self, name, properties):
        self.by_name[name] = self
        self.name = name
        self.properties = properties

    @classmethod
    def get_from_label(cls, label):
        if ("(" not in label) or (")" not in label):
            return None
        start = label.find("(")
        end = label.find(")")
        name = label[start + 1:end]
        return cls.by_name.get(name, None)

    @classmethod
    def create(cls, name, properties):
        _ = cls(name, properties)

    @classmethod
    def load_from_file(cls, filename):
        try:
            with open(filename,'rb') as f:
                mat = dill.load(f)
                cls.by_name[mat.name] = mat
                return mat.name
        except:
            pass

    @classmethod
    def load_from_zipfile(cls, filename):
        with zipfile.ZipFile(filename) as z:
            for matfile in z.namelist():
                with z.open(matfile) as f:
                    try:
                        mat = dill.load(f)
                        cls.by_name[mat.name] = mat
                    except:
                        pass

    def change_of_direction(self, ray, normal_vector):
        """
        Computes the karma of the material


        """
        pass
        # Compute new direction (ray.current_material)


class VolumeMaterial(Material, object):
    def __init__(self, name, properties):
        """
        Initializes a Volume Material. The properties parameter must be a dict with the physical properties
        describing the material. At least, the following must be provided:
        'index_of_refraction': index of refraction of the material, as a function of its wavelength, only real part.
        'extinction_coefficient': imaginary part of the index of refraction of the material as a function of its wavelength.
        """
        super(VolumeMaterial, self).__init__(name, properties)
        self.kind = 'Volume'

    def change_of_direction(self, ray, normal_vector):
        wavelength = ray.wavelength
        if 'Matrix_reflectance_thin_film' in ray.current_medium.properties: # the ray is traveling inside thin film material and sure refraction
            ray.current_medium = self
            return OpticalState(ray.polarization_vectors[-1], ray.directions[-1], Phenomenon.REFRACTION)
        if 'Matrix_reflectance_thin_film' in self.properties: # the ray impacts on thin film material
            n1 = ray.materials[-1].properties['index_of_refraction'](ray.wavelength)
            if 'extinction_coefficient' in ray.current_medium.properties:
                n1 = n1 + 1j * ray.materials[-1].properties['extinction_coefficient'](ray.wavelength)
            n_front = self.properties['index_of_refraction_front'](ray.wavelength)
            n_back = self.properties['index_of_refraction_back'](ray.wavelength)
            if n_front == ray.materials[-1].properties['index_of_refraction'](ray.wavelength): # impacts on front material
                n2 = n_back
            else:
                n2 = n_front
            energy_absorbed_thin_film, optical_state = calculate_state_thin_film(ray.directions[-1], normal_vector, n1, n2,
                                                                         ray.polarization_vectors[-1],
                                                                         self.properties, ray.wavelength)
            ray.energy = ray.energy * (1.0 - energy_absorbed_thin_film)
            if optical_state.phenomenon == Phenomenon.REFRACTION:
                ray.current_medium = self		
            return optical_state
        else:
            n1 = ray.current_medium.properties['index_of_refraction'](wavelength)
            if 'extinction_coefficient' in ray.current_medium.properties:
                n1 = n1 + 1j * ray.current_medium.properties['extinction_coefficient'](wavelength)
            n2 = self.properties['index_of_refraction'](wavelength)
            if 'extinction_coefficient' in self.properties:
                n2 = n2 + 1j * self.properties['extinction_coefficient'](wavelength)
            optical_state = refraction(ray.directions[-1], normal_vector, n1, n2,
                                                                ray.polarization_vectors[-1])
            if optical_state.phenomenon == Phenomenon.REFRACTION:
                ray.current_medium = self
            return optical_state


def create_simple_volume_material(name, index_of_refraction, attenuation_coefficient=None):
    VolumeMaterial.create(name, {'index_of_refraction': constant_function(index_of_refraction),
                                 'attenuation_coefficient': constant_function(attenuation_coefficient)})


def create_wavelength_volume_material(name, file_index_of_refraction):
    # file_index_of_refraction with three columns: wavelenth in nm, real(index of refraction), imaginary(index of refraction)
    data_refraction = np.loadtxt(file_index_of_refraction, usecols=(0, 1, 2))
    wavelength_values = data_refraction[:, 0]
    n_values = data_refraction[:, 1]
    k_values = data_refraction[:, 2]
    VolumeMaterial.create(name, {'index_of_refraction': tabulated_function(wavelength_values, n_values),
                                 'extinction_coefficient': tabulated_function(wavelength_values, k_values)})

def create_PV_material(name, file_index_of_refraction):
    # file_index_of_refraction with three columns: wavelenth in nm, real(index of refraction), imaginary(index of refraction)
    data_refraction = np.loadtxt(file_index_of_refraction, usecols=(0, 1, 2))
    wavelength_values = data_refraction[:, 0]
    n_values = data_refraction[:, 1]
    k_values = data_refraction[:, 2]
    VolumeMaterial.create(name, {'index_of_refraction': tabulated_function(wavelength_values, n_values),
                                 'extinction_coefficient': tabulated_function(wavelength_values, k_values),
                                 'PV_material' : True})

def create_polarized_thin_film(name, file_thin_film, file_front, file_back):
    # thin film material calculated by TMM method, six columns:
	# wavelenth in nm, angle in deg., reflectance s-polarized (perpendicular), reflectance p-polarized (parallel),  transmittance s-polarized, transmittance p-polarized
    # the values in coating_material should be in the corresponding order columns
    data = np.loadtxt(file_thin_film)
    data_reflectance = data[:,[0,1,2,3]]
    data_transmittance = data[:,[0,1,4,5]]
    if file_front is not 'Vacuum':
        data_refraction_front = np.loadtxt(file_front, usecols=(0, 1, 2))
        wavelength_values_front = data_refraction_front[:, 0]
        n_values_front = data_refraction_front[:, 1]
        k_values_front = data_refraction_front[:, 2]
        index_of_refraction_front = tabulated_function(wavelength_values_front, n_values_front)
        extinction_coefficient_front = tabulated_function(wavelength_values_front, k_values_front)
    else:
        index_of_refraction_front = constant_function(1.0)
        extinction_coefficient_front = constant_function(0.0)
    if file_back is not 'Vacuum':
        data_refraction_back = np.loadtxt(file_back, usecols=(0, 1, 2))
        wavelength_values_back = data_refraction_back[:, 0]
        n_values_back = data_refraction_back[:, 1]
        k_values_back = data_refraction_back[:, 2]
        index_of_refraction_back = tabulated_function(wavelength_values_back, n_values_back)
        extinction_coefficient_back = tabulated_function(wavelength_values_back, k_values_back)
    else:
        index_of_refraction_back = constant_function(1.0)
        extinction_coefficient_back = constant_function(0.0)
    VolumeMaterial.create(name, {'Matrix_reflectance_thin_film': matrix_reflectance(data_reflectance),
                                 'Matrix_transmittance_thin_film': matrix_reflectance(data_transmittance),
                                 'index_of_refraction_front': index_of_refraction_front,
                                 'extinction_coefficient_front': extinction_coefficient_front,
                                 'index_of_refraction_back': index_of_refraction_back,
                                 'extinction_coefficient_back': extinction_coefficient_back,
                                 'thin_film' : True})

create_simple_volume_material("Vacuum", 1.0, 0.0)
# noinspection PyNoneFunctionAssignment
vacuum_medium = VolumeMaterial.by_name["Vacuum"]


class SurfaceMaterial(Material, object):
    def __init__(self, name, properties_front, properties_back=None):
        """
        Initializes a Surface Material. The properties parameter must be a dict with the physical properties
        describing the material. At least, the following must be provided:
        'probability_of_reflexion': probability that a photon gets reflected, as a function of its wavelength.
        'probability_of_absortion': probability that a photon gets absorbed, as a function of its wavelength.
        'probability_of_transmitance': probability that a photon passes through the material, as a function of its wavelength.
        """
        super(SurfaceMaterial, self).__init__(name, properties_front)
        self.properties_front = properties_front
        if properties_back is None:
            self.properties_back = properties_front
        else:
            self.properties_back = properties_back
        self.kind = 'Surface'

    @classmethod
    def create(cls, name, properties_front, properties_back=None):
        _ = cls(name, properties_front, properties_back)

    def decide_phenomenon(self, ray, normal_vector, properties, nearby_material):
        phenomena = [Phenomenon.REFLEXION, Phenomenon.ABSORTION, Phenomenon.TRANSMITANCE]
        polarization_vector = ray.polarization_vectors[-1]
        perpendicular_polarized = False
        probabilities = None
        if 'TW_model' in properties:
            b_constant = properties['b_constant']
            c_constant = properties['c_constant']
            absortion_ratio = TW_absorptance_ratio(normal_vector, b_constant, c_constant, ray.directions[-1])
            absortion = properties['probability_of_absortion'](ray.properties['wavelength']) * absortion_ratio
            por = 1.0 - absortion
            probabilities = [por, absortion, 0] # Here I assume no transmitance
        if 'Matrix_polarized_reflectance_coating' in properties:  # polarized_coating_layer
            n1 = ray.current_medium.properties['index_of_refraction'](ray.wavelength)
            if 'extinction_coefficient' in ray.current_medium.properties:
                n1 = ray.current_medium.properties['index_of_refraction'](ray.wavelength) + 1j * ray.current_medium.properties['extinction_coefficient'](ray.wavelength)
            n2 = nearby_material.properties['index_of_refraction'](ray.wavelength)
            if 'extinction_coefficient' in nearby_material.properties:
                n2 = nearby_material.properties['index_of_refraction'](ray.wavelength) + 1j * nearby_material.properties['extinction_coefficient'](ray.wavelength)		
            results = calculate_probabilities_polarizaton_coating(ray.directions[-1], normal_vector, n1, n2,
                                                                         ray.polarization_vectors[-1],
                                                                         properties, ray.wavelength)
            probabilities = [results[0], results[1], results[2]]
            polarization_vector = results[3]
            perpendicular_polarized = results[4] 
        if 'metallic_material' in properties:
            n1 = ray.current_medium.properties['index_of_refraction'](ray.wavelength)
            if 'extinction_coefficient' in ray.current_medium.properties:
                n1 = ray.current_medium.properties['index_of_refraction'](ray.wavelength) + 1j * ray.current_medium.properties['extinction_coefficient'](ray.wavelength)
            n2 = properties['index_of_refraction'](ray.wavelength)
            if 'extinction_coefficient' in properties:
                n2 = properties['index_of_refraction'](ray.wavelength) + 1j * properties['extinction_coefficient'](ray.wavelength)
            results = calculate_reflexion_metallic(ray.directions[-1], normal_vector, n1, n2, polarization_vector)
            probabilities = [results[0], results[1], results[2]]
            polarization_vector = results[3]
            perpendicular_polarized = results[4]

        if probabilities is None:
            try:
                por = properties['probability_of_reflexion'](ray.properties['wavelength'])
            except:
                por = 1.0
            try:
                poa = properties['probability_of_absortion'](ray.properties['wavelength'])
            except:
                poa = 1 - por
            try:
                pot = properties['probability_of_transmitance'](ray.properties['wavelength'])
            except:
                pot = 0.0

            probabilities = [por, poa, pot]
        phenomenon = np.random.choice(phenomena, 1, p=probabilities)[0]
        return phenomenon, polarization_vector, perpendicular_polarized

    def change_of_direction_by_absortion(self, ray, normal_vector, properties):
        if properties['energy_collector']:
            return OpticalState(Base.Vector(0.0, 0.0, 0.0), Base.Vector(0.0, 0.0, 0.0), Phenomenon.GOT_ABSORBED)
        else:
            return OpticalState(Base.Vector(0.0, 0.0, 0.0), Base.Vector(0.0, 0.0, 0.0), Phenomenon.ABSORTION)

    def change_of_direction_by_reflexion(self, ray, normal_vector, properties, polarization_vector_calculated_before):
        if 'specular_material' in properties:
            state = reflexion(ray.directions[-1], normal_vector, ray.polarization_vectors[-1], polarization_vector_calculated_before)
            try:
                sigma_1 = properties['sigma_1']               
                try:
                    sigma_2 = properties['sigma_2']
                    try:
                        k = properties['k']
                        return double_gaussian_dispersion(normal_vector, state, sigma_1, sigma_2, k)
                    except:
                        k = 0.5
                        return double_gaussian_dispersion(normal_vector, state, sigma_1, sigma_2, k)
                except:
                    pass
                return single_gaussian_dispersion(normal_vector, state, sigma_1)                 
            except:
                pass
            return state
        if 'lambertian_material' in properties:
            state = lambertian_reflexion(ray.directions[-1], normal_vector)
            return state
			

    def change_of_direction_by_transmitance(self, ray, normal_vector, nearby_material, perpendicular_polarized):
        n1 = ray.current_medium.properties['index_of_refraction'](ray.wavelength)
        if 'extinction_coefficient' in ray.current_medium.properties:
            n1 = ray.current_medium.properties['index_of_refraction'](ray.wavelength) + 1j * ray.current_medium.properties['extinction_coefficient'](ray.wavelength)
        n2 = nearby_material.properties['index_of_refraction'](ray.wavelength)
        if 'extinction_coefficient' in nearby_material.properties:
            n2 = nearby_material.properties['index_of_refraction'](ray.wavelength) + 1j * nearby_material.properties['extinction_coefficient'](ray.wavelength)
        if n1 == n2: # transparent_simple_layer
            state = OpticalState(ray.polarization_vectors[-1], ray.directions[-1], Phenomenon.REFRACTION)
        else:
            state = shure_refraction(ray.directions[-1], normal_vector, n1, n2,
                                                                ray.polarization_vectors[-1],
                                                                perpendicular_polarized)
        ray.current_medium = nearby_material
        return state


    def change_of_direction(self, ray, normal_vector, nearby_material):

        if ray.directions[-1].dot(normal_vector) < 0:  # Ray intercepted on the frontside of the surface
            properties = self.properties_front
        else:  # Ray intercepted on the backside of the surface
            properties = self.properties_back

        results = self.decide_phenomenon(ray, normal_vector, properties, nearby_material)
        phenomenon = results[0]
        if ray.polarization_vectors[-1] == results[1]: # polarization_vector not calculated
            polarization_vector_calculated_before = False
        else: # polarization_vector calculated before
            polarization_vector_calculated_before = True
        ray.polarization_vectors[-1] = results[1]
        perpendicular_polarized = results[2] # True or False

        if phenomenon == Phenomenon.REFLEXION:
            return self.change_of_direction_by_reflexion(ray, normal_vector, properties, polarization_vector_calculated_before)
        elif phenomenon == Phenomenon.ABSORTION:
            return self.change_of_direction_by_absortion(ray, normal_vector, properties)
        elif phenomenon == Phenomenon.TRANSMITANCE:
            return self.change_of_direction_by_transmitance(ray, normal_vector, nearby_material, perpendicular_polarized)


def create_opaque_simple_layer(name):
    SurfaceMaterial.create(name, {'probability_of_reflexion': constant_function(0.0),
                                  'probability_of_absortion': constant_function(1.0),
                                  'probability_of_transmitance': constant_function(0.0),
                                  'energy_collector': False})


def create_transparent_simple_layer(name, pot):
    SurfaceMaterial.create(name, {'probability_of_reflexion': constant_function(1 - pot),
                                  'probability_of_absortion': constant_function(0),
                                  'probability_of_transmitance': constant_function(pot),
                                  'energy_collector': False,
                                  'specular_material': True})

								  
def create_absorber_simple_layer(name, poa):
    SurfaceMaterial.create(name, {'probability_of_reflexion': constant_function(1 - poa),
                                  'probability_of_absortion': constant_function(poa),
                                  'probability_of_transmitance': constant_function(0),
                                  'energy_collector': True,
                                  'specular_material': True})
								  
								  
def create_absorber_lambertian_layer(name, poa):
    SurfaceMaterial.create(name, {'probability_of_reflexion': constant_function(1 - poa),
                                  'probability_of_absortion': constant_function(poa),
                                  'probability_of_transmitance': constant_function(0),
                                  'energy_collector': True,
                                  'lambertian_material': True})
								  

def create_absorber_TW_model_layer(name, poa, b_constant, c_constant):
    SurfaceMaterial.create(name, {'probability_of_reflexion': constant_function(1 - poa),
                                  'probability_of_absortion': constant_function(poa),
                                  'probability_of_transmitance': constant_function(0),
                                  'energy_collector': True,
                                  'lambertian_material': True,
                                  'TW_model': True,
                                  'b_constant': b_constant,
                                  'c_constant': c_constant})


def create_reflector_specular_layer(name, por, sigma_1 = None, sigma_2 = None, k = None):
    SurfaceMaterial.create(name, {'probability_of_reflexion': constant_function(por),
                                  'probability_of_absortion': constant_function(1 - por),
                                  'probability_of_transmitance': constant_function(0),
                                  'energy_collector': False,
                                  'specular_material': True,
                                  'sigma_1': sigma_1,
                                  'sigma_2': sigma_2,
                                  'k': k})


def create_reflector_lambertian_layer(name, por):
    SurfaceMaterial.create(name, {'probability_of_reflexion': constant_function(por),
                                  'probability_of_absortion': constant_function(1 - por),
                                  'probability_of_transmitance': constant_function(0),
                                  'energy_collector': False,
                                  'lambertian_material': True})
							  

def create_metallic_specular_layer(name, file_index_of_refraction, sigma_1 = None, sigma_2 = None, k = None):
    # file_index_of_refraction with three columns: wavelenth in nm, real(index of refraction), imaginary(index of refraction)
    data_refraction = np.loadtxt(file_index_of_refraction, usecols=(0, 1, 2))
    wavelength_values = data_refraction[:, 0]
    n_values = data_refraction[:, 1]
    k_values = data_refraction[:, 2]
    SurfaceMaterial.create(name, {'index_of_refraction': tabulated_function(wavelength_values, n_values),
                                 'extinction_coefficient': tabulated_function(wavelength_values, k_values),
                                  'energy_collector': False,
                                  'specular_material': True,
                                  'metallic_material': True,
                                  'sigma_1': sigma_1,
                                  'sigma_2': sigma_2,
                                  'k': k})
								  
def create_metallic_lambertian_layer(name, file_index_of_refraction):
    # file_index_of_refraction with three columns: wavelenth in nm, real(index of refraction), imaginary(index of refraction)
    data_refraction = np.loadtxt(file_index_of_refraction, usecols=(0, 1, 2))
    wavelength_values = data_refraction[:, 0]
    n_values = data_refraction[:, 1]
    k_values = data_refraction[:, 2]
    SurfaceMaterial.create(name, {'index_of_refraction': tabulated_function(wavelength_values, n_values),
                                 'extinction_coefficient': tabulated_function(wavelength_values, k_values),
                                  'energy_collector': False,
                                  'lambertian_material': True,
                                  'metallic_material': True})
								  
								  
def create_polarized_coating_reflector_layer(name, coating_file, sigma_1 = None, sigma_2 = None, k = None):
    # coating_material with four columns: wavelenth in nm, angle in deg., reflectance s-polarized (perpendicular), reflectance p-polarized (parallel)
    # the values in coating_material should be in the corresponding order columns
    data_material = np.loadtxt(coating_file, usecols=(0,1,2,3))
    SurfaceMaterial.create(name, {'Matrix_polarized_reflectance_coating': matrix_reflectance(data_material),
                                  'probability_of_absortion': constant_function(0),
                                  'energy_collector': False,
                                  'specular_material': True,
                                  'transparent_material': False,
                                  'sigma_1': sigma_1,
                                  'sigma_2': sigma_2,
                                  'k': k})

								  
def create_polarized_coating_transparent_layer(name, coating_file):
    # coating_material with four columns: wavelenth in nm, angle in deg., reflectance s-polarized (perpendicular), reflectance p-polarized (parallel)
    # the values in coating_material should be in the corresponding order columns
    data_material = np.loadtxt(coating_file, usecols=(0,1,2,3))
    SurfaceMaterial.create(name, {'Matrix_polarized_reflectance_coating': matrix_reflectance(data_material),
                                  'probability_of_absortion': constant_function(0),
                                  'energy_collector': False,
                                  'specular_material': False,
                                  'transparent_material': True})

								   
def create_polarized_coating_absorber_layer(name, coating_file):
    # coating_material with four columns: wavelenth in nm, angle in deg., reflectance s-polarized (perpendicular), reflectance p-polarized (parallel)
    # the values in coating_material should be in the corresponding order columns
    data_material = np.loadtxt(coating_file, usecols=(0,1,2,3))
    SurfaceMaterial.create(name, {'Matrix_polarized_reflectance_coating': matrix_reflectance(data_material),
                                  'probability_of_transmitance': constant_function(0),
                                  'energy_collector': True,
                                  'specular_material': False,
                                  'transparent_material': False,
                                  'lambertian_material': True})								  

								  
def create_two_layers_material(name, layer_front, layer_back):
    SurfaceMaterial.create(name, Material.by_name[layer_front].properties,
	                       Material.by_name[layer_back].properties)


# endregion

class Scene:
    """
    Class used to define the Scene. It encodes all the objects 
    that interfere with light rays. 
    """

    def __init__(self, objects):
        self.faces = []  # All the faces in the Scene
        self.solids = []  # All the solids in the Scene
        self.materials = {}  # Assign the materials to objects
        # self.sum_energy = 0.0 # Energy of the Scene
        self.epsilon = 40 * 10.0 ** (-6.0)  # Tolerance for solid containment # 10 nm.
        self.boundbox = None

        for obj in objects:
            # noinspection PyNoneFunctionAssignment
            material = Material.get_from_label(obj.Label)
            if not material:
                continue
            solids = obj.Shape.Solids
            faces = obj.Shape.Faces
            if solids:  # Object is a solid
                for solid in solids:
                    self.materials[solid] = material
            else:  # Object is made of faces
                for face in faces:
                    self.materials[face] = material
            self.solids.extend(solids)
            self.faces.extend(faces)
            # self.materials[obj] = material ## Cal?
            if not self.boundbox:
                self.boundbox = obj.Shape.BoundBox
            else:
                self.boundbox.add(obj.Shape.BoundBox)

        self.diameter = self.boundbox.DiagonalLength
        self.remove_duplicate_faces()

    def remove_duplicate_faces(self):
        faces_with_material = [face for face in self.faces if
                               face in self.materials]
        faces_no_material = [face for face in self.faces if
                             face not in self.materials]
        if not faces_no_material:
            return
        complex_no_material = faces_no_material[0]
        for face in faces_no_material[1:]:
            complex_no_material = complex_no_material.fuse(face)
        for face in faces_with_material:
            complex_no_material = complex_no_material.cut(face)
        self.faces = faces_with_material
        self.faces.extend(complex_no_material.Faces)

    def solid_at_point(self, point):
        """        
        Returns the solid that a point is inside.
        """
        for solid in self.solids:
            if solid.isInside(point, self.epsilon, False):
                return solid
        return None

    def face_at_point(self, point):
        """        
        Returns the face that a point is inside.
        """
        for face in self.faces:
            if face.isInside(point, self.epsilon, True):
                return face
        return None


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
                # noinspection PyUnusedLocal
                min_area = area
                best_origin = p + v1 * minx + v2 * miny
                best_v1 = v1
                best_v2 = v2
                best_length1 = length1
                best_length2 = length2
            # noinspection PyUnboundLocalVariable
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


class Experiment:
    """
    Sets up and runs and experiment in a given scene with a given light source.
    If show_in_doc is given, the emitting region is drawn in the FreeCAD active document.
    If show_in_doc is given, the rays could be drawn in the FreeCAD active document (using the run function).
    """

    def __init__(self, scene, light_source, number_of_rays, show_in_doc=None):
        self.scene = scene
        self.light_source = light_source
        if show_in_doc:
            self.light_source.emitting_region.add_to_document(show_in_doc)
        self.number_of_rays = number_of_rays
        self.wavelengths = []
        self.captured_energy_Th = 0
        self.captured_energy_PV = 0
        self.Th_energy = []
        self.Th_wavelength = []
        self.PV_energy = []
        self.PV_wavelength = []
        self.PV_values = []
        self.points_absorber_Th = []
        random_congruential(time.time()) # TODO: change location

    def run(self, show_in_doc=None):
        for i in np.arange(0,self.number_of_rays,1):
            ray = self.light_source.emit_ray()
            ray.run()
            self.wavelengths.append(ray.wavelength)
            if show_in_doc:
                ray.add_to_document(show_in_doc)
            if ray.got_absorbed:
                self.captured_energy_Th += ray.energy
                self.Th_energy.append(ray.energy)
                self.Th_wavelength.append(ray.wavelength)
                self.points_absorber_Th.append((ray.energy, 
                                             ray.points[-1].x, ray.points[-1].y, ray.points[-1].z,
                                             ray.points[-2].x, ray.points[-2].y, ray.points[-2].z,
                                             ray.normals[-1].x, ray.normals[-1].y, ray.normals[-1].z))
            else:
                self.Th_energy.append(0.0)
                self.Th_wavelength.append(ray.wavelength)
            if ray.in_PV:
                PV_energy_absorbed = np.sum(ray.PV_absorbed)
                self.captured_energy_PV += PV_energy_absorbed 
                self.PV_energy.append(PV_energy_absorbed)
                self.PV_wavelength.append(ray.wavelength)
                length = len(ray.PV_values)
                if length > 0:
                    for z in np.arange(0,length,1):
                        self.PV_values.append(ray.PV_values[z])
            else:
                self.PV_energy.append(0.0)
                self.PV_wavelength.append(ray.wavelength)
