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


# ---
# Helper functions for reflexion and refraction
# ---
def reflexion(incident, normal):
    """
    Implementation of law of reflexion
    """
    # We assume all vectors are normalized and the normal.dot(incident) < 0
    return incident - normal * 2.0 * normal.dot(incident), "Reflexion"


def lambertian_reflexion(normal):
    """
    Implementation of lambertian reflection for diffusely reflecting surface
    """
    # We assume the normal is normalized and the normal.dot(incident) < 0
    dot =  - 1.0
    while (dot < 0.0): 
        random_vector = Base.Vector(myrandom()-0.5, myrandom()-0.5, myrandom()-0.5)
        random_vector.normalize()
        dot = normal.dot(random_vector)
    return random_vector, "Reflexion"
	
	
def single_gaussian_dispersion(normal, reflected, sigma_1):
    """
    Implementation of single gaussian dispersion based on the ideal reflected direction
    """
    rad = np.pi/180.0
    u = myrandom()
    theta = (-2. * sigma_1 ** 2. * np.log(u)) ** 0.5 / 1000.0 / rad
    v = reflected[0]
    axis_1 = normal.cross(reflected[0])
    rotation_1 = Base.Rotation(axis_1,theta)
    new_v1 = rotation_1.multVec(v)
    u = myrandom()
    phi = 360. * u
    axis_2 = v
    rotation_2 = Base.Rotation(axis_2,phi)
    new_v2 = rotation_2.multVec(new_v1)    
    return new_v2, "Reflexion"    


def double_gaussian_dispersion(normal, reflected, sigma_1, sigma_2, k):
    """
    Implementation of double gaussian dispersion based on the ideal reflected direction
    """
    rad = np.pi/180.0
    k_ran = myrandom()
    u = myrandom()
    if k_ran < k :
        theta = (-2. * sigma_1 ** 2. * np.log(u)) ** 0.5 / 1000.0 / rad
    else:
	    theta = (-2. * sigma_2 ** 2. * np.log(u)) ** 0.5 / 1000.0 / rad
    v = reflected[0]
    axis_1 = normal.cross(reflected[0])
    rotation_1 = Base.Rotation(axis_1,theta)
    new_v1 = rotation_1.multVec(v)
    u = myrandom()
    phi = 360. * u
    axis_2 = v
    rotation_2 = Base.Rotation(axis_2,phi)
    new_v2 = rotation_2.multVec(new_v1)    
    return new_v2, "Reflexion"


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
    incidence_angle = np.arccos (normal.dot(incident) * (-1.0))
    incidence_angle_deg = incidence_angle * 180.0 / np.pi
    if incidence_angle_deg < 80.0:
        absortion_ratio = 1.0 - b_constant * (1.0 / np.cos(incidence_angle) - 1.0) ** c_constant
    else:
        y0 = 1.0 - b_constant * (1.0 / np.cos(80.0 * np.pi / 180.0) - 1.0) ** c_constant
        m = y0 / 10.0
        absortion_ratio = y0 - m * (incidence_angle_deg - 80.0)
    return absortion_ratio
	

	# noinspection PyUnresolvedReferences,PyUnresolvedReferences
def refraction(incident, normal, n1, n2, perpendicular_polarized):
    """
    Implementation of Snell's law of refraction
    """
    # https://en.wikipedia.org/wiki/Snell's_law#Vector_form
    mynormal = normal * 1.0
    if mynormal.dot(incident) > 0:
        # noinspection PyAugmentAssignment
        mynormal = mynormal * (-1.0)
    r = n1 / n2
    c1 = - mynormal.dot(incident) # cos (incidence_angle)
    c2sq = 1.0 - r * r * (1.0 - c1 * c1) # cos (refracted_angle) ** 2
    if c2sq < 0: # total internal reflection
        return reflexion(incident, normal)
    c2 = c2sq ** 0.5 # cos (refracted_angle)
    # https://en.wikipedia.org/wiki/Fresnel_equations # Fresnel equations
    if perpendicular_polarized:
        R = ((n1 * c1 - n2 * c2) / (n1 * c1 + n2 * c2)) ** 2. # reflectance for s-polarized (perpendicular) light
    else:
        R = ((n1 * c2 - n2 * c1) / (n1 * c2 + n2 * c1)) ** 2. # reflectance for p-polarized (parallel) light
    if myrandom() < R:
        return reflexion(incident, normal)
    else:		
        return incident * r + mynormal * (r * c1 - c2), "Refraction"


def polar_to_cartesian(phi, theta):
    """
    Convert polar coordinates (given in degrees) to cartesian
    :param phi:
    :param theta:
    :return:
    """
    rad = np.arccos(-1.0) / 180.0
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
    seed = np.remainder(random_congruential.previous * a + c , m)
    random_congruential.previous = seed
    return seed * rm


# ---
# Define the random algorithm
# ---
myrandom = random_congruential
#myrandom = random.random
	
	
# ---
# Helper function for dispersions and polarization vector  
# ---	
def dispersion_from_main_direction(main_direction, theta, phi):
    v = main_direction
    v_p = Base.Vector(v[1],-v[0],0)
    if v_p == Base.Vector(0,0,0): # to avoid null vector at mynormal and incident parallel vectors
        v_p = Base.Vector(1,0,0)
    v_p.normalize()
    rotation_1 = Base.Rotation(v_p,theta)
    new_v1 = rotation_1.multVec(v)
    axis_2 = main_direction
    rotation_2 = Base.Rotation(axis_2,phi)
    new_v2 = rotation_2.multVec(new_v1)
    return new_v2    


def dispersion_polarization(main_direction,polarization_vector,theta,phi):
    v = main_direction
    v_p = Base.Vector(v[1],-v[0],0)
    if v_p == Base.Vector(0,0,0): # to avoid null vector at mynormal and incident parallel vectors
        v_p = Base.Vector(1,0,0)
    v_p.normalize()
    rotation_1 = Base.Rotation(v_p,theta)
    new_v1 = rotation_1.multVec(polarization_vector)
    axis_2 = main_direction
    rotation_2 = Base.Rotation(axis_2,phi)
    new_v2 = rotation_2.multVec(new_v1)	    
    return new_v2  	


def random_polarization(direction):
    v_p = Base.Vector(direction[1],-direction[0],0)
    if v_p == Base.Vector(0,0,0):
        v_p = Base.Vector(1,0,0)
    v_p.normalize()
    phi = 360.0 * myrandom()
    rotation = Base.Rotation(direction,phi)
    new_v_p = rotation.multVec(v_p)	
    return new_v_p   
	

# ---
# Classes for materials
# ---

# region Materials

class Material(object):
    """
    Class used to represent materials and their physical properties
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

    def change_of_direction(self, ray, normal_vector):
        pass
        # Compute new direction (ray.current_material)


class VolumeMaterial(Material, object):
    def __init__(self, name, properties):
        """
        Initializes a Volume Material. The properties parameter must be a dict with the physical properties
        describing the material. At least, the following must be provided:
        'index_of_refraction': index of refraction of the material, as a function of its wavelength.
        'extinction_coefficient': extinction coefficient of the material in m-1
        """
        super(VolumeMaterial, self).__init__(name, properties)
        self.kind = 'Volume'

    def change_of_direction(self, ray, normal_vector):
        wavelength = ray.wavelength
        n1 = ray.current_medium.properties['index_of_refraction'](wavelength)
        n2 = self.properties['index_of_refraction'](wavelength)
        direction, phenomenon = refraction(ray.directions[-1], normal_vector, n1, n2, ray.perpendicular_polarized)
        if phenomenon == "Refraction":
            ray.current_medium = self
        return direction, phenomenon
        pass


def create_simple_volume_material(name, index_of_refraction, extinction_coefficient = None):
    VolumeMaterial.create(name, {'index_of_refraction': constant_function(index_of_refraction),
                                 'extinction_coefficient': constant_function(extinction_coefficient)})


create_simple_volume_material("Vacuum", 1.0,0.0)
# noinspection PyNoneFunctionAssignment
vacuum_medium = VolumeMaterial.by_name["Vacuum"]


class SurfaceMaterial(Material, object):
    def __init__(self, name, properties_front, properties_back = None):
        """
        Initializes a Surface Material. The properties parameter must be a dict with the physical properties
        describing the material. At least, the following must be provided:
        'probability_of_reflexion': probability that a photon gets reflected, as a function of its wavelength.
        'probability_of_absortion': probability that a photon gets absorbed, as a function of its wavelength.
        'transmitance': probability that a photon passes through the material, as a function of its wavelength.
        """
        super(SurfaceMaterial, self).__init__(name, properties_front)
        self.properties_front = properties_front
        if properties_back is None:
            self.properties_back = properties_front
        else:
            self.properties_back = properties_back
        self.kind = 'Surface'

    @classmethod
    def create(cls, name, properties_front, properties_back = None):
        _ = cls(name, properties_front, properties_back)

    def change_of_direction(self, ray, normal_vector):
        # TODO: Implementar bé. El que hi ha és provisional
        phenomena = ["Reflexion", "Absortion"]
        if ray.directions[-1].dot(normal_vector) < 0 : # Ray intercepted on the frontside of the surface
            normal = normal_vector
            properties = self.properties_front
        else: # Ray intercepted on the backside of the surface
            normal = normal_vector * (-1.0)
            properties = self.properties_back
        if 'TW_model' in properties:
            b_constant = properties['b_constant']
            c_constant = properties['c_constant'] 
            absortion_ratio = TW_absorptance_ratio(normal, b_constant, c_constant, ray.directions[-1])
            abs = properties['probability_of_absortion'](ray.properties['wavelength']) * absortion_ratio
            por = 1.0 - abs
            probabilities = [por, abs]
        else:
		    probabilities = [properties['probability_of_reflexion'](ray.properties['wavelength']),
                             properties['probability_of_absortion'](ray.properties['wavelength'])]
        phenomenon = np.random.choice(phenomena, 1, p = probabilities)[0]
        if phenomenon == 'Reflexion':
            if 'specular_material' in properties:
                reflected = reflexion(ray.directions[-1], normal)
                if 'sigma_1' in properties:
                    sigma_1 = properties['sigma_1']
                    if 'sigma_2' in properties:
                        sigma_2 = properties['sigma_2']
                        k = properties['k']
                        return double_gaussian_dispersion(normal, reflected, sigma_1, sigma_2, k)
                    return single_gaussian_dispersion(normal, reflected, sigma_1)						   
                return reflected
            if 'lambertian_material' in properties:
                reflected = lambertian_reflexion(normal)
                return reflected
            if 'dispersion_factor' in properties:
                reflected = reflected + myrandom()
                return reflected
        if phenomenon == "Absortion":
            if properties['energy_collector']:
                return Base.Vector(0.0, 0.0, 0.0), "Got_Absorbed"
            else:
                return Base.Vector(0.0, 0.0, 0.0), "Absortion"
        pass

    def scatter_direction(self, ray, direction):
        # TODO: pensar
        pass


def create_opaque_simple_material(name, por):
    SurfaceMaterial.create(name, {'probability_of_reflexion': constant_function(por),
                                  'probability_of_absortion': constant_function(1 - por),
                                  'transmitance': constant_function(0),
                                  'energy_collector': False,
                                  'specular_material': True}, None)


def create_transparent_simple_material(name, por):
    SurfaceMaterial.create(name, {'probability_of_reflexion': constant_function(por),
                                  'probability_of_absortion': constant_function(0),
                                  'transmitance': constant_function(1 - por)})

								  
def create_absorber_simple_material(name, por):
    SurfaceMaterial.create(name, {'probability_of_reflexion': constant_function(por),
                                  'probability_of_absortion': constant_function(1 - por),
                                  'transmitance': constant_function(0),
                                  'energy_collector': True,
                                  'lambertian_material': True})

								  
def simple_reflector_twolayers(name, por_front, por_back):
    SurfaceMaterial.create(name, {'probability_of_reflexion': constant_function(por_front),
                                  'probability_of_absortion': constant_function(1 - por_front),
                                  'transmitance': constant_function(0),
                                  'specular_material': True,
                                  'energy_collector': False},
                                 {'probability_of_reflexion': constant_function(por_back),
                                  'probability_of_absortion': constant_function(1 - por_back),
                                  'transmitance': constant_function(0),
                                  'specular_material': True,
                                  'energy_collector': False})
								  
								  
def create_absorber_lambertian_layer(name, por):
    SurfaceMaterial.create(name, {'probability_of_reflexion': constant_function(por),
                                  'probability_of_absortion': constant_function(1 - por),
                                  'transmitance': constant_function(0),
                                  'energy_collector': True,
                                  'lambertian_material': True})


def create_absorber_TW_model_layer(name, por, b_constant, c_constant):
    SurfaceMaterial.create(name, {'probability_of_reflexion': constant_function(por),
                                  'probability_of_absortion': constant_function(1 - por),
                                  'transmitance': constant_function(0),
                                  'energy_collector': True,
                                  'lambertian_material': True,
                                  'TW_model': True,
                                  'b_constant': b_constant,
                                  'c_constant': c_constant})
								  
								  
def create_reflector_specular_layer(name, por):
    SurfaceMaterial.create(name, {'probability_of_reflexion': constant_function(por),
                                  'probability_of_absortion': constant_function(1 - por),
                                  'transmitance': constant_function(0),
                                  'energy_collector': False,
                                  'specular_material': True})
								  

def create_reflector_onegaussian_layer(name, por, sigma):
    SurfaceMaterial.create(name, {'probability_of_reflexion': constant_function(por),
                                  'probability_of_absortion': constant_function(1 - por),
                                  'transmitance': constant_function(0),
                                  'energy_collector': False,
                                  'specular_material': True,
                                  'sigma_1': sigma})

								  
def create_reflector_twogaussian_layer(name, por, sigma_1, sigma_2, k):
    SurfaceMaterial.create(name, {'probability_of_reflexion': constant_function(por),
                                  'probability_of_absortion': constant_function(1 - por),
                                  'transmitance': constant_function(0),
                                  'energy_collector': False,
                                  'specular_material': True,
                                  'sigma_1': sigma_1,
                                  'sigma_2': sigma_2,
                                  'k' : k})

								  
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
        self.epsilon = 0.00001  # Tolerance for solid containment
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
            if solid.isInside(point, self.epsilon, True):
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
    that the projection of the scene on the plane is enclosed in the rectangle
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
            return best_origin, best_v1, best_v2, best_length1, best_length2

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


class SunWindowBuie(SunWindow,object):
    """
    Class inherit from SunWindow, used to define a sun window with the sun profile based on Buie Equations.
    The lenght of the window is 10% greater than the SunWindow class.
	
    A source of direct light from the sun according to the Buie model:
    
    Buie, D., 2005. Corrigendum to "The effective size of the solar cone for 
    solar concentrating systems" [Solar Energy 74 (2003) 417-427]. 
    Solar Energy 79, 568-570. 5.2
    """

		
    def __init__(self, scene, direction, CircumSolarRatio):
        super(SunWindowBuie,self).__init__(scene, direction)
        self.length1 *= 1.1
        self.length2 *= 1.1
        self.aperture = self.length1 * self.length2
        self.origin = self.origin - self.v1 * self.length1 * 0.05 - self.v2 * self.length2 * 0.05
        self.CSR = CircumSolarRatio
        self.SD = 4.65
        self.SS = 43.6
        self.aa1 = SunWindowBuie.a1(self.CSR,4.65)
        self.aa2 = SunWindowBuie.a2(self.CSR,self.SD,self.SS)
        self.CDF_Disk_Region = SunWindowBuie.CDF_disk_region(self.aa1,self.SD)
        self.main_direction = direction

		
    def random_point(self):
        return (self.origin + self.v1 * self.length1 * myrandom() +
                self.v2 * self.length2 * myrandom())

				
    def random_direction(self):
        ran_1 = myrandom()
        if ran_1 < 1.0 - self.CSR:
            th_u = SunWindowBuie.th_solar_disk_region(ran_1,self.aa1,self.SD)/1000.0*180.0/np.pi
        else:
            th_u = SunWindowBuie.th_circumsolar_region(ran_1,self.CSR,self.SD,self.SS,self.aa2)/1000.0*180.0/np.pi 
        v = self.main_direction
        rotation_1 = Base.Rotation(self.v1,th_u)
        new_v1 = rotation_1.multVec(v)
        ran_2 = myrandom()
        phi = 360. * ran_2
        axis_2 = v
        rotation_2 = Base.Rotation(axis_2,phi)
        new_v2 = rotation_2.multVec(new_v1)	    
        return new_v2

				
    @staticmethod
    def solar_disk__density_distribution(angle):
        return 2.0 *np.pi * np.cos(0.326*angle) / np.cos(0.305*angle) * angle
		
		
    @staticmethod
    def circumsolar__density_distribution(angle,CSR):
        gamma = 2.2 * np.log(0.52 * CSR) * CSR **(0.43) - 0.1
        kappa = 0.9 * np.log(13.5 * CSR) * CSR ** (-0.3)
        return 2.0 * np.pi * np.exp(kappa) * angle ** (gamma + 1.0)

    @staticmethod		
    def a1(CSR,SD):
        """ Parameter a1 needed for the normalization of the probability distribution in thedisk region"""    
        f = SunWindowBuie.solar_disk__density_distribution
        th = np.arange(0.0,SD,0.001)
        f_th = np.vectorize(f)
        aa1 = ( 1.0 - CSR ) / np.trapz(f_th(th), dx = 0.001)		
        return aa1
		
    @staticmethod		
    def a2(CSR,SD,SS):
        """ Parameter a2 needed for the normalization of the probability distribution in the circumsolar region"""    
        f = SunWindowBuie.circumsolar__density_distribution
        th = np.arange(SD,SS,0.001)
        f_th = np.vectorize(f)
        aa2 = CSR / np.trapz(f_th(th,CSR), dx = 0.001)		
        return aa2
		
    @staticmethod	
    def CDF_disk_region(a1,SD):
        """ Cumulative Distribution Function in the solar disk region"""
        f = SunWindowBuie.solar_disk__density_distribution
        th = np.arange(0.0,SD,0.001)
        f_th = np.vectorize(f)
        cumulative = np.cumsum(f_th(th))
        CDF = (th, a1 * cumulative / 1000.0)			
        return CDF
		
    @staticmethod
    def th_solar_disk_region(u,a1,SD):
        """ Random angle based on the probability distribution in the circumsolar region"""
        CDF = SunWindowBuie.CDF_disk_region(a1,SD)
        idx = (np.abs(CDF[1] - u)).argmin()
        return  CDF[0][idx]	
    
    @staticmethod
    def th_circumsolar_region(u,CSR,SD,SS,aa2):
        """ Random angle based on the CDF in the circumsolar region"""
        gamma = 2.2 * np.log(0.52 * CSR) * CSR **(0.43) - 0.1
        kappa = 0.9 * np.log(13.5 * CSR) * CSR ** (-0.3)
        u_csr = (u - 1.0) + CSR # Since CSR-CDF starts at zero
        f1 = u_csr  * (gamma + 2.0) / (aa2 * 2 * np.pi * np.exp(kappa))
        f2 = SD ** (gamma + 2.0)
        exp = (1.0 / (gamma + 2.0))
        th_u = np.power(f1 + f2, exp)
        return  th_u
        		
		
class Ray:
    """
    Class used to model a sun ray and its polarization vector. It keeps information of the path it 
    describes as it interacts with the scene and its energy.
    """

    def __init__(self, scene, origin, direction, properties):
        self.scene = scene
        self.points = [origin]
        self.directions = [direction]
        self.materials = [vacuum_medium]
        self.properties = properties
        self.wavelength = properties['wavelength']
        self.energy = properties['energy']
        self.polarization_vectors = [properties['polarization_vector']]
        self.finished = False
        self.got_absorbed = False
        self.current_medium = vacuum_medium

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

    def next_direction(self, face):
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
            polarization_vector, direction, phenomenon = material.change_of_direction(self, normal)# TODO polarization_vector 
        else:
            # face is not active
            point_plus_delta = current_point + current_direction * self.scene.epsilon
            next_solid = self.scene.solid_at_point(point_plus_delta)
            nearby_material = self.scene.materials.get(next_solid, vacuum_medium)
            polarization_vector, direction, phenomenon = nearby_material.change_of_direction(self, normal)# TODO polarization_vector
        next_material = None
        if phenomenon == 'Refraction':
            # noinspection PyUnboundLocalVariable
            next_material = nearby_material
        elif phenomenon == 'Reflexion':
            next_material = current_material
        elif phenomenon == 'Absortion':
            next_material = None
        elif phenomenon == 'Got_Absorbed': # it is needed? Review
            next_material = None
        return polarization_vector, direction, next_material, phenomenon

    def update_energy(self):
        # TODO: @Ramon
        point_1 = self.points[-1]
        point_2 = self.points[-2]
        middle_point = point_1.add(point_2) * 0.5
        actual_solid = self.scene.solid_at_point(middle_point)
        if actual_solid:
            if 'extinction_coefficient' in self.scene.materials[actual_solid].properties:
                k = self.scene.materials[actual_solid].properties['extinction_coefficient'](self.wavelength)
                if k > 0:
                    d = point_1.distanceToPoint(point_2)
                    self.energy = self.energy * np.exp(- k * d / 1000.0)

    def run(self, max_hops=20):
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
            polarization_vector, vector, material, phenomenon = self.next_direction(face) # TODO polarization_vector
            self.directions.append(vector)
            self.materials.append(material)
            self.polarization_vectors.append(polarization_vector)
            self.update_energy()
            if phenomenon == 'Absortion':
                self.finished = True
            if phenomenon == 'Got_Absorbed':
                self.got_absorbed = True
                self.finished = True

    def add_to_document(self, doc):
        lshape_wire = Part.makePolygon(self.points)
        my_shape_ray = doc.addObject("Part::Feature", "Ray")
        my_shape_ray.Shape = lshape_wire


class LightSource:
    """
    Sets up a light source with a given scene, a given emitting region and a given light spectrum.
    The emitting region must provide the main direction.
    Light spectrum could be: a constant value (for a single wavelength in nanometers), or a spectrum distribution.
    The distribution (dispersion) for the main direction is provided in "direction_distribution".
    The polarization_vector is a Base.Vector for polarized light. If is not given unpolarized light is generated.
    """

    def __init__(self, scene, emitting_region, light_spectrum, initial_energy, direction_distribution = None,
                 polarization_vector = None):
        self.scene = scene
        self.emitting_region = emitting_region
        self.light_spectrum = light_spectrum
        self.initial_energy = initial_energy
        self.direction_distribution = direction_distribution
        self.polarization_vector = polarization_vector

    def emit_ray(self):
        point = self.emitting_region.random_point()
        main_direction = self.emitting_region.main_direction # emitting main direction
        direction = main_direction
        if self.direction_distribution is not None: # main direction has a distribution
            theta = self.direction_distribution.angle_distribution(myrandom())
            phi = 360.0 * myrandom()
            direction = dispersion_from_main_direction(main_direction, theta, phi)
            if self.polarization_vector: # single polarization vector is active
                polarization_vector = dispersion_polarization(main_direction,self.polarization_vector,theta,phi)
        if self.polarization_vector is None: # unpolarization is active 
            polarization_vector = random_polarization(direction) # random polarization from light direction		  
        if callable(self.light_spectrum):
            wavelength = pick_random(self.light_spectrum) # TODO light spectrum is active (nanometers)
        else:
            wavelength = self.light_spectrum # experiment with a single wavelength (nanometers)
        ray = Ray(self.scene,point,direction,{'wavelength':wavelength,
                                         'energy':self.initial_energy,
                                         'polarization_vector': polarization_vector})
        return ray		

		
class Experiment:
    """
    Sets up and runs and experiment in a given scene with a given light source.
    If show_in_doc is given, the emitting region is drawn in the FreeCAD active document.
    If show_in_doc is given, the rays could be drawn in the FreeCAD active document (using the run function).
    """

    def __init__(self, scene, light_source, number_of_rays, show_in_doc = None):
        self.scene = scene
        self.light_source = light_source
        if show_in_doc:
            self.light_source.emitting_region.add_to_document(show_in_doc)
        self.number_of_rays = number_of_rays
        self.captured_energy = 0
        random_congruential(time.time()) # TODO: change location

    def run(self, show_in_doc=None):
        for _ in xrange(self.number_of_rays):
            ray = self.light_source.emit_ray()
            ray.run()
            if show_in_doc:
                ray.add_to_document(show_in_doc)
            if ray.got_absorbed:
                self.captured_energy += ray.energy
