import json
import numpy as np
import zipfile
import dill
from .logging_unit import logger
from autologging import traced
from .optics import *
from .math import *


class NumpyEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return json.JSONEncoder.default(self, obj)


# region Materials

# ---
# Classes for materials
# ---

def plain_properties_to_properties(plain_properties):
    properties = {}
    for key in plain_properties:
        plain_property = plain_properties[key]
        prop_type = plain_property['type']
        prop_value = plain_property['value']
        if prop_type == 'scalar':
            properties[key] = prop_value
        if prop_type == 'constant':
            properties[key] = constant_function(prop_value)
        if prop_type == 'tabulated':
            properties[key] = tabulated_function(np.array(prop_value[0]), np.array(prop_value[1]))
        if prop_type == 'matrix':
            properties[key] = matrix_reflectance(np.array(prop_value))
    properties['plain_properties'] = plain_properties
    return properties


def properties_to_plain_properties(properties):
    return properties.get('plain_properties', None)


@traced(logger)
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

    def __init__(self, name, properties={}):
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
    def load_from_json_fileobject(cls, f):
        info = json.load(f)
        if type(info).__name__ == 'dict':
            info = [info]
        for mat_spec in info:
            kind = mat_spec['kind']
            name = mat_spec['name']
            if kind == 'Volume':
                mat = VolumeMaterial(name, {})
                plain_properties = mat_spec['plain_properties']
                properties = plain_properties_to_properties(plain_properties)
                mat.properties = properties
            if kind == 'Surface':
                mat = SurfaceMaterial(name, {})
                plain_properties_back = mat_spec['plain_properties_back']
                mat.properties_back = plain_properties_to_properties(plain_properties_back)
                plain_properties_front = mat_spec['plain_properties_front']
                mat.properties_front = plain_properties_to_properties(plain_properties_front)
        return name

    @classmethod
    def load_from_json_file(cls, filename):
        try:
            with open(filename, 'rb') as f:
                return cls.load_from_json_fileobject(f)
        except:
            logger.exception("error in processing file %s", filename)

    @classmethod
    def load_from_json_zip(cls, filename):
        try:
            with zipfile.ZipFile(filename) as z:
                for matfile in z.namelist():
                    with z.open(matfile) as f:
                        cls.load_from_json_fileobject(f)
        except:
            logger.exception("error in processing file %s", filename)

    @classmethod
    def load_from_file(cls, filename):
        try:
            with open(filename, 'rb') as f:
                mat = dill.load(f)
                cls.by_name[mat.name] = mat
                return mat.name
        except:
            logger.exception("error in processing file %s", filename)

    @classmethod
    def load_from_zipfile(cls, filename):
        with zipfile.ZipFile(filename) as z:
            for matfile in z.namelist():
                with z.open(matfile) as f:
                    try:
                        mat = dill.load(f)
                        cls.by_name[mat.name] = mat
                    except:
                        logger.exception("error in processing file %s", matfile)

    def change_of_direction(self, ray, normal_vector):
        """
        Computes the karma of the material


        """
        pass
        # Compute new direction (ray.current_material)

    def save_to_json_file(self, filename):
        with open(filename, 'w') as f:
            f.write(self.to_json())


@traced(logger)
class VolumeMaterial(Material, object):
    def __init__(self, name, properties={}):
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
        if 'Matrix_reflectance_thin_film' in ray.current_medium.properties:  # the ray is traveling inside thin film material and sure refraction
            ray.current_medium = self
            return OpticalState(ray.polarization_vectors[-1], ray.directions[-1], Phenomenon.REFRACTION)
        if 'Matrix_reflectance_thin_film' in self.properties:  # the ray impacts on thin film material
            n1 = ray.materials[-1].properties['index_of_refraction'](ray.wavelength)
            if 'extinction_coefficient' in ray.current_medium.properties:
                n1 = n1 + 1j * ray.materials[-1].properties['extinction_coefficient'](ray.wavelength)
            n_front = self.properties['index_of_refraction_front'](ray.wavelength)
            n_back = self.properties['index_of_refraction_back'](ray.wavelength)
            if n_front == ray.materials[-1].properties['index_of_refraction'](
                    ray.wavelength):  # impacts on front material
                n2 = n_back
            else:
                n2 = n_front
            energy_absorbed_thin_film, optical_state = calculate_state_thin_film(ray.directions[-1], normal_vector, n1,
                                                                                 n2,
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

    def to_json(self):
        return json.dumps(
            {
                'name': self.name,
                'kind': self.kind,
                'plain_properties': self.properties.get('plain_properties', None)
            }, cls=NumpyEncoder, indent=4
        )


# def create_simple_volume_material(name, index_of_refraction, attenuation_coefficient=None):
#     VolumeMaterial.create(name, {'index_of_refraction': constant_function(index_of_refraction),
#                                  'attenuation_coefficient': constant_function(attenuation_coefficient)})

def create_simple_volume_material(name, index_of_refraction, attenuation_coefficient=None):
    plain_properties = {
        'index_of_refraction': {
            'type': 'constant',
            'value': index_of_refraction
        },
        'attenuation_coefficient': {
            'type': 'constant',
            'value': attenuation_coefficient
        }
    }
    material = VolumeMaterial(name, {})
    material.properties = plain_properties_to_properties(plain_properties)
    return material


# def create_wavelength_volume_material(name, file_index_of_refraction):
#     # file_index_of_refraction with three columns: wavelenth in nm, real(index of refraction), imaginary(index of refraction)
#     data_refraction = np.loadtxt(file_index_of_refraction, usecols=(0, 1, 2))
#     wavelength_values = data_refraction[:, 0]
#     n_values = data_refraction[:, 1]
#     k_values = data_refraction[:, 2]
#     VolumeMaterial.create(name, {'index_of_refraction': tabulated_function(wavelength_values, n_values),
#                                  'extinction_coefficient': tabulated_function(wavelength_values, k_values)})

def create_wavelength_volume_material(name, file_index_of_refraction):
    # file_index_of_refraction with three columns: wavelenth in nm, real(index of refraction), imaginary(index of refraction)
    data_refraction = np.loadtxt(file_index_of_refraction, usecols=(0, 1, 2))
    wavelength_values = data_refraction[:, 0]
    n_values = data_refraction[:, 1]
    k_values = data_refraction[:, 2]
    plain_properties = {
        'index_of_refraction': {
            'type': 'tabulated',
            'value': [wavelength_values, n_values]
        },
        'extinction_coefficient': {
            'type': 'tabulated',
            'value': [wavelength_values, k_values]
        }
    }
    material = VolumeMaterial(name)
    material.properties = plain_properties_to_properties(plain_properties)
    return material


# def create_PV_material(name, file_index_of_refraction):
#     # file_index_of_refraction with three columns: wavelenth in nm, real(index of refraction), imaginary(index of refraction)
#     data_refraction = np.loadtxt(file_index_of_refraction, usecols=(0, 1, 2))
#     wavelength_values = data_refraction[:, 0]
#     n_values = data_refraction[:, 1]
#     k_values = data_refraction[:, 2]
#     VolumeMaterial.create(name, {'index_of_refraction': tabulated_function(wavelength_values, n_values),
#                                  'extinction_coefficient': tabulated_function(wavelength_values, k_values),
#                                  'PV_material' : True})

def create_PV_material(name, file_index_of_refraction):
    # file_index_of_refraction with three columns: wavelenth in nm, real(index of refraction), imaginary(index of refraction)
    data_refraction = np.loadtxt(file_index_of_refraction, usecols=(0, 1, 2))
    wavelength_values = data_refraction[:, 0]
    n_values = data_refraction[:, 1]
    k_values = data_refraction[:, 2]
    plain_properties = {
        'index_of_refraction': {
            'type': 'tabulated',
            'value': [wavelength_values, n_values]
        },
        'extinction_coefficient': {
            'type': 'tabulated',
            'value': [wavelength_values, k_values]
        },
        'PV_material': {
            'type': 'scalar',
            'value': True
        }
    }
    material = VolumeMaterial(name)
    material.properties = plain_properties_to_properties(plain_properties)
    return material


# def create_polarized_thin_film(name, file_thin_film, file_front, file_back):
#     # thin film material calculated by TMM method, six columns:
# 	# wavelenth in nm, angle in deg., reflectance s-polarized (perpendicular), reflectance p-polarized (parallel),  transmittance s-polarized, transmittance p-polarized
#     # the values in coating_material should be in the corresponding order columns
#     data = np.loadtxt(file_thin_film)
#     data_reflectance = data[:,[0,1,2,3]]
#     data_transmittance = data[:,[0,1,4,5]]
#     if file_front is not 'Vacuum':
#         data_refraction_front = np.loadtxt(file_front, usecols=(0, 1, 2))
#         wavelength_values_front = data_refraction_front[:, 0]
#         n_values_front = data_refraction_front[:, 1]
#         k_values_front = data_refraction_front[:, 2]
#         index_of_refraction_front = tabulated_function(wavelength_values_front, n_values_front)
#         extinction_coefficient_front = tabulated_function(wavelength_values_front, k_values_front)
#     else:
#         index_of_refraction_front = constant_function(1.0)
#         extinction_coefficient_front = constant_function(0.0)
#     if file_back is not 'Vacuum':
#         data_refraction_back = np.loadtxt(file_back, usecols=(0, 1, 2))
#         wavelength_values_back = data_refraction_back[:, 0]
#         n_values_back = data_refraction_back[:, 1]
#         k_values_back = data_refraction_back[:, 2]
#         index_of_refraction_back = tabulated_function(wavelength_values_back, n_values_back)
#         extinction_coefficient_back = tabulated_function(wavelength_values_back, k_values_back)
#     else:
#         index_of_refraction_back = constant_function(1.0)
#         extinction_coefficient_back = constant_function(0.0)
#     VolumeMaterial.create(name, {'Matrix_reflectance_thin_film': matrix_reflectance(data_reflectance),
#                                  'Matrix_transmittance_thin_film': matrix_reflectance(data_transmittance),
#                                  'index_of_refraction_front': index_of_refraction_front,
#                                  'extinction_coefficient_front': extinction_coefficient_front,
#                                  'index_of_refraction_back': index_of_refraction_back,
#                                  'extinction_coefficient_back': extinction_coefficient_back,
#                                  'thin_film' : True})

def create_polarized_thin_film(name, file_thin_film, file_front, file_back):
    # thin film material calculated by TMM method, six columns:
    # wavelenth in nm, angle in deg., reflectance s-polarized (perpendicular), reflectance p-polarized (parallel),  transmittance s-polarized, transmittance p-polarized
    # the values in coating_material should be in the corresponding order columns
    data = np.loadtxt(file_thin_film)
    data_reflectance = data[:, [0, 1, 2, 3]]
    data_transmittance = data[:, [0, 1, 4, 5]]
    if file_front is not 'Vacuum':
        data_refraction_front = np.loadtxt(file_front, usecols=(0, 1, 2))
        wavelength_values_front = data_refraction_front[:, 0]
        n_values_front = data_refraction_front[:, 1]
        k_values_front = data_refraction_front[:, 2]
        index_of_refraction_front = {
            'type': 'tabulated',
            'value': [wavelength_values_front, n_values_front]
        }
        extinction_coefficient_front = {
            'type': 'tabulated',
            'value': [wavelength_values_front, k_values_front]
        }
        # index_of_refraction_front = tabulated_function(wavelength_values_front, n_values_front)
        # extinction_coefficient_front = tabulated_function(wavelength_values_front, k_values_front)
    else:
        index_of_refraction_front = {
            'type': 'constant',
            'value': 1.0
        }
        extinction_coefficient_front = {
            'type': 'constant',
            'value': 0.0
        }
        # index_of_refraction_front = constant_function(1.0)
        # extinction_coefficient_front = constant_function(0.0)
    if file_back is not 'Vacuum':
        data_refraction_back = np.loadtxt(file_back, usecols=(0, 1, 2))
        wavelength_values_back = data_refraction_back[:, 0]
        n_values_back = data_refraction_back[:, 1]
        k_values_back = data_refraction_back[:, 2]
        index_of_refraction_back = {
            'type': 'tabulated',
            'value': [wavelength_values_back, n_values_back]
        }
        extinction_coefficient_back = {
            'type': 'tabulated',
            'value': [wavelength_values_back, k_values_back]
        }
        # index_of_refraction_back = tabulated_function(wavelength_values_back, n_values_back)
        # extinction_coefficient_back = tabulated_function(wavelength_values_back, k_values_back)
    else:
        index_of_refraction_back = {
            'type': 'constant',
            'value': 1.0
        }
        extinction_coefficient_back = {
            'type': 'constant',
            'value': 0.0
        }
        # index_of_refraction_back = constant_function(1.0)
        # extinction_coefficient_back = constant_function(0.0)
    plain_properties = {
        'Matrix_reflectance_thin_film': {
            'type': 'matrix',
            'value': data_reflectance
        },
        'Matrix_transmittance_thin_film': {
            'type': 'matrix',
            'value': data_transmittance
        },
        'index_of_refraction_front': index_of_refraction_front,
        'extinction_coefficient_front': extinction_coefficient_front,
        'index_of_refraction_back': index_of_refraction_back,
        'extinction_coefficient_back': extinction_coefficient_back,
        'thin_film': {
            'type': 'scalar',
            'value': True
        },
    }
    material = VolumeMaterial(name)
    material.properties = plain_properties_to_properties(plain_properties)
    return material


create_simple_volume_material("Vacuum", 1.0, 0.0)
# noinspection PyNoneFunctionAssignment
vacuum_medium = VolumeMaterial.by_name["Vacuum"]


@traced(logger)
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
            probabilities = [por, absortion, 0]  # Here I assume no transmitance
        if 'Matrix_polarized_reflectance_coating' in properties:  # polarized_coating_layer
            n1 = ray.current_medium.properties['index_of_refraction'](ray.wavelength)
            if 'extinction_coefficient' in ray.current_medium.properties:
                n1 = ray.current_medium.properties['index_of_refraction'](ray.wavelength) + 1j * \
                     ray.current_medium.properties['extinction_coefficient'](ray.wavelength)
            n2 = nearby_material.properties['index_of_refraction'](ray.wavelength)
            if 'extinction_coefficient' in nearby_material.properties:
                n2 = nearby_material.properties['index_of_refraction'](ray.wavelength) + 1j * \
                     nearby_material.properties['extinction_coefficient'](ray.wavelength)
            results = calculate_probabilities_polarizaton_coating(ray.directions[-1], normal_vector, n1, n2,
                                                                  ray.polarization_vectors[-1],
                                                                  properties, ray.wavelength)
            probabilities = [results[0], results[1], results[2]]
            polarization_vector = results[3]
            perpendicular_polarized = results[4]
        if 'metallic_material' in properties:
            n1 = ray.current_medium.properties['index_of_refraction'](ray.wavelength)
            if 'extinction_coefficient' in ray.current_medium.properties:
                n1 = ray.current_medium.properties['index_of_refraction'](ray.wavelength) + 1j * \
                     ray.current_medium.properties['extinction_coefficient'](ray.wavelength)
            n2 = properties['index_of_refraction'](ray.wavelength)
            if 'extinction_coefficient' in properties:
                n2 = properties['index_of_refraction'](ray.wavelength) + 1j * properties['extinction_coefficient'](
                    ray.wavelength)
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
            state = reflexion(ray.directions[-1], normal_vector, ray.polarization_vectors[-1],
                              polarization_vector_calculated_before)
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
            n1 = ray.current_medium.properties['index_of_refraction'](ray.wavelength) + 1j * \
                 ray.current_medium.properties['extinction_coefficient'](ray.wavelength)
        n2 = nearby_material.properties['index_of_refraction'](ray.wavelength)
        if 'extinction_coefficient' in nearby_material.properties:
            n2 = nearby_material.properties['index_of_refraction'](ray.wavelength) + 1j * nearby_material.properties[
                'extinction_coefficient'](ray.wavelength)
        if n1 == n2:  # transparent_simple_layer
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
        if ray.polarization_vectors[-1] == results[1]:  # polarization_vector not calculated
            polarization_vector_calculated_before = False
        else:  # polarization_vector calculated before
            polarization_vector_calculated_before = True
        ray.polarization_vectors[-1] = results[1]
        perpendicular_polarized = results[2]  # True or False

        if phenomenon == Phenomenon.REFLEXION:
            return self.change_of_direction_by_reflexion(ray, normal_vector, properties,
                                                         polarization_vector_calculated_before)
        elif phenomenon == Phenomenon.ABSORTION:
            return self.change_of_direction_by_absortion(ray, normal_vector, properties)
        elif phenomenon == Phenomenon.TRANSMITANCE:
            return self.change_of_direction_by_transmitance(ray, normal_vector, nearby_material,
                                                            perpendicular_polarized)

    def to_json(self):
        return json.dumps(
            {
                'name': self.name,
                'kind': self.kind,
                'plain_properties_back': self.properties_back.get('plain_properties', None),
                'plain_properties_front': self.properties_front.get('plain_properties', None),
            }, cls=NumpyEncoder
        )


# def create_opaque_simple_layer(name):
#     SurfaceMaterial.create(name, {'probability_of_reflexion': constant_function(0.0),
#                                   'probability_of_absortion': constant_function(1.0),
#                                   'probability_of_transmitance': constant_function(0.0),
#                                   'energy_collector': False})

def create_opaque_simple_layer(name):
    plain_properties = {
        'probability_of_reflexion': {
            'type': 'constant',
            'value': 0.0
        },
        'probability_of_absortion': {
            'type': 'constant',
            'value': 1.0
        },
        'probability_of_transmitance': {
            'type': 'constant',
            'value': 0.0
        },
        'energy_collector': {
            'type': 'scalar',
            'value': False
        }
    }
    properties = plain_properties_to_properties(plain_properties)
    material = SurfaceMaterial(name, properties, properties)
    return material


# def create_transparent_simple_layer(name, pot):
#     SurfaceMaterial.create(name, {'probability_of_reflexion': constant_function(1 - pot),
#                                   'probability_of_absortion': constant_function(0),
#                                   'probability_of_transmitance': constant_function(pot),
#                                   'energy_collector': False,
#                                   'specular_material': True})

def create_transparent_simple_layer(name, pot):
    plain_properties = {
        'probability_of_reflexion': {
            'type': 'constant',
            'value': 1 - pot
        },
        'probability_of_absortion': {
            'type': 'constant',
            'value': 0.0
        },
        'probability_of_transmitance': {
            'type': 'constant',
            'value': pot
        },
        'energy_collector': {
            'type': 'scalar',
            'value': False
        },
        'specular_material': {
            'type': 'scalar',
            'value': True
        }
    }
    properties = plain_properties_to_properties(plain_properties)
    material = SurfaceMaterial(name, properties, properties)
    return material


# def create_absorber_simple_layer(name, poa):
#     SurfaceMaterial.create(name, {'probability_of_reflexion': constant_function(1 - poa),
#                                   'probability_of_absortion': constant_function(poa),
#                                   'probability_of_transmitance': constant_function(0),
#                                   'energy_collector': True,
#                                   'specular_material': True})

def create_absorber_simple_layer(name, poa):
    plain_properties = {
        'probability_of_reflexion': {
            'type': 'constant',
            'value': 1 - poa
        },
        'probability_of_absortion': {
            'type': 'constant',
            'value': poa
        },
        'probability_of_transmitance': {
            'type': 'constant',
            'value': 0.0
        },
        'energy_collector': {
            'type': 'scalar',
            'value': True
        },
        'specular_material': {
            'type': 'scalar',
            'value': True
        }
    }
    properties = plain_properties_to_properties(plain_properties)
    material = SurfaceMaterial(name, properties, properties)
    return material


# def create_absorber_lambertian_layer(name, poa):
#     SurfaceMaterial.create(name, {'probability_of_reflexion': constant_function(1 - poa),
#                                   'probability_of_absortion': constant_function(poa),
#                                   'probability_of_transmitance': constant_function(0),
#                                   'energy_collector': True,
#                                   'lambertian_material': True})

def create_absorber_lambertian_layer(name, poa):
    plain_properties = {
        'probability_of_reflexion': {
            'type': 'constant',
            'value': 1 - poa
        },
        'probability_of_absortion': {
            'type': 'constant',
            'value': poa
        },
        'probability_of_transmitance': {
            'type': 'constant',
            'value': 0.0
        },
        'energy_collector': {
            'type': 'scalar',
            'value': True
        },
        'lambertian_material': {
            'type': 'scalar',
            'value': True
        }
    }
    properties = plain_properties_to_properties(plain_properties)
    material = SurfaceMaterial(name, properties, properties)
    return material


# def create_absorber_TW_model_layer(name, poa, b_constant, c_constant):
#     SurfaceMaterial.create(name, {'probability_of_reflexion': constant_function(1 - poa),
#                                   'probability_of_absortion': constant_function(poa),
#                                   'probability_of_transmitance': constant_function(0),
#                                   'energy_collector': True,
#                                   'lambertian_material': True,
#                                   'TW_model': True,
#                                   'b_constant': b_constant,
#                                   'c_constant': c_constant})

def create_absorber_TW_model_layer(name, poa, b_constant, c_constant):
    plain_properties = {
        'probability_of_reflexion': {
            'type': 'constant',
            'value': 1 - poa
        },
        'probability_of_absortion': {
            'type': 'constant',
            'value': poa
        },
        'probability_of_transmitance': {
            'type': 'constant',
            'value': 0.0
        },
        'energy_collector': {
            'type': 'scalar',
            'value': True
        },
        'lambertian_material': {
            'type': 'scalar',
            'value': True
        },
        'TW_model': {
            'type': 'scalar',
            'value': True
        },
        'b_constant': {
            'type': 'scalar',
            'value': b_constant
        },
        'c_constant': {
            'type': 'scalar',
            'value': c_constant
        }
    }
    properties = plain_properties_to_properties(plain_properties)
    material = SurfaceMaterial(name, properties, properties)
    return material


# def create_reflector_specular_layer(name, por, sigma_1 = None, sigma_2 = None, k = None):
#     SurfaceMaterial.create(name, {'probability_of_reflexion': constant_function(por),
#                                   'probability_of_absortion': constant_function(1 - por),
#                                   'probability_of_transmitance': constant_function(0),
#                                   'energy_collector': False,
#                                   'specular_material': True,
#                                   'sigma_1': sigma_1,
#                                   'sigma_2': sigma_2,
#                                   'k': k})

def create_reflector_specular_layer(name, por, sigma_1=None, sigma_2=None, k=None):
    plain_properties = {
        'probability_of_reflexion': {
            'type': 'constant',
            'value': por
        },
        'probability_of_absortion': {
            'type': 'constant',
            'value': 1 - por
        },
        'probability_of_transmitance': {
            'type': 'constant',
            'value': 0.0
        },
        'energy_collector': {
            'type': 'scalar',
            'value': False
        },
        'specular_material': {
            'type': 'scalar',
            'value': True
        },
        'sigma_1': {
            'type': 'scalar',
            'value': sigma_1
        },
        'sigma_2': {
            'type': 'scalar',
            'value': sigma_2
        },
        'k': {
            'type': 'scalar',
            'value': k
        }
    }
    properties = plain_properties_to_properties(plain_properties)
    material = SurfaceMaterial(name, properties, properties)
    return material


# def create_reflector_lambertian_layer(name, por):
#     SurfaceMaterial.create(name, {'probability_of_reflexion': constant_function(por),
#                                   'probability_of_absortion': constant_function(1 - por),
#                                   'probability_of_transmitance': constant_function(0),
#                                   'energy_collector': False,
#                                   'lambertian_material': True})

def create_reflector_lambertian_layer(name, por):
    plain_properties = {
        'probability_of_reflexion': {
            'type': 'constant',
            'value': por
        },
        'probability_of_absortion': {
            'type': 'constant',
            'value': 1 - por
        },
        'probability_of_transmitance': {
            'type': 'constant',
            'value': 0.0
        },
        'energy_collector': {
            'type': 'scalar',
            'value': False
        },
        'lambertian_material': {
            'type': 'scalar',
            'value': True
        },
    }
    properties = plain_properties_to_properties(plain_properties)
    material = SurfaceMaterial(name, properties, properties)
    return material


# def create_metallic_specular_layer(name, file_index_of_refraction, sigma_1 = None, sigma_2 = None, k = None):
#     # file_index_of_refraction with three columns: wavelenth in nm, real(index of refraction), imaginary(index of refraction)
#     data_refraction = np.loadtxt(file_index_of_refraction, usecols=(0, 1, 2))
#     wavelength_values = data_refraction[:, 0]
#     n_values = data_refraction[:, 1]
#     k_values = data_refraction[:, 2]
#     SurfaceMaterial.create(name, {'index_of_refraction': tabulated_function(wavelength_values, n_values),
#                                  'extinction_coefficient': tabulated_function(wavelength_values, k_values),
#                                   'energy_collector': False,
#                                   'specular_material': True,
#                                   'metallic_material': True,
#                                   'sigma_1': sigma_1,
#                                   'sigma_2': sigma_2,
#                                   'k': k})

def create_metallic_specular_layer(name, file_index_of_refraction, sigma_1=None, sigma_2=None, k=None):
    # file_index_of_refraction with three columns: wavelenth in nm, real(index of refraction), imaginary(index of refraction)
    data_refraction = np.loadtxt(file_index_of_refraction, usecols=(0, 1, 2))
    wavelength_values = data_refraction[:, 0]
    n_values = data_refraction[:, 1]
    k_values = data_refraction[:, 2]
    plain_properties = {
        'index_of_refraction': {
            'type': 'tabulated',
            'value': [wavelength_values, n_values]
        },
        'extinction_coefficient': {
            'type': 'tabulated',
            'value': [wavelength_values, k_values]
        },
        'energy_collector': {
            'type': 'scalar',
            'value': False
        },
        'specular_material': {
            'type': 'scalar',
            'value': True
        },
        'metallic_material': {
            'type': 'scalar',
            'value': True
        },
        'sigma_1': {
            'type': 'scalar',
            'value': sigma_1
        },
        'sigma_2': {
            'type': 'scalar',
            'value': sigma_2
        },
        'k': {
            'type': 'scalar',
            'value': k
        }
    }
    properties = plain_properties_to_properties(plain_properties)
    material = SurfaceMaterial(name, properties, properties)
    return material


# def create_metallic_lambertian_layer(name, file_index_of_refraction):
#     # file_index_of_refraction with three columns: wavelenth in nm, real(index of refraction), imaginary(index of refraction)
#     data_refraction = np.loadtxt(file_index_of_refraction, usecols=(0, 1, 2))
#     wavelength_values = data_refraction[:, 0]
#     n_values = data_refraction[:, 1]
#     k_values = data_refraction[:, 2]
#     SurfaceMaterial.create(name, {'index_of_refraction': tabulated_function(wavelength_values, n_values),
#                                  'extinction_coefficient': tabulated_function(wavelength_values, k_values),
#                                   'energy_collector': False,
#                                   'lambertian_material': True,
#                                   'metallic_material': True})

def create_metallic_lambertian_layer(name, file_index_of_refraction):
    # file_index_of_refraction with three columns: wavelenth in nm, real(index of refraction), imaginary(index of refraction)
    data_refraction = np.loadtxt(file_index_of_refraction, usecols=(0, 1, 2))
    wavelength_values = data_refraction[:, 0]
    n_values = data_refraction[:, 1]
    k_values = data_refraction[:, 2]
    plain_properties = {
        'index_of_refraction': {
            'type': 'tabulated',
            'value': [wavelength_values, n_values]
        },
        'extinction_coefficient': {
            'type': 'tabulated',
            'value': [wavelength_values, k_values]
        },
        'energy_collector': {
            'type': 'scalar',
            'value': False
        },
        'lambertian_material': {
            'type': 'scalar',
            'value': True
        },
        'metallic_material': {
            'type': 'scalar',
            'value': True
        },
    }
    properties = plain_properties_to_properties(plain_properties)
    material = SurfaceMaterial(name, properties, properties)
    return material


# def create_polarized_coating_reflector_layer(name, coating_file, sigma_1 = None, sigma_2 = None, k = None):
#     # coating_material with four columns: wavelenth in nm, angle in deg., reflectance s-polarized (perpendicular), reflectance p-polarized (parallel)
#     # the values in coating_material should be in the corresponding order columns
#     data_material = np.loadtxt(coating_file, usecols=(0,1,2,3))
#     SurfaceMaterial.create(name, {'Matrix_polarized_reflectance_coating': matrix_reflectance(data_material),
#                                   'probability_of_absortion': constant_function(0),
#                                   'energy_collector': False,
#                                   'specular_material': True,
#                                   'transparent_material': False,
#                                   'sigma_1': sigma_1,
#                                   'sigma_2': sigma_2,
#                                   'k': k})

def create_polarized_coating_reflector_layer(name, coating_file, sigma_1=None, sigma_2=None, k=None):
    # coating_material with four columns: wavelenth in nm, angle in deg., reflectance s-polarized (perpendicular), reflectance p-polarized (parallel)
    # the values in coating_material should be in the corresponding order columns
    data_material = np.loadtxt(coating_file, usecols=(0, 1, 2, 3))
    plain_properties = {
        'Matrix_polarized_reflectance_coating': {
            'type': 'matrix',
            'value': data_material
        },
        'probability_of_absortion': {
            'type': 'constant',
            'value': 0
        },
        'energy_collector': {
            'type': 'scalar',
            'value': False
        },
        'specular_material': {
            'type': 'scalar',
            'value': True
        },
        'transparent_material': {
            'type': 'scalar',
            'value': False
        },
        'sigma_1': {
            'type': 'scalar',
            'value': sigma_1
        },
        'sigma_2': {
            'type': 'scalar',
            'value': sigma_2
        },
        'k': {
            'type': 'scalar',
            'value': k
        }
    }
    properties = plain_properties_to_properties(plain_properties)
    material = SurfaceMaterial(name, properties, properties)
    return material


# def create_polarized_coating_transparent_layer(name, coating_file):
#     # coating_material with four columns: wavelenth in nm, angle in deg., reflectance s-polarized (perpendicular), reflectance p-polarized (parallel)
#     # the values in coating_material should be in the corresponding order columns
#     data_material = np.loadtxt(coating_file, usecols=(0,1,2,3))
#     SurfaceMaterial.create(name, {'Matrix_polarized_reflectance_coating': matrix_reflectance(data_material),
#                                   'probability_of_absortion': constant_function(0),
#                                   'energy_collector': False,
#                                   'specular_material': False,
#                                   'transparent_material': True})

def create_polarized_coating_transparent_layer(name, coating_file):
    # coating_material with four columns: wavelenth in nm, angle in deg., reflectance s-polarized (perpendicular), reflectance p-polarized (parallel)
    # the values in coating_material should be in the corresponding order columns
    data_material = np.loadtxt(coating_file, usecols=(0, 1, 2, 3))
    plain_properties = {
        'Matrix_polarized_reflectance_coating': {
            'type': 'matrix',
            'value': data_material
        },
        'probability_of_absortion': {
            'type': 'constant',
            'value': 0
        },
        'energy_collector': {
            'type': 'scalar',
            'value': False
        },
        'specular_material': {
            'type': 'scalar',
            'value': False
        },
        'transparent_material': {
            'type': 'scalar',
            'value': True
        },
    }
    properties = plain_properties_to_properties(plain_properties)
    material = SurfaceMaterial(name, properties, properties)
    return material


# def create_polarized_coating_absorber_layer(name, coating_file):
#     # coating_material with four columns: wavelenth in nm, angle in deg., reflectance s-polarized (perpendicular), reflectance p-polarized (parallel)
#     # the values in coating_material should be in the corresponding order columns
#     data_material = np.loadtxt(coating_file, usecols=(0,1,2,3))
#     SurfaceMaterial.create(name, {'Matrix_polarized_reflectance_coating': matrix_reflectance(data_material),
#                                   'probability_of_transmitance': constant_function(0),
#                                   'energy_collector': True,
#                                   'specular_material': False,
#                                   'transparent_material': False,
#                                   'lambertian_material': True})

def create_polarized_coating_absorber_layer(name, coating_file):
    # coating_material with four columns: wavelenth in nm, angle in deg., reflectance s-polarized (perpendicular), reflectance p-polarized (parallel)
    # the values in coating_material should be in the corresponding order columns
    data_material = np.loadtxt(coating_file, usecols=(0, 1, 2, 3))
    plain_properties = {
        'Matrix_polarized_reflectance_coating': {
            'type': 'matrix',
            'value': data_material
        },
        'probability_of_transmitance': {
            'type': 'constant',
            'value': 0
        },
        'energy_collector': {
            'type': 'scalar',
            'value': True
        },
        'specular_material': {
            'type': 'scalar',
            'value': False
        },
        'transparent_material': {
            'type': 'scalar',
            'value': False
        },
        'lambertian_material': {
            'type': 'scalar',
            'value': True
        },
    }
    properties = plain_properties_to_properties(plain_properties)
    material = SurfaceMaterial(name, properties, properties)
    return material


def create_two_layers_material(name, layer_front, layer_back):
    SurfaceMaterial.create(name, Material.by_name[layer_front].properties_front,
                           Material.by_name[layer_back].properties_back)

# endregion
