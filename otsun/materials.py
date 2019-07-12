"""Module otsun.materials for treating materials

The module relies on a basic class `Material` with two subclasses
`VolumeMaterial` and `SurfaceMaterial`, and several subclasses of them
for specific materials.
"""


import json
import zipfile
from FreeCAD import Base
from .optics import Phenomenon, OpticalState, reflexion, refraction, matrix_reflectance,\
    calculate_reflectance, simple_polarization_reflexion, simple_polarization_refraction, \
    simple_reflexion, shure_refraction, lambertian_reflexion
from .math import arccos, parallel_orthogonal_components, rad_to_deg, myrandom, normalize,\
    constant_function, correct_normal, tabulated_function
from numpy import sqrt
import numpy as np
from autologging import traced
from logging_unit import logger


class NumpyEncoder(json.JSONEncoder):
    """Wrapper to dump numpy arrays as json"""
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return json.JSONEncoder.default(self, obj)


@traced(logger)
class Material(object):
    """
    Class used to represent materials and their physical properties

    The `properties` dictionary holds the physical properties of the material.
    Its contents are specific to the kind of material.

    Parameters
    ----------
    name : str
        Name of the material
    properties : dict
        Dictionary with physical properties of the material

    Attributes
    ----------
    classname : str
        String identifying the class
    """

    by_name = {}
    """
    Dict that associates the name of each created material with the material itself
    """


    def __init__(self, name, properties=None):
        self.by_name[name] = self
        self.name = name
        self.classname = ""
        if properties is None:
            properties = {}
        self.properties = properties

    @staticmethod
    def plain_properties_to_properties(plain_properties):
        """
        Converts properties of a material in plain format (json) to internal format

        Parameters
        ----------
        plain_properties : dict

        Returns
        -------
        dict
        """
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
                properties[key] = tabulated_function(
                    np.array(prop_value[0]), np.array(prop_value[1]))
            if prop_type == 'matrix':
                properties[key] = matrix_reflectance(np.array(prop_value))
        properties['plain_properties'] = plain_properties
        return properties

    @staticmethod
    def properties_to_plain_properties(properties):
        """
        Converts properties of a material in internal format to plain (json ready) format

        Since the plain properties are stored in the internal format,
        no need for conversions

        Parameters
        ----------
        properties : dict

        Returns
        -------
        dict
        """
        return properties.get('plain_properties', None)

    @classmethod
    def get_from_label(cls, label):
        """
        Returns the material given its label

        Given a `label` of an object (from a FreeCad document) of the form
        XXX(MatYYY), returns the material whose name is MatYYY

        Parameters
        ----------
        label : str

        Returns
        -------
        Material
        """
        if ("(" not in label) or (")" not in label):
            return None
        start = label.find("(")
        end = label.find(")")
        name = label[start + 1:end]
        return cls.by_name.get(name, None)

    @classmethod
    def create(cls, name, properties):
        """Wrapper to create a material"""
        _ = cls(name, properties)

    @classmethod
    def load_from_json(cls, info):
        if type(info).__name__ == 'dict':
            info = [info]
        names = []
        for mat_spec in info:
            classname = mat_spec['classname']
            logger.debug(classname)
            the_class = globals()[classname]
            name = mat_spec['name']
            names.append(name)
            if issubclass(the_class, TwoLayerMaterial):
                name_front_layer = mat_spec['name_front_layer']
                name_back_layer = mat_spec['name_back_layer']
                mat = TwoLayerMaterial(name, name_front_layer, name_back_layer)
            else:
                plain_properties = mat_spec['plain_properties']
                properties = the_class.plain_properties_to_properties(plain_properties)
                mat = Material(name, properties)
            mat.__class__ = the_class
        if len(names) == 1:
            return names[0]
        else:
            return names

    @classmethod
    def load_from_json_fileobject(cls, f):
        """
        Load materials from a json fileobject

        If the file contains a single dict, then it means that it contains a
        single material. Otherwise it contains an array, where each entry is a dict
        representing a material.

        Parameters
        ----------
        f : BinaryIO
            File object

        Returns
        -------
        str
            String with the name of the last material imported from the file
        """
        info = json.load(f)
        return cls.load_from_json(info)

    @classmethod
    def load_from_json_file(cls, filename):
        """
        Load materials from a json file

        Parameters
        ----------
        filename : str
            Name of the file

        Returns
        -------
        str
            String with the name of the last material imported from the file
        """
        try:
            with open(filename, 'rb') as f:
                return cls.load_from_json_fileobject(f)
        except IOError:
            logger.exception("error in processing file %s", filename)

    @classmethod
    def load_from_json_zip(cls, filename):
        """
        Load all materials from a zip file
        Parameters
        ----------
        filename : str
            Name of the file
        """
        try:
            with zipfile.ZipFile(filename) as z:
                for matfile in z.namelist():
                    try:
                        with z.open(matfile) as f:
                            cls.load_from_json_fileobject(f)
                    except:
                        logger.exception("File %s in zip contains errors", matfile)
        except IOError:
            logger.exception("error in processing file %s", filename)


    def get_n(self, wavelength):
        """
        Returns the (complex) refractive index at a certain wavelength

        Parameters
        ----------
        wavelength : float

        Returns
        -------
            complex
        """
        n = self.properties['index_of_refraction'](wavelength)
        if 'extinction_coefficient' in self.properties:
            kappaf = self.properties['extinction_coefficient']
            kappa = kappaf(wavelength)
            return n + 1j * kappa
        else:
            return n


    def change_of_optical_state(self, *args):
        """
        Computes how a ray behaves when interacting with the material.
        MUST be subclassed

        The material where the ray is actually located is held in
        ray.current_material.

        Parameters
        ----------
        args
            Variable length argument list

        Returns
        -------
        OpticalState
        """
        pass

    def to_json(self):
        """Converts material to json. MUST be subclassed"""
        return ""

    def save_to_json_file(self, filename):
        """
        Save material to json file

        Parameters
        ----------
        filename : str
            Name of the file
        """
        with open(filename, 'w') as f:
            f.write(self.to_json())

    @staticmethod
    def all_to_json():
        """
        Convert all materials to json

        Returns
        -------
        list of str
        """
        materials = Material.by_name.values()
        simple_mats = [material for material in materials if
                       not isinstance(material, TwoLayerMaterial)]
        composite_mats = [material for material in materials if
                          isinstance(material, TwoLayerMaterial)]
        all_mats = simple_mats + composite_mats
        return [mat.to_json() for mat in all_mats]

    @staticmethod
    def save_all_to_json_file(filename):
        """
        Saves all materials to text file

        Parameters
        ----------
        filename
        """
        with open(filename, 'w') as f:
            f.write('[\n')
            f.write(',\n'.join(Material.all_to_json()))
            f.write('\n]')

@traced(logger)
class VolumeMaterial(Material):
    """
    Subclass of Volume for materials with volume

    The `properties` parameter must be a dict with the physical properties
    describing the material. At least, the following must be provided:
    'index_of_refraction': float (index of refraction of the material, as a function
    of its wavelength, only real part)
    'extinction_coefficient': float (imaginary part of the index of refraction
    of the material as a function of its wavelength)
    """
    def __init__(self, name, properties=None):
        super(VolumeMaterial, self).__init__(name, properties)

    def change_of_optical_state(self, ray, normal_vector):
        """
        Compute the change of optical state

        Computes the new optical state when `ray` hits the material
        at a point with given `normal_vector`

        # TODO: Change name
        Parameters
        ----------
        ray : Ray
        normal_vector : Base.Vector

        Returns
        -------
        OpticalState
        """
        wavelength = ray.wavelength
        if isinstance(ray.current_medium(), PolarizedThinFilm):
            return OpticalState(ray.current_polarization(),
                                ray.current_direction(), Phenomenon.REFRACTION, self)
        else:
            n1 = ray.current_medium().get_n(wavelength)
            n2 = self.get_n(wavelength)
            optical_state = refraction(
                ray.current_direction(), normal_vector, n1, n2,
                ray.current_polarization())
            if optical_state.phenomenon == Phenomenon.REFRACTION:
                optical_state.material = self
            else:
                optical_state.material = ray.current_medium()
            return optical_state

    def to_json(self):
        return json.dumps(
            {
                'name': self.name,
                'classname': self.__class__.__name__,
                'plain_properties': self.properties.get(
                    'plain_properties', None)
            }, cls=NumpyEncoder, indent=4
        )


class SimpleVolumeMaterial(VolumeMaterial):
    """
    # TODO: Document
    """
    def __init__(self, name, index_of_refraction, attenuation_coefficient=None):
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
        super(SimpleVolumeMaterial, self).__init__(name, {})
        self.properties = Material.plain_properties_to_properties(plain_properties)


class WavelengthVolumeMaterial(VolumeMaterial):
    """
    # TODO: Document
    """
    def __init__(self, name, file_index_of_refraction):
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
        super(WavelengthVolumeMaterial, self).__init__(name)
        self.properties = Material.plain_properties_to_properties(plain_properties)


class PVMaterial(VolumeMaterial):
    """
    # TODO: Document
    """
    def __init__(self, name, file_index_of_refraction):
        # file_index_of_refraction with three columns:
        # wavelenth in nm,
        # real(index of refraction),
        # imaginary(index of refraction)
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
        super(PVMaterial, self).__init__(name)
        self.properties = Material.plain_properties_to_properties(plain_properties)

    def get_PV_data(self, ray, energy_before):
        """

        Parameters
        ----------
        ray : Ray
            Ray that has passed through the PV material
        energy_before : float
            Energy of the ray before passing through the PV material

        Returns
        -------
        float, tuple of floats
        """
        alpha = self.properties['extinction_coefficient'](
            ray.wavelength) * 4 * np.pi / (ray.wavelength / 1E6)  # mm-1
        # TODO: @Ramon: Revisar angle (sona raro)
        angle_incident = arccos(
            - ray.last_normal.dot(ray.current_direction())) * 180.0 / np.pi
        point_1 = ray.points[-1]
        point_2 = ray.points[-2]
        return (energy_before - ray.energy,
                (point_2.x, point_2.y, point_2.z, point_1.x, point_1.y, point_1.z,
                 energy_before, ray.energy, ray.wavelength, alpha, angle_incident)
                )

class PolarizedThinFilm(VolumeMaterial):
    """
    # TODO: Document
    """
    def __init__(self, name, file_thin_film, file_front, file_back):
        # thin film material calculated by TMM method, six columns:
        # wavelenth in nm, angle in deg.,
        # reflectance s-polarized (perpendicular),
        # reflectance p-polarized (parallel),  transmittance s-polarized,
        # transmittance p-polarized
        # the values in coating_material should be in the corresponding
        # order columns
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
        else:
            index_of_refraction_front = {
                'type': 'constant',
                'value': 1.0
            }
            extinction_coefficient_front = {
                'type': 'constant',
                'value': 0.0
            }
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
        else:
            index_of_refraction_back = {
                'type': 'constant',
                'value': 1.0
            }
            extinction_coefficient_back = {
                'type': 'constant',
                'value': 0.0
            }
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
        super(PolarizedThinFilm, self).__init__(name)
        self.properties = Material.plain_properties_to_properties(plain_properties)

    @staticmethod
    def calculate_state_thin_film(incident, normal_vector, n1, n2, polarization_vector,
                                  properties, wavelength):
        """
        # TODO: Document

        Parameters
        ----------
        incident
        normal_vector
        n1
        n2
        polarization_vector
        properties
        wavelength

        Returns
        -------

        """
        # TODO: document
        # returns optical state of the ray in thin film material
        normal = correct_normal(normal_vector, incident)
        backside = False
        if normal != normal_vector:
            backside = True
        r = n1 / n2
        c1 = - normal.dot(incident)
        # cos (incident_angle)
        c2sq = 1.0 - r * r * (1.0 - c1 * c1)
        # cos (refracted_angle) ** 2
        if c2sq.real < 0:
            # total internal reflection
            state = reflexion(incident, normal, polarization_vector)
            return (0.0, state)
            # no energy is abosrbed in the thinfilm
        c2 = sqrt(c2sq)
        # cos (refracted_angle)
        parallel_v, perpendicular_v, normal_parallel_plane = \
            parallel_orthogonal_components(polarization_vector, incident, normal)
        # parallel and perpendicular components of polarization vector
        # and orthogonal vector of the parallel plane
        ref_per = perpendicular_v.Length ** 2.0 / polarization_vector.Length ** 2.0
        # weight of perpendicular component: 0 < ref_per < 1
        if backside:
            # Ray intercepted on the backside of the transparent surface
            inc_angle = rad_to_deg(arccos(c2.real))
        else:
            inc_angle = rad_to_deg(arccos(c1))
        reflectance_matrix = properties['Matrix_reflectance_thin_film']
        r_matrix = reflectance_matrix(inc_angle, wavelength)
        # reflectance dependent of incidence angle and wavelength
        # We decide the polarization projection onto the parallel / perpendicular plane
        if myrandom() < ref_per:
            reflectance = calculate_reflectance(r_matrix, inc_angle, wavelength)[0]
            # reflectance for s-polarized (perpendicular) light
            perpendicular_polarized = True
            polarization_vector = normalize(perpendicular_v)
        else:
            reflectance = calculate_reflectance(r_matrix, inc_angle, wavelength)[1]
            # reflectance for p-polarized (parallel) light
            perpendicular_polarized = False
            polarization_vector = normalize(parallel_v)
        if myrandom() < reflectance:
            # ray reflected
            reflected_direction = simple_reflexion(incident, normal).normalize()
            if not perpendicular_polarized:
                # reflexion changes the parallel component of incident polarization
                polarization_vector = simple_polarization_reflexion(
                    incident, normal, normal_parallel_plane, polarization_vector)
            return 0.0, OpticalState(polarization_vector, reflected_direction, Phenomenon.REFLEXION)
        else:
            # ray refracted: computing the refracted direction and energy absorbed in the thinfilm
            transmittance_matrix = properties['Matrix_transmittance_thin_film']
            t_matrix = transmittance_matrix(inc_angle, wavelength)
            # transmittance dependent of incidence angle and wavelength
            if perpendicular_polarized:
                transmittance = calculate_reflectance(t_matrix, inc_angle, wavelength)[0]
            else:
                transmittance = calculate_reflectance(t_matrix, inc_angle, wavelength)[1]
            factor_energy_absorbed = (1 - reflectance - transmittance) / (1 - reflectance)
            refracted_direction = incident * r.real + \
                                  normal * (r.real * c1 - c2.real)
            refracted_direction.normalize()
            if not perpendicular_polarized:
                # refraction changes the parallel component of incident polarization
                polarization_vector = \
                    simple_polarization_refraction(
                        incident, normal, normal_parallel_plane, c2, polarization_vector)
            return (factor_energy_absorbed,
                    OpticalState(polarization_vector, refracted_direction, Phenomenon.REFRACTION))

    def change_of_optical_state(self, ray, normal_vector):
        """
        # TODO: Document
        Parameters
        ----------
        ray
        normal_vector

        Returns
        -------

        """
        # the ray impacts on thin film material
        n1 = ray.current_medium().get_n(ray.wavelength)
        n_front = self.properties['index_of_refraction_front'](ray.wavelength)
        k_front = self.properties['extinction_coefficient_front'](ray.wavelength)
        n_back = self.properties['index_of_refraction_back'](ray.wavelength)
        k_back = self.properties['extinction_coefficient_back'](ray.wavelength)
        properties = self.properties
        if n1 == n_front + 1j * k_front:
            n2 = n_back + 1j * k_back
        else:
            n2 = n_front + 1j * k_front
        factor_energy_absorbed, optical_state = (
            self.calculate_state_thin_film(
                ray.current_direction(), normal_vector, n1, n2,
                ray.current_polarization(),
                properties, ray.wavelength))
        optical_state.extra_data['factor_energy_absorbed'] = \
            factor_energy_absorbed
        if optical_state.phenomenon == Phenomenon.REFRACTION:
            optical_state.material = self
        else:
            optical_state.material = ray.current_medium()
        return optical_state


vacuum_medium = SimpleVolumeMaterial("Vacuum", 1.0, 0.0)


@traced(logger)
class SurfaceMaterial(Material):
    """
    # TODO: Document
    """
    def __init__(self, name, properties):
        """
        Initializes a Surface Material. The properties parameter must be a
        dict with the physical properties
        describing the material. At least, the following must be provided:
        'probability_of_reflexion': probability that a photon gets reflected,
        as a function of its wavelength.
        'probability_of_absortion': probability that a photon gets absorbed,
        as a function of its wavelength.
        'probability_of_transmitance': probability that a photon passes through
        the material, as a function of its wavelength.
        """
        super(SurfaceMaterial, self).__init__(name, properties)
        self.properties = properties

    @classmethod
    def create(cls, name, properties):
        _ = cls(name, properties)


    def compute_probabilities_and_polarizations(self, ray, normal_vector, nearby_material):
        """
        # TODO: Document

        Parameters
        ----------
        ray
        normal_vector
        nearby_material

        Returns
        -------

        """
        properties = self.properties
        try:
            por = properties['probability_of_reflexion'](ray.wavelength)
        except KeyError:
            por = 1.0
        try:
            poa = properties['probability_of_absortion'](ray.wavelength)
        except KeyError:
            poa = 1 - por
        try:
            pot = properties['probability_of_transmitance'](ray.wavelength)
        except KeyError:
            pot = 0.0

        return [por, poa, pot], ray.current_polarization(), False


    def decide_phenomenon(self, ray, normal_vector, nearby_material):
        """
        # TODO: Document

        Parameters
        ----------
        ray
        normal_vector
        nearby_material

        Returns
        -------

        """
        phenomena = [
            Phenomenon.REFLEXION,
            Phenomenon.ABSORPTION,
            Phenomenon.TRANSMITTANCE]

        probabilities, polarization_vector, perpendicular_polarized = \
            self.compute_probabilities_and_polarizations(ray, normal_vector,
                                                         nearby_material)
        phenomenon = np.random.choice(phenomena, 1, p=probabilities)[0]
        return phenomenon, polarization_vector, perpendicular_polarized

    def change_of_optical_state_by_absortion(self, ray, normal_vector):
        """
        # TODO: Document

        Parameters
        ----------
        ray
        normal_vector

        Returns
        -------

        """
        properties = self.properties
        if properties['energy_collector']:
            return (OpticalState(Base.Vector(0.0, 0.0, 0.0),
                                 Base.Vector(0.0, 0.0, 0.0),
                                 Phenomenon.ENERGY_ABSORBED,
                                 self))
        else:
            return (OpticalState(Base.Vector(0.0, 0.0, 0.0),
                                 Base.Vector(0.0, 0.0, 0.0),
                                 Phenomenon.ABSORPTION,
                                 self))

    def change_of_optical_state_by_reflexion(self, ray, normal_vector,
                                             polarization_vector_calculated_before):
        """
        # TODO @Ramon: Mirar l'embull del polarization_vector_calculated_before
        # TODO: Document

        Parameters
        ----------
        ray
        normal_vector
        polarization_vector_calculated_before

        Returns
        -------

        """
        properties = self.properties
        if 'specular_material' in properties:
            state = reflexion(ray.current_direction(), normal_vector,
                              ray.current_polarization(),
                              polarization_vector_calculated_before)
            state.material = ray.current_medium()
            state.apply_dispersion(properties, normal_vector)
            return state
        if 'lambertian_material' in properties:
            state = lambertian_reflexion(ray.current_direction(), normal_vector)
            state.material = ray.current_medium()
            return state

    def change_of_optical_state_by_transmitance(self, ray, normal_vector,
                                                nearby_material,
                                                perpendicular_polarized):
        """
        # TODO: Document

        Parameters
        ----------
        ray
        normal_vector
        nearby_material
        perpendicular_polarized

        Returns
        -------

        """
        n1 = ray.current_medium().get_n(ray.wavelength)
        n2 = nearby_material.get_n(ray.wavelength)
        if n1 == n2:  # transparent_simple_layer
            state = OpticalState(ray.current_polarization(),
                                 ray.current_direction(), Phenomenon.REFRACTION, nearby_material)
        else:
            state = shure_refraction(ray.current_direction(), normal_vector, n1, n2,
                                     ray.current_polarization(),
                                     perpendicular_polarized)
            state.material = nearby_material
        return state

    def change_of_optical_state(self, ray, normal_vector, nearby_material):
        """
        # TODO: Document

        Parameters
        ----------
        ray
        normal_vector
        nearby_material

        Returns
        -------

        """
        material = self
        results = material.decide_phenomenon(ray, normal_vector, nearby_material)
        phenomenon = results[0]
        if ray.current_polarization() == results[1]:
            # polarization_vector not calculated
            polarization_vector_calculated_before = False
        else:
            # polarization_vector calculated before
            polarization_vector_calculated_before = True
        ## ray.polarization_vectors[-1] = results[1]
        # TODO: Caution!!!! previous line!
        perpendicular_polarized = results[2]  # True or False
        if phenomenon == Phenomenon.REFLEXION:
            state = material.change_of_optical_state_by_reflexion(
                ray, normal_vector,
                polarization_vector_calculated_before)
            return state
        elif phenomenon == Phenomenon.ABSORPTION:
            state = material.change_of_optical_state_by_absortion(
                ray, normal_vector)
            return state
        elif phenomenon == Phenomenon.TRANSMITTANCE:
            state = material.change_of_optical_state_by_transmitance(
                ray, normal_vector, nearby_material,
                perpendicular_polarized)
            return state

    @classmethod
    def from_plain_properties(cls, name, plain_properties):
        """
        # TODO: Document

        Parameters
        ----------
        name : str
        plain_properties : dict

        Returns
        -------
            Material
        """
        properties = Material.plain_properties_to_properties(plain_properties)
        material = cls(name, properties)
        return material

    def to_json(self):
        return json.dumps(
            {
                'name': self.name,
                'classname': self.__class__.__name__,
                'plain_properties': self.properties.get(
                    'plain_properties', None),
            }, cls=NumpyEncoder
        )


@traced(logger)
class OpaqueSimpleLayer(SurfaceMaterial):
    """
    # TODO: Document
    """
    def __init__(self, name):
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
        properties = Material.plain_properties_to_properties(plain_properties)
        super(OpaqueSimpleLayer, self).__init__(name, properties)

    def change_of_optical_state(self, ray, normal_vector, nearby_material):
        # absorption in opaque material: the ray is killed
        return OpticalState(Base.Vector(0.0, 0.0, 0.0),
                            Base.Vector(0.0, 0.0, 0.0),
                            Phenomenon.ABSORPTION, self)


@traced(logger)
class TransparentSimpleLayer(SurfaceMaterial):
    """
    # TODO: Document
    """
    def __init__(self, name, pot):
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
        properties = Material.plain_properties_to_properties(plain_properties)
        super(TransparentSimpleLayer, self).__init__(name, properties)

    def change_of_optical_state(self, ray, normal_vector, nearby_material):
        properties = self.properties
        reflectance = properties['probability_of_reflexion'](ray.wavelength)
        if myrandom() < reflectance:
            polarization_vector = ray.current_polarization()
            incident = ray.current_direction()
            state = reflexion(incident, normal_vector, polarization_vector, False)
            state.material = ray.current_medium()
            state.apply_dispersion(properties, normal_vector)
            return state
        else:
            # refraction in transparent layer
            n1 = ray.current_medium().get_n(ray.wavelength)
            n2 = nearby_material.get_n(ray.wavelength)
            if n1 == n2:  # transparent_simple_layer
                state = OpticalState(ray.current_polarization(),
                                     ray.current_direction(),
                                     Phenomenon.REFRACTION,
                                     nearby_material)
            else:
                state = shure_refraction(ray.current_direction(), normal_vector,
                                         n1, n2,ray.current_polarization())
                state.material = nearby_material
            return state


@traced(logger)
class AbsorberSimpleLayer(SurfaceMaterial):
    """
    # TODO: Document
    """
    def __init__(self, name, poa):
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
        properties = Material.plain_properties_to_properties(plain_properties)
        super(AbsorberSimpleLayer, self).__init__(name, properties)

    def change_of_optical_state(self, ray, normal_vector, nearby_material):
        properties = self.properties
        absortion = properties['probability_of_absortion'](ray.wavelength)
        reflectance = 1.0 - absortion
        if myrandom() < reflectance:
            state = reflexion(ray.current_direction(),
                              normal_vector, ray.current_polarization(), False)
            state.material = ray.current_medium()
            return state
        else:
            # absorption in absorber material: the energy is abosrbed
            return OpticalState(Base.Vector(0.0, 0.0, 0.0),
                                Base.Vector(0.0, 0.0, 0.0),
                                Phenomenon.ENERGY_ABSORBED, self)


@traced(logger)
class AbsorberLambertianLayer(SurfaceMaterial):
    """
    # TODO: Document
    """
    def __init__(self, name, poa):
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
        properties = Material.plain_properties_to_properties(plain_properties)
        super(AbsorberLambertianLayer, self).__init__(name, properties)

    def change_of_optical_state(self, ray, normal_vector, nearby_material):
        properties = self.properties
        absortion = properties['probability_of_absortion'](ray.wavelength)
        reflectance = 1.0 - absortion
        if myrandom() < reflectance:
            polarization_vector = ray.current_polarization()
            state = lambertian_reflexion(ray.current_direction(), normal_vector)
            state.material = ray.current_medium()
            return state
        else:
            # absorptiontion in absorber material: the ray is killed
            return OpticalState(Base.Vector(0.0, 0.0, 0.0),
                                Base.Vector(0.0, 0.0, 0.0),
                                Phenomenon.ENERGY_ABSORBED, self)

@traced(logger)
class AbsorberTWModelLayer(SurfaceMaterial):
    """
    # TODO: Document
    """
    def __init__(self, name, poa, b_constant, c_constant):
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
        properties = Material.plain_properties_to_properties(plain_properties)
        super(AbsorberTWModelLayer, self).__init__(name, properties)

    @staticmethod
    def tw_absorptance_ratio(normal_vector, b_constant, c_constant, incident):
        """Angular Solar Absorptance model for selective absorber material.

        Given by the formula 1 - b * (1/cos - 1) ** c, based on:
        Tesfamichael, T., and Wackelgard, E., 2000, "Angular Solar Absorptance and
        Incident Angle Modifier of Selective Absorbers for Solar Thermal Collectors,"
        Sol. Energy, 68, pp. 335-341.

        Parameters
        ----------
        normal : Base.Vector
            normal vector of the surface at the point of incidence
        b_constant : float
        c_constant : float
        incident : Base.Vector
            direction vector of the incident ray

        Returns
        -------
        float

        """
        # We assume the normal is normalized.
        normal = correct_normal(normal_vector, incident)
        c1 = - normal.dot(incident)
        inc_angle = rad_to_deg(arccos(c1))
        # incidence angle
        if inc_angle < 80.0:
            absorption_ratio = 1.0 - b_constant * abs((1.0 / c1 - 1.0)) ** c_constant
        else:
            y0 = 1.0 - b_constant * (1.0 / np.cos(80.0 * np.pi / 180.0) - 1.0) ** c_constant
            m = y0 / 10.0
            absorption_ratio = y0 - m * (inc_angle - 80.0)
        return absorption_ratio

    def change_of_optical_state(self, ray, normal_vector, nearby_material):
        properties = self.properties
        b_constant = properties['b_constant']
        c_constant = properties['c_constant']
        absortion_ratio = self.tw_absorptance_ratio(
            normal_vector, b_constant, c_constant, ray.current_direction())
        absortion = properties['probability_of_absortion'](ray.wavelength) * absortion_ratio
        reflectance = 1.0 - absortion
        if myrandom() < reflectance:
            polarization_vector = ray.current_polarization()
            state = reflexion(ray.current_direction(), normal_vector, polarization_vector, False)
            state.material = ray.current_medium()
            return state
        else:
            # absorptiontion in absorber material: the ray is killed
            return OpticalState(Base.Vector(0.0, 0.0, 0.0),
                                Base.Vector(0.0, 0.0, 0.0),
                                Phenomenon.ENERGY_ABSORBED, self)

@traced(logger)
class ReflectorSpecularLayer(SurfaceMaterial):
    """
    # TODO: Document
    """
    def __init__(self, name, por, sigma_1=None, sigma_2=None, k=None):
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
        properties = Material.plain_properties_to_properties(plain_properties)
        super(ReflectorSpecularLayer, self).__init__(name, properties)

    def change_of_optical_state(self, ray, normal_vector, nearby_material):
        properties = self.properties
        reflectance = properties['probability_of_reflexion'](ray.wavelength)
        if myrandom() < reflectance:
            polarization_vector = ray.current_polarization()
            incident = ray.current_direction()
            state = reflexion(incident, normal_vector, polarization_vector, False)
            state.material = ray.current_medium()
            state.apply_dispersion(properties, normal_vector)
            return state
        else:
            # refraction in metallic layer: the ray is killed
            return OpticalState(Base.Vector(0.0, 0.0, 0.0),
                                Base.Vector(0.0, 0.0, 0.0),
                                Phenomenon.ABSORPTION, self)

@traced(logger)
class ReflectorLambertianLayer(SurfaceMaterial):
    """
    # TODO: Document
    """
    def __init__(self, name, por):
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
        properties = Material.plain_properties_to_properties(plain_properties)
        super(ReflectorLambertianLayer, self).__init__(name, properties)

    def change_of_optical_state(self, ray, normal_vector, nearby_material):
        properties = self.properties
        reflectance = properties['probability_of_reflexion'](ray.wavelength)
        if myrandom() < reflectance:
            incident = ray.current_direction()
            state = lambertian_reflexion(incident, normal_vector)
            state.material = ray.current_medium()
            return state
        else:
            # refraction in metallic layer: the ray is killed
            return OpticalState(Base.Vector(0.0, 0.0, 0.0),
                                Base.Vector(0.0, 0.0, 0.0),
                                Phenomenon.ABSORPTION, self)


@traced(logger)
class MetallicLayer(SurfaceMaterial):
    """
    # TODO: Document
    """
    # TODO: material obsolet o be que sigui el pare de
    # MetallicSpecularLayer i MetallicLambertianLayer


@traced(logger)
class MetallicSpecularLayer(MetallicLayer):
    """
    # TODO: Document
    """
    def __init__(self, name, file_index_of_refraction,
                 sigma_1=None, sigma_2=None, k=None):
        # file_index_of_refraction with three columns: wavelenth in nm,
        # real(index of refraction), imaginary(index of refraction)
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
        properties = Material.plain_properties_to_properties(plain_properties)
        super(MetallicSpecularLayer, self).__init__(name, properties)

    def change_of_optical_state(self, ray, normal_vector, nearby_material):
        properties = self.properties
        polarization_vector = ray.current_polarization()
        n1 = ray.current_medium().get_n(ray.wavelength)
        n2 = self.get_n(ray.wavelength)
        incident = ray.current_direction()
        state = refraction(incident, normal_vector, n1, n2, polarization_vector)
        if state.phenomenon == Phenomenon.REFLEXION:
            state.material = ray.current_medium()
            state.apply_dispersion(properties, normal_vector)
            return state
        if state.phenomenon == Phenomenon.REFRACTION:
            # refraction in metallic layer: the ray is killed
            return OpticalState(Base.Vector(0.0, 0.0, 0.0),
                                Base.Vector(0.0, 0.0, 0.0),
                                Phenomenon.ABSORPTION, self)


@traced(logger)
class MetallicLambertianLayer(MetallicLayer):
    """
    # TODO: Document
    """
    def __init__(self, name, file_index_of_refraction):
        # file_index_of_refraction with three columns:
        # wavelenth in nm, real(index of refraction),
        # imaginary(index of refraction)
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
        properties = Material.plain_properties_to_properties(plain_properties)
        super(MetallicLambertianLayer, self).__init__(name, properties)

    def change_of_optical_state(self, ray, normal_vector, nearby_material):
        polarization_vector = ray.current_polarization()
        n1 = ray.current_medium().get_n(ray.wavelength)
        n2 = self.get_n(ray.wavelength)
        incident = ray.current_direction()
        state = refraction(incident, normal_vector, n1, n2, polarization_vector, True)
        if state.phenomenon == Phenomenon.REFLEXION:
            state.material = ray.current_medium()
            return state
        if state.phenomenon == Phenomenon.REFRACTION:
            # refraction in metallic layer: the ray is killed
            return OpticalState(Base.Vector(0.0, 0.0, 0.0),
                                Base.Vector(0.0, 0.0, 0.0),
                                Phenomenon.ABSORPTION, self)


@traced(logger)
class PolarizedCoatingLayer(SurfaceMaterial):
    """
    # TODO: Document
    """
    def __init__(self, *args):
        super(PolarizedCoatingLayer, self).__init__(*args)



@traced(logger)
class PolarizedCoatingReflectorLayer(PolarizedCoatingLayer):
    """
    # TODO: Document
    """
    def __init__(self, name, coating_file, sigma_1=None, sigma_2=None, k=None):
        # coating_material with four columns: wavelenth in nm,
        # angle in deg., reflectance s-polarized (perpendicular),
        # reflectance p-polarized (parallel)
        # the values in coating_material should be in the corresponding
        # order columns
        data_material = np.loadtxt(coating_file, usecols=(0, 1, 2, 3))
        plain_properties = {
            'Matrix_reflectance_coating': {
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
        properties = Material.plain_properties_to_properties(plain_properties)
        super(PolarizedCoatingReflectorLayer, self).__init__(name, properties)

    def change_of_optical_state(self, ray, normal_vector, nearby_material):
        polarization_vector = ray.current_polarization()
        incident = ray.current_direction()
        wavelength = ray.wavelength
        properties = self.properties
        normal = correct_normal(normal_vector, incident)
        c1 = - normal.dot(incident)
        inc_angle = rad_to_deg(arccos(c1))
        # incidence angle
        parallel_v, perpendicular_v, normal_parallel_plane = \
            parallel_orthogonal_components(polarization_vector, incident, normal)
        ref_per = perpendicular_v.Length ** 2.0 / polarization_vector.Length ** 2.0
        reflectance_matrix = properties['Matrix_reflectance_coating']
        r_matrix = reflectance_matrix(inc_angle, wavelength)
        # reflectance dependent of incidence angle and wavelength
        # We decide the polarization projection onto the parallel / perpendicular plane
        if myrandom() < ref_per:
            reflectance = calculate_reflectance(r_matrix, inc_angle, wavelength)[0]
            # reflectance for s-polarized (perpendicular) light
            perpendicular_polarized = True
            polarization_vector = normalize(perpendicular_v)
        else:
            reflectance = calculate_reflectance(r_matrix, inc_angle, wavelength)[1]
            # reflectance for p-polarized (parallel) light
            perpendicular_polarized = False
            polarization_vector = normalize(parallel_v)
        if myrandom() < reflectance:
            # ray reflected
            reflected = simple_reflexion(incident, normal).normalize()
            if not perpendicular_polarized:
                # reflexion changes the parallel component of incident polarization
                polarization_vector = simple_polarization_reflexion(
                    incident, normal, normal_parallel_plane, polarization_vector)
            state = OpticalState(polarization_vector, reflected, Phenomenon.REFLEXION, self)
            state.material = ray.current_medium()
            state.apply_dispersion(properties, normal_vector)
            return state
        else:
            # ray is killed in the coating reflector
            return OpticalState(Base.Vector(0.0, 0.0, 0.0),
                                Base.Vector(0.0, 0.0, 0.0),
                                Phenomenon.ABSORPTION)


@traced(logger)
class PolarizedCoatingTransparentLayer(PolarizedCoatingLayer):
    """
    # TODO: Document
    """
    def __init__(self, name, coating_file):
        # coatingmaterial calculated by TMM method, six columns:
        # wavelength in nm, angle in deg.,
        # reflectance s-polarized (perpendicular),
        # reflectance p-polarized (parallel),  transmittance s-polarized,
        # transmittance p-polarized
        # the values in coating_material should be in the corresponding
        # order columns
        data = np.loadtxt(coating_file)
        data_reflectance = data[:, [0, 1, 2, 3]]
        data_transmittance = data[:, [0, 1, 4, 5]]
        plain_properties = {
            'Matrix_reflectance_coating': {
                'type': 'matrix',
                'value': data_reflectance
            },
            'Matrix_transmittance_coating': {
                'type': 'matrix',
                'value': data_transmittance
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
        properties = Material.plain_properties_to_properties(plain_properties)
        super(PolarizedCoatingTransparentLayer, self).__init__(name, properties)

    def change_of_optical_state(self, ray, normal_vector, nearby_material):
        polarization_vector = ray.current_polarization()
        incident = ray.current_direction()
        wavelength = ray.wavelength
        properties = self.properties
        normal = correct_normal(normal_vector, incident)
        n1 = ray.current_medium().get_n(ray.wavelength)
        n2 = nearby_material.get_n(ray.wavelength)
        r = n1 / n2
        c1 = - normal.dot(incident)
        # cos (incident_angle)
        c2sq = 1.0 - r * r * (1.0 - c1 * c1)
        # cos (refracted_angle) ** 2
        if c2sq.real < 0:
            # total internal reflection
            state = reflexion(incident, normal, polarization_vector)
            state.material = ray.current_medium()
            return state
        c2 = sqrt(c2sq)
        # cos (refracted_angle)
        parallel_v, perpendicular_v, normal_parallel_plane = \
            parallel_orthogonal_components(polarization_vector, incident, normal)
        # parallel and perpendicular components of polarization vector
        # and orthogonal vector of the parallel plane
        ref_per = perpendicular_v.Length ** 2.0 / polarization_vector.Length ** 2.0
        # weight of perpendicular component: 0 < ref_per < 1
        inc_angle = rad_to_deg(arccos(c1))
        reflectance_matrix = properties['Matrix_reflectance_coating']
        r_matrix = reflectance_matrix(inc_angle, wavelength)
        # reflectance dependent of incidence angle and wavelength
        # We decide the polarization projection onto the parallel / perpendicular plane
        if myrandom() < ref_per:
            reflectance = calculate_reflectance(r_matrix, inc_angle, wavelength)[0]
            # reflectance for s-polarized (perpendicular) light
            perpendicular_polarized = True
            polarization_vector = normalize(perpendicular_v)
        else:
            reflectance = calculate_reflectance(r_matrix, inc_angle, wavelength)[1]
            # reflectance for p-polarized (parallel) light
            perpendicular_polarized = False
            polarization_vector = normalize(parallel_v)
        if myrandom() < reflectance:
            # ray reflected
            reflected_direction = simple_reflexion(incident, normal)
            if not perpendicular_polarized:
                # reflexion changes the parallel component of incident polarization
                polarization_vector = simple_polarization_reflexion(
                    incident, normal, normal_parallel_plane, polarization_vector)
            state = OpticalState(polarization_vector, reflected_direction, Phenomenon.REFLEXION)
            state.material = ray.current_medium()
            return state
        else:
            # ray refracted: computing the refracted direction and energy absorbed in coating
            transmittance_matrix = properties['Matrix_transmittance_coating']
            t_matrix = transmittance_matrix(inc_angle, wavelength)
            # transmittance dependent of incidence angle and wavelength
            if perpendicular_polarized:
                transmittance = calculate_reflectance(t_matrix, inc_angle, wavelength)[0]
            else:
                transmittance = calculate_reflectance(t_matrix, inc_angle, wavelength)[1]
            factor_energy_absorbed = (1 - reflectance - transmittance) / (1 - reflectance)
            refracted_direction = incident * r.real + \
                                  normal * (r.real * c1 - c2.real)
            refracted_direction.normalize()
            if not perpendicular_polarized:
                # refraction changes the parallel component of incident polarization
                polarization_vector = \
                    simple_polarization_refraction(
                        incident, normal, normal_parallel_plane, c2, polarization_vector)
            optical_state = OpticalState(polarization_vector, refracted_direction, Phenomenon.REFRACTION, nearby_material)
            optical_state.extra_data['factor_energy_absorbed'] = \
                factor_energy_absorbed
            return optical_state

@traced(logger)
class PolarizedCoatingAbsorberLayer(PolarizedCoatingLayer):
    """
    # TODO: Document
    """
    def __init__(self, name, coating_file):
        # coating_material with four columns: wavelenth in nm, angle in deg.,
        # reflectance s-polarized (perpendicular),
        # reflectance p-polarized (parallel)
        # the values in coating_material should be in the corresponding order
        # columns
        data_material = np.loadtxt(coating_file, usecols=(0, 1, 2, 3))
        plain_properties = {
            'Matrix_reflectance_coating': {
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
                'value': True
            },
            'transparent_material': {
                'type': 'scalar',
                'value': False
            },
            'lambertian_material': {
                'type': 'scalar',
                'value': False
            },
        }
        properties = Material.plain_properties_to_properties(plain_properties)
        super(PolarizedCoatingAbsorberLayer, self).__init__(name, properties)

    def change_of_optical_state(self, ray, normal_vector, nearby_material):
        polarization_vector = ray.current_polarization()
        incident = ray.current_direction()
        wavelength = ray.wavelength
        properties = self.properties
        normal = correct_normal(normal_vector, incident)
        c1 = - normal.dot(incident)
        inc_angle = rad_to_deg(arccos(c1))
        # incidence angle
        parallel_v, perpendicular_v, normal_parallel_plane = \
            parallel_orthogonal_components(polarization_vector, incident, normal)
        ref_per = perpendicular_v.Length ** 2.0 / polarization_vector.Length ** 2.0
        reflectance_matrix = properties['Matrix_reflectance_coating']
        r_matrix = reflectance_matrix(inc_angle, wavelength)
        # reflectance dependent of incidence angle and wavelength
        # We decide the polarization projection onto the parallel / perpendicular plane
        if myrandom() < ref_per:
            reflectance = calculate_reflectance(r_matrix, inc_angle, wavelength)[0]
            # reflectance for s-polarized (perpendicular) light
            perpendicular_polarized = True
            polarization_vector = normalize(perpendicular_v)
        else:
            reflectance = calculate_reflectance(r_matrix, inc_angle, wavelength)[1]
            # reflectance for p-polarized (parallel) light
            perpendicular_polarized = False
            polarization_vector = normalize(parallel_v)
        if myrandom() < reflectance:
            # ray reflected
            reflected = simple_reflexion(incident, normal).normalize()
            if not perpendicular_polarized:
                # reflexion changes the parallel component of incident polarization
                polarization_vector = simple_polarization_reflexion(
                    incident, normal, normal_parallel_plane, polarization_vector)
            state = OpticalState(polarization_vector, reflected, Phenomenon.REFLEXION, self)
            state.material = ray.current_medium()
            return state
        else:
            # ray refracted: ray energy will be absorbed in the coating absorber
            state = OpticalState(Base.Vector(0.0, 0.0, 0.0),
                                 Base.Vector(0.0, 0.0, 0.0),
                                 Phenomenon.ENERGY_ABSORBED, self)
            return state


@traced(logger)
class TwoLayerMaterial(Material):
    """
    # TODO: Document
    """
    def __init__(self, name, name_front_layer, name_back_layer):
        super(TwoLayerMaterial, self).__init__(name, {})
        self.name_front_layer = name_front_layer
        self.name_back_layer = name_back_layer
        self.front_material = Material.by_name[name_front_layer]
        self.back_material = Material.by_name[name_back_layer]

    def to_json(self):
        return json.dumps(
            {
                'name': self.name,
                'classname': 'TwoLayerMaterial',
                'name_front_layer': self.name_front_layer,
                'name_back_layer': self.name_back_layer
            }, cls=NumpyEncoder
        )

    def change_of_optical_state(self, ray, normal_vector, nearby_material):
        if ray.current_direction().dot(normal_vector) < 0:
            # Ray intercepted on the frontside of the surface
            material = self.front_material
        else:
            # Ray intercepted on the backside of the surface
            material = self.back_material
        return material.change_of_optical_state(ray, normal_vector, nearby_material)
