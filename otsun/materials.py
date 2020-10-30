"""Module otsun.materials for treating materials

The module relies on a basic class `Material` with two subclasses
`VolumeMaterial` and `SurfaceMaterial`, and several subclasses of them
for specific materials.
"""

import json
import zipfile
from FreeCAD import Base
from .optics import Phenomenon, OpticalState, reflection, refraction, matrix_reflectance, \
    calculate_reflectance, simple_polarization_reflection, simple_polarization_refraction, \
    simple_reflection, shure_refraction, lambertian_reflection
from .math import arccos, parallel_orthogonal_components, rad_to_deg, myrandom, normalize, \
    constant_function, correct_normal, tabulated_function
from numpy import sqrt
import numpy as np
from autologging import traced
from .logging_unit import logger


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
        XXX(MatYYY,OTHER_INFO), returns the material whose name is MatYYY

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
        string = label[start + 1:end]
        name = string.split(',')[0]
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
    Subclass of Material for materials with volume

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
                                ray.current_direction(), Phenomenon.REFRACTION, self)  # TODO: Set solid
        else:
            n1 = ray.current_medium().get_n(wavelength)
            n2 = self.get_n(wavelength)
            optical_state = refraction(ray.current_direction(), normal_vector, n1, n2, ray.current_polarization())
            if optical_state.phenomenon == Phenomenon.REFRACTION:
                optical_state.material = self  # TODO: Set solid
            else:
                optical_state.material = ray.current_medium()  # TODO: Set solid
            return optical_state

    def to_json(self):
        """
        Dumps the material in json format.
        """
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
    Subclass of `VolumeMaterial` for those materials with constant properties.
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
    Subclass of `VolumeMaterial` for materials with tabulated index of refraction.
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
    Subclass of `VolumeMaterial` for photovoltaic materials.
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
        Computes the photovoltaic data stored in a ray.

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
                 energy_before, ray.energy, ray.wavelength, alpha, angle_incident, self.name)
                )


class PolarizedThinFilm(VolumeMaterial):
    """
    Subclass of `VolumeMaterial` for polarized thin film materials.
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
        Helper function for the computation of the optical state once the ray has passed through the film
        """
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
            state = reflection(incident, normal, polarization_vector)
            return 0.0, state
            # no energy is abosrbed in the thinfilm
        c2 = sqrt(c2sq)
        # cos (refracted_angle)
        if c2.real > 1:
            # avoiding invalid solutions
            c2 = 1        
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
            reflected_direction = simple_reflection(incident, normal).normalize()
            if not perpendicular_polarized:
                # reflection changes the parallel component of incident polarization
                polarization_vector = simple_polarization_reflection(
                    incident, normal, normal_parallel_plane, polarization_vector)
            return 0.0, OpticalState(polarization_vector, reflected_direction,
                                     Phenomenon.REFLEXION)  # TODO: Set solid
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
            refracted_direction = incident * r.real + normal * (r.real * c1 - c2.real)
            refracted_direction.normalize()
            if not perpendicular_polarized:
                # refraction changes the parallel component of incident polarization
                polarization_vector = \
                    simple_polarization_refraction(
                        incident, normal, normal_parallel_plane, c2, polarization_vector)
            return (factor_energy_absorbed,
                    OpticalState(polarization_vector, refracted_direction,
                                 Phenomenon.REFRACTION))  # TODO: Set solid

    def change_of_optical_state(self, ray, normal_vector):
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
            optical_state.material = self  # TODO: Set solid
        else:
            optical_state.material = ray.current_medium()  # TODO: Set solid
        return optical_state


vacuum_medium = SimpleVolumeMaterial("Vacuum", 1.0, 0.0)

@traced(logger)
class SurfaceMaterial(Material):
    """
    Subclass of Material for 2-dimensional materials (without volume)

    The `properties` parameter must be a dict with the physical properties
    describing the material. At least, the following must be provided:
    'probability_of_reflection': probability that a photon gets reflected,
    as a function of its wavelength.
    'probability_of_absorption': probability that a photon gets absorbed,
    as a function of its wavelength.
    'probability_of_transmittance': probability that a photon passes through
    the material, as a function of its wavelength.
    """

    def __init__(self, name, properties):
        super(SurfaceMaterial, self).__init__(name, properties)
        self.properties = properties

    @classmethod
    def create(cls, name, properties):
        _ = cls(name, properties)

    def compute_probabilities(self, ray):
        """
        Computes the tuple of probabilities that a ray hitting the surface gets reflected, absorbed or transmitted
        """
        properties = self.properties
        try:
            por = properties['probability_of_reflection'](ray.wavelength)
        except KeyError:
            por = 1.0
        try:
            poa = properties['probability_of_absorption'](ray.wavelength)
        except KeyError:
            poa = 1 - por
        try:
            pot = properties['probability_of_transmittance'](ray.wavelength)
        except KeyError:
            pot = 0.0

        return [por, poa, pot]

    def decide_phenomenon(self, ray):
        """
        Decides which phenomenon will take place when a ray hits the surface.
        """
        phenomena = [
            Phenomenon.REFLEXION,
            Phenomenon.ABSORPTION,
            Phenomenon.TRANSMITTANCE]

        probabilities = self.compute_probabilities(ray)
        phenomenon = np.random.choice(phenomena, 1, p=probabilities)[0]
        return phenomenon

    def change_of_optical_state(self, ray, normal_vector, nearby_material):
        phenomenon = self.decide_phenomenon(ray)
        properties = self.properties
        if phenomenon == Phenomenon.REFLEXION:
            polarization_vector = ray.current_polarization()
            incident = ray.current_direction()
            if self.properties.get('lambertian_material', False):
                state = lambertian_reflection(ray.current_direction(), normal_vector)
            else:
                state = reflection(incident, normal_vector, polarization_vector, False)
            state.material = ray.current_medium()  # TODO: Set solid
            state.apply_dispersion(properties, normal_vector)
            return state
        if phenomenon == Phenomenon.ABSORPTION:
            if self.properties.get('thermal_material', False):
                return (OpticalState(Base.Vector(0.0, 0.0, 0.0),
                                     Base.Vector(0.0, 0.0, 0.0),
                                     Phenomenon.ENERGY_ABSORBED,
                                     self))  # TODO: Set solid
            else:
                return (OpticalState(Base.Vector(0.0, 0.0, 0.0),
                                     Base.Vector(0.0, 0.0, 0.0),
                                     Phenomenon.ABSORPTION,
                                     self))  # TODO: Set solid
        if phenomenon == Phenomenon.TRANSMITTANCE:
            # refraction in transparent layer
            n1 = ray.current_medium().get_n(ray.wavelength)
            n2 = nearby_material.get_n(ray.wavelength)
            if n1 == n2:  # transparent_simple_layer
                state = OpticalState(ray.current_polarization(),
                                     ray.current_direction(),
                                     Phenomenon.REFRACTION,
                                     nearby_material)  # TODO: Set solid
            else:
                state = shure_refraction(ray.current_direction(), normal_vector,
                                         n1, n2, ray.current_polarization(),
                                         self.properties.get('lambertian_material', False))
                state.material = nearby_material  # TODO: Set solid
            return state

    @classmethod
    def from_plain_properties(cls, name, plain_properties):
        """
        Loads a material from its properties stored in a simple dictionary

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
    Subclass of `SurfaceMaterial` for completely opaque layers.
    """

    def __init__(self, name):
        plain_properties = {
            'probability_of_reflection': {
                'type': 'constant',
                'value': 0.0
            },
            'probability_of_absorption': {
                'type': 'constant',
                'value': 1.0
            },
            'probability_of_transmittance': {
                'type': 'constant',
                'value': 0.0
            },
            'thermal_material': {
                'type': 'scalar',
                'value': False
            }
        }
        properties = Material.plain_properties_to_properties(plain_properties)
        super(OpaqueSimpleLayer, self).__init__(name, properties)


@traced(logger)
class TransparentSimpleLayer(SurfaceMaterial):
    """
    Subclass of `SurfaceMaterial` for transparent layers.
    """

    def __init__(self, name, pot):
        plain_properties = {
            'probability_of_reflection': {
                'type': 'constant',
                'value': 1.0 - pot
            },
            'probability_of_absorption': {
                'type': 'constant',
                'value': 0.0
            },
            'probability_of_transmittance': {
                'type': 'constant',
                'value': pot
            },
            'thermal_material': {
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


@traced(logger)
class AbsorberSimpleLayer(SurfaceMaterial):
    """
    Subclass of `SurfaceMaterial` for absorber layers with behaviour independent of the wavelength.
    """

    def __init__(self, name, poa):
        plain_properties = {
            'probability_of_reflection': {
                'type': 'constant',
                'value': 1.0 - poa
            },
            'probability_of_absorption': {
                'type': 'constant',
                'value': poa
            },
            'probability_of_transmittance': {
                'type': 'constant',
                'value': 0.0
            },
            'thermal_material': {
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


@traced(logger)
class AbsorberLambertianLayer(SurfaceMaterial):
    """
    Subclass of `SurfaceMaterial` for absorber layers with lambertian behaviour when reflecting rays.
    """

    def __init__(self, name, poa):
        plain_properties = {
            'probability_of_reflection': {
                'type': 'constant',
                'value': 1.0 - poa
            },
            'probability_of_absorption': {
                'type': 'constant',
                'value': poa
            },
            'probability_of_transmittance': {
                'type': 'constant',
                'value': 0.0
            },
            'thermal_material': {
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


@traced(logger)
class ReflectorSpecularLayer(SurfaceMaterial):
    """
    Subclass of `SurfaceMaterial` for reflector specular layers.
    """

    def __init__(self, name, por, sigma_1=None, sigma_2=None, k=None):
        plain_properties = {
            'probability_of_reflection': {
                'type': 'constant',
                'value': por
            },
            'probability_of_absorption': {
                'type': 'constant',
                'value': 1.0 - por
            },
            'probability_of_transmittance': {
                'type': 'constant',
                'value': 0.0
            },
            'thermal_material': {
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


@traced(logger)
class ReflectorLambertianLayer(SurfaceMaterial):
    """
    Subclass of `SurfaceMaterial` for reflector layers with lambertian behaviour.
    """

    def __init__(self, name, por):
        plain_properties = {
            'probability_of_reflection': {
                'type': 'constant',
                'value': por
            },
            'probability_of_absorption': {
                'type': 'constant',
                'value': 1.0 - por
            },
            'probability_of_transmittance': {
                'type': 'constant',
                'value': 0.0
            },
            'thermal_material': {
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


@traced(logger)
class AbsorberTWModelLayer(SurfaceMaterial):
    """
    Subclass of `SurfaceMaterial` for materials using Tesfamichael-Wackelgard model.
    """

    def __init__(self, name, poa, b_constant, c_constant):
        plain_properties = {
            'probability_of_reflection': {
                'type': 'constant',
                'value': 1.0 - poa
            },
            'probability_of_absorption': {
                'type': 'constant',
                'value': poa
            },
            'probability_of_transmittance': {
                'type': 'constant',
                'value': 0.0
            },
            'thermal_material': {
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
        normal_vector : Base.Vector
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
        absorption_ratio = self.tw_absorptance_ratio(
            normal_vector, b_constant, c_constant, ray.current_direction())
        absorption = properties['probability_of_absorption'](ray.wavelength) * absorption_ratio
        reflectance = 1.0 - absorption
        if myrandom() < reflectance:
            polarization_vector = ray.current_polarization()
            state = reflection(ray.current_direction(), normal_vector, polarization_vector, False)
            state.material = ray.current_medium()  # TODO: Set solid
            return state
        else:
            # absorption in absorber material: the ray is killed
            return OpticalState(Base.Vector(0.0, 0.0, 0.0),
                                Base.Vector(0.0, 0.0, 0.0),
                                Phenomenon.ENERGY_ABSORBED, self)  # TODO: Set solid


@traced(logger)
class MetallicLayer(SurfaceMaterial):
    """
    Subclass of `SurfaceMaterial` for metallic layers.
    """


@traced(logger)
class MetallicSpecularLayer(MetallicLayer):
    """
    Subclass of `MetallicLayer` for metallic layers with specular properties.
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
            'index_of_refraction': {  # TODO: @Ramon Crec que aixo no s'empra per res
                'type': 'tabulated',
                'value': [wavelength_values, n_values]
            },
            'extinction_coefficient': {
                'type': 'tabulated',
                'value': [wavelength_values, k_values]
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
            state.material = ray.current_medium()  # TODO: Set solid
            state.apply_dispersion(properties, normal_vector)
            return state
        if state.phenomenon == Phenomenon.REFRACTION:
            # refraction in metallic layer: the ray is killed
            return OpticalState(Base.Vector(0.0, 0.0, 0.0),
                                Base.Vector(0.0, 0.0, 0.0),
                                Phenomenon.ABSORPTION, self)  # TODO: Set solid


@traced(logger)
class MetallicLambertianLayer(MetallicLayer):
    """
    Subclass of `MetallicLayer` for metallic layers with lambertian behaviour.
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
            state.material = ray.current_medium()  # TODO: Set solid
            return state
        if state.phenomenon == Phenomenon.REFRACTION:
            # refraction in metallic layer: the ray is killed
            return OpticalState(Base.Vector(0.0, 0.0, 0.0),
                                Base.Vector(0.0, 0.0, 0.0),
                                Phenomenon.ABSORPTION, self)  # TODO: Set solid


@traced(logger)
class PolarizedCoatingLayer(SurfaceMaterial):
    """
    Subclass of `SurfaceMaterial` for polarized coating layers
    """

    def __init__(self, *args):
        super(PolarizedCoatingLayer, self).__init__(*args)

    def precompute_change_of_optical_state(self, ray, normal_vector):
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
            reflected = simple_reflection(incident, normal).normalize()
            if not perpendicular_polarized:
                # reflection changes the parallel component of incident polarization
                polarization_vector = simple_polarization_reflection(
                    incident, normal, normal_parallel_plane, polarization_vector)
            state = OpticalState(polarization_vector, reflected, Phenomenon.REFLEXION, self)
            state.material = ray.current_medium()  # TODO: Set solid
            state.apply_dispersion(properties, normal_vector)
            return state
        else:
            return (inc_angle, incident, perpendicular_polarized,
                    reflectance, normal_parallel_plane)


@traced(logger)
class PolarizedCoatingReflectorLayer(PolarizedCoatingLayer):
    """
    Subclass of `PolarizedCoatingLayer` for reflector layers.
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
        new_state = self.precompute_change_of_optical_state(ray, normal_vector)
        if isinstance(new_state, OpticalState):
            return new_state
        else:
            # ray is killed in the coating reflector
            return OpticalState(Base.Vector(0.0, 0.0, 0.0),
                                Base.Vector(0.0, 0.0, 0.0),
                                Phenomenon.ABSORPTION)  # TODO: Set solid


@traced(logger)
class PolarizedCoatingAbsorberLayer(PolarizedCoatingLayer):
    """
    Subclass of `PolarizedCoatingLayer` for absorber layers.
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
            'probability_of_transmittance': {
                'type': 'constant',
                'value': 0
            },
        }
        properties = Material.plain_properties_to_properties(plain_properties)
        super(PolarizedCoatingAbsorberLayer, self).__init__(name, properties)

    def change_of_optical_state(self, ray, normal_vector, nearby_material):
        new_state = self.precompute_change_of_optical_state(ray, normal_vector)
        if isinstance(new_state, OpticalState):
            return new_state
        else:
            # ray refracted: ray energy will be absorbed in the coating absorber
            state = OpticalState(Base.Vector(0.0, 0.0, 0.0),
                                 Base.Vector(0.0, 0.0, 0.0),
                                 Phenomenon.ENERGY_ABSORBED, self)  # TODO: Set solid
            return state


@traced(logger)
class PolarizedCoatingTransparentLayer(PolarizedCoatingLayer):
    """
    Subclass of `PolarizedCoatingLayer` for transparent layers.
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
            state = reflection(incident, normal, polarization_vector)
            state.material = ray.current_medium()  # TODO: Set solid
            return state
        c2 = sqrt(c2sq)
        if c2.real > 1:
            # avoiding invalid solutions
            c2 = 1
        new_state = self.precompute_change_of_optical_state(ray, normal_vector)
        if isinstance(new_state, OpticalState):
            return new_state
        else:
            (inc_angle, incident, perpendicular_polarized,
             reflectance, normal_parallel_plane) = new_state
            # ray refracted: computing the refracted direction and energy absorbed in coating
            transmittance_matrix = properties['Matrix_transmittance_coating']
            t_matrix = transmittance_matrix(inc_angle, wavelength)
            # transmittance dependent of incidence angle and wavelength
            if perpendicular_polarized:
                transmittance = calculate_reflectance(t_matrix, inc_angle, wavelength)[0]
            else:
                transmittance = calculate_reflectance(t_matrix, inc_angle, wavelength)[1]
            factor_energy_absorbed = (1 - reflectance - transmittance) / (1 - reflectance)
            refracted_direction = incident * r.real + normal * (r.real * c1 - c2.real)
            refracted_direction.normalize()
            if not perpendicular_polarized:
                # refraction changes the parallel component of incident polarization
                polarization_vector = \
                    simple_polarization_refraction(
                        incident, normal, normal_parallel_plane, c2, polarization_vector)
            optical_state = OpticalState(polarization_vector, refracted_direction,
                                         Phenomenon.REFRACTION, nearby_material)  # TODO: Set solid
            optical_state.extra_data['factor_energy_absorbed'] = \
                factor_energy_absorbed
            return optical_state


@traced(logger)
class TwoLayerMaterial(Material):
    """
    Subclass of `Material` for surface materials formed by two layers (back and front)
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
