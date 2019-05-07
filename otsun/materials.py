"""Module otsun.materials for treating materials

The module relies on a basic class `Material` with two subclasses
`VolumeMaterial` and `SurfaceMaterial`, and several subclasses of them
for specific materials.
"""


import json
import zipfile
import dill
from .optics import *
from .math import *
#from .source import Ray
#from .optics import OpticalState


class NumpyEncoder(json.JSONEncoder):
    """Wrapper to dump numpy arrays as json"""
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return json.JSONEncoder.default(self, obj)


# ---
# Classes for materials
# ---


@traced(logger)
class Material(object):
    """
    Class used to represent materials and their physical properties

    The `properties` dictionary holds the physical properties of the material.
    Its contents are specific to the kind of material.

    Class Variables
    ---------------
    by_name : dict
        associates each material name to the material itself

    Attributes
    ----------
    name : str
        Holds the name of the material
    kind : str
        Holds the type of the material ("Surface" or "Volume")
    properties : dict
        Dictionary with physical properties of the material
    """
    by_name = {}

    def __init__(self, name, properties=None):
        self.by_name[name] = self
        self.name = name
        self.kind = ""
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
    def load_from_json_fileobject(cls, f):
        """
        Load materials from a json fileobject

        If the file contains a single dict, then it means that it contains a
        single material. Otherwise it contains an array, where each entry is a dict
        representing a material.

        Parameters
        ----------
        f : file
            File object

        Returns
        -------
        str
            String with the name of the last material imported from the file
        """
        info = json.load(f)
        if type(info).__name__ == 'dict':
            info = [info]
        name = ""
        for mat_spec in info:
            classname = mat_spec['classname']
            logger.debug(classname)
            kind = mat_spec['kind']
            name = mat_spec['name']
            if kind == 'Volume':
                mat = VolumeMaterial(name, {})
                plain_properties = mat_spec['plain_properties']
                properties = Material.plain_properties_to_properties(plain_properties)
                mat.properties = properties
                the_class = globals()[classname]
                mat.__class__ = the_class
            if kind == 'Surface':
                if classname == "TwoLayerMaterial":
                    name_front_layer = mat_spec['name_front_layer']
                    name_back_layer = mat_spec['name_back_layer']
                    mat = TwoLayerMaterial(name, name_front_layer, name_back_layer)
                else:
                    mat = SurfaceMaterial(name, {})
                    plain_properties = mat_spec['plain_properties']
                    properties = Material.plain_properties_to_properties(plain_properties)
                    mat.properties = properties
                    the_class = globals()[classname]
                    mat.__class__ = the_class
        return name

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

        Returns
        -------
        """
        try:
            with zipfile.ZipFile(filename) as z:
                for matfile in z.namelist():
                    with z.open(matfile) as f:
                        cls.load_from_json_fileobject(f)
        except IOError:
            logger.exception("error in processing file %s", filename)

    @classmethod
    def load_from_file(cls, filename):
        """
        Load materials from binary format

        ..deprecated::
            The binary format will be deprecated soon... move to json

        Parameters
        ----------
        filename : str
            Name of the file

        Returns
        -------
        str
            Name of the imported material
        """
        try:
            with open(filename, 'rb') as f:
                mat = dill.load(f)
                cls.by_name[mat.name] = mat
                return mat.name
        except IOError:
            logger.exception("error in processing file %s", filename)

    @classmethod
    def load_from_zipfile(cls, filename):
        """
        Load materials from zipfile of files in binary format

        ..deprecated::
            The binary format will be deprecated soon... move to json

        Parameters
        ----------
        filename : str
            Name of the zip file

        Returns
        -------
        """
        with zipfile.ZipFile(filename) as z:
            for matfile in z.namelist():
                with z.open(matfile) as f:
                    try:
                        mat = dill.load(f)
                        cls.by_name[mat.name] = mat
                    except IOError:
                        logger.exception("error in processing file %s", matfile)

    def change_of_direction(self, ray, normal_vector, *args):
        """
        Computes how a ray behaves when interacting with the material.
        MUST be subclassed

        The material where the ray is actually located is held in
        ray.current_material.

        Parameters
        ----------
        ray : Ray
        normal_vector : Vector
        args
            Variable length argument list

        Returns
        -------
        OpticalState
        """
        pass

    def to_json(self):
        """Converts material to json. MUST be subclassed"""
        pass

    def save_to_json_file(self, filename):
        """
        Save material to json file

        Parameters
        ----------
        filename : str
            Name of the file

        Returns
        -------
        """
        with open(filename, 'w') as f:
            f.write(self.to_json())

    @staticmethod
    def all_to_json():
        materials = Material.by_name.values()
        simple_mats = [material for material in materials if
                       not isinstance(material, TwoLayerMaterial)]
        composite_mats = [material for material in materials if
                          isinstance(material, TwoLayerMaterial)]
        all_mats = simple_mats + composite_mats
        return [mat.to_json() for mat in all_mats]

    @staticmethod
    def save_all_to_json_file(filename):
        with open(filename, 'w') as f:
            f.write('[\n')
            f.write(',\n'.join(Material.all_to_json()))
            f.write('\n]')

@traced(logger)
class VolumeMaterial(Material):
    """
    Subclass of Volume for materials with volume

    TODO: comment on properties

    """
    def __init__(self, name, properties={}):
        """
        Initializes a Volume Material.

        The `properties` parameter must be a dict with the physical properties
        describing the material. At least, the following must be provided:
        'index_of_refraction': index of refraction of the material, as a
        function of its wavelength, only real part.
        'extinction_coefficient': imaginary part of the index of refraction
        of the material as a function of its wavelength.

        Parameters
        ----------
        name : str
            Name of the material
        properties : dict
            Properties of the material
        """
        super(VolumeMaterial, self).__init__(name, properties)
        self.kind = 'Volume'

    def change_of_direction(self, ray, normal_vector):
        wavelength = ray.wavelength
        if 'Matrix_reflectance_thin_film' in ray.current_medium.properties:
            # the ray is traveling inside thin film material and
            # then it will certainly refract
            ray.current_medium = self
            return OpticalState(ray.polarization_vectors[-1],
                                ray.directions[-1], Phenomenon.REFRACTION)
        if 'Matrix_reflectance_thin_film' in self.properties:
            # the ray impacts on thin film material
            n1 = ray.materials[-1].properties['index_of_refraction'](
                ray.wavelength)
            if 'extinction_coefficient' in ray.current_medium.properties:
                n1 += 1j * \
                      ray.materials[-1].properties['extinction_coefficient'](
                          ray.wavelength)
            n_front = self.properties['index_of_refraction_front'](ray.wavelength)
            n_back = self.properties['index_of_refraction_back'](ray.wavelength)
            if n_front == ray.materials[-1].properties['index_of_refraction'](
                    ray.wavelength):
                # impacts on front material
                n2 = n_back
            else:
                n2 = n_front
            energy_absorbed_thin_film, optical_state = (
                calculate_state_thin_film(
                    ray.directions[-1], normal_vector, n1,
                    n2,
                    ray.polarization_vectors[-1],
                    self.properties, ray.wavelength))
            ray.energy = ray.energy * (1.0 - energy_absorbed_thin_film)
            if optical_state.phenomenon == Phenomenon.REFRACTION:
                ray.current_medium = self
            return optical_state
        else:
            n1 = ray.current_medium.properties['index_of_refraction'](
                wavelength)
            if 'extinction_coefficient' in ray.current_medium.properties:
                n1 += 1j * \
                      ray.current_medium.properties['extinction_coefficient'](
                          wavelength)
            n2 = self.properties['index_of_refraction'](wavelength)
            if 'extinction_coefficient' in self.properties:
                n2 += 1j * self.properties['extinction_coefficient'](wavelength)
            optical_state = refraction(
                ray.directions[-1], normal_vector, n1, n2,
                ray.polarization_vectors[-1])
            if optical_state.phenomenon == Phenomenon.REFRACTION:
                ray.current_medium = self
            return optical_state

    def to_json(self):
        return json.dumps(
            {
                'name': self.name,
                'kind': self.kind,
                'classname': self.__class__.__name__,
                'plain_properties': self.properties.get(
                    'plain_properties', None)
            }, cls=NumpyEncoder, indent=4
        )


class SimpleVolumeMaterial(VolumeMaterial):
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
        super(SimpleVolumeMaterial,self).__init__(name, {})
        self.properties = Material.plain_properties_to_properties(plain_properties)


class WavelengthVolumeMaterial(VolumeMaterial):
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
        super(WavelengthVolumeMaterial,self).__init__(name)
        self.properties = Material.plain_properties_to_properties(plain_properties)


class PVMaterial(VolumeMaterial):
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
        super(PVMaterial,self).__init__(name)
        self.properties = Material.plain_properties_to_properties(plain_properties)


class PolarizedThinFilm(VolumeMaterial):
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
        super(PolarizedThinFilm,self).__init__(name)
        self.properties = Material.plain_properties_to_properties(plain_properties)


vacuum_medium = SimpleVolumeMaterial("Vacuum", 1.0, 0.0)


@traced(logger)
class SurfaceMaterial(Material):
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
        self.kind = 'Surface'

    @classmethod
    def create(cls, name, properties):
        _ = cls(name, properties)

    def decide_phenomenon(self, ray, normal_vector, nearby_material):
        # TODO : Humanize
        properties = self.properties
        phenomena = [
            Phenomenon.REFLEXION,
            Phenomenon.ABSORPTION,
            Phenomenon.TRANSMITTANCE]
        polarization_vector = ray.polarization_vectors[-1]
        perpendicular_polarized = False
        probabilities = None
        if 'TW_model' in properties:
            b_constant = properties['b_constant']
            c_constant = properties['c_constant']
            absortion_ratio = tw_absorptance_ratio(
                normal_vector, b_constant, c_constant, ray.directions[-1])
            absortion = properties['probability_of_absortion'](
                ray.properties['wavelength']) * absortion_ratio
            por = 1.0 - absortion
            probabilities = [por, absortion, 0]  # Here I assume no transmitance
        if 'Matrix_polarized_reflectance_coating' in properties:
            # polarized_coating_layer
            n1 = ray.current_medium.properties['index_of_refraction'](
                ray.wavelength)
            if 'extinction_coefficient' in ray.current_medium.properties:
                n1 = ray.current_medium.properties['index_of_refraction'](
                    ray.wavelength) + 1j * \
                     ray.current_medium.properties['extinction_coefficient'](
                         ray.wavelength)
            n2 = nearby_material.properties['index_of_refraction'](
                ray.wavelength)
            if 'extinction_coefficient' in nearby_material.properties:
                n2 = nearby_material.properties['index_of_refraction'](
                    ray.wavelength) + 1j * \
                     nearby_material.properties['extinction_coefficient'](
                         ray.wavelength)
            results = calculate_probabilities_polarizaton_coating(
                ray.directions[-1], normal_vector, n1, n2,
                ray.polarization_vectors[-1],
                properties, ray.wavelength)
            probabilities = [results[0], results[1], results[2]]
            polarization_vector = results[3]
            perpendicular_polarized = results[4]
        if 'metallic_material' in properties:
            n1 = ray.current_medium.properties['index_of_refraction'](
                ray.wavelength)
            if 'extinction_coefficient' in ray.current_medium.properties:
                n1 = ray.current_medium.properties['index_of_refraction'](
                    ray.wavelength) + 1j * \
                     ray.current_medium.properties['extinction_coefficient'](
                         ray.wavelength)
            n2 = properties['index_of_refraction'](ray.wavelength)
            if 'extinction_coefficient' in properties:
                n2 = properties['index_of_refraction'](ray.wavelength) + \
                     1j * properties['extinction_coefficient'](ray.wavelength)
            results = calculate_reflexion_metallic(
                ray.directions[-1], normal_vector, n1, n2, polarization_vector)
            probabilities = [results[0], results[1], results[2]]
            polarization_vector = results[3]
            perpendicular_polarized = results[4]

        if probabilities is None:
            try:
                por = properties['probability_of_reflexion'](
                    ray.properties['wavelength'])
            except KeyError:
                por = 1.0
            try:
                poa = properties['probability_of_absortion'](
                    ray.properties['wavelength'])
            except KeyError:
                poa = 1 - por
            try:
                pot = properties['probability_of_transmitance'](
                    ray.properties['wavelength'])
            except KeyError:
                pot = 0.0

            probabilities = [por, poa, pot]
        phenomenon = np.random.choice(phenomena, 1, p=probabilities)[0]
        return phenomenon, polarization_vector, perpendicular_polarized

    def change_of_direction_by_absortion(self, ray, normal_vector):
        properties = self.properties
        if properties['energy_collector']:
            return (OpticalState(Base.Vector(0.0, 0.0, 0.0),
                                 Base.Vector(0.0, 0.0, 0.0),
                                 Phenomenon.GOT_ABSORBED))
        else:
            return (OpticalState(Base.Vector(0.0, 0.0, 0.0),
                                 Base.Vector(0.0, 0.0, 0.0),
                                 Phenomenon.ABSORPTION))

    def change_of_direction_by_reflexion(self, ray, normal_vector,
                                         polarization_vector_calculated_before):
        properties = self.properties
        if 'specular_material' in properties:
            state = reflexion(ray.directions[-1], normal_vector,
                              ray.polarization_vectors[-1],
                              polarization_vector_calculated_before)
            if properties.get('sigma_1',None):
                sigma_1 = properties['sigma_1']
                if properties.get('sigma_2',None):
                    sigma_2 = properties['sigma_2']
                    k = properties.get('k', None) or 0.5
                    return double_gaussian_dispersion(
                        normal_vector, state, sigma_1, sigma_2, k)
                return single_gaussian_dispersion(normal_vector, state, sigma_1)
            return state
        if 'lambertian_material' in properties:
            state = lambertian_reflexion(ray.directions[-1], normal_vector)
            return state

    def change_of_direction_by_transmitance(self, ray, normal_vector,
                                            nearby_material,
                                            perpendicular_polarized):
        cur_med_prop = ray.current_medium.properties
        n1 = cur_med_prop['index_of_refraction'](ray.wavelength)
        if 'extinction_coefficient' in cur_med_prop:
            n1 = cur_med_prop['index_of_refraction'](ray.wavelength) + \
                 1j * cur_med_prop['extinction_coefficient'](ray.wavelength)
        nea_med_prop = nearby_material.properties
        n2 = nea_med_prop['index_of_refraction'](ray.wavelength)
        if 'extinction_coefficient' in nea_med_prop:
            n2 = nea_med_prop['index_of_refraction'](ray.wavelength) +\
                 1j * nea_med_prop['extinction_coefficient'](ray.wavelength)
        if n1 == n2:  # transparent_simple_layer
            state = OpticalState(ray.polarization_vectors[-1],
                                 ray.directions[-1], Phenomenon.REFRACTION)
        else:
            state = shure_refraction(ray.directions[-1], normal_vector, n1, n2,
                                     ray.polarization_vectors[-1],
                                     perpendicular_polarized)
        ray.current_medium = nearby_material
        return state

    def change_of_direction(self, ray, normal_vector, nearby_material):
        if isinstance(self, TwoLayerMaterial):
            if ray.directions[-1].dot(normal_vector) < 0:
                # Ray intercepted on the frontside of the surface
                material = self.front_material
            else:  # Ray intercepted on the backside of the surface
                material = self.back_material
        else:
            material = self

        results = material.decide_phenomenon(ray, normal_vector, nearby_material)
        phenomenon = results[0]
        if ray.polarization_vectors[-1] == results[1]:
            # polarization_vector not calculated
            polarization_vector_calculated_before = False
        else:
            # polarization_vector calculated before
            polarization_vector_calculated_before = True
        ray.polarization_vectors[-1] = results[1]
        perpendicular_polarized = results[2]  # True or False
        if phenomenon == Phenomenon.REFLEXION:
            return material.change_of_direction_by_reflexion(
                ray, normal_vector,
                polarization_vector_calculated_before)
        elif phenomenon == Phenomenon.ABSORPTION:
            return material.change_of_direction_by_absortion(
                ray, normal_vector)
        elif phenomenon == Phenomenon.TRANSMITTANCE:
            return material.change_of_direction_by_transmitance(
                ray, normal_vector, nearby_material,
                perpendicular_polarized)

    @classmethod
    def from_plain_properties(cls, plain_properties):
        properties = Material.plain_properties_to_properties(plain_properties)
        material = cls(properties,properties)
        return material

    def to_json(self):
        return json.dumps(
            {
                'name': self.name,
                'kind': self.kind,
                'classname': self.__class__.__name__,
                'plain_properties': self.properties.get(
                    'plain_properties', None),
            }, cls=NumpyEncoder
        )


@traced(logger)
class OpaqueSimpleLayer(SurfaceMaterial):
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


@traced(logger)
class TransparentSimpleLayer(SurfaceMaterial):
    def __init__(self,name,pot):
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


@traced(logger)
class AbsorberSimpleLayer(SurfaceMaterial):
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


@traced(logger)
class AbsorberLambertianLayer(SurfaceMaterial):
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
        super(AbsorberLambertianLayer,self).__init__(name, properties)


@traced(logger)
class AbsorberTWModelLayer(SurfaceMaterial):
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
        properties = Material.plain_properties_to_properties(plain_properties)
        super(AbsorberTWModelLayer,self).__init__(name, properties)


@traced(logger)
class ReflectorSpecularLayer(SurfaceMaterial):
    def __init__(self,name, por, sigma_1=None, sigma_2=None, k=None):
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
        super(ReflectorSpecularLayer,self).__init__(name, properties)


@traced(logger)
class ReflectorLambertianLayer(SurfaceMaterial):
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
        super(ReflectorLambertianLayer,self).__init__(name, properties)


@traced(logger)
class MetallicSpecularLayer(SurfaceMaterial):
    def __init__(self, name,file_index_of_refraction,
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
        super(MetallicSpecularLayer,self).__init__(name, properties)


@traced(logger)
class MetallicLambertianLayer(SurfaceMaterial):
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
        super(MetallicLambertianLayer,self).__init__(name, properties)


@traced(logger)
class PolarizedCoatingReflectorLayer(SurfaceMaterial):
    def __init__(self, name, coating_file, sigma_1=None, sigma_2=None, k=None):
        # coating_material with four columns: wavelenth in nm,
        # angle in deg., reflectance s-polarized (perpendicular),
        # reflectance p-polarized (parallel)
        # the values in coating_material should be in the corresponding
        # order columns
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
        properties = Material.plain_properties_to_properties(plain_properties)
        super(PolarizedCoatingReflectorLayer,self).__init__(name, properties)


@traced(logger)
class PolarizedCoatingTransparentLayer(SurfaceMaterial):
    def __init__(self, name, coating_file):
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
        properties = Material.plain_properties_to_properties(plain_properties)
        super(PolarizedCoatingTransparentLayer,self).__init__(name, properties)


@traced(logger)
class PolarizedCoatingAbsorberLayer(SurfaceMaterial):
    def __init__(self, name, coating_file):
        # coating_material with four columns: wavelenth in nm, angle in deg.,
        # reflectance s-polarized (perpendicular),
        # reflectance p-polarized (parallel)
        # the values in coating_material should be in the corresponding order
        # columns
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
        properties = Material.plain_properties_to_properties(plain_properties)
        super(PolarizedCoatingAbsorberLayer,self).__init__(name, properties)


@traced(logger)
class TwoLayerMaterial(SurfaceMaterial):
    def __init__(self, name, name_front_layer, name_back_layer):
        super(SurfaceMaterial, self).__init__(name, {})
        self.name_front_layer = name_front_layer
        self.name_back_layer = name_back_layer
        self.front_material = Material.by_name[name_front_layer]
        self.back_material = Material.by_name[name_back_layer]
        self.kind = 'Surface'

    def to_json(self):
        return json.dumps(
            {
                'name': self.name,
                'kind': self.kind,
                'classname': 'TwoLayerMaterial',
                'name_front_layer': self.name_front_layer,
                'name_back_layer': self.name_back_layer
            }, cls=NumpyEncoder
        )
