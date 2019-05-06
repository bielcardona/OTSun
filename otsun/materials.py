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
        if properties is None:
            properties = {}
        self.properties = properties

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

    @staticmethod
    def plain_properties_to_properties(plain_properties):
        """
        Converts properties of a material in plain format (dict from json) to the internal format

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
                properties[key] = tabulated_function(np.array(prop_value[0]), np.array(prop_value[1]))
            if prop_type == 'matrix':
                properties[key] = SurfaceMaterial.matrix_reflectance(np.array(prop_value))
        properties['plain_properties'] = plain_properties
        return properties

    @staticmethod
    def properties_to_plain_properties(properties):
        """
        Converts properties of a material in internal format to plain (json ready) format

        Since the plain properties are stored in the internal format, no need for conversions

        Parameters
        ----------
        properties : dict

        Returns
        -------
        dict
        """
        return properties.get('plain_properties', None)

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
                mat = SurfaceMaterial(name, {})
                plain_properties_back = mat_spec['plain_properties_back']
                mat.properties_back = Material.plain_properties_to_properties(plain_properties_back)
                plain_properties_front = mat_spec['plain_properties_front']
                mat.properties_front = Material.plain_properties_to_properties(plain_properties_front)
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
        Computes how a ray behaves when interacting with the material. MUST be subclassed

        The material where the ray is actually located is held in ray.current_material.

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
        return ""

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


@traced(logger)
class VolumeMaterial(Material, object):
    """
    Subclass of Volume for materials with volume

    TODO: comment on properties

    """
    def __init__(self, name, properties=None):
        """
        Initializes a Volume Material.

        The `properties` parameter must be a dict with the physical properties
        describing the material. At least, the following must be provided:
        'index_of_refraction': index of refraction of the material, as a function of its wavelength, only real part.
        'extinction_coefficient': imaginary part of the index of refraction of the material as a function of its wavelength.

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
            # the ray is traveling inside thin film material and then it will certainly refract
            ray.current_medium = self
            return OpticalState(ray.polarization_vectors[-1], ray.directions[-1], Phenomenon.REFRACTION)
        # The next code has been moved to subclass
        # if 'Matrix_reflectance_thin_film' in self.properties:
        #     # the ray impacts on thin film material
        #     n1 = ray.materials[-1].properties['index_of_refraction'](ray.wavelength)
        #     if 'extinction_coefficient' in ray.current_medium.properties:
        #         n1 = n1 + 1j * ray.materials[-1].properties['extinction_coefficient'](ray.wavelength)
        #     n_front = self.properties['index_of_refraction_front'](ray.wavelength)
        #     n_back = self.properties['index_of_refraction_back'](ray.wavelength)
        #     if n_front == ray.materials[-1].properties['index_of_refraction'](ray.wavelength):
        #         # impacts on front material
        #         n2 = n_back
        #     else:
        #         n2 = n_front
        #     energy_absorbed_thin_film, optical_state = calculate_state_thin_film(ray.directions[-1], normal_vector, n1,
        #                                                                          n2,
        #                                                                          ray.polarization_vectors[-1],
        #                                                                          self.properties, ray.wavelength)
        #     ray.energy = ray.energy * (1.0 - energy_absorbed_thin_film)
        #     if optical_state.phenomenon == Phenomenon.REFRACTION:
        #         ray.current_medium = self
        #     return optical_state
        # else:
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
                'classname': self.__class__.__name__,
                'plain_properties': self.properties.get('plain_properties', None)
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
        super(PVMaterial,self).__init__(name)
        self.properties = Material.plain_properties_to_properties(plain_properties)


class PolarizedThinFilm(VolumeMaterial):
    def __init__(self, name, file_thin_film, file_front, file_back):
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

    @staticmethod
    def calculate_state_thin_film(incident, normal, n1, n2, polarization_vector, properties, wavelength):
        """

        Parameters
        ----------
        incident
        normal
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
        mynormal = normal * 1.0
        backside = False
        if mynormal.dot(incident) > 0:  # Ray intercepted on the backside of the surface
            mynormal = mynormal * (-1.0)
            backside = True
        r = n1.real / n2.real
        c1 = - mynormal.dot(incident)  # cos (incidence_angle)
        c2sq = 1.0 - r * r * (1.0 - c1 * c1)  # cos (refracted_angle) ** 2
        if c2sq.real < 0:  # total internal reflection
            return 0.0, reflexion(incident, normal, polarization_vector)
        c2 = sqrt(c2sq)  # cos (refracted_angle)

        normal_parallel_plane = incident.cross(mynormal)  # normal vector of the parallel plane
        if normal_parallel_plane == Base.Vector(0, 0,
                                                0):  # to avoid null vector at mynormal and incident parallel vectors
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

        if backside:  # Ray intercepted on the backside of the transparent surface
            angle = np.arccos(c2.real) * 180.0 / np.pi
        else:
            angle = np.arccos(c1) * 180.0 / np.pi
        reflectance_matrix = properties['Matrix_reflectance_thin_film']
        r_matrix = reflectance_matrix(angle, wavelength)
        if myrandom() < ref_per:
            r = SurfaceMaterial.calculate_reflectance(r_matrix, angle, wavelength)[
                0]  # reflectance for s-polarized (perpendicular) light
            perpendicular_polarized = True
            polarization_vector = perpendicular_v.normalize()
        else:
            angle = np.arccos(c1) * 180.0 / np.pi
            r = SurfaceMaterial.calculate_reflectance(r_matrix, angle, wavelength)[1]  # reflectance for p-polarized (parallel) light
            polarization_vector = parallel_v.normalize()
        if myrandom() < r:  # ray reflected
            return 0.0, reflexion(incident, normal, polarization_vector)
        else:
            transmittance_matrix = properties['Matrix_transmittance_thin_film']
            t_matrix = transmittance_matrix(angle, wavelength)
            if perpendicular_polarized:
                t = SurfaceMaterial.calculate_reflectance(t_matrix, angle, wavelength)[0]
            else:
                t = SurfaceMaterial.calculate_reflectance(t_matrix, angle, wavelength)[1]
            energy_absorbed_thin_film = (1 - r - t) / (1 - r)
            refracted_direction = incident * r.real + mynormal * (r.real * c1.real - c2.real)
            return energy_absorbed_thin_film, OpticalState(polarization_vector, refracted_direction,
                                                           Phenomenon.REFRACTION)

    def change_of_direction(self, ray, normal_vector):
        #here we assume 'Matrix_reflectance_thin_film' in self.properties:
        # the ray impacts on thin film material
        n1 = ray.materials[-1].properties['index_of_refraction'](ray.wavelength)
        if 'extinction_coefficient' in ray.current_medium.properties:
            n1 = n1 + 1j * ray.materials[-1].properties['extinction_coefficient'](ray.wavelength)
        n_front = self.properties['index_of_refraction_front'](ray.wavelength)
        n_back = self.properties['index_of_refraction_back'](ray.wavelength)
        if n_front == ray.materials[-1].properties['index_of_refraction'](ray.wavelength):
            # impacts on front material
            n2 = n_back
        else:
            n2 = n_front
        energy_absorbed_thin_film, optical_state = self.calculate_state_thin_film(
            ray.directions[-1], normal_vector, n1,
            n2,
            ray.polarization_vectors[-1],
            self.properties, ray.wavelength)
        ray.energy = ray.energy * (1.0 - energy_absorbed_thin_film)
        if optical_state.phenomenon == Phenomenon.REFRACTION:
            ray.current_medium = self
        return optical_state


vacuum_medium = SimpleVolumeMaterial("Vacuum", 1.0, 0.0)


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

    @staticmethod
    def calculate_reflexion_metallic(incident, normal, n1, n2, polarization_vector):
        """Implementation of Fresnel equations for metallic materials

        Parameters
        ----------
        incident
        normal
        n1
        n2
        polarization_vector

        Returns
        -------

        """
        # TODO: document
        mynormal = normal * 1.0
        if mynormal.dot(incident) > 0:  # Ray intercepted on the backside of the surface
            # noinspection PyAugmentAssignment
            mynormal = mynormal * (-1.0)
        r = n1 / n2
        c1 = - mynormal.dot(incident)  # cos (incidence_angle)
        c2 = sqrt(1.0 - r * r * (1.0 - c1 * c1))  # cos (refracted_angle)

        normal_parallel_plane = incident.cross(mynormal)  # normal vector of the parallel plane
        if normal_parallel_plane == Base.Vector(0, 0,
                                                0):  # to avoid null vector at mynormal and incident parallel vectors
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
            r = a * a.conjugate()  # reflectance for s-polarized (perpendicular) light
            perpendicular_polarized = True
            polarization_vector = perpendicular_v.normalize()
        else:
            a = (n1 * c2 - n2 * c1) / (n1 * c2 + n2 * c1)
            r = a * a.conjugate()  # reflectance for p-polarized (parallel) light
            polarization_vector = parallel_v.normalize()
        if myrandom() < r.real:  # ray reflected
            return 1, 0, 0, polarization_vector, perpendicular_polarized, True
        else:  # ray refracted
            return 0, 1, 0, polarization_vector, perpendicular_polarized, True

    # ---
    # Helper function for reflectance depending on the wavelength for coatings layers.
    # ---

    @staticmethod
    def matrix_reflectance(data_material):
        """
        data_material: wavelenth in nm, angle in deg., reflectance s-polarized (perpendicular), reflectance p-polarized (parallel)
        the values should be in the corresponding order columns with constants steps
    	We generate a matrix for the interpolation depending on the x=angle and y=wavelength values (lambda values)

        Parameters
        ----------
        data_material

        Returns
        -------

        """
        # TODO: document
        steps = np.diff(data_material[:, 1])
        try:
            a = min(steps[steps > 0.0])
        except ValueError:
            a = 1E-4
        steps = np.diff(data_material[:, 0])
        try:
            w = min(steps[steps > 0.0])
        except ValueError:
            w = 1E-4
        return lambda x, y: [row for row in data_material if
                             (x + a > row[1] > x - a) and (y + w > row[0] > y - w)]

    @staticmethod
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
        if len(r_matrix) == 1:  # interpolation is not needed
            return r_matrix[0, 2], r_matrix[0, 3]

        if len(r_matrix) == 2:  # simple interpolation
            if r_matrix[0, 0] == r_matrix[1, 0]:  # identical wavelengths; hence angle interpolation (column 1)
                x_values = [r_matrix[0, 1], r_matrix[1, 1]]
                y_values = [r_matrix[0, 2], r_matrix[1, 2]]
                r_per = np.interp(angle, x_values, y_values)

                x_values = [r_matrix[0, 1], r_matrix[1, 1]]
                y_values = [r_matrix[0, 3], r_matrix[1, 3]]
                r_par = np.interp(angle, x_values, y_values)
                return r_per, r_par

            if r_matrix[0, 1] == r_matrix[1, 1]:  # identical angles; hence wavelength interpolation (column 0)
                x_values = [r_matrix[0, 0], r_matrix[1, 0]]
                y_values = [r_matrix[0, 2], r_matrix[1, 2]]
                r_per = np.interp(wavelength, x_values, y_values)

                x_values = [r_matrix[0, 0], r_matrix[1, 0]]
                y_values = [r_matrix[0, 3], r_matrix[1, 3]]
                r_par = np.interp(wavelength, x_values, y_values)
                return r_per, r_par

        if len(r_matrix) == 4:  # double interpolation
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

    @staticmethod
    def calculate_probabilities_polarizaton_coating(incident, normal, n1, n2, polarization_vector, properties,
                                                    wavelength):
        """

        Parameters
        ----------
        incident
        normal
        n1
        n2
        polarization_vector
        properties
        wavelength

        Returns
        -------

        """
        # Only used if Matrix_polarized_reflectance_coating in properties
        # (PolarizedCoatingReflectorLayer, PolarizedCoatingTransparentLayer, PolarizedCoatingAbsorberLayer)
        # TODO: document
        # returns probability of Reflexion, probability of Absortion, probability of Transmitance, polarization_vector
        mynormal = normal * 1.0
        backside = False
        if mynormal.dot(incident) > 0:  # Ray intercepted on the backside of the surface
            # noinspection PyAugmentAssignment
            mynormal = mynormal * (-1.0)
            backside = True
        r = n1 / n2
        c1 = - mynormal.dot(incident)  # cos (incidence_angle)
        c2sq = 1.0 - r * r * (1.0 - c1 * c1)  # cos (refracted_angle) ** 2
        if properties['transparent_material']:  # transparent coating
            if c2sq.real < 0:  # total internal reflection
                return reflexion(incident, normal, polarization_vector)
        c2 = sqrt(c2sq)  # cos (refracted_angle)
        normal_parallel_plane = incident.cross(mynormal)  # normal vector of the parallel plane
        if normal_parallel_plane == Base.Vector(0, 0,
                                                0):  # to avoid null vector at mynormal and incident parallel vectors
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

        if backside == True and properties['transparent_material']:  # Ray intercepted on the backside of the surface
            angle = np.arccos(c2.real) * 180.0 / np.pi
        else:
            angle = np.arccos(c1) * 180.0 / np.pi
        reflectance_matrix = properties['Matrix_polarized_reflectance_coating']
        r_matrix = reflectance_matrix(angle, wavelength)
        if myrandom() < ref_per:
            r = SurfaceMaterial.calculate_reflectance(r_matrix, angle, wavelength)[
                0]  # reflectance for s-polarized (perpendicular) light
            perpendicular_polarized = True
            polarization_vector = perpendicular_v.normalize()
        else:
            angle = np.arccos(c1) * 180.0 / np.pi
            r = SurfaceMaterial.calculate_reflectance(r_matrix, angle, wavelength)[1]  # reflectance for p-polarized (parallel) light
            polarization_vector = parallel_v.normalize()
        if myrandom() < r:  # ray reflected
            return 1, 0, 0, polarization_vector, perpendicular_polarized
        else:  # ray refracted or absorbed
            if properties['energy_collector']:  # absorber coating
                return 0, 1, 0, polarization_vector, perpendicular_polarized
            if properties['specular_material']:  # reflector coating
                return 0, 1, 0, polarization_vector, perpendicular_polarized
            if properties['transparent_material']:  # transparent coating
                return 0, 0, 1, polarization_vector, perpendicular_polarized

    @staticmethod
    def compute_polarizations(ray, normal_vector, properties, nearby_material):
        polarization_vector = ray.polarization_vectors[-1]
        perpendicular_polarized = False
        if 'Matrix_polarized_reflectance_coating' in properties:  # polarized_coating_layer
            n1 = ray.current_medium.properties['index_of_refraction'](ray.wavelength)
            if 'extinction_coefficient' in ray.current_medium.properties:
                n1 = ray.current_medium.properties['index_of_refraction'](ray.wavelength) + 1j * \
                     ray.current_medium.properties['extinction_coefficient'](ray.wavelength)
            n2 = nearby_material.properties['index_of_refraction'](ray.wavelength)
            if 'extinction_coefficient' in nearby_material.properties:
                n2 = nearby_material.properties['index_of_refraction'](ray.wavelength) + 1j * \
                     nearby_material.properties['extinction_coefficient'](ray.wavelength)
            results = SurfaceMaterial.calculate_probabilities_polarizaton_coating(
                ray.directions[-1], normal_vector, n1, n2,
                ray.polarization_vectors[-1],
                properties, ray.wavelength)
            polarization_vector = results[3]
            perpendicular_polarized = results[4]
        elif 'metallic_material' in properties:
            n1 = ray.current_medium.properties['index_of_refraction'](ray.wavelength)
            if 'extinction_coefficient' in ray.current_medium.properties:
                n1 = ray.current_medium.properties['index_of_refraction'](ray.wavelength) + 1j * \
                     ray.current_medium.properties['extinction_coefficient'](ray.wavelength)
            n2 = properties['index_of_refraction'](ray.wavelength)
            if 'extinction_coefficient' in properties:
                n2 = properties['index_of_refraction'](ray.wavelength) + 1j * properties['extinction_coefficient'](
                    ray.wavelength)
            results = SurfaceMaterial.calculate_reflexion_metallic(ray.directions[-1], normal_vector, n1, n2, polarization_vector)
            polarization_vector = results[3]
            perpendicular_polarized = results[4]

        else:
            pass

        return polarization_vector, perpendicular_polarized


    @staticmethod
    def compute_probabilities(ray, normal_vector, properties, nearby_material):
        polarization_vector = ray.polarization_vectors[-1]
        perpendicular_polarized = False
        probabilities = None
        # if 'TW_model' in properties:
        #     b_constant = properties['b_constant']
        #     c_constant = properties['c_constant']
        #     absortion_ratio = AbsorberTWModelLayer.tw_absorptance_ratio(normal_vector, b_constant, c_constant, ray.directions[-1])
        #     absortion = properties['probability_of_absortion'](ray.properties['wavelength']) * absortion_ratio
        #     por = 1.0 - absortion
        #     probabilities = [por, absortion, 0]  # Here I assume no transmitance
        if 'Matrix_polarized_reflectance_coating' in properties:  # polarized_coating_layer
            n1 = ray.current_medium.properties['index_of_refraction'](ray.wavelength)
            if 'extinction_coefficient' in ray.current_medium.properties:
                n1 = ray.current_medium.properties['index_of_refraction'](ray.wavelength) + 1j * \
                     ray.current_medium.properties['extinction_coefficient'](ray.wavelength)
            n2 = nearby_material.properties['index_of_refraction'](ray.wavelength)
            if 'extinction_coefficient' in nearby_material.properties:
                n2 = nearby_material.properties['index_of_refraction'](ray.wavelength) + 1j * \
                     nearby_material.properties['extinction_coefficient'](ray.wavelength)
            results = SurfaceMaterial.calculate_probabilities_polarizaton_coating(
                ray.directions[-1], normal_vector, n1, n2,
                ray.polarization_vectors[-1],
                properties, ray.wavelength)
            probabilities = [results[0], results[1], results[2]]
        if 'metallic_material' in properties:
            n1 = ray.current_medium.properties['index_of_refraction'](ray.wavelength)
            if 'extinction_coefficient' in ray.current_medium.properties:
                n1 = ray.current_medium.properties['index_of_refraction'](ray.wavelength) + 1j * \
                     ray.current_medium.properties['extinction_coefficient'](ray.wavelength)
            n2 = properties['index_of_refraction'](ray.wavelength)
            if 'extinction_coefficient' in properties:
                n2 = properties['index_of_refraction'](ray.wavelength) + 1j * properties['extinction_coefficient'](
                    ray.wavelength)
            results = SurfaceMaterial.calculate_reflexion_metallic(ray.directions[-1], normal_vector, n1, n2, polarization_vector)
            probabilities = [results[0], results[1], results[2]]

        if probabilities is None:
            try:
                por = properties['probability_of_reflexion'](ray.properties['wavelength'])
            except KeyError:
                por = 1.0
            try:
                poa = properties['probability_of_absortion'](ray.properties['wavelength'])
            except KeyError:
                poa = 1 - por
            try:
                pot = properties['probability_of_transmitance'](ray.properties['wavelength'])
            except KeyError:
                pot = 0.0

            probabilities = [por, poa, pot]
        return probabilities

    @staticmethod
    def decide_phenomenon(ray, normal_vector, properties, nearby_material):
        # TODO : Humanize
        phenomena = [Phenomenon.REFLEXION, Phenomenon.ABSORPTION, Phenomenon.TRANSMITTANCE]
        probabilities = SurfaceMaterial.compute_probabilities(
            ray, normal_vector, properties, nearby_material)
        polarization_vector, perpendicular_polarized = SurfaceMaterial.compute_polarizations(
            ray, normal_vector, properties, nearby_material)

        phenomenon = np.random.choice(phenomena, 1, p=probabilities)[0]
        return phenomenon, polarization_vector, perpendicular_polarized

    @staticmethod
    def change_of_direction_by_absortion(ray, normal_vector, properties):
        if properties['energy_collector']:
            return OpticalState(Base.Vector(0.0, 0.0, 0.0), Base.Vector(0.0, 0.0, 0.0), Phenomenon.GOT_ABSORBED)
        else:
            return OpticalState(Base.Vector(0.0, 0.0, 0.0), Base.Vector(0.0, 0.0, 0.0), Phenomenon.ABSORPTION)

    @staticmethod
    def change_of_direction_by_reflexion(ray, normal_vector, properties, polarization_vector_calculated_before):
        if 'specular_material' in properties:
            state = reflexion(ray.directions[-1], normal_vector, ray.polarization_vectors[-1],
                              polarization_vector_calculated_before)
            if properties.get('sigma_1',None):
                sigma_1 = properties['sigma_1']
                if properties.get('sigma_2',None):
                    sigma_2 = properties['sigma_2']
                    k = properties.get('k', None) or 0.5
                    return double_gaussian_dispersion(normal_vector, state, sigma_1, sigma_2, k)
                return single_gaussian_dispersion(normal_vector, state, sigma_1)
            return state
        if 'lambertian_material' in properties:
            state = lambertian_reflexion(ray.directions[-1], normal_vector)
            return state

    @staticmethod
    def change_of_direction_by_transmitance(ray, normal_vector, nearby_material, perpendicular_polarized):
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
        elif phenomenon == Phenomenon.ABSORPTION:
            return self.change_of_direction_by_absortion(ray, normal_vector, properties)
        elif phenomenon == Phenomenon.TRANSMITTANCE:
            return self.change_of_direction_by_transmitance(ray, normal_vector, nearby_material,
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
                'plain_properties_back': self.properties_back.get('plain_properties', None),
                'plain_properties_front': self.properties_front.get('plain_properties', None),
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
        super(OpaqueSimpleLayer, self).__init__(name, properties, properties)


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
        super(TransparentSimpleLayer, self).__init__(name, properties, properties)


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
        super(AbsorberSimpleLayer, self).__init__(name, properties, properties)


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
        super(AbsorberLambertianLayer,self).__init__(name, properties, properties)


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
        super(AbsorberTWModelLayer,self).__init__(name, properties, properties)

    @staticmethod
    def tw_absorptance_ratio(normal, b_constant, c_constant, incident):
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
        # Only used if TW_model in properties (AbsorberTWModelLayer)
        # We assume the normal is normalized.
        my_normal = normal * 1.0
        if my_normal.dot(incident) > 0:  # Ray intercepted on the backside of the surface
            my_normal = my_normal * (-1.0)
        incidence_angle = np.arccos(my_normal.dot(incident) * (-1.0))
        incidence_angle_deg = incidence_angle * 180.0 / np.pi
        if incidence_angle_deg < 80.0:
            absorption_ratio = 1.0 - b_constant * (1.0 / np.cos(incidence_angle) - 1.0) ** c_constant
        else:
            y0 = 1.0 - b_constant * (1.0 / np.cos(80.0 * np.pi / 180.0) - 1.0) ** c_constant
            m = y0 / 10.0
            absorption_ratio = y0 - m * (incidence_angle_deg - 80.0)
        return absorption_ratio

    @staticmethod
    def compute_probabilities(ray, normal_vector, properties, nearby_material):
        b_constant = properties['b_constant']
        c_constant = properties['c_constant']
        absortion_ratio = AbsorberTWModelLayer.tw_absorptance_ratio(
            normal_vector, b_constant, c_constant,ray.directions[-1])
        absortion = properties['probability_of_absortion'](
            ray.properties['wavelength']) * absortion_ratio
        por = 1.0 - absortion
        return [por, absortion, 0]  # Here I assume no transmitance


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
        super(ReflectorSpecularLayer,self).__init__(name, properties, properties)


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
        super(ReflectorLambertianLayer,self).__init__(name, properties, properties)


@traced(logger)
class MetallicLayer(SurfaceMaterial):
    def __init__(self, *args):
        super(MetallicLayer, self).__init__(*args)


@traced(logger)
class MetallicSpecularLayer(MetallicLayer):
    def __init__(self, name,file_index_of_refraction, sigma_1=None, sigma_2=None, k=None):
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
        properties = Material.plain_properties_to_properties(plain_properties)
        super(MetallicSpecularLayer,self).__init__(name, properties, properties)


@traced(logger)
class MetallicLambertianLayer(MetallicLayer):
    def __init__(self, name, file_index_of_refraction):
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
        properties = Material.plain_properties_to_properties(plain_properties)
        super(MetallicLambertianLayer,self).__init__(name, properties, properties)


@traced(logger)
class PolarizedCoatingReflectorLayer(SurfaceMaterial):
    def __init__(self, name, coating_file, sigma_1=None, sigma_2=None, k=None):
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
        properties = Material.plain_properties_to_properties(plain_properties)
        super(PolarizedCoatingReflectorLayer,self).__init__(name, properties, properties)


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
        super(PolarizedCoatingTransparentLayer,self).__init__(name, properties, properties)


@traced(logger)
class PolarizedCoatingAbsorberLayer(SurfaceMaterial):
    def __init__(self, name, coating_file):
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
        properties = Material.plain_properties_to_properties(plain_properties)
        super(PolarizedCoatingAbsorberLayer,self).__init__(name, properties, properties)


@traced(logger)
class TwoLayerMaterial(SurfaceMaterial):
    def __init__(self, name, layer_front, layer_back):
        super(TwoLayerMaterial,self).__init__(name,
                                              Material.by_name[layer_front].properties_front,
                                              Material.by_name[layer_back].properties_back)

