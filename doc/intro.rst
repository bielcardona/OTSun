Introduction
============

The python package OTSun (imported as ``otsun``) is composed by
different modules, some of them collecting helper functions that
implement both mathematical and optical methods, and other ones built
around each of the main classes that define the functionality of the
package. To ease the reading of this section, where we shall be dealing
with classes and instances of those classes, whenever a class (say
``Scene`` or ``LightSource``) is considered, the downcase version of
their identifiers (``scene`` and ``light_source`` in our example) will
indicate instances of those classes. Except for those classes defined by
FreeCAD (like ``Part.Face`` or ``Base.Vector``), we refer to the
corresponding section in this manuscript for the definition and
initialization options for each of these classes. For instance, later
we explain how to create an ``experiment`` (an instance of the
``Experiment`` class), giving an ``scene`` and a ``light_source``, which
are instances of ``Scene`` and ``LightSource``, described in their respective
sections.

A typical use of the package involves the creation of an ``experiment``,
which specifies the solar optics experiment to be run, defined by
certain data that includes a ``scene`` and a ``light_source``. The
``scene`` holds the data of all the different objects that interact with
light rays, included in a FreeCAD document,
and where each of them has a ``material`` associated describing its
optical behaviour. Eventually, some objects may have movements,
implemented by a ``multi_tracking``, and in such a case these elements
will be moved to maximize the absorbed energy. When the experiment is
run, the ``light_source`` creates ``ray``\ s, which interact with the
scene until they either leave the scene or are absorbed, and in this
last case the collected energy (among other data) is stored for future
analysis.

We comment now the main classes that have been implemented, together
with its basic functionality. The complete documentation of the API can
be found at https://otsun.readthedocs.io.

.. _experiment-class:

The ``Experiment`` class
------------------------

An ``experiment`` is initialized giving the parameters that define it:
An ``scene`` and a ``light_source`` that describe the physical
environment where the experiment takes place, and the number of rays
that have to be simulated. The execution of the experiment is launched
with the ``experiment.run()`` method, and once it is finished, the
information that has been recollected is found in instance variables
like ``experiment.captured_energy_Th`` and
``experiment.captured_energy_PV``, that give the overall thermal and
photovoltaic (respectively) energy that has been collected by the active
materials found in the scene.

.. _scene-class:

The ``Scene`` class
-------------------

Instances of the class ``Scene`` hold the data used to describe the
physical objects present in an experiment, stored in three main
variables, ``faces``, ``solids`` and ``materials``. Each object in the
array ``faces`` (resp. ``solids``) is a ``Part.Face`` (resp.
``Part.Solid``) object of FreeCAD that represents a surface (resp.
volume) that can affect the propagation of a ray incident with it. The
dictionary ``materials`` assigns to each face or solid a ``material``
that describes its optical properties.

Such a ``scene`` is initialized with an array ``objects``, all whose
elements are instances of ``Part.Feature``, and typically they are all
the objects included in a FreeCAD document. The assignation of materials
to each object is done by looking at its label. Namely, if an object
``obj`` in a FreeCAD document has a label of the form
“Label\ ``(mat_name)``”, then the assigned material
``scene.materials[obj]`` will be the ``Material`` instance ``mat`` such
that ``mat.name`` is ``mat_name``.

For instance, the file ``ParabolicTrough.FCStd`` in the public
repository https://github.com/bielcardona/OTSunSuppMat contains a model prepared to analyze a parabolic trough
collector. The parabolic mirror (in blue) is made by extruding a
parabolic segment, and its label is ``Parabolic_reflector(Mir1)``. It
means that when imported with OTSun, it will have an associated material
whose parameter ``name`` is ``Mir1``. In turn, this material has to be
properly defined so that it behaves as a mirror. Other
elements in this model are the central cylindrical surface (in red),
labeled ``Cylindrical_absorber(Abs1)``, and its covering (in green),
labeled ``Tube_glass(Glass1)``; hence, materials named ``Abs1`` and
``Glass1`` have to be defined as an absorber surface and as a
transparent volume material, respectively.

.. _lightsource-class:

The ``LightSource`` class
-------------------------

Instances of the class ``LightSource`` are used to model the source of
rays in an experiment. There are many parameters that define its
behaviour, like its ``emitting_region``, describing the physical
location of the source of the rays to be emitted, and its
``light_spectrum`` and ``direction_distribution``, describing
respectively the distribution of wavelengths and directions of the rays
to be emitted.

The parameter ``light_spectrum`` can either be a constant, meaning that
all rays will be emitted with the same specified wavelength (given in
nanometers), or a cumulative distribution function (CDF)
:math:`F(\lambda)` which is defined by interpolation on the discrete
values :math:`(\lambda_i,F(\lambda_i))` stored in an array
:math:`((\lambda_1,\lambda_2,\ldots),(F(\lambda_1),F(\lambda_2)\ldots)`.

The ``emitting_region`` has to be an instance of any class that
implements the method ``random_point()``, which returns a random point
from where a ray will be emitted, and has an attribute
``main_direction``, giving the direction of the emitted ray. For
convenience, the class ``SunWindow`` implements such an emitting region
as a plane rectangle :math:`\Pi` in the space, orthogonal to a fixed
direction :math:`u`, and such that all the objects in the scene are
contained in the rectangular semi-prism
:math:`\{\Pi+\xi u\mid \xi\ge 0\}`.

The parameter ``direction_distribution`` can either be ``None`` (meaning
that the emitted rays are emitted in the main direction) or an instance
of a class that implements the method ``angle_distribution()``, giving a
random angle (in degrees) of deviation for the emitted ray with respect
to the main direction :math:`u`. For convenience, the class
``BuieDistribution`` implements such deviation according to the Buie
distribution, determined by its circumsolar
ratio (CSR), which is a parameter of the class.

The ``Ray`` class
-----------------

Instances of the class ``Ray`` model light rays, which are emitted by a
``light_source``. A ``ray`` is initialized giving its initial
``optical_state``, as well as the ``scene`` where it will travel.
Instances of ``OpticalState`` gather some relevant information of a
light ray at a given moment, like its ``direction``, ``polarization``
and ``material``, giving, respectively, the direction and polarization
vector of the ray, and the material (optical medium) where it is
traveling. When the method ``ray.run()`` is called, the propagation of
the ray inside the scene starts to be simulated. A simplified version of
the iteration process is:

#. Find the closest intersection of ``ray`` with objects in ``scene``.

#. If no intersection is found, the ``ray`` is lost and the simulation
   is finished.

#. If the first intersection is with an object having a determined
   ``material``, then the method ``material.change_of_optical_state()``
   is called (with different parameters that determine how the ray hits
   the material), which decides if the ray is reflected or refracted
   (and gives the next optical state) or that the ray has been absorbed
   by some active optical element.

#. If the ray has been reflected or refracted, go to step 1. Otherwise,
   the simulation is finished.

.. _material-class:

The ``Material`` class
----------------------

The ``Material`` class is the most complex of all the classes
implemented in OTSun, since there are many kinds of materials, and their
optical properties need to be explicitly defined. There are two main
subclasses, ``SurfaceMaterial`` and ``VolumeMaterial``, corresponding,
respectively, to materials that can be assumed to be two-dimensional
(like first surface mirrors and selective absorbers) or not (like
glasses, second surface mirrors, PV active materials, thin films, …).
Any material has an important property, ``material.name``, indicating
how it will be called when identifying objects in a scene. The physical properties of a ``material`` are encoded in
``material.properties``, a dictionary whose contents depend on the kind
of material.

Any user willing to use his own materials in his experiments can
subclass ``SurfaceMaterial`` or ``VolumeMaterial`` to adapt the contents
of ``material.properties``, which implement the specific properties of
the materials. The user must override the method
``material.change_of_optical_state()`` to implement the computation of
how the interaction with the material changes the optical state
(direction, polarization, etc.) of a ray.

Additionally, since it is interesting to store externally the properties
of materials, the method ``material.to_json()`` and the class method
``SubclassedMaterial.load_from_json(info)`` should be implemented. The
first one must convert any information stored in ``material.properties``
into a serializable dictionary, and the second one must use this
dictionary to reconstruct the ``material.properties`` dictionary.

The ``VolumeMaterial`` class
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Instances of ``VolumeMaterial`` represent the optical properties of
physical objects whose depth is not negligible, like glasses or PV
active materials, where the ray energy attenuation is determined by the
Beer–Lambert law. In this case, the method
``material.change_of_optical_state()`` is generically implemented using
the law of reflection, Snell’s law, and Fresnel’s equations, but
any user could subclass it and implement some other optical behaviour of
the material.

Some subclasses of this class are provided, so that materials appearing
usually in the field of solar collectors can be used without further
implementation. For example:

-  ``SimpleVolumeMaterial``, representing a material with constant
   optical parameters (refraction index and absorption coefficient,
   given in :math:`\textrm{mm}^{-1}`).

-  ``WavelengthVolumeMaterial``, where the index of refraction is
   complex (:math:`\tilde n =n - i\kappa`) and depends on the wavelength
   of the ray. These values are computed by interpolation from data
   given in tabulated form with rows
   :math:`(\lambda, n(\lambda),\kappa(\lambda))`. Note that the
   imaginary part of the refractive index is the so called the
   extinction coefficient, and the absorption coefficient is calculated
   as :math:`\alpha = 4 \pi \kappa / \lambda`. The wavelengths are given
   in nanometers.

-  ``PolarizedThinFilm``, which represents a thin layer, such as an
   optical coating, where the thickness and light coherence (that
   enables interference) can not be considered as negligible in the
   simulation. The data values are given in tabulated form with rows
   :math:`(\lambda, \theta, R_s(\lambda,\theta), R_p(\lambda,\theta), T_s(\lambda,\theta), T_p(\lambda,\theta))`,
   where :math:`\theta` is the incidence angle, :math:`R` and :math:`T`
   denote the power reflection and transmission coefficients
   respectively, and sub-indexes :math:`s` and :math:`p` denote
   respectively the perpendicular and parallel ray polarization.
   Wavelengths are given in nanometers and incidence angles in degrees.
   We remark that it is precisely in this case where the ray equations
   are complemented by the so-called fully-coherent medium transfer
   matrix formalism (TMM).

-  ``PVMaterial``, which represents the active material in photovoltaic
   cells such as semiconductors or any other material with that
   functionality. The photo-absorption in such materials is characterized by
   their extinction coefficient. The values of the index of refraction
   :math:`(\tilde n =n - i\kappa)`, which depends on the light
   wavelength, are given in tabulated form as in the
   ``WavelengthVolumeMaterial`` case.

The ``SurfaceMaterial`` class
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Any ``surface_material`` represents a two-dimensional physical object,
in the sense that its third dimension is negligible, or simply that its
behaviour does not depend on it. Examples of these objects are front
surface mirrors, selective absorbers, metallic coatings, …. In a first
approximation, the interaction of a ray with such a material can result
in a reflection, an absorption or a transmittance, each with a given
probability that may depend on the wavelength of the ray and are stored
in the dictionary ``p=material.properties``. Hence,
``material.change_of_optical_state()`` generically implements these
different phenomena. This behaviour is also affected by other properties
of the material, like the booleans:

-  ``p['lambertian_material']``, indicating that, in the case of
   reflection, the direction of the reflected ray should be a random
   vector, instead of that computed using the law of reflection.

-  ``p['thermal_material']``, indicating that, in case of absorption,
   the energy is absorbed and processed, instead of lost in the
   material.

Some more specific materials are provided by subclassing
``SurfaceMaterial`` and overriding the ``change_of_optical_state()``
method. Some examples of these specific materials are:

-  ``AbsorberTWModelLayer``, represents a thermal absorber where its
   absorption depends on the incidence angle, :math:`\theta`, according
   to
   :math:`\alpha =\alpha_{0}  \{ 1-b (\frac{1}{\cos \theta} -1  )^c  \}`. The following
   data values are given: :math:`\alpha_{0}, {b}, {c}`. In this case,
   the boolean property ``p['thermal_material']`` is ``True``.

-  ``MetallicSpecularLayer``, represents a metal surface, such as the
   silver coating in second surface mirrors. Fresnel equations are
   considered and its characterization is defined by the complex index
   of refraction :math:`(\tilde n =n - i\kappa)` depending on the light
   wavelength. The data values are given in tabulated form like in the
   ``WavelengthVolumeMaterial`` case.

-  ``MetallicLambertianLayer``, represents a metal surface where Fresnel
   equations are considered, but if the ray is reflected, a total
   diffuse reflection model with Lambertian scattering is used. In this
   material, the boolean property ``p['lambertian_material']`` is
   ``True``. Also, its characterization is defined by the complex index
   of refraction :math:`(\tilde n =n - i\kappa)` depending on the light
   wavelength. The data values are given in tabulated form like in the
   ``WavelengthVolumeMaterial`` case.

-  ``PolarizedCoatingLayer``, and its subclasses
   ``PolarizedCoatingReflectorLayer``,
   ``PolarizedCoatingTransparentLayer``,
   ``PolarizedCoatingAbsorberLayer``, that represent thin layers such as
   optical coatings. The difference with the ``PolarizedThinFilm`` is
   that the thickness of such material is negligible. The data values
   are given as in the ``PolarizedThinFilm`` case. Depending on the role
   of the material, three cases are defined: reflector (no light
   transmission is possible), transparent (reflection, absorption and
   transmission are possible), and thermal absorber material (the
   boolean property ``p['thermal_material']`` is ``True`` and no light
   transmission is possible). In each case, the parameters are given
   analogously to the case of ``PolarizedThinFilm``.

The ``MultiTracking`` class
---------------------------

The class ``MultiTracking`` is designed to implement movements of the
active elements in a ``scene`` so that the rays emitted by a given
``light_source`` tend to be focused on a target (in case that the
attribute ``target`` is set to a point) or tend to return it to the
source (in case that the attribute is not set). That is,
``MultiTracking`` can be used either to orient the solar collector to
the sun or to direct rays to a target, as happens with the segment
mirrors of a Linear Fresnel Collector (LFR) or the heliostats in solar
power tower plants.

Movements of elements are implemented by the helper class ``Joint``, and
its subclasses ``CentralJoint`` and ``AxialJoint``. The former
implements rotations around a given point in space (that is, with two
degrees of freedom), while in the latter the rotations are around an
axis (and hence with a single degree of freedom). Each kind of joint can
be easily represented by a geometrical object in FreeCAD, either by a
``Vertex`` or an ``Edge`` with two points.

To describe the movement of a concrete element in the ``scene``, one
needs to associate to this object a ``joint``, but since the goal is to
direct the rays to a specified region, one also needs to specify the
corresponding *principal vector*. Here, by the *principal vector*, we
mean the direction that best approaches the normal of the mobile
element. When ``multi_tracking.target`` is not set, the element will be
moved so that this vector points to the source; otherwise, the movement
will be computed so that a solar ray reflected on the plane normal to
the *principal vector* and passing through the ``joint`` hits the point
stored in ``multi_tracking.target``.

We associate objects in the scene to joints using the following
convention: Instead of giving to the object under
consideration a label of the form “Label\ ``(mat_name)``”, where
``mat_name`` is the identifier of the ``material`` of the object, we use
a label of the form “Label\ ``(mat_name,joint_name,normal_name)``” or
“Label\ ``(mat_name,joint_name,normal_name,target_name)``”, where
``joint_name`` is the label of the FreeCAD object that describes the
joint (i.e. either a ``Vertex`` or a ``Edge``), ``normal_name`` is the
label of the FreeCAD ``Edge`` whose direction is the *principal vector*
of the optical element, and ``target_name`` (if present) is the label of
the FreeCAD object acting as target.

A ``multi_tracking`` is created by giving the ``scene`` (which includes
the elements that describe the joints, together with their principal
vectors and targets, if needed) and the ``light_source``, a
``Base.Vector`` giving the main direction of the sun rays. Once it is
created, the method ``target_tracking.make_movements()`` transforms the
scene, rotating conveniently the elements, so that the scene behaves as
explained above.

