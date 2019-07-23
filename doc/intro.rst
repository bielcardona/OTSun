Introduction
============

The python package ``otsun`` is composed by different modules, some of
them collecting helper functions that implement both mathematical and
optical methods, and other ones built around each of the main classes
that define the functionality of the package. To ease the reading of
this section, where we shall be dealing with classes and instances of
those classes, whenever a class (say ``Scene`` or ``LightSource``) is
considered, the downcase version of their identifiers (``scene`` and
``light_source`` in our example) will indicate instances of those
variables.

A typical use of the package involves the creation of an ``experiment``,
which specifies the solar experiment to be run, defined by certain data
that includes an ``scene`` and a ``light_source``. The ``scene`` holds
the data of all the different objects that interact with light rays,
included in a FreeCAD document, and where
each of them has associated a ``material`` describing its optical
behaviour. When the experiment is run, the ``light_source`` creates
``ray``\ s, which interact with the scene until they either leave the
scene or are absorbed, and in this last case the collected energy (among
other data) is stored for future analysis.

We comment now the main classes that have been implemented, together
with its basic functionallity. The complete documentation can be found
in XXXXXX.

The ``Experiment`` class
------------------------

An ``experiment`` is initialized giving the parameters that define it:
An ``scene`` and a ``light_source`` that describe the physical
environment where the experiment takes place, and the number of rays
that has to be simulated. The execution of the experiment is launched
with the ``experiment.run()`` method, and once it is finished, the
information that has been recollected is found in instance variables
like ``experiment.captured_energy_Th`` and
``experiment.captured_energy_PV``, that give the overall thermal and
photovoltaic (respectively) energy that has been collected by active
materials found in the scene.

The ``Scene`` class
-------------------

Instances of the class ``Scene`` hold the data used to describe the
physical objects present in an experiment, stored in three main
variables, ``faces``, ``solids`` and ``materials``. Each object in the
array ``faces`` (resp. ``solids``) is a ``Part.Face`` (resp.
``Part.Solid``) object of FreeCAD that represents a surface (resp.
volume) that can affect the propagation of a ray incident with it. The
dictionary ``materials`` assigns to each face or solid a ``material``
that describes its properties.

Such an ``scene`` is initialized with an array ``objects``, all whose
elements are instances of ``Part.Feature``, and typically they are all
the objects included in a FreeCAD document. The assignation of materials
to each element is done by looking at its label. Namely, if an object
``obj`` in a FreeCAD document has a label of the form
“Label\ ``(mat_name)``”, then the assigned material
``scene.materials[obj]`` will be the ``Material`` instance ``mat`` such
that ``mat.name`` is ``mat_name``.

The ``LightSource`` class
-------------------------

Instances of the class ``LightSource`` are used to model the sources of
rays in an experiment. There are many parameters that define its
behaviour, like its ``emitting_region``, describing the physical
location of the sources of the rays to be emitted, and its
``light_spectrum`` and ``direction_distribution``, describing
respectively the frequencies and directions of the rays to be emitted.

The parameter ``light_spectrum`` can either be a constant, meaning that
all rays will be emitted with the same specified frequency, or a
probability distribution with CDF :math:`F(\omega)` represented by an
array of the form
:math:`((\omega_1,\omega_2,\dots),(F(\omega_1),F(\omega_2)\dots)`.

The ``emitting_region`` has to be an instance of any class that
implements the method ``random_point()``, which returns a random point
from where a ray will be emitted, and has an attribute
``main_direction``, giving the direction of the emitted ray. For
convenience, the class ``SunWindow`` implements such an emitting region
as a rectangle in space :math:`\Pi` orthogonal to a fixed direction
:math:`\vec u`, and such that all the objects in the scene are contained
in the rectangular semi-prism
:math:`\{\Pi+\lambda\vec u\mid \lambda\ge 0\}`.

The parameter ``direction_distribution`` can either be ``None`` or an
instance of a class that implements the method ``angle_distribution()``,
giving a random angle of deviation for the emitted ray with respect to
the main direction. For convenience, the class ``BuieDistribution``
implements such deviation according to a Buie distribution, determined
by its circumsolar ratio, which is a parameter of the class.

The ``Ray`` class
-----------------

Instances of the class ``Ray`` model light rays, which are emitted by a
``light_source``. A ``ray`` is initialized giving its initial
``optical_state``, as well as the ``scene`` where it will travel.
Instances of ``OpticalState`` gather some relevant information of a
light ray at a given moment, like its ``direction``, ``polarization``
and ``material``, giving respectively the direction and polarization
vectors, and the material where it is traveling. When the method
``ray.run()`` is called, the propagation of the ray inside the scene
starts to be simulated. A simplified version of the iteration process
is:

#. Find the closest intersection of ``ray`` with objects in ``scene``.

#. If no intersection is found, the ``ray`` is lost and the simulation
   is finished.

#. If the first intersection is with an object having a determined
   ``material``, then the method ``material.change_of_optical_state()``
   is called (with different parameters that determine how the ray hits
   the material), which decides if the ray is reflected or refracted
   (and gives the next optical state) or that the ray has been absorbed
   by some active optical element.

#. If the ray has been reflected or refracted, go to step 1.

The ``Material`` class
----------------------

The ``Material`` class is the most complex of all the classes
implemented in ``otsun``, since there are many different kinds of
materials, and their optical properties need to be explicitly stated.
There are two main subclasses, ``SurfaceMaterial`` and
``VolumeMaterial``, corresponding respectively to materials that can be
assumed to be two-dimensional (like mirrors and absorbers) or not (like
glasses, XXXX). Any material has an important property,
``material.name``, indicating how it will be called when identifying
objects in a scene. The physical properties of a ``material`` are
encoded in ``material.properties``, a dictionary whose contents depend
on the kind of material.

Any user willing to use his own materials in his experiments can
subclass ``SurfaceMaterial`` or ``VolumeMaterial`` to adapt the contents
of ``material.properties``, that implement the specific properties of
the materials. The user must subclass the method
``material.change_of_optical_state()`` in order to implement the
computation of how the interaction with the material changes the optical
state (direction, polarization, etc.) of the ray. Also, since it is
interesting to store externally the properties of materials, the method
``material.to_json()`` and the class method
``SubclassedMaterial.load_from_json(info)`` should be implemented. The
first one must convert any information stored in ``material.properties``
to a serializable dictionary, and the second one must use this
dictionary to reconstruct the ``material.properties`` dictionary.

The ``VolumeMaterial`` class
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Instances of ``VolumeMaterial`` represent optical properties of physical
objects whose depth is not negligible, like glasses or XXXXXX. In this
case the method ``material.change_of_optical_state()`` is generically
implemented using Snell’s laws of refraction, but any user could
subclass it and implement some other optical behaviour of the material.
Some subclasses of this class are provided, so that materials appearing
usually in the field of solar collectors can be used without further
implementations.

Some of this materials are:

-  ``SimpleVolumeMaterial``, representing a material with constant index
   of refraction and attenuation coefficient.

-  ``WavelengthVolumeMaterial``, where the indices of refraction an
   attenuation depend on the wavelength of the ray, which are given in
   tabulated form.

-  ``PolarizedThinFilm``, XXXXXXXXX

The ``SurfaceMaterial`` class
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Any ``surface_material`` represents a two-dimensional physical object,
in the sense that its third dimension is negligible, or simply that its
behaviour does not depend on it. Examples of these objects are mirrors,
optical collectors, XXXXX. In a first approximation, the interaction of
a ray with such a material can result in a reflection, an absorption or
a transmittance, each with a given probability that may depend on the
wavelength of the ray and are stored in the dictionary
``p=material.properties``. Hence, ``material.change_of_optical_state()``
generically implements these different phenomena. These behaviour is
also affected by other properties of the material, like the booleans:

-  ``p['lambertian_material']``, indicating that, in case of reflection,
   the direction of the reflected ray should be a random vector, instead
   of that computed using the law of reflection.

-  ``p['energy_collector']``, indicating that, in case of an absorption,
   the energy is absorbed and processed, instead of lost in the
   material.

Some more specific materials are provided by subclassing
``VolumeMaterial`` and overriding the ``change_of_optical_state()``
method. Some examples of these specific materials are:

-  ``AbsorberTWModelLayer``, XXXX

-  ``PolarizedCoatingLayer``, and its subclasses
   ``PolarizedCoatingReflectorLayer``,
   ``PolarizedCoatingTransparenLayer``,
   ``PolarizedCoatingAbsorberLayer``,

