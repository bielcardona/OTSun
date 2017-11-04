import sys
sys.path.append("/usr/lib/freecad")
sys.path.append("/usr/lib/freecad/lib")
sys.path.append("/vagrant")
# General modules

import FreeCAD
import Part
from FreeCAD import Base
import math
import numpy as np
import os
import raytrace
import time
import multiprocessing

# Alias for active Document
doc = FreeCAD.newDocument()
# Alias for gui active Document
# gui = Gui.activeDocument()

# adds a Part Plane object type to the document and assigns the parameters
doc.addObject("Part::Plane", "Plane_1")
doc.Plane_1.Length = 100.00
doc.Plane_1.Width = 100.00
doc.Plane_1.Placement = Base.Placement(Base.Vector(-50.00, -50.00, -0.10), Base.Rotation(0.00, 0.00, 0.00, 1.00))
doc.Plane_1.Label = 'Plane_1(Abs1)'
# to change the color of Plane_1
# gui.getObject("Plane_1").ShapeColor = (1.0,0.0,0.0)

L = 4.0

# adds a Part Plane object type to the document and assigns the parameters
doc.addObject("Part::Box", "Box_1")
doc.Box_1.Length = 100.00 * L
doc.Box_1.Width = 100.00 * 1
doc.Box_1.Height = 4.00
doc.Box_1.Placement = Base.Placement(Base.Vector(-50.00 * L, -50.00 * 1, 0.0), App.Rotation(App.Vector(0, 0, 1), 0))
doc.Box_1.Label = 'Box_1(Glass1)'
# to make transparent the glass box
# gui.getObject("Box_1").Transparency = 80

# adds a Part Face object to the document and assigns the parameters based on polygon geometry
p1 = Base.Vector(50.0 * L, -50.0 * 1, 0.0)
p2 = Base.Vector(50.0 * L, 50.0 * 1, 0.0)
p3 = Base.Vector(50.0 * L, 50.0 * 1, 4.0)
p4 = Base.Vector(50.0 * L, -50.0 * 1, 4.0)
wire = Part.makePolygon([p1, p2, p3, p4], True)
face = Part.Face(wire)
face_Shape = doc.addObject("Part::Face", "Plane_2")
face_Shape.Shape = face
doc.Plane_2.Label = 'Plane_1(Opa1)'
# gui.getObject("Plane_2").ShapeColor = (1.0,1.0,1.0)

FreeCAD.ActiveDocument.recompute()
# Gui.activeDocument().activeView().viewAxonometric()
# Gui.SendMsgToActiveView("ViewFit")

raytrace.create_opaque_simple_material("opa1", 0.0)
raytrace.create_absorber_lambertian_layer('lamb1', 0.0)
raytrace.create_two_layers_material("Abs1", "lamb1", "opa1")
raytrace.create_two_layers_material("Opa1", "opa1", "opa1")
# raytrace.create_simple_volume_material("Glass1", 1.5)
file_BK7 = './BK7.txt'
raytrace.create_wavelength_volume_material("Glass1", file_BK7)
# file_AR = 'D:Ramon_2015/RECERCA/RETOS-2015/Tareas/Proves-FreeCAD-2/materials-PV/AR.txt'
# raytrace.create_polarized_coating_transparent_layer("file_AR", file_AR, file_BK7)
# raytrace.create_two_layers_material("AR1","file_AR","file_AR")

sel = doc.Objects
current_scene = raytrace.Scene(sel)
phi_ini = 0
phi_end = 0.3
phi_end = phi_end + 0.00001
phi_step = 0.1
theta_ini = 0
theta_end = 0
theta_end = theta_end + 0.00001
theta_step = 2
number_photons = 800
aperture_collector = 100. * 100.
data_file_spectrum = './ASTMG173-direct.txt'
light_spectrum = raytrace.create_CDF_from_PDF(data_file_spectrum)
outfile = open('./kk.txt', 'w')


def make_run(ph, th):
    main_direction = raytrace.polar_to_cartesian(ph, th) * -1.0  # Sun direction vector
    emitting_region = raytrace.SunWindow(current_scene, main_direction)
    l_s = raytrace.LightSource(current_scene, emitting_region, light_spectrum, 1.0)
    exp = raytrace.Experiment(current_scene, l_s, number_photons)
    exp.run()
    efficiency = (exp.captured_energy / aperture_collector) / (
        exp.number_of_rays / exp.light_source.emitting_region.aperture)
    t1 = time.time()
    print ("%s %s %s %s" % (ph, th, efficiency, t1 - t0))
    outfile.write("%s %s %s %s\n" % (ph, th, efficiency, t1 - t0))


t0 = time.time()

print "Sequential run"

for ph in np.arange(phi_ini, phi_end, phi_step):
    for th in np.arange(theta_ini, theta_end, theta_step):
        # p = multiprocessing.Process(target=make_run, args=(ph,th,))
        # jobs.append(p)
        # p.start()
        make_run(ph, th)

print "Time: %s" % (time.time() - t0)

t0 = time.time()

print "Parallel run"

print "CPUs %d" % multiprocessing.cpu_count()
jobs = []
for ph in np.arange(phi_ini, phi_end, phi_step):
    for th in np.arange(theta_ini, theta_end, theta_step):
        p = multiprocessing.Process(target=make_run, args=(ph, th,))
        jobs.append(p)
        p.start()
        # make_run(ph, th)
for p in jobs:
    p.join()
print "Time: %s" % (time.time() - t0)

outfile.close()
