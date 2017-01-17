import raytrace
import FreeCAD
from FreeCAD import Base

raytrace.create_opaque_simple_material("Mir1", 0.885)
raytrace.create_opaque_simple_material("Abs1", 1-0.841)
raytrace.create_simple_volume_material("Glass1", 1.5)

current_doc = FreeCAD.activeDocument()
sel = current_doc.Objects
current_scene = raytrace.Scene(sel)
exp = raytrace.Experiment(current_scene, Base.Vector(1, 1, 1), 100, 1.0, 1.0, current_doc)
exp.run(current_doc)
