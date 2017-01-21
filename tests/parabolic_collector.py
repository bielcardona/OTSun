# General modules

import FreeCAD
import Part
from FreeCAD import Base
import raytrace

doc = FreeCAD.newDocument()

p = Part.Parabola()
p.Focal = 647.0
edge_p = Part.Edge(p, -1845.0 / 2.0, 1845.0 / 2.0)
edge_p.rotate(Base.Vector(0.0, 0.0, 0.0), Base.Vector(0.0, 1.0, 0.0), -90.0)
edge_p.rotate(Base.Vector(0.0, 0, 0.0), Base.Vector(0.0, 0.0, 1.0), +90.0)
edge_p.translate(Base.Vector(0.0, 0.0, 0.0))
my_shape = doc.addObject("Part::Feature", "Parabola")
my_shape.Shape = edge_p

extruded_parabola = doc.addObject("Part::Extrusion", "Extruded_parabola")
extruded_parabola.Label = "Extruded_parabola(Mir1)"
extruded_parabola.Base = doc.Parabola
extruded_parabola.Dir = (0.0, 10347.0, 0.0)
extruded_parabola.Solid = (False)
extruded_parabola.TaperAngle = (0.0)

circle_abs = Part.makeCircle(34.0 / 2.0, Base.Vector(0.0, 0.0, 647.0), Base.Vector(0.0, 1.0, 0.0))
circle_abs_Part = doc.addObject("Part::Feature", "circle_abs_Part")
circle_abs_Part.Shape = circle_abs

abs_circle_abs_Extrude = doc.addObject("Part::Extrusion", "Abs_circle_abs_Extrude")
abs_circle_abs_Extrude.Label = "Abs_circle_abs_Extrude(Abs1)"
# define the extrusion of the parabola
abs_circle_abs_Extrude.Base = FreeCAD.ActiveDocument.circle_abs_Part
abs_circle_abs_Extrude.Dir = (0.0, 10347.0, 0.0)
abs_circle_abs_Extrude.Solid = (False)
abs_circle_abs_Extrude.TaperAngle = (0.0)

cylinder_out = doc.addObject("Part::Cylinder", "cylinder_out")
cylinder_out.Radius = 56.0 / 2.0
cylinder_out.Height = 10347.0
cylinder_out.Angle = 360.0
cylinder_out.Placement = FreeCAD.Placement(Base.Vector(0.0, 0.0, 647.0),
                                           FreeCAD.Rotation(Base.Vector(1, 0.0, 0.0), -90.0))

cylinder_in = doc.addObject("Part::Cylinder", "cylinder_in")
cylinder_in.Radius = 56.0 / 2.0 - 2.5
cylinder_in.Height = 10347.0 + 2.0
cylinder_in.Angle = 360.0
cylinder_in.Placement = FreeCAD.Placement(Base.Vector(0.0, 0.0 - 1.0, 647.0),
                                          FreeCAD.Rotation(Base.Vector(1.0, 0.0, 0.0), -90.0))

tube_glass = doc.addObject("Part::Cut", "tube_glass")
tube_glass.Label = "tube_glass(Glass1)"
tube_glass.Base = doc.cylinder_out
tube_glass.Tool = doc.cylinder_in

doc.recompute()

raytrace.create_opaque_simple_material("Mir1", 0.885)
raytrace.create_opaque_simple_material("Abs1", 1 - 0.841)
raytrace.create_simple_volume_material("Glass1", 1.5)

# current_doc = FreeCAD.activeDocument()
sel = doc.Objects
current_scene = raytrace.Scene(sel)
exp = raytrace.Experiment(current_scene, Base.Vector(0, 0, -1), 100, 1.0, 1.0, None)
exp.run(None)
print exp.captured_energy
