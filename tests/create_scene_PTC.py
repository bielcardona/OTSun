# General modules

import Part
from FreeCAD import Base
import math
import numpy as np
import os

# Alias for active Document
doc = App.activeDocument()
# Alias for gui active Document
gui = Gui.activeDocument()

# create a Part object that is a Parabola in the XY plane (the parabola is infinite).
p = Part.Parabola()
# define de Focal distance in the X axe
p.Focal = 647.0
# create a Part object defined by the edge of the parabola with the limits in the Y axe
edge_p = Part.Edge(p, -1845.0/2.0, 1845.0/2.0)
# adds a Part object type to the document and assigns the shape representation of the edge_p
MyShape = doc.addObject("Part::Feature","Parabola")
MyShape.Shape = edge_p

# transformating the edge_p
edge_p.rotate(App.Vector(0.0,0.0,0.0),App.Vector(0.0,1.0,0.0),-90.0)
edge_p.rotate(App.Vector(0.0,0,0.0),App.Vector(0.0,0.0,1.0),+90.0)
edge_p.translate(App.Base.Vector(0.0,0.0,0.0))
MyShape.Shape = edge_p

#create an Extrusion object to extrude the edge_p
Extruded_parabola = doc.addObject("Part::Extrusion","Extruded_parabola")
doc.Extruded_parabola.Label = "Extruded_parabola(Mir1)"
# define the extrusion for the parabola
doc.Extruded_parabola.Base = doc.Parabola
doc.Extruded_parabola.Dir = (0.0,10347.0,0.0)
doc.Extruded_parabola.Solid = (False)
doc.Extruded_parabola.TaperAngle = (0.0)

# to recalculate the whole document
FreeCAD.ActiveDocument.recompute()

# to change the color of the parabola
gui.getObject("Extruded_parabola").ShapeColor = (0.0,0.67,1.00)

# create a generic circle for the absorber
circle_abs = Part.makeCircle(34.0/2.0, Base.Vector(0.0,0.0,647.0), Base.Vector(0.0,1.0,0.0))
#create a FreeCAD object with Part attributes
circle_abs_Part = doc.addObject("Part::Feature","circle_abs_Part")
# asign the circle to circle_abs_Part
circle_abs_Part.Shape = circle_abs

#create a Extrusion object for the extrusion
Abs_circle_abs_Extrude = doc.addObject("Part::Extrusion","Abs_circle_abs_Extrude")
doc.Abs_circle_abs_Extrude.Label = "Abs_circle_abs_Extrude(Abs1)"
# define the extrusion of the parabola
doc.Abs_circle_abs_Extrude.Base = FreeCAD.ActiveDocument.circle_abs_Part
doc.Abs_circle_abs_Extrude.Dir = (0.0,10347.0,0.0)
doc.Abs_circle_abs_Extrude.Solid = (False)
doc.Abs_circle_abs_Extrude.TaperAngle = (0.0)

# to recalculate the whole document
FreeCAD.ActiveDocument.recompute()

# to change the color of the absorber tube
gui.getObject("Abs_circle_abs_Extrude").ShapeColor = (1.0,0.0,0.0)

#create a FreeCAD object with Cylinder attributes
cylinder_out = doc.addObject("Part::Cylinder","cylinder_out")
cylinder_out.Radius = 56.0/2.0
cylinder_out.Height = 10347.0
cylinder_out.Angle = 360.0
cylinder_out.Placement = App.Placement(App.Vector(0.0,0.0,647.0),App.Rotation(App.Vector(1,0.0,0.0),-90.0))

#create a FreeCAD object with Cylinder attributes
cylinder_in = doc.addObject("Part::Cylinder","cylinder_in")
cylinder_in.Radius = 56.0/2.0 - 2.5
cylinder_in.Height = 10347.0 + 2.0
cylinder_in.Angle = 360.0
cylinder_in.Placement = App.Placement(App.Vector(0.0,0.0-1.0,647.0),App.Rotation(App.Vector(1.0,0.0,0.0),-90.0))

#create a FreeCAD object with Part Cut attributes
tube_glass = doc.addObject("Part::Cut","tube_glass")
doc.tube_glass.Label = "tube_glass(Glass1)"
doc.tube_glass.Base = doc.cylinder_out
doc.tube_glass.Tool = doc.cylinder_in
gui.hide("cylinder_out")
gui.hide("cylinder_in")
FreeCAD.ActiveDocument.recompute()

# to make transparent the glass tube
gui.getObject("tube_glass").Transparency = 80

FreeCAD.ActiveDocument.recompute()
Gui.activeDocument().activeView().viewAxonometric()
Gui.SendMsgToActiveView("ViewFit")

# execfile("H:Ramon_2015/RECERCA/RETOS-2015/Tareas/Proves-FreeCAD-1/create_scene_PTC.py")
# execfile("H:Ramon_2015/RECERCA/RETOS-2015/Tareas/Proves-FreeCAD-1/raytrace.py")
