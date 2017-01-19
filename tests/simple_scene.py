import Part
from FreeCAD import Base

doc = App.activeDocument()

# We put an sphere of glass
esfera = Part.makeSphere(5, Base.Vector(0,0,10),Base.Vector(0,0,10),-90,90,360)
shapeobj = doc.addObject("Part::Feature","")
shapeobj.Label = "Esfera(Glass1)"
shapeobj.Shape = esfera
viewobj = shapeobj.ViewObject
viewobj.Transparency = 80

# We put a cube of glass
cub = Part.makeBox(10,10,10)
cub.translate(Base.Vector(10,0,0))
shapeobj = doc.addObject("Part::Feature","")
shapeobj.Label = "Cub(Glass1)"
shapeobj.Shape = cub
viewobj = shapeobj.ViewObject
viewobj.Transparency = 80

# We put an absorber on a face of the cube (and paint it red)
cara1 = cub.Faces[1]
shapeobj = doc.addObject("Part::Feature","")
shapeobj.Label = "Cara1(Abs1)"
shapeobj.Shape = cara1
viewobj = shapeobj.ViewObject
viewobj.ShapeColor = (0.8, 0.0, 0.0, 0.0)

# We put a mirror on a face of the cube (and paint it blue)
cara5 = cub.Faces[5]
shapeobj = doc.addObject("Part::Feature","")
shapeobj.Label = "Cara5(Mir1)"
shapeobj.Shape = cara5
viewobj = shapeobj.ViewObject
viewobj.ShapeColor = (0.0, 0.0, 0.8, 0.0)


FreeCAD.ActiveDocument.recompute()
Gui.activeDocument().activeView().viewAxonometric()
Gui.SendMsgToActiveView("ViewFit")

