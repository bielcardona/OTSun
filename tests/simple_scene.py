import Part
from FreeCAD import Base

doc = App.activeDocument()

esfera = Part.makeSphere(5, Base.Vector(0,0,10),Base.Vector(0,0,10),-90,90,360)
shapeobj = doc.addObject("Part::Feature","")
shapeobj.Label = "Esfera(Glass1)"
shapeobj.Shape = esfera

cub = Part.makeBox(10,10,10)
cub.translate(Base.Vector(10,0,0))
shapeobj = doc.addObject("Part::Feature","")
shapeobj.Label = "Cub(Glass1)"
shapeobj.Shape = cub

cara1 = cub.Faces[1]
shapeobj = doc.addObject("Part::Feature","")
shapeobj.Label = "Cara1(Abs1)"
shapeobj.Shape = cara1

cara5 = cub.Faces[5]
shapeobj = doc.addObject("Part::Feature","")
shapeobj.Label = "Cara5(Mir1)"
shapeobj.Shape = cara5


FreeCAD.ActiveDocument.recompute()
Gui.activeDocument().activeView().viewAxonometric()
Gui.SendMsgToActiveView("ViewFit")

