import raytrace
import FreeCAD
from FreeCAD import Base
import Part

raytrace.create_opaque_simple_material("Mir1", 0.885)
raytrace.create_opaque_simple_material("Abs1", 1 - 0.841)
raytrace.create_simple_volume_material("Glass1", 1.5)

doc = FreeCAD.newDocument()

# We put an sphere of glass
esfera = Part.makeSphere(5, Base.Vector(0, 0, 10), Base.Vector(0, 0, 10), -90, 90, 360)
shapeobj = doc.addObject("Part::Feature", "")
shapeobj.Label = "Esfera(Glass1)"
shapeobj.Shape = esfera

# We put a cube of glass
cub = Part.makeBox(10, 10, 10)
cub.translate(Base.Vector(10, 0, 0))
shapeobj = doc.addObject("Part::Feature", "")
shapeobj.Label = "Cub(Glass1)"
shapeobj.Shape = cub

# We put an absorber on a face of the cube (and paint it red)
cara1 = cub.Faces[1]
shapeobj = doc.addObject("Part::Feature", "")
shapeobj.Label = "Cara1(Abs1)"
shapeobj.Shape = cara1

# We put a mirror on a face of the cube (and paint it blue)
cara5 = cub.Faces[5]
shapeobj = doc.addObject("Part::Feature", "")
shapeobj.Label = "Cara5(Mir1)"
shapeobj.Shape = cara5

doc.recompute()

sel = doc.Objects
current_scene = raytrace.Scene(sel)
exp = raytrace.Experiment(current_scene, Base.Vector(1, 1, 1), 100, 1.0, 1.0, None)
exp.run(None)
print exp.captured_energy
