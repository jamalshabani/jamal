from firedrake import *
import numpy as np

# Import gmesh
mesh = Mesh("trajectory.msh")

V = FunctionSpace(mesh, 'CG', 1)
VV = VectorFunctionSpace(mesh, 'CG', 1, dim = 2)

u = Function(VV, name = "Time displacement")
uu = Function(VV, name = "Displacement")
uuu = Function(VV, name = "Product")
u_x = Function(V)
u_y = Function(V)

bcs = DirichletBC(VV, Constant((0, 0)), 7)

t = 0
vtk = File("haha/haha.pvd")

for n in range(11):
     u_star = Constant((cos(pi * t - pi/2), sin(pi * t - pi/2)))
     u.interpolate(u_star)

     uu_star = Constant((8, 6))
     uu.interpolate(uu_star)

     u_array = u.vector().array().reshape(16068, 2)
     uu_array = uu.vector().array().reshape(16068, 2)

     uuu_array = u_array * uu_array
     print(u_array)
     print(uu_array)
     print("Product = ", uuu_array)


     uuu.vector()[:] = uuu_array

     vtk.write(u, uu,uuu, time = t)
     t = t + 0.1