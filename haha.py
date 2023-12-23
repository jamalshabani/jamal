from firedrake import *

# Import gmesh
mesh = Mesh("trajectory.msh")

V = FunctionSpace(mesh, 'CG', 1)
VV = VectorFunctionSpace(mesh, 'CG', 1, dim = 2)

u = Function(VV, name = "Time displacement")
uu = Function(VV, name = "Displacement")
uuu = Function(VV, name = "Product")
u_x = Function(V)
u_y = Function(V)

t = 0
vtk = File("haha/haha.pvd")
for n in range(11):
     u_star = Constant((2*cos(pi * t - pi/2), 2*sin(pi * t - pi/2) + 1))
     uu_star = Constant((8,6))

     u.interpolate(u_star)
     uu.interpolate(uu_star)

     uuu = inner(u,uu)
     
     print(type(uuu))

     vtk.write(u, uu, time = t)
     t = t + 0.1