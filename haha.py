from firedrake import *

# Import gmesh
mesh = Mesh("trajectory.msh")

V = FunctionSpace(mesh, 'CG', 1)
VV = VectorFunctionSpace(mesh, 'CG', 1, dim = 2)

u = Function(VV)
u_y = Function(V)

t = 0
vtk = File("haha/haha.pvd")
for n in range(11):
     u_star = Constant((2*cos(pi * t - pi/2), 2*sin(pi * t - pi/2) + 1))
     print(u_star, t)
     u.interpolate(u_star)
     u_y.interpolate(u.sub(0))
     vtk.write(u, u_y time = t)
     t = t + 0.1