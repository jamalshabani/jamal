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

bcs = DirichletBC(VV, Constant((0, 0)), 7)

t = 0
vtk = File("haha/haha.pvd")
for n in range(11):
     u_star = Constant((2*cos(pi * t - pi/2), 2*sin(pi * t - pi/2) + 1))

     u.interpolate(u_star).apply(bcs)

     vtk.write(u, time = t)
     t = t + 0.1