from fenics import *
from fenics_adjoint import *

mesh = UnitSquareMesh(50, 50)
V = FunctionSpace(mesh, "CG", 1)

def u_target(V, t):
     data = Expression("16*x[0]*(x[0]-1)*x[1]*(x[1]-1)*sin(pi*t)", t = t, degree = 4)
     u_tar = Function(V)
     u_tar.interpolate(data, V)
     return u_tar