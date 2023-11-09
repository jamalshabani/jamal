from firedrake import *
from firedrake_adjoint import *

N = 64
mesh = UnitSquareMesh(N, N)

V = FunctionSpace(mesh, "CG", 1)
W = FunctionSpace(mesh, "R", 0)

u = Function(V, name="state")
v = TestFunction(V)

c = Constant(1)
f = Function(V)
f.assign(1.0)

u = Function(V)
v = TestFunction(V)
bc = DirichletBC(V, Constant(1), "on_boundary")

F = inner(grad(u), grad(v))*dx - f**2*v*dx
solve(F == 0, u, bc)

# J = assemble(c**2*u*dx)
J = Functional(c**2*u*dx)
Jhat = ReducedFunctional(J, Control(c))
Jhat.derivative()