from fenics import *
from fenics_adjoint import *

# Define the domain [-1, 1] x [-1, 1]
n = 200
mesh = UnitSquareMesh(n, n)

# Define the function spaces
V = FunctionSpace(mesh, "CG", 1)
W = FunctionSpace(mesh, "DG", 0)

# Define the functions
u = Function(V, name = "Solution")
m = Function(W, name = "Control")
v = TestFunction(V)


# Solve the forward model to create the tape
F = (inner(grad(u), grad(v)) - m*v)*dx

def boundary(x, on_boundary):
     return on_boundary
bc = DirichletBC(V, Constant(0.0), boundary)

solve(F == 0, u, bcs = bc)

# Define the functional of interest
x, y = SpatialCoordinate(mesh)
u_d = exp(-1/(1-x*x)-1/(1-y*y))
ud = project(u_d, V)
J = assemble((0.5*inner(u - u_d, u - u_d))*dx)

# Define the reduced functional
J_hat = ReducedFunctional(J, SteadyParameter(m))
# Solve the optimisation problem
m_opt = minimize(J_hat, method = "L-BFGS-B", bounds = (0.0, 0.5),options = {"gtol": 1e-16, "ftol": 1e-16})
File("problem/u.pvd") << u
File("problem/ud.pvd") << ud
File("problem/m.pvd") << m_opt