from fenics import *
from fenics_adjoint import *

# Define the domain [-1, 1] x [-1, 1]
n = 200
mesh = RectangleMesh(-1, -1, 1, 1, n, n)

# Define the function spaces
V = FunctionSpace(mesh, "CG", degree = 1)
W = FunctionSpace(mesh, "DG", degree = 0)

# Define the functions
u = Function(V, name = "Solution")
m = Function(W, name = "Control")
v = TestFunction(V)


# Solve the forward model to create the tape
F = (inner(grad(u), grad(v)) - m*v)*dx
bc = DirichletBC(V, 0.0, "on_boundary")
solve(F == 0, u, bc)

# Define the functional of interest
u_d = exp(-1/(1-x[0]*x[0])-1/(1-x[1]*x[1]))
J = Functional((0.5*inner(u - u_d, u - u_d))*dx)

# Define the reduced functional
J_hat = ReducedFunctional(J, Control(m))
# Solve the optimisation problem
m_opt = minimize(J_hat, method = "L-BFGS-B", bounds = (0.0, 0.5),options = {"gtol": 1e-16, "ftol": 1e-16})