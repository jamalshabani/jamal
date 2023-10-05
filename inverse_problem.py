from fenics import *
from fenics_adjoint import *

# Next, we define the expressions for observational data :math:`d` and the
# viscosity :math:`\nu`.

data = Expression("16*x[0]*(x[0]-1)*x[1]*(x[1]-1)*sin(pi*t)", t=0, degree=4)
nu = Constant(1e-5)

mesh = UnitSquareMesh(8, 8)
V = FunctionSpace(mesh, "CG", 1)

# ... and time:

dt = Constant(0.1)
T = 1

