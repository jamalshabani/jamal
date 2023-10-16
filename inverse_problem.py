def parse():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-vs', '--volume_s', type = float, default = 0.3, help = 'Volume percentage for structural material')
    parser.add_argument('-vr', '--volume_r', type = float, default = 0.3, help = 'Volume percentage for responsive material')
    parser.add_argument('-n', '--maxit', type = int, default = 2000, help = 'Number of iterations')
    parser.add_argument('-k', '--kappa', type = float, default = 5.0e-4, help = 'Weight of Modica-Mortola')
    parser.add_argument('-e', '--epsilon', type = float, default = 4.0e-3, help = 'Phase-field regularization parameter')
    parser.add_argument('-o', '--output', type = str, default = 'output', help = 'Output folder')
    parser.add_argument('-es', '--esmodulus', type = float, default = 1.0e1, help = 'Elastic Modulus for structural material')
    parser.add_argument('-er', '--ermodulus', type = float, default = 1.0e-1, help = 'Elastic Modulus for responsive material')
    parser.add_argument('-p', '--power_p', type = float, default = 2.0, help = 'Power for elasticity interpolation')
    parser.add_argument('-q', '--power_q', type = float, default = 2.0, help = 'Power for multiple-well function')
    parser.add_argument('-s', '--steamy', type = float, default = 1.0, help = 'Initial stimulus')

    options = parser.parse_args()
    return options

options = parse()

from fenics import *
from fenics_adjoint import *
from collections import OrderedDict

# Next, we define the expressions for observational data :math:`d` and the
# viscosity :math:`\nu`.

data = Expression("16*x[0]*(x[0]-1)*x[1]*(x[1]-1)*sin(pi*t)", t=0, degree=4)
nu = Constant(1e-5)

mesh = UnitSquareMesh(8, 8)
V = FunctionSpace(mesh, "CG", 1)

# ... and time:

dt = Constant(0.1)
T = 2

# We are considering a time-distributed forcing as control. In the next step,
# we create one control function for each timestep in the model, and store all
# controls in a dictionary that maps timestep to control function:

ctrls = OrderedDict()
u_sol = OrderedDict()
u_tar = OrderedDict()

vtkfile = File("problem/u_tar.pvd")
t = float(dt)
n = 1
while n <= T*10:
    data.t = t
    ctrls[n] = Function(V)
    u_sol[n] = Function(V)
    u_tar[n] = interpolate(data, V)
    t += float(dt)
    n = n + 1

for key, value in u_tar.items():
    print(key, value)

for n in range(20):
    vtkfile.write(u_tar[n+1], n+1)



# The following function implements a heat equation solver in FEniCS,
# and constructs the first functional term.

def solve_heat(ctrls):
    u = TrialFunction(V)
    v = TestFunction(V)

    f = Function(V, name="source")
    u_0 = Function(V, name="solution")
    uu = Function(V, name="solution")
    d = Function(V, name="data")

    F = ( (u - u_0)/dt*v + nu*inner(grad(u), grad(v)) - f*v)*dx
    a, L = lhs(F), rhs(F)
    bc = DirichletBC(V, 0, "on_boundary")

    t = float(dt)

    j = 0.5*float(dt)*assemble((u_0 - d)**2*dx)

    while t <= T:
        # Update source term from control array
        f.assign(ctrls[t])

        # Update data function
        data.t = t
        d.assign(interpolate(data, V))

        # Solve PDE
        solve(a == L, u_0, bc)
        solve(a == L, uu, bc)
        File('problem/solution.pvd') << (uu, t)
        File("problem/u_tar.pvd") << (u_tar[t], t)
        u_0.assign(uu)

        # Implement a trapezoidal rule 
        if t > T - float(dt):
           weight = 0.5
        else:
           weight = 1

        j += weight*float(dt)*assemble((u_0 - d)**2*dx)

        # Update time
        t += float(dt)

    return u_0, d, j

u, d, j = solve_heat(ctrls)

# With this preparation steps, we are now ready to define the functional.
# First we discretise the regularisation term
#
# .. math::
#             \frac{\alpha}{2} \int_0^T \int_\Omega \dot f^2 \textrm{d} \Omega \text{d}t
#
# Note, that :math:`f` is a piecewise linear function in time over the time intervals :math:`K = [(0, \delta t), (\delta t, 2 \delta t), \dots, (T-\delta
# t, T)]`. Thus, we can write the integral as a sum over all intervals
#
# .. math::
#             \frac{\alpha}{2} \sum_{a_k, b_k \in K} \int_{a_k}^{b_k} \int_\Omega \dot f(t)^2 \textrm{d} \Omega\text{d}t
#
# Discretising the time-derivative yields:
#
# .. math::
#             \frac{\alpha}{2} \sum_K \int_{a_k}^{b_k}
#             \int_\Omega \left(\frac{f(b_k)-
#             f(a_k)}{b_k-a_k}\right)^2\textrm{d}\Omega \\
#             = \frac{\alpha}{2} \sum_K (b_k-a_k)^{-1}
#             \int_\Omega \left(f(b_k)- f(a_k)\right)^2\textrm{d}\Omega
#
#
# In code this is translates to:

alpha = Constant(1e-1)
regularisation = alpha/2*sum([1/dt*(fb-fa)**2*dx for fb, fa in
    zip(list(ctrls.values())[1:], list(ctrls.values())[:-1])])

# We add the regularisation term to the first functional term and define define the controls:

J = j + assemble(regularisation)
m = [Control(c) for c in ctrls.values()]

# Finally, we define the reduced functional and solve the optimisation problem:

rf = ReducedFunctional(J, m)
opt_ctrls = minimize(rf, options={"maxiter": 50})

for i in range(len(opt_ctrls)):
    File("problem/controls.pvd") << opt_ctrls[i]
print(type(opt_ctrls[10]))


# If we solve this optimisation problem with varying :math:`\alpha` parameters,
# we observe that we get different behaviour in the controls: the higher the
# alpha value, the "smoother" the control function becomes. The following plots
# show the optimised control evaluated at the middle point :math:`(0.5, 0.5)`
# over time for different :math:`\alpha` values:

# .. image:: control_alpha=0.0001.png
#     :scale: 45
#     :align: left
# .. image:: control_alpha=0.001.png
#     :scale: 45
#     :align: right
# .. image:: control_alpha=0.01.png
#     :scale: 45
#     :align: left
# .. image:: control_alpha=0.1.png
#     :scale: 45
#     :align: right
