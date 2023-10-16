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

    options = parser.parse_args()
    return options

options = parse()

from fenics import *
from fenics_adjoint import *
from collections import OrderedDict

# Next, we define the expressions for observational data :math:`d` and the
# viscosity :math:`\nu`.

# Import gmesh
xml_mesh = "mesh/leastsquare.xml"
mesh = Mesh(xml_mesh)

subdomains = MeshFunction('size_t', mesh, "mesh/leastsquare_physical_region.xml")
boundaries = MeshFunction('size_t', mesh, "mesh/leastsquare_facet_region.xml")

ds = Measure('ds', domain = mesh, subdomain_data = boundaries)
dx = Measure('dx', domain = mesh, subdomain_data = subdomains)

Id = Identity(mesh.geometric_dimension()) #Identity tensor

V = FunctionSpace(mesh, "CG", 1)
VV = VectorFunctionSpace(mesh, "CG", 1)

# Create initial designs
rhoi = Function(V, name = "Materual density")
rhos = Function(V, name = "Structural material")  # Structural material 1(Blue)
rhor = Function(V, name = "Responsive material")  # Responsive material 2(Red)

# Create initial design
rhos.interpolate(Constant(options.volume_s))
rhor.interpolate(Constant(options.volume_r))
###### End Initial Design #####

# Define the constants and parameters used in the problem
kappa = Constant(options.kappa)
delta = Constant(1.0e-3)
epsilon = Constant(options.epsilon)
kappa_d_e = Constant(kappa / epsilon)
kappa_m_e = Constant(kappa * epsilon)

# Define the boundary/traction force
f = Constant((0.0, -1.0))
#u_star = Constant((0.0, t))

# Young's modulus of the beam and poisson ratio
E_v = Constant(delta)
E_s = Constant(options.esmodulus)
E_r = Constant(options.ermodulus)
nu = Constant(0.3) #nu poisson ratio

mu_v = E_v/(2 * (1 + nu))
lambda_v = (E_v * nu)/((1 + nu) * (1 - 2 * nu))

mu_s = E_s/(2 * (1 + nu))
lambda_s = (E_s * nu)/((1 + nu) * (1 - 2 * nu))

mu_r = E_r/(2 * (1 + nu))
lambda_r = (E_r * nu)/((1 + nu) * (1 - 2 * nu))

def v_v(rhos, rhor):
    return 1 - rhos - rhor

def v_s(rhos):
    return rhos

def v_r(rhor):
    return rhor

def k(rhos, rhor):
    return delta + (1 - delta) * (rhos + rhor)

# Define h_v(rho)=rho_v^(p)
def h_v(rhos, rhor):
    return pow((1 - rhos - rhor), options.power_p)

# Define h_s(rho)=rho_s^(p)
def h_s(rhos):
    return pow(rhos, options.power_p)

# Define h_r(rho)=rho_r^(p)
def h_r(rhor):
    return pow(rhor, options.power_p)

# Define the double-well potential function
# W(x, y) = x^q * (1 - x)^q
def W(x):
    return pow(x, options.power_q) * pow((1 - x), options.power_q)

def WW(x, y):
	return pow((x + y), options.power_q) * pow((1 - x), options.power_q) * pow((1 - y), options.power_q)

# Define strain tensor epsilon(u)
def epsilon(u):
    return 0.5 * (grad(u) + grad(u).T)

# Define the residual stressescc
def sigma_A(A, Id):
	return lambda_r * tr(A) * Id + 2 * mu_r * A

# Define the stress tensor sigma_v(u) for void
def sigma_v(u, Id):
    return lambda_v * tr(epsilon(u)) * Id + 2 * mu_v * epsilon(u)

# Define the stress tensor sigma_s(u) for structural material
def sigma_s(u, Id):
    return lambda_s * tr(epsilon(u)) * Id + 2 * mu_s * epsilon(u)

# Define the stress tensor sigma_r(u) for responsive material
def sigma_r(u, Id):
    return lambda_r * tr(epsilon(u)) * Id + 2 * mu_r * epsilon(u)

dt = Constant(0.1)
T = 1

ctrls = OrderedDict()
t = float(dt)
while t <= T:
    ctrls[t] = Function(V)
    t += float(dt)

# The following function implements a heat equation solver in FEniCS,
# and constructs the first functional term.

def solve_pdes(ctrls):
    s = TrialFunction(V)
    w = TestFunction(V)

    # Define test function and beam displacement
    v = TestFunction(VV)
    u = Function(VV, name = "Displacement")

    g = Function(V, name="source")
    s_0 = Function(V, name="solution")

    F = ( (s - s_0)/dt*w + nu*inner(grad(s), grad(w)) - g*w)*dx
    a, L = lhs(F), rhs(F)
    # The left side of the beam is clamped

    bcs = DirichletBC(VV, Constant((0, 0)), boundaries, 7)
    bc = DirichletBC(V, Constant(0.0), "on_boundary")

    # Define the weak form for forward PDE
    a_forward_v = h_v(rhos, rhor) * inner(sigma_v(u, Id), epsilon(v)) * dx
    a_forward_s = h_s(rhos) * inner(sigma_s(u, Id), epsilon(v)) * dx
    a_forward_r = h_r(rhor) * inner(sigma_r(u, Id), epsilon(v)) * dx
    a_forward = a_forward_v + a_forward_s + a_forward_r

    L_forward = inner(f, v) * ds(8) + s_0 * h_r(rhor) * inner(sigma_A(Id, Id), epsilon(v)) * dx
    R_fwd = a_forward - L_forward

    t = float(dt)
    u_star = Constant((0.0, 0.0))

    j = 0.5*float(dt)*assemble((u - u_star)**2*dx)

    while t <= T:
        # Update source term from control array
        g.assign(ctrls[t])
        u_star = Constant((0.0, t))

        # Solve PDEs
        solve(a == L, s_0, bcs = bc)
        
        solve(R_fwd == 0, u, bcs = bcs)

        j += 0.5*float(dt)*assemble((u - u_star)**2*dx)
        

        # Update time
        t += float(dt)

    return s_0, u, j

u, d, j = solve_pdes(ctrls)

alpha = Constant(1e-1)
regularisation = alpha/2*sum([1/dt*(gb - ga)**2*dx for gb, ga in
    zip(list(ctrls.values())[1:], list(ctrls.values())[:-1])])

# We add the regularisation term to the first functional term and define define the controls:
# Define the Modica-Mortola functional
func1 = kappa_d_e * WW(rhos, rhor) * dx

func2_sub1 = inner(grad(v_v(rhos, rhor)), grad(v_v(rhos, rhor))) * dx
func2_sub2 = inner(grad(v_s(rhos)), grad(v_s(rhos))) * dx
func2_sub3 = inner(grad(v_r(rhor)), grad(v_r(rhor))) * dx
func2 = kappa_m_e * (func2_sub1 + func2_sub2 + func2_sub3)
P = func1 + func2

J = j + assemble(regularisation) + assemble(P)
ms = Control(rhos)
mr = Control(rhor)
m = [Control(c) for c in ctrls.values()].append(ms).append(mr)
print(type(mc))
print(len(mc))

# Finally, we define the reduced functional and solve the optimisation problem:

Jhat = ReducedFunctional(J, [ms, mr, m])

# Define box contraints
lm = 0.0
um = 1.0

boxconstraints = [(lm, um) for i in range(12)]
print(len(boxconstraints))

volumes_constraint = UFLInequalityConstraint((options.volume_s - rhos)*dx, [ms, mr])
volumer_constraint = UFLInequalityConstraint((options.volume_r - rhor)*dx, [ms, mr])

problem = MinimizationProblem(Jhat, bounds = boxconstraints, constraints = [volumes_constraint, volumer_constraint])

parameters = {"acceptable_tol": 1.0e-5, "maximum_iterations": options.maxit}
solver = IPOPTSolver(problem, parameters = parameters)
rhos_opt, rhor_opt, opt_ctrls = solver.solve()

#opt_ctrls = minimize(rf, options={"maxiter": 50})

