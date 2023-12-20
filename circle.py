def parse():
	import argparse
	parser = argparse.ArgumentParser()
	parser.add_argument('-tao_type', '--tao_type', type = str, default = 'bncg', help = 'TAO algorithm type')
	parser.add_argument('-tao_max_funcs', '--tao_max_funcs', type = int, default = 10000, help = 'TAO maximum functions evaluations')
	parser.add_argument('-tao_monitor', '--tao_monitor', action = 'store_true', help = 'TAO monitor')
	parser.add_argument('-tao_ls_monitor', '--tao_ls_monitor', action = 'store_true', help = 'TAO line search monitor')
	parser.add_argument('-ls', '--lagrange_s', type = float, default = 0.02, help = 'Lagrange multiplier for structural material')
	parser.add_argument('-lr', '--lagrange_r', type = float, default = 4.5, help = 'Lagrange multiplier for responsive material')
	parser.add_argument('-tao_ls_type', '--tao_ls_type', type = str, default = 'more-thuente', help = "TAO line search")
	parser.add_argument('-tao_view', '--tao_view', action = 'store_true', help = "View convergence details")
	parser.add_argument('-tao_max_it', '--tao_max_it', type = int, default = 1000, help = 'Number of TAO iterations')
	parser.add_argument('-vs', '--volume_s', type = float, default = 0.3, help = 'Volume percentage for structural material')
	parser.add_argument('-vr', '--volume_r', type = float, default = 0.3, help = 'Volume percentage for responsive material')
	parser.add_argument('-k', '--kappa', type = float, default = 1.0e-3, help = 'Weight of Modica-Mortola')
	parser.add_argument('-e', '--epsilon', type = float, default = 3.0e-3, help = 'Phase-field regularization parameter')
	parser.add_argument('-o', '--output', type = str, default = 'test1', help = 'Output folder')
	parser.add_argument('-m', '--mesh', type = str, default = 'trajectory.msh', help = 'Dimensions of meshed beam')
	parser.add_argument('-es', '--esmodulus', type = float, default = 1.0e-1, help = 'Elastic Modulus for structural material')
	parser.add_argument('-er', '--ermodulus', type = float, default = 1.0e1, help = 'Elastic Modulus for responsive material')
	parser.add_argument('-p', '--power_p', type = float, default = 2.0, help = 'Power for elasticity interpolation')
	parser.add_argument('-q', '--power_q', type = float, default = 1.0, help = 'Power for multiple-well function')
	parser.add_argument('-s', '--steamy', type = float, default = 1.0, help = 'Initial stimulus')
	parser.add_argument('-g', '--heatsource', type = float, default = 0.0, help = 'Heat source')
	options = parser.parse_args()
	return options

options = parse()

from firedrake import *
from petsc4py import PETSc
import time
import numpy as np

start = time.time()

# Import gmesh
mesh = Mesh(options.mesh)
Id = Identity(mesh.geometric_dimension()) #Identity tensor

# Define the function spaces
V = FunctionSpace(mesh, 'CG', 1)
VV = VectorFunctionSpace(mesh, 'CG', 1, dim = 2)
VVV = VectorFunctionSpace(mesh, 'CG', 1, dim = 3)

# Create initial design
###### Begin Initial Design #####
mesh_coordinates = mesh.coordinates.dat.data[:]
M = len(mesh_coordinates)

rho =  Function(VVV, name = "Design variable")
rho_i = Function(V, name = "Material density")
rhos = Function(V, name = "Structural material")  # Structural material 1(Blue)
rhor = Function(V, name = "Responsive material")  # Responsive material 2(Red)
g = Function(V, name = "Heat source")

x, y = SpatialCoordinate(mesh)
rhos.interpolate(Constant(options.volume_s))
rhos.interpolate(Constant(1.0), mesh.measure_set("cell", 4))

rhor.interpolate(Constant(options.volume_r))
rhor.interpolate(Constant(0.0), mesh.measure_set("cell", 4))

g.interpolate(Constant(options.heatsource))

rho = as_vector([rhos, rhor, g])
rho = interpolate(rho, VVV)
###### End Initial Design #####

# Define the constant parameter used in the problem
kappa = Constant(options.kappa)
print(kappa)
lagrange_r = Constant(options.lagrange_r)
lagrange_s = Constant(options.lagrange_s)

# Total volume of the domain |omega|
omega = assemble(interpolate(Constant(1.0), V) * dx)

delta = Constant(1.0e-3)
epsilon = Constant(options.epsilon)
kappa_d_e = Constant(kappa / epsilon)
kappa_m_e = Constant(kappa * epsilon)

kappa_d_ge = Constant(kappa / epsilon / 10)
kappa_m_ge = Constant(kappa * epsilon / 10)

# Define the traction force and predescribed displacement
def u_starx(t):
     return 2 * cos(pi * t - pi/2)

def u_stary(t):
     return 2 * sin(pi * t - pi/2) + 1

def u_staryff(t):
	if (0 <= t < 0.25):
		return 4 * t, -1
		
	if (0.25 <= t < 0.75):
		return 2 - 4 * t, 1
		
	if (0.75 <= t <= 1):
		return 4 * t - 4, -1

# f = Constant((0, -1.0))

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

def v_v(rho):
	return 1 - rho.sub(0) - rho.sub(1)

def v_s(rho):
	return rho.sub(0)

def v_r(rho):
	return rho.sub(1)

def v_g(rho):
	return rho.sub(2)

# Define h_v(rho)=rho_v^(p)
def h_v(rho):
	return pow((1 - rho.sub(0) - rho.sub(1)), options.power_p)

# Define h_s(rho)=rho_s^(p)
def h_s(rho):
	return pow(rho.sub(0), options.power_p)

# Define h_r(rho)=rho_r^(p)
def h_r(rho):
	return pow(rho.sub(1), options.power_p)

# Define l(rho)
def k(rho):
    return delta + (1 - delta) * (rho.sub(0) + rho.sub(1))
    # return 1


def g(rho):
	return pow(rho.sub(2), options.power_q)

def WW(rho):
	return pow(rho.sub(2), options.power_p) * pow((1 - rho.sub(2)), options.power_p)

# Define the double-well potential function
# W(x, y) = (x + y)^q * (1 - x)^q * (1 - y)^q
def W(rho):
	return pow((rho.sub(0) + rho.sub(1)), options.power_q) * pow((1 - rho.sub(0)), options.power_q) * pow((1 - rho.sub(1)), options.power_q)

# Define strain tensor epsilon(u)
def epsilon(u):
	return 0.5 * (grad(u) + grad(u).T)

# Define the residual stresses
def sigma_A(A, Id):
	return lambda_r * tr(A) * Id + 2 * mu_r * A

# Define the stress tensor sigma_v(u) for void
def sigma_v(u, Id):
	return lambda_v * div(u) * Id + 2 * mu_v * epsilon(u)

# Define the stress tensor sigma_s(u) for structural material
def sigma_s(u, Id):
	return lambda_s * div(u) * Id + 2 * mu_s * epsilon(u)

# Define the stress tensor sigma_r(u) for responsive material
def sigma_r(u, Id):
	return lambda_r * div(u) * Id + 2 * mu_r * epsilon(u)

# Define test function and beam displacement
v = TestFunction(VV)
u = Function(VV, name = "Displacement")
us = Function(VV, name = "Displacement")
p = Function(VV, name = "Adjoint variable")

w = TestFunction(V)
s = Function(V, name = "Stimulus")
q = Function(V, name = "Adjoint variable heat")

# The left side of the beam is clamped
bcs = DirichletBC(VV, Constant((0, 0)), 7)
bcss = DirichletBC(V, Constant(0), 7)

T = 1.0            # Final time
num_steps = 20     # Number of time steps
dt = T / num_steps # Time step size

# Define the objective function
# J = 0.5 * inner(u - u_star, u - u_star) * dx(4)

func1 = kappa_d_e * W(rho) * dx
func2_sub1 = inner(grad(v_v(rho)), grad(v_v(rho))) * dx
func2_sub2 = inner(grad(v_s(rho)), grad(v_s(rho))) * dx
func2_sub3 = inner(grad(v_r(rho)), grad(v_r(rho))) * dx
func2 = kappa_m_e * (func2_sub1 + func2_sub2 + func2_sub3)
P = func1 + func2

func3 = kappa_d_ge * WW(rho) * dx
func4 = kappa_m_ge * inner(grad(v_g(rho)), grad(v_g(rho))) * dx
PP = func3 + func4

func5 = lagrange_s * v_s(rho) * dx
func6 = lagrange_r * v_r(rho) * dx

func7 = pow(v_v(rho), 2) * pow(g(rho), 2) * dx
func8 = pow(v_s(rho), 2) * pow(g(rho), 2) * dx

# Objective function + Modica-Mortola functional + Volume penalties
# JJ = J + P + PP +  func5 + func6 + func7 + func8

# Define weak form for the heat conduction
# Solve for "s"
# Define initial value for stimulus
s_0 = interpolate(Constant(0.0), V)
R_heat_forward = s * w * dx + dt * k(rho) * inner(grad(s), grad(w)) * dx - (s_0 + dt * g(rho)) * w * dx

# Define adjoint weak form for the heat conduction
# Solve for "q"
# Define final value for stimulus
q_n = interpolate(Constant(0.0), V)
a_heat_adjoint = -(q - q_n)/dt * w * dx + k(rho) * inner(grad(q), grad(w)) * dx
L_heat_adjoint = w * h_r(rho) * inner(Id, epsilon(p)) * dx
R_heat_adjoint = a_heat_adjoint - L_heat_adjoint

a_heat_adjoint = k(rho) * inner(grad(w), grad(q)) * dx
L_heat_adjoint = w * h_r(rho) * inner(Id, epsilon(p)) * dx
R_heat_adjoint = a_heat_adjoint - L_heat_adjoint

# Define the weak form for forward PDE
# Solve for "u"
a_forward_v = h_v(rho) * inner(sigma_v(u, Id), epsilon(v)) * dx
a_forward_s = h_s(rho) * inner(sigma_s(u, Id), epsilon(v)) * dx
a_forward_r = h_r(rho) * inner(sigma_r(u, Id), epsilon(v)) * dx
a_forward = a_forward_v + a_forward_s + a_forward_r

L_forward_s = s * h_r(rho) * inner(Id, epsilon(v)) * dx
R_fwd_s = a_forward - L_forward_s


# Define the weak form for adjoint PDE
# Solve for "p"
a_adjoint_v = h_v(rho) * inner(sigma_v(v, Id), epsilon(p)) * dx
a_adjoint_s = h_s(rho) * inner(sigma_s(v, Id), epsilon(p)) * dx
a_adjoint_r = h_r(rho) * inner(sigma_r(v, Id), epsilon(p)) * dx
a_adjoint = a_adjoint_v + a_adjoint_s + a_adjoint_r


# Define the Lagrangian
a_lagrange_v = h_v(rho) * inner(sigma_v(u, Id), epsilon(p)) * dx
a_lagrange_s = h_s(rho) * inner(sigma_s(u, Id), epsilon(p)) * dx
a_lagrange_r = h_r(rho) * inner(sigma_r(u, Id), epsilon(p)) * dx
a_lagrange   = a_lagrange_v + a_lagrange_s + a_lagrange_r


R_heat_lagrange = s * q * dx + dt * inner(grad(s), grad(q)) * dx - (s_0 + dt * g(rho)) * q * dx
a_heat_lagrange = k(rho) * inner(grad(s), grad(q)) * dx
L_heat_lagrange = inner(g(rho), q) * dx
R_heat_lagrange = a_heat_lagrange - L_heat_lagrange

#L = JJ - R_lagrange - R_heat_lagrange


# Beam .pvd file for saving designs
# beam = File(options.output + '/beam.pvd')

dJdrhos = Function(V, name = "Gradient s")
dJdrhor = Function(V, name = "Gradient r")
dJdrhog = Function(V, name = "Gradient g")

rho_res = Function(V, name = "Responsive")
rho_str = Function(V, name = "Structural")
rho_g = Function(V, name = "Heat source")
ave_u = Function(V, name = "Average displacement")

stimulus = Function(V, name = "Stimulus")

N = M * 3
indexs = []
indexr = []
indexg = []

for i in range(N):
	if (i%3) == 0:
		indexs.append(i)
	if (i%3) == 1:
		indexr.append(i)
	if (i%3) == 2:
		indexg.append(i)

u_array = np.empty([num_steps], dtype=object)
s_array = np.empty([num_steps], dtype=object)

def FormObjectiveGradient(tao, x, G):
	Obj = 0

	# Print volume fraction of structural material
	volume_s = assemble(v_s(rho) * dx)/omega
	print("The volume fraction(Vs) is {}".format(volume_s))

	# Print volume fraction of responsive material
	volume_r = assemble(v_r(rho) * dx)/omega
	print("The volume fraction(Vr) is {}".format(volume_r))
	print(" ")

	i = tao.getIterationNumber()
	t = 0

	if (i % 5) == 0:
		rho_i.interpolate(rho.sub(1) - rho.sub(0))
		stimulus.interpolate(s)
		rho_str.interpolate(rho.sub(0))
		rho_res.interpolate(rho.sub(1))

		beam = File(options.output + '/iteration-{}/beam.pvd'.format(i))
		vtkfile = File(options.output + '/iteration-{}/ustimulus.pvd'.format(i))

		for n in range(num_steps):
			t += dt

			rho_g.interpolate(sin(2 * pi * t) * rho.sub(2))

			R_heat_forward2 = s * w * dx + dt * k(rho) * inner(grad(s), grad(w)) * dx - (s_0 + dt * rho_g) * w * dx
			solve(R_heat_forward2 == 0, s, bcs = bcss)
			s_0.assign(s)


			# Step 2: Solve forward PDE
			solve(R_fwd_s == 0, u, bcs = bcs)
			# ave_u.interpolate(assemble(u * dx(4)))
			beam.write(rho_i, rho_str, rho_res, rho_g, s, u, time = t)
			vtkfile.write(s, u, ave_u, time = t)
	
	t = 0
	for n in range(num_steps):
		# Update time
		m = n / 10 +  dt
		# print(t)

		# if (i%5) == 0:

		# 	rho_i.interpolate(rho.sub(1) - rho.sub(0))
		# 	stimulus.interpolate(s)
		# 	rho_str.interpolate(rho.sub(0))
		# 	rho_res.interpolate(rho.sub(1))
		# 	rho_g.interpolate(rho.sub(2))
		# 	solve(R_fwd_s == 0, u, bcs = bcs)
		# 	beam.write(rho_i, stimulus, rho_str, rho_res, rho_g, u)

		with rho.dat.vec as rho_vec:
			rho_vec.set(0.0)
			rho_vec.axpy(1.0, x)

		u_stary, ff = u_staryff(t)

		u_star = Constant((0, u_stary))

		f = Constant((0, ff))

		# print(t, u_star, f)


		L_adjoint = inner(u - u_star, v) * dx(4)
		R_adj = a_adjoint - L_adjoint

		L_forward = inner(f, v) * ds(8) + h_r(rho) * inner(s*Id, epsilon(v)) * dx
		R_fwd = a_forward - L_forward

		L_lagrange = inner(f, p) * ds(8) + h_r(rho) * inner(s*Id, epsilon(p)) * dx
		R_lagrange = a_lagrange - L_lagrange

		# Step 1: Solve heat conduction
		solve(R_heat_forward == 0, s, bcs = bcss)
		s_0.assign(s)

		# Step 2: Solve forward PDE
		solve(R_fwd == 0, u, bcs = bcs)

		# Step 3: Solve adjoint PDE
		solve(R_adj == 0, p, bcs = bcs)

		# Step 4: Solve asjoint Heat
		solve(R_heat_adjoint == 0, q, bcs = bcss)
		q_n.assign(q)

		Obj = Obj + 0.5 * float(dt) * inner(u - u_star, u - u_star) * dx(4)

		t += dt

		# Evaluate the objective function
		# objective_value = assemble(J)
		# print("The value of objective function is {}".format(objective_value))
	JJ = Obj + P + PP +  func5 + func6 + func7 + func8
	L = JJ - R_lagrange - R_heat_lagrange

	# Compute gradiet w.r.t rhos and rhor and s
	dJdrhos.interpolate(10 * assemble(derivative(L, rho.sub(0))).riesz_representation(riesz_map="l2"))
	dJdrhos.interpolate(Constant(0.0), mesh.measure_set("cell", 4))

	dJdrhor.interpolate(10 *  assemble(derivative(L, rho.sub(1))).riesz_representation(riesz_map="l2"))
	dJdrhor.interpolate(Constant(0.0), mesh.measure_set("cell", 4))

	dJdrhog.interpolate(10 * assemble(derivative(L, rho.sub(2))).riesz_representation(riesz_map="l2"))

	G.setValues(indexs, dJdrhos.vector().array())
	G.setValues(indexr, dJdrhor.vector().array())
	G.setValues(indexg, dJdrhog.vector().array())

	f_val = assemble(L)
	return f_val

# Setting lower and upper bounds
lb = as_vector((0, 0, -1))
ub = as_vector((1, 1, 1))
lb = interpolate(lb, VVV)
ub = interpolate(ub, VVV)

with lb.dat.vec as lb_vec:
	rho_lb = lb_vec

with ub.dat.vec as ub_vec:
	rho_ub = ub_vec

# Setting TAO solver
tao = PETSc.TAO().create(PETSc.COMM_SELF)
tao.setType('bncg')
tao.setObjectiveGradient(FormObjectiveGradient, None)
tao.setVariableBounds(rho_lb, rho_ub)
tao.setFromOptions()

# Initial design guess
with rho.dat.vec as rho_vec:
	x = rho_vec.copy()

# Solve the optimization problem
tao.solve(x)
tao.destroy()

end = time.time()
print("\nExecution time (in seconds):", (end - start))
