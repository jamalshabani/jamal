from firedrake import *

T = 2.0            # final time
num_steps = 50     # number of time steps
dt = T / num_steps # time step size

# Create mesh and define function space
nx = ny = 30
mesh = UnitSquareMesh(nx, ny)
V = FunctionSpace(mesh, 'CG', 1)

# Define boundary condition
bc = DirichletBC(V, Constant(0), [1,2,3,4])

# Define initial value
u_n = interpolate(Constant(1.0), V)

# Define variational problem
u = Function(V)
v = TestFunction(V)
f = Constant(1.0)

F = u*v*dx + dt*inner(grad(u), grad(v))*dx - (u_n + dt*f)*v*dx

# Create VTK file for saving solution
vtkfile = File('test/solution.pvd')

# Time-stepping
u = Function(V)
t = 0
for n in range(num_steps):

    # Update current time
    t += dt

    # Compute solution
    solve(F == 0, u, bcs = bc)

    # Save to file and plot solution
    vtkfile << (u, t)

    # Update previous solution
    u_n.assign(u)