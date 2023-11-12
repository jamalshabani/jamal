from firedrake import *
from petsc4py import PETSc

T = 2.0            # final time
num_steps = 10    # number of time steps
dt = T / num_steps # time step size

# Create mesh and define function space
# Import gmesh
mesh = Mesh(options.mesh)
Id = Identity(mesh.geometric_dimension()) #Identity tensor

# Define the function spaces
V = FunctionSpace(mesh, 'CG', 1)
VV = VectorFunctionSpace(mesh, 'CG', 1, dim = 2)

f = Function(VV)


f = Constant((1,1))
f = interpolate(f, VV)

File("test1/f.pvd").write(f)

print(list(range(11)))

# Define boundary condition
bc = DirichletBC(V, Constant(0), [1,2,3,4])

# Define initial value
u_n = interpolate(Constant(1.0), V)

# Define variational problem
u = Function(V)
v = TestFunction(V)
g = interpolate(Constant(1.0), V)

#F = (u - u_n)/dt*v*dx + inner(grad(u), grad(v))*dx - g*v*dx
F = u*v*dx + dt*inner(grad(u), grad(v))*dx - (u_n+dt*g)*v*dx
# Create VTK file for saving solution
vtkfile = File('test1/solution.pvd')

# Time-stepping
t = 0
for n in range(num_steps):

    # Update current time
    t += dt

    # Compute solution
    solve(F == 0, u, bcs = bc)

    # Save to file and plot solution
    vtkfile.write(u, time = t)

    # Update previous solution
    u_n.assign(u)