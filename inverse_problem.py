from firedrake import *
from petsc4py import PETSc
import time
import numpy as np

start = time.time()

mesh = UnitSquareMesh(100, 100)

File("problem/mesh.pvd").write(mesh)