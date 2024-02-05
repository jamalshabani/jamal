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