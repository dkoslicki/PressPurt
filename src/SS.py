# this will implement equation 3.40 and 3.42 in the manuscript
# volume of sign sensitivity region when perturbing multiple entries
import numpy as np
import sympy as sp
import scipy.stats as st
from NumSwitch import interval_of_stability


def intervals(aij, x=0.01):
	"""
	This function defines an interval of size int_size around aij without changing the sign of aij
	:param aij: the (i,j) entry of a matrix (scalar)
	:param int_size: the size of the desired interval
	:return: an (ordered) list defining the endpoints of the interval
	"""
	if aij == 0:
		lower_bound = 0
		upper_bound = 0
	elif aij > 0:
		if x <= 2*aij:
			lower_bound = -x/2.
		elif x > 2*aij:
			lower_bound = -aij
		#lower_bound = np.piecewise(aij, [x <= 2*aij, x > 2*aij], [-x/2., -aij])
		if x <= 2*aij:
			upper_bound = x/2.
		elif x > 2*aij:
			upper_bound = x - aij
		#upper_bound = np.piecewise(aij, [x <= 2*aij, x >= 2*aij], [x/2., x - aij])
	elif aij < 0:
		if x <= 2*np.abs(aij):
			lower_bound = -x/2.
		elif x > 2*np.abs(aij):
			lower_bound = -x + np.abs(aij)
		#lower_bound = np.piecewise(aij, [x <= 2*np.abs(aij), x >= 2*np.abs(aij)], [-x/2., -x + np.abs(aij)])
		if x <= 2*np.abs(aij):
			upper_bound = x/2.
		elif x > 2*np.abs(aij):
			upper_bound = np.abs(aij)
		#upper_bound = np.piecewise(aij, [x <= 2*np.abs(aij), x >= 2*np.abs(aij)], [x/2., np.abs(aij)])
	return [lower_bound, upper_bound]


# TODO: will also want a function that gives me 1 or 0 depending on stability/non-stability

#entries_to_perturb = np.ones((4,4))
entries_to_perturb = np.array([[1,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]])
A = sp.Matrix(np.array([[-0.237, -1, 0, 0], [0.1, -0.015, -1, -1], [0, 0.1, -0.015, -1], [0, .045, 0.1, -0.015]]))
Ainv = sp.Matrix(A.inv())

# get the variables we are going to perturb
pert_locations_i, pert_locations_j = np.where(entries_to_perturb)
symbol_string = ""
for i, j in zip(pert_locations_i, pert_locations_j):
	symbol_string += "eps_%d_%d" % (i, j)
if len(pert_locations_i) == 1:
	symbol_tup = [sp.symbols(symbol_string)]
else:
	symbol_tup = list(sp.symbols(symbol_string))

# create the matrix with the symbolic perturbation values
B_temp = sp.Matrix(np.zeros(A.shape))
iter = 0
for i, j in zip(pert_locations_i, pert_locations_j):
	B_temp[i, j] = symbol_tup[iter]
	iter += 1

# form the symbolic matrix (A+B)^(-1)./A^(-1)
B = sp.Matrix(B_temp)
AplusB = A + B
AplusBinv = AplusB.inv()
AplusBinvDivAinv = sp.Matrix(np.zeros(A.shape))
for i in range(A.shape[0]):
	for j in range(A.shape[1]):
		AplusBinvDivAinv[i, j] = AplusBinv[i, j] / Ainv[i, j]

# lambdafy the symbolic quantity
AplusBinvDivAinveval = sp.lambdify(symbol_tup, AplusBinvDivAinv, "numpy")

# if any of the entries of AplusBinvDivAinveval are < 0, then a sign switch has occurred





###############################################################
# Tests
assert np.allclose(intervals(4, 1), [-0.5, 0.5])
assert np.allclose(intervals(4, 100), [-4, 96])
assert np.allclose(intervals(-4, 1), [-0.5, 0.5])
assert np.allclose(intervals(-4, 100), [-96, 4])