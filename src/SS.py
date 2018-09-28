# this will implement equation 3.40 and 3.42 in the manuscript
# volume of sign sensitivity region when perturbing multiple entries
import numpy as np
import sympy as sp
import scipy.stats as st
from NumSwitch import interval_of_stability
#entries_to_perturb = np.ones((4,4))
entries_to_perturb = np.array([[1,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]])
Aigp = np.array([[-0.237, -1, 0, 0], [0.1, -0.015, -1, -1], [0, 0.1, -0.015, -1], [0, .045, 0.1, -0.015]])
Aigpinv = np.linalg.inv(Aigp)
pert_locations_i, pert_locations_j = np.where(entries_to_perturb)
symbol_string = ""
for i,j in zip(pert_locations_i, pert_locations_j):
	symbol_string += "eps_%d_%d" % (i, j)
symbol_tup = sp.symbols(symbol_string)

# for an initial approach, we will take a uniform distribution with a specified length, that don't change the sign

A = Aigp
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



###############################################################
# Tests
assert np.allclose(intervals(4, 1), [-0.5, 0.5])
assert np.allclose(intervals(4, 100), [-4, 96])
assert np.allclose(intervals(-4, 1), [-0.5, 0.5])
assert np.allclose(intervals(-4, 100), [-96, 4])