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


# if any of the entries of AplusBinvDivAinveval are < 0, then a sign switch has occurred
def exists_switch(eps_dict, AplusBinvDivAinvEval):
	"""
	Takes in a dictionary of with keys eps symbols string (use symbol.name), values the values they are to be evaluated at.
	Returns 1 if a switch has occurred, 0 otherwise
	:param eps_dict: dictionary {eps_symbols: eps_values}
	:return: 0 or 1
	"""
	AplusBinvDivAinvEvalulated = AplusBinvDivAinvEval(**eps_dict)
	(m, n) = AplusBinvDivAinvEvalulated.shape
	switch = 0
	for i in range(m):
		for j in range(n):
			if AplusBinvDivAinvEvalulated[i, j] < 0:
				switch = 1
				break
		if switch == 1:
			break
	return switch


def is_stable(A):
	[s, _] = np.linalg.eig(A)
	if all([np.real(i) < 0 for i in s]):
		return 1
	else:
		return 0

#entries_to_perturb = np.ones((4,4))
#entries_to_perturb = np.array([[1,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]])
entries_to_perturb = np.array([[1,1,0,0],[1,1,1,1],[0,1,1,1],[0,1,1,1]])  # TODO: programmatically make it perturb all non-zero values
A = sp.Matrix(np.array([[-0.237, -1, 0, 0], [0.1, -0.015, -1, -1], [0, 0.1, -0.015, -1], [0, .045, 0.1, -0.015]]))
Ainv = sp.Matrix(A.inv())

# get the variables we are going to perturb
pert_locations_i, pert_locations_j = np.where(entries_to_perturb)
symbol_string = ""
for i, j in zip(pert_locations_i, pert_locations_j):
	symbol_string += "eps_%d_%d " % (i, j)
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

# component-wise division
AplusBinvDivAinv = sp.Matrix(np.zeros(A.shape))
for i in range(A.shape[0]):
	for j in range(A.shape[1]):
		AplusBinvDivAinv[i, j] = AplusBinv[i, j] / Ainv[i, j]

# lambdafy the symbolic quantity
AplusBinvDivAinvEval = sp.lambdify(symbol_tup, AplusBinvDivAinv, "numpy")
AplusBEval = sp.lambdify(symbol_tup, AplusB, "numpy")


num_iterates = 20000
interval_length = 0.015
switch_count = 0
is_stable_count = 0

# initialize the dictionaries of values to pass
eps_dicts = []
for iterate in range(num_iterates):
	eps_dicts.append(dict())

# for each one of the symbols, sample from the appropriate distribution
for symbol in symbol_tup:
	symbol_name = symbol.name
	i = eval(symbol_name.split('_')[1])
	j = eval(symbol_name.split('_')[2])
	interval = intervals(A[i, j], interval_length)
	dist = st.uniform(interval[0], interval[1])
	vals = dist.rvs(num_iterates)
	iter = 0
	for eps_dict in eps_dicts:
		eps_dict[symbol_name] = vals[iter]
		iter += 1

# check for sign switches and stability
for eps_dict in eps_dicts:
	if is_stable(AplusBEval(**eps_dict)):
		is_stable_count += 1
	if exists_switch(eps_dict, AplusBinvDivAinvEval) and is_stable(AplusBEval(**eps_dict)):
		switch_count += 1
print(switch_count / float(num_iterates))
print(switch_count / float(is_stable_count))








###############################################################
# Tests
assert np.allclose(intervals(4, 1), [-0.5, 0.5])
assert np.allclose(intervals(4, 100), [-4, 96])
assert np.allclose(intervals(-4, 1), [-0.5, 0.5])
assert np.allclose(intervals(-4, 100), [-96, 4])