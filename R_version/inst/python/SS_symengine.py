# this will implement equation 3.40 and 3.42 in the manuscript
# volume of sign sensitivity region when perturbing multiple entries
import numpy as np
#import sympy as sp
import symengine as sp
import scipy.stats as st
from NumSwitch import interval_of_stability
import timeit


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
def exists_switch(eps_dict, AplusBinvDivAinvEval, symbol_tup, n):
	"""
	Takes in a dictionary of with keys eps symbols string (use symbol.name), values the values they are to be evaluated at.
	Returns 1 if a switch has occurred, 0 otherwise
	:param eps_dict: dictionary {eps_symbols: eps_values}
	:return: 0 or 1
	"""
	val_list = []
	for var in symbol_tup:
		val_list.append(eps_dict[var.name])
	AplusBinvDivAinvEvalulated = AplusBinvDivAinvEval(*val_list)[0].reshape((n,n))
	#(m, n) = AplusBinvDivAinvEvalulated.shape
	switch = 0
	for i in range(n):
		for j in range(n):
			if AplusBinvDivAinvEvalulated[i, j] < 0:
				switch = 1
				break
		if switch == 1:
			break
	return switch


def check_switch_matrix(eps_dict, AplusBinvDivAinvEval):
	"""
	Takes in a dictionary of with keys eps symbols string (use symbol.name), values the values they are to be evaluated at.
	Returns binary matrix telling if a switch occurred in that entry
	:param eps_dict: dictionary {eps_symbols: eps_values}
	:return: binary numpy array
	"""
	AplusBinvDivAinvEvalulated = AplusBinvDivAinvEval(**eps_dict)
	(m, n) = AplusBinvDivAinvEvalulated.shape
	switch = np.zeros((m, n))
	for i in range(m):
		for j in range(n):
			if AplusBinvDivAinvEvalulated[i, j] < 0:
				switch[i, j] = 1
	return switch


def is_stable(A):
	"""
	Check if the input matrix is asymptotically stable
	:param A: input matrix
	:return: Bool (1 iff asymptotically stable)
	"""
	[s, _] = np.linalg.eig(A)
	if all([np.real(i) < 0 for i in s]):
		return 1
	else:
		return 0


def get_entries_to_perturb(A):
	"""
	Returns a binary matrix indicating which entries to perturb
	:param A: input matrix (numpy or sympy matrix)
	:return: binary matrix of same dimensions as A
	"""
	m, n = A.shape
	res = np.zeros((m, n))
	for i in range(m):
		for j in range(n):
			if A[i, j] != 0:
				res[i, j] = 1
	return res



def SS(A, num_iterates=10000, interval_length=0.01):
	"""
	Computes equation 3.42: the volume of the number of perturbations that cause a sign switch in some part of the matrix
	:param A:  input matrix, numpy array
	:return: Scalar (percent of perturbations that caused some sign switch)
	"""
	if not is_stable(A):
		raise Exception(
			"The input matrix is not stable itself (one or more eigenvalues have non-negative real part). Cannot continue analysis.")
	t0 = timeit.default_timer()
	entries_to_perturb = get_entries_to_perturb(A)
	Ainv = sp.Matrix(np.linalg.inv(A).tolist())
	(m,n) = A.shape
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
	B_temp = np.zeros(A.shape).tolist() #sp.Matrix(np.zeros(A.shape).tolist())
	iter = 0
	for i, j in zip(pert_locations_i, pert_locations_j):
		B_temp[i][j] = symbol_tup[iter]
		iter += 1
	#B_temp = sp.Matrix(B_temp)
	t1 = timeit.default_timer()
	print("preprocess time: %f" %(t1 - t0))
	t0 = timeit.default_timer()
	# form the symbolic matrix (A+B)^(-1)./A^(-1)
	print(B_temp)
	B = sp.Matrix(B_temp)
	print(B)
	AplusB = sp.Matrix(A.tolist()) + B
	AplusBinv = AplusB.inv()
	#AplusBinv = AplusB.inv(method='LU', try_block_diag=True)
	#AplusBinv = AplusB.inv(method='ADJ', try_block_diag=True)
	t1 = timeit.default_timer()
	print("AplusBinv time: %f" %(t1 - t0))
	t0 = timeit.default_timer()
	# component-wise division
	AplusBinvDivAinv = sp.Matrix(np.zeros(A.shape).tolist())
	for i in range(A.shape[0]):
		for j in range(A.shape[1]):
			AplusBinvDivAinv[i, j] = AplusBinv[i, j] / float(Ainv[i, j])
	t1 = timeit.default_timer()
	print("AplusBinvDivAinv time: %f" % (t1 - t0))
	t0 = timeit.default_timer()
	# lambdafy the symbolic quantity
	AplusBinvDivAinvEval = sp.lambdify(symbol_tup, AplusBinvDivAinv, "numpy")
	AplusBEval = sp.lambdify(symbol_tup, AplusB, "numpy")
	t1 = timeit.default_timer()
	print("lambdify time: %f" % (t1 - t0))

	#num_iterates = 10000
	#interval_length = 0.01
	switch_count = 0
	is_stable_count = 0

	# initialize the dictionaries of values to pass
	eps_dicts = []
	for iterate in range(num_iterates):
		eps_dicts.append(dict())
	t0 = timeit.default_timer()
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
	t1 = timeit.default_timer()
	print("Sample time: %f" % (t1 - t0))
	# check for sign switches and stability
	t0 = timeit.default_timer()
	for eps_dict in eps_dicts:
		val_list = []
		for var in symbol_tup:
			val_list.append(eps_dict[var.name])
		stab_indicator = is_stable(AplusBEval(*val_list)[0].reshape(A.shape))
		if stab_indicator:
			is_stable_count += 1
		if exists_switch(eps_dict, AplusBinvDivAinvEval, symbol_tup, n) and stab_indicator:
			switch_count += 1
	t1 = timeit.default_timer()
	print("Actual eval time: %f" % (t1 - t0))
	#print(switch_count)
	#print(num_iterates)
	#print(is_stable_count)
	return switch_count / float(num_iterates)  # this is how it was done in the paper/mathematica
	#return switch_count / float(is_stable_count)  # I think this is the correct way to do it








###############################################################
# Tests
def tests():
	"""
	Run all the tests
	:return: None
	"""
	assert np.allclose(intervals(4, 1), [-0.5, 0.5])
	assert np.allclose(intervals(4, 100), [-4, 96])
	assert np.allclose(intervals(-4, 1), [-0.5, 0.5])
	assert np.allclose(intervals(-4, 100), [-96, 4])

	# SS function tests
	# IGP
	A = np.array([[-0.237, -1, 0, 0], [0.1, -0.015, -1, -1], [0, 0.1, -0.015, -1], [0, .045, 0.1, -0.015]])
	ss = SS(A, num_iterates=5000, interval_length=0.01)
	assert abs(ss - 0.37) < 0.1

	# Tri-diagonal
	A = np.array([[-0.237, -1, 0, 0], [0.1, -0.015, -1, 0], [0, 0.1, -0.015, -1], [0, 0, 0.1, -0.015]])
	ss = SS(A, num_iterates=5000, interval_length=0.01)
	assert abs(ss - 0.0) < 0.01

	# Other
	A = np.array([[-0.337, -1, 0, 0], [0.1, -0.015, -1, -1], [0, 0.1, -0.015, -1], [0, .045, 0.1, -0.015]])
	try:
		ss = SS(A, num_iterates=5000, interval_length=0.01)
	except:
		pass  # it should throw an error, since this matrix is not stable to begin with

	# Other 2
	A = np.array([[-0.237, -1, 0, 0], [0.1, -0.015, -1, -1], [0, 0.1, -0.015, -1], [0, .005, 0.1, -0.015]])
	ss = SS(A, num_iterates=5000, interval_length=0.01)
	assert abs(ss - 0.34) < 0.1

	A = np.array([[-0.237, -1, 0, 0], [0.1, 0.015, -1, 0], [0, 0.1, -0.015, -1], [0, 0, 0.1, -0.015]])
	ss = SS(A, num_iterates=5000, interval_length=0.01)
	assert abs(ss - 0.99) < 0.01
	# TODO: this give 99% switch, whereas mathematica gives it only 50%, something weird is going on!!!!!
	# TODO: figured out this is an issue with Mathematica, see NaiveSS.py for corroboration that this is correct

###############################################################
# Numerical accuracy of AplusBinvDivAinvEval
#eps_vals = list(np.ones(len(symbol_tup)))
#eps_dict = dict(zip([i.name for i in symbol_tup], eps_vals))
#print(check_switch_matrix(eps_dict, AplusBinvDivAinvEval))
#print(AplusBinvDivAinvEval(**eps_dict))
# python
#[[ 0.10876355 -0.         -0.         -0.        ]
# [ 0.37958414 -0.11109425  0.          0.        ]
# [ 0.97667782 -0.28584778 -0.21597541 -0.        ]
# [ 0.6636542  -0.19423404 -9.52164049 -0.10449546]]
# Mathematica
#0.108764	0.	0.	0.
#0.379584	-0.111094	-1.10792*10^-17	2.58702*10^-18
#0.976678	-0.285848	-0.215975	0.
#0.663654	-0.194234	-9.52164	-0.104495

# TODO: there are differences here between mathematica and sympy, due to numerical issues
# example: Mathematica returns -1.10792*10^-17, yet sympy gives a 0
