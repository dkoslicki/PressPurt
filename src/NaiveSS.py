# For sanity purposes, let's implement the super naive version to determine sign stability
import numpy as np
import scipy.stats as st
#from multiprocessing import Pool  # Much faster without dummy (threading)
import multiprocessing
from pathos.multiprocessing import ProcessingPool as Pool
import scipy

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


def exists_switch(Ainv, Apinv):
	"""
	Takes in a dictionary of with keys eps symbols string (use symbol.name), values the values they are to be evaluated at.
	Returns 1 if a switch has occurred, 0 otherwise

	:param eps_dict: dictionary {eps_symbols: eps_values}
	:return: 0 or 1
	"""
	(m, n) = Ainv.shape
	switch = 0
	for i in range(m):
		for j in range(n):
			if np.sign(Ainv[i, j]) != np.sign(Apinv[i, j]):
				switch = 1
				break
		if switch == 1:
			break
	return switch

# other tri diag
#A = np.array([[-0.237, -1, 0, 0], [0.1, 0.015, -1, 0], [0, 0.1, -0.015, -1], [0, 0, 0.1, -0.015]])

# IGP
#A = np.array([[-0.237, -1, 0, 0], [0.1, -0.015, -1, -1], [0, 0.1, -0.015, -1], [0, .045, 0.1, -0.015]])

# Tri
#A = np.array([[-0.237, -1, 0, 0], [0.1, -0.015, -1, 0], [0, 0.1, -0.015, -1], [0, 0, 0.1, -0.015]])

#A = np.array([[-0.237, -1, 0, 0], [0.1, -0.015, -1, -1], [0, 0.1, -0.015, -1], [0, .005, 0.1, -0.015]])

#A = np.array([[-0.237, -1, 0, 0], [0.1, 0.015, -1, 0], [0, 0.1, -0.015, -1], [0, 0, 0.1, -0.015]])


def naive_SS(A, num_iterates, interval_length, num_threads):
	"""
	Computes the sign sensitivity (expected number of perturbations that lead to a sign switch in the inverse Jacobian) when perturbing multiple entries via Monte Carlo Sampling.

	:param A: input matrix
	:param num_iterates: number of Monte Carlo samples to make
	:param interval_length: length of interval to draw samples from
	:param num_threads: number of threads to use in the multithreading
	:return: float
	"""
	Ainv = np.linalg.inv(A)
	m, n = A.shape

	# get the pertubation intervals
	intervals_array = np.zeros((m,n,2))
	for i in range(m):
		for j in range(n):
			intervals_array[i, j, :] = intervals(A[i, j], interval_length)

	# Array of perturbation values
	pert_array = np.zeros((m, n, num_iterates))
	for i in range(m):
		for j in range(n):
			interval = intervals_array[i, j, :]
			dist = st.uniform(interval[0], interval[1])
			vals = dist.rvs(num_iterates)
			pert_array[i, j, :] = vals

	# helper function defining the work to do for one iterate
	def helper(it):
		stable_q = False
		switch_q = False
		Ap = A + pert_array[:, :, it]  # form the perturbation
		if is_stable(Ap):
			stable_q = True
			#Apinv = np.linalg.inv(Ap)
			Apinv = scipy.linalg.inv(Ap, overwrite_a=True, check_finite=False)
			is_switch = exists_switch(Ainv, Apinv)
			if is_switch:
				switch_q = True
		return stable_q, switch_q

	# do the work
	pool = Pool(processes=num_threads)
	res = pool.map(helper, range(num_iterates))

	# collect the results
	stable_counter = 0
	switch_counter = 0
	for stable_q, switch_q in res:
		if stable_q:
			stable_counter += 1
			if switch_q:
				switch_counter += 1
	# TODO: note, it might be nice to return stable_counter too, as this will indicate
	# how close the matrix is to being unstable (stable_counter small means many perturbations lead to instability)
	if stable_counter == 0:
		raise Exception("No perturbations lead to a stable matrix. Decrease the interval_length (-l) and/or increase the number of iterates. (-n)")
	return switch_counter / float(stable_counter)

	# old, serial way
	#stable_counter = 0
	#switch_counter = 0
	#for it in range(num_iterates):
	#	Ap = np.array(A)
	#	for i in range(m):
	#		for j in range(n):
	#			Ap[i, j] += pert_array[i, j, it]
	#	Apinv = np.linalg.inv(Ap)
	#	is_switch = exists_switch(Ainv, Apinv)
	#	if is_stable(Ap):
	#		stable_counter += 1
	#		if is_switch:
	#			switch_counter += 1
	#return switch_counter/float(stable_counter)

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
	ss = naive_SS(A, num_iterates=5000, interval_length=0.01)
	assert abs(ss - 0.37) < 0.1

	# Tri-diagonal
	A = np.array([[-0.237, -1, 0, 0], [0.1, -0.015, -1, 0], [0, 0.1, -0.015, -1], [0, 0, 0.1, -0.015]])
	ss = naive_SS(A, num_iterates=5000, interval_length=0.01)
	assert abs(ss - 0.0) < 0.01

	# Other
	A = np.array([[-0.337, -1, 0, 0], [0.1, -0.015, -1, -1], [0, 0.1, -0.015, -1], [0, .045, 0.1, -0.015]])
	try:
		ss = naive_SS(A, num_iterates=5000, interval_length=0.01)
	except:
		pass  # it should throw an error, since this matrix is not stable to begin with

	# Other 2
	A = np.array([[-0.237, -1, 0, 0], [0.1, -0.015, -1, -1], [0, 0.1, -0.015, -1], [0, .005, 0.1, -0.015]])
	ss = naive_SS(A, num_iterates=5000, interval_length=0.01)
	assert abs(ss - 0.34) < 0.1

	A = np.array([[-0.237, -1, 0, 0], [0.1, 0.015, -1, 0], [0, 0.1, -0.015, -1], [0, 0, 0.1, -0.015]])
	ss = naive_SS(A, num_iterates=5000, interval_length=0.01)
	assert abs(ss - 0.99) < 0.01
# TODO: this give 99% switch, whereas mathematica gives it only 50%, something weird is going on!!!!!
# TODO: figured out this is an issue with Mathematica, see NaiveSS.py for corroboration that this is correct



