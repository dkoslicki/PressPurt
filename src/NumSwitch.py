# This script will implement equations 3.4, 3.5, and 3.6
import numpy as np
import scipy.integrate as integrate
import scipy.stats as st
from scipy import special

def ind_switch(Ainv, eps, i, j, k, l):
	"""
	This function implements equation 3.4, telling if a sign switch has occured
	:param Ainv: input matrix, inverted, numpy array
	:param eps: perturbation size (scalar)
	:param i: index
	:param j: index
	:param k: index
	:param l: index
	:return: 1 or 0
	"""
	n = Ainv.shape[0]
	if any([index>=n for index in [i, j, k, l]]):
		raise Exception("Matrix is only %dx%d, invalid choice of subscripts: %d, %d, %d, %d" % (n, n, i, j, k, l))
	sherm_ratio = eps*Ainv[i, k]*Ainv[l, j]/(Ainv[i,j]*(1+eps*Ainv[l, k]))
	if 1 - sherm_ratio < 0:
		return 1
	elif 1 - sherm_ratio > 0:
		return 0
	else:
		raise Exception("Invalid sherman morrison ratio: %f" % sherm_ratio)


def NS(Ainv, eps, k, l):
	"""
	This function implements equation 3.5: gives the number of sign switches
	:param Ainv: inverse input matrix, numpy array
	:param eps: perturbation size (scalar)
	:param k: index
	:param l: inex
	:return: Natural number
	"""
	n = Ainv.shape[0]
	if any([index >= n for index in [k, l]]):
		raise Exception("Matrix is only %dx%d, invalid choice of subscripts: %d, %d." % (n, n, k, l))
	ns_sum = 0
	for i in range(n):
		for j in range(n):
			ns_sum += ind_switch(Ainv, eps, i, j, k, l)
	return ns_sum


def interval_of_stability(A, Ainv, k, l, max_bound=10, num_sample=100):
	"""
	This function will return the interval of stability when perturbing just the (k,l) entry of A
	:param A: numpy array
	:param Ainv: inverse of A
	:param k: index
	:param l: index
	:param max_bound: some are always stable on one side, only go out this far), scalar > 0
	:param num_sample: number of points (natural number) to sample to determine stability
	:return: (a:scalar, b:scalar) the largest interval where A[k,l]+eps is stable
	"""
	# Input matrix needs to be stable
	[s, _] = np.linalg.eig(A)
	if not all([np.real(i) < 0 for i in s]):
		raise Exception("The input matrix is not stable itself (one or more eigenvalues have non-negative real part). Cannot continue analysis.")

	# Find an initial region of stability. Use contrapositive of theorem 3.1
	if A[k, l] > 0 and Ainv[l, k] < 0:  # pos neg
		initial_interval = (-A[k, l], -1/float(Ainv[l, k]))
	if A[k, l] < 0 and Ainv[l, k] > 0:  # neg pos
		initial_interval = (-1 / float(Ainv[l, k]), -A[k, l])
	if A[k, l] < 0 and Ainv[l, k] < 0:  # neg neg
		initial_interval = (-max_bound, min([-1 / float(Ainv[l, k]), -A[k, l]]))
	if A[k, l] > 0 and Ainv[l, k] > 0:  # pos pos
		initial_interval = (max([-1 / float(Ainv[l, k]), -A[k, l]]), max_bound)
	# now need to sample to determine the actual region of stability
	# we'll basically do a grid search and look for the real parts of the eigenvalues begin negative
	to_sample = np.linspace(initial_interval[0], initial_interval[1], num_sample)
	#return to_sample
	zero_matrix = np.zeros(A.shape)  # matrix we will be perturbing by
	eig_values = []
	for pert_val in to_sample:
		zero_matrix[k, l] = pert_val  # update the perturb value
		[s, _] = np.linalg.eig(A + zero_matrix)
		eig_values.append(max([np.real(i) for i in s]))  # look at the largest eigenvalue
	#return eig_values
	# now, return the largest interval where all the eigenvalues have negative real part
	zero_loc = int(np.argmin(np.abs(to_sample)))
	upper_bound = initial_interval[1]
	for i in range(len(to_sample) - zero_loc):
		if eig_values[zero_loc + i] >= 0:
			upper_bound = to_sample[zero_loc + i - 1]
			break
	lower_bound = initial_interval[0]
	for i in range(zero_loc):
		if eig_values[zero_loc - i] >= 0:
			lower_bound = to_sample[zero_loc - i + 1]
			break
	return (lower_bound, upper_bound)


def interval_of_stability_crawl(A, Ainv, k, l, max_bound=10, step_size=.0001):
	"""
	This function will return the interval of stability when perturbing just the (k,l) entry of A
	Difference with this implementation is that it decides to crawl out from the origin, looking for instability
	:param A: numpy array
	:param Ainv: inverse of A
	:param k: index
	:param l: index
	:param max_bound: some are always stable on one side, only go out this far), scalar > 0
	:param step_size: step size for perturbation values
	:return: (a:scalar, b:scalar) the largest interval where A[k,l]+eps is stable
	"""
	# Input matrix needs to be stable
	[s, _] = np.linalg.eig(A)
	if not all([np.real(i) < 0 for i in s]):
		raise Exception("The input matrix is not stable itself (one or more eigenvalues have non-negative real part). Cannot continue analysis.")

	# Find an initial region of stability. Use contrapositive of theorem 3.1
	if A[k, l] > 0 and Ainv[l, k] < 0:  # pos neg
		initial_interval = (-A[k, l], -1/float(Ainv[l, k]))
	if A[k, l] < 0 and Ainv[l, k] > 0:  # neg pos
		initial_interval = (-1 / float(Ainv[l, k]), -A[k, l])
	if A[k, l] < 0 and Ainv[l, k] < 0:  # neg neg
		initial_interval = (-max_bound, min([-1 / float(Ainv[l, k]), -A[k, l]]))
	if A[k, l] > 0 and Ainv[l, k] > 0:  # pos pos
		initial_interval = (max([-1 / float(Ainv[l, k]), -A[k, l]]), max_bound)
	# now need to sample to determine the actual region of stability
	# we'll basically do a grid search and look for the real parts of the eigenvalues begin negative
	# crawl to the right
	zero_matrix = np.zeros(A.shape)
	eps = step_size
	all_neg_eigs = True
	while all_neg_eigs is True:
		zero_matrix[k, l] = eps
		[s, _] = np.linalg.eig(A + zero_matrix)
		max_eig_val = max([np.real(i) for i in s])
		if max_eig_val >= 0 or eps >= initial_interval[1]:
			eps -= step_size
			all_neg_eigs = False
		else:
			eps += step_size
	upper_bound = eps
	# crawl to the right
	zero_matrix = np.zeros(A.shape)
	eps = -step_size
	all_neg_eigs = True
	while all_neg_eigs is True:
		zero_matrix[k, l] = eps
		[s, _] = np.linalg.eig(A + zero_matrix)
		max_eig_val = max([np.real(i) for i in s])
		if max_eig_val >= 0 or eps <= initial_interval[0]:
			eps += step_size
			all_neg_eigs = False
		else:
			eps -= step_size
	lower_bound = eps
	return (lower_bound, upper_bound)


def exp_num_switch(A, Ainv, k, l, dist=None):
	"""
	This implements equation 3.6: the expected number of sign switches
	:param A: The original input matrix (numpy array)
	:param Ainv: The inverse of A (pre-computed), numpy array
	:param k: index
	:param l: index
	:param dist: A probability density function: function of a single argument
	:return:
	"""
	if A[k, l] == 0:
		raise Exception("You can only perturb non-zero entries: A[%d, %d] is zero." % (k, l))
	# use a default normal distribution
	# TODO: make sure this stays in the region of stability
	if not dist:
		pass# TODO: put custom pdf in here
	return



# test cases
# test matrix
Atri = np.array([[-0.237, -1, 0, 0], [0.1, -0.015, -1, 0], [0, 0.1, -0.015, -1], [0, 0, 0.1, -0.015]])
Aigp = np.array([[-0.237, -1, 0, 0], [0.1, -0.015, -1, -1], [0, 0.1, -0.015, -1], [0, .045, 0.1, -0.015]])
Atriinv = np.linalg.inv(Atri)
Aigpinv = np.linalg.inv(Aigp)
assert NS(Aigpinv, -.2, 3, 2) == 9
assert NS(Aigpinv, .2, 3, 2) == 1
assert NS(Aigpinv, -.5, 1, 2) == 16
assert NS(Aigpinv, .5, 1, 2) == 1

interval = interval_of_stability(Aigp, Aigpinv, 0, 0, num_sample=1000)
assert np.abs(interval[0] - -0.08298659074229853) < .001
assert np.abs(interval[1] - 0.10901820241962716) < .001

interval = interval_of_stability_crawl(Aigp, Aigpinv, 0, 0, step_size=.0001)
assert np.abs(interval[0] - -0.08298659074229853) < .001
assert np.abs(interval[1] - 0.10901820241962716) < .001



