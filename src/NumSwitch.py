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
	if any([index >= n for index in [i, j, k, l]]):
		raise Exception("Matrix is only %dx%d, invalid choice of subscripts: %d, %d, %d, %d" % (n, n, i, j, k, l))
	if 1+eps*Ainv[l, k] == 0:
		raise Exception("Division by zero in sherman ratio")
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


def interval_of_stability(A, Ainv, k, l, max_bound=10, num_sample=1000):
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
	(to_sample, step_size) = np.linspace(initial_interval[0], initial_interval[1], num_sample, retstep=True)
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
	return (lower_bound + step_size, upper_bound - step_size)


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
	# TODO: this is too slow when it's stable unbounded on one end
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
	return (lower_bound + step_size, upper_bound - step_size)


def exp_num_switch(A, Ainv, k, l, num_sample=1000, dist=None, interval=None):
	"""
	This implements equation 3.6: the expected number of sign switches
	:param A: The original input matrix (numpy array)
	:param Ainv: The inverse of A (pre-computed), numpy array
	:param k: index
	:param l: index
	:param dist: A probability density function: function of a single argument
	:param interval: the interval of stability (optional, can be pre-computed with interval_of_stability())
	:return:
	"""
	m, n = A.shape
	if A[k, l] == 0:
		raise Exception("You can only perturb non-zero entries: A[%d, %d] is zero." % (k, l))
	# get the region of stability
	if not interval:
		interval = interval_of_stability(A, Ainv, k, l, num_sample=num_sample)
	# NS is a step function, so find the respective values and intervals
	(to_sample, step_size) = np.linspace(interval[0], interval[1], num_sample, retstep=True)
	ns_values = [NS(Ainv, eps, k, l) for eps in to_sample]
	to_compute_args = []
	initial_loc = to_sample[0]
	initial_val = ns_values[0]
	for i in range(1, len(ns_values)):
		if ns_values[i] != initial_val:
			if i < len(ns_values) - 1:
				terminal_loc = to_sample[i - 1] + step_size/float(2)  # mid-point like rule
			else:
				terminal_loc = to_sample[i - 1]
			terminal_val = ns_values[i - 1]
			to_append = (initial_val, (initial_loc, terminal_loc))  # format is (ns_value, (interval_start, interval_end))
			to_compute_args.append(to_append)
			initial_loc = to_sample[i] - step_size/float(2)
			initial_val = ns_values[i]
	# make sure to add the last value
	to_compute_args.append((initial_val, (initial_loc, to_sample[-1])))
	#return to_compute_args
	# get the right distribution function
	if not dist:
		my_std = (interval[1] - interval[0])/2.
		my_mean = 0
		a, b = (interval[0] - my_mean) / my_std, (interval[1] - my_mean) / my_std
		dist = st.truncnorm(a, b, 0, my_std)
	# Now do the integration
	exp_value = 0
	for (value, (start, end)) in to_compute_args:
		exp_value += value*(dist.cdf(end) - dist.cdf(start))
	return exp_value/float(m*n)

# This has been checked against Mathematica
# a, b = (-.5 - 0) / .75, (1 - 0) / .75
#dist = st.truncnorm(a, b, 0, .75)
#x = np.linspace(-.6, 1.1, 100)
#fig, ax = plt.subplots(1, 1)
#ax.plot(x, dist.pdf(x),'r-',label='norm pdf')
#fig.show()
# dist.pdf(0)





#######################################################################
# test cases
# test matrix
Atri = np.array([[-0.237, -1, 0, 0], [0.1, -0.015, -1, 0], [0, 0.1, -0.015, -1], [0, 0, 0.1, -0.015]])
Aigp = np.array([[-0.237, -1, 0, 0], [0.1, -0.015, -1, -1], [0, 0.1, -0.015, -1], [0, .045, 0.1, -0.015]])
Atriinv = np.linalg.inv(Atri)
Aigpinv = np.linalg.inv(Aigp)

# num switch tests
assert NS(Aigpinv, -.2, 3, 2) == 9
assert NS(Aigpinv, .2, 3, 2) == 1
assert NS(Aigpinv, -.5, 1, 2) == 16
assert NS(Aigpinv, .5, 1, 2) == 1

# interval of stability tests
interval = interval_of_stability(Aigp, Aigpinv, 0, 0, num_sample=1000)
assert np.abs(interval[0] - -0.08298659074229853) < .001
assert np.abs(interval[1] - 0.10901820241962716) < .001
interval = interval_of_stability_crawl(Aigp, Aigpinv, 0, 0, step_size=.0001)
assert np.abs(interval[0] - -0.08298659074229853) < .001
assert np.abs(interval[1] - 0.10901820241962716) < .001
interval = interval_of_stability(Aigp, Aigpinv, 1, 0, num_sample=1000)
assert np.abs(interval[0] - -0.025934401526958355) < .02
assert np.abs(interval[1] - 10) < .02

# This test case takes too long
#interval = interval_of_stability_crawl(Aigp, Aigpinv, 1, 0, step_size=.0001)
#assert np.abs(interval[0] - -0.025934401526958355) < .001
#assert np.abs(interval[1] - 10) < .001

# expected num switch tests
exp_value = exp_num_switch(Aigp, Aigpinv, 3, 2)
assert np.abs(exp_value - 0.081) < 0.01
exp_value = exp_num_switch(Aigp, Aigpinv, 1, 2)
assert np.abs(exp_value - 0.048) < 0.01
exp_value = exp_num_switch(Aigp, Aigpinv, 0, 0)
assert np.abs(exp_value - 0.032) < 0.01
exp_value = exp_num_switch(Aigp, Aigpinv, 1, 1)
assert np.abs(exp_value - 0.12) < 0.01
exp_value = exp_num_switch(Aigp, Aigpinv, 2, 2)
assert np.abs(exp_value - 0.16) < 0.01
exp_value = exp_num_switch(Aigp, Aigpinv, 3, 3)
assert np.abs(exp_value - 0.17) < 0.01
exp_value = exp_num_switch(Atri, Atriinv, 0, 0)
assert np.abs(exp_value - 0) < 0.01