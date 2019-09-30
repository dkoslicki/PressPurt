# This script will implement equations 3.4, 3.5, and 3.6
import numpy as np
import scipy.stats as st
import warnings
from scipy.optimize import bisect, brentq
import os
import pandas as pd
import numpy

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
		raise Exception("Division by zero in sherman ratio. This perturbation caused the matrix to be non-invertible.")
	if Ainv[i, j] != 0:  # if it's non-zero use the pre-divided version
		sherm_ratio = eps*Ainv[i, k]*Ainv[l, j]/(Ainv[i,j]*(1+eps*Ainv[l, k]))
		if 1 - sherm_ratio < 0:
			return 1
		elif 1 - sherm_ratio > 0:
			return 0
		else:
			raise Exception("Invalid sherman morrison ratio: %f" % sherm_ratio)
	else:  # in the case where the inverse is zero, use the non-divided version
		if eps*Ainv[i, k]*Ainv[l, j] != 0:  # count any change from 0 as a switch
			return 1
		else:  # else you still have a zero here, so no switch
			return 0


def NS(Ainv, eps, k, l):
	"""
	This function implements equation 3.5: gives the number of sign switches
	when perturbing the k, l entry

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


def largest_root(A, k, l, eps):
	"""
	Returns the largest real part of the eigenvalues of A when perturbing the (k,l) entry by eps.

	:param A: numpy matrix
	:param k: int (row)
	:param l: int (column)
	:param eps: float (perturbation value)
	:return:
	"""
	zero_matrix = np.zeros(A.shape)
	zero_matrix[k, l] = eps
	[s, _] = np.linalg.eig(A + zero_matrix)
	return max([np.real(i) for i in s])


def interval_of_stability(A, Ainv, k, l, max_bound=10):
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

	if A[k, l] == 0:
		warnings.warn("You are attempting to perturb a zero entry: A[%d, %d] == 0." % (k, l))

	# Find an initial region of stability. Use contrapositive of theorem 3.1
	if A[k, l] > 0 and Ainv[l, k] < 0:  # pos neg
		initial_interval = (-A[k, l], -1/float(Ainv[l, k]))
	elif A[k, l] < 0 and Ainv[l, k] > 0:  # neg pos
		initial_interval = (-1 / float(Ainv[l, k]), -A[k, l])
	elif A[k, l] < 0 and Ainv[l, k] < 0:  # neg neg
		initial_interval = (-max_bound, min([-1 / float(Ainv[l, k]), -A[k, l]]))
	elif A[k, l] > 0 and Ainv[l, k] > 0:  # pos pos
		initial_interval = (max([-1 / float(Ainv[l, k]), -A[k, l]]), max_bound)
	else:
		initial_interval = (-max_bound, max_bound)
	initial_interval = list(initial_interval)
	if initial_interval[1] > max_bound:
		initial_interval[1] = max_bound
	if initial_interval[0] < -max_bound:
		initial_interval[0] = -max_bound
	initial_interval = tuple(initial_interval)
		#num_sample *= 2*max_bound  # if you fall in this case, better crank up the number of samples
	# now need to sample to determine the actual region of stability

	# brentq root finding, basically binary search
	if largest_root(A, k, l, initial_interval[1]) < 0:
		upper_bound = initial_interval[1]
	else:
		upper_bound = brentq(lambda x: largest_root(A, k, l, x), 0, initial_interval[1])
	if largest_root(A, k, l, initial_interval[0]) < 0:  # it's monotonic
		lower_bound = initial_interval[0]
	else:
		lower_bound = brentq(lambda x: largest_root(A, k, l, x), initial_interval[0], 0)
	small_num = .00000001
	if (lower_bound + small_num) > (upper_bound - small_num):
		return (lower_bound, upper_bound)
	else:
		return (lower_bound + small_num, upper_bound - small_num)
	# bisection method, depreciated since bisect is so much faster.
	# we'll basically do a grid search and look for the real parts of the eigenvalues begin negative
	#(to_sample, step_size) = np.linspace(initial_interval[0], initial_interval[1], num_sample, retstep=True)
	#zero_matrix = np.zeros(A.shape)  # matrix we will be perturbing by
	#eig_values = []
	#for pert_val in to_sample:
	#	zero_matrix[k, l] = pert_val  # update the perturb value
	#	[s, _] = np.linalg.eig(A + zero_matrix)
	#	eig_values.append(max([np.real(i) for i in s]))  # look at the largest eigenvalue
	# now, return the largest interval where all the eigenvalues have negative real part
	#zero_loc = int(np.argmin(np.abs(to_sample)))
	#upper_bound = initial_interval[1]
	#for i in range(len(to_sample) - zero_loc):
	#	if eig_values[zero_loc + i] >= 0:
	#		upper_bound = to_sample[zero_loc + i - 1]
	#		break
	#lower_bound = initial_interval[0]
	#for i in range(zero_loc):
	#	if eig_values[zero_loc - i] >= 0:
	#		lower_bound = to_sample[zero_loc - i + 1]
	#		break
	# adjust them, only if it doesn't screw up the order of the interval (for very small in magnitude
	# lower and upper bounds)
	#if upper_bound > step_size:
	#	upper_bound -= step_size
	#if np.abs(lower_bound) > step_size:
	#	lower_bound += step_size
	#return (lower_bound, upper_bound)


def exp_num_switch(A, Ainv, k, l, num_sample=1000, dist=None, interval=None):
	"""
	This implements equation 3.6: the expected number of sign switches

	:param A: The original input matrix (numpy array)
	:param Ainv: The inverse of A (pre-computed), numpy array
	:param k: index
	:param l: index
	:param num_sample: Number of samples to take to determine the interval of stability and expected value
	:param dist: A probability density function: function of a single argument
	:param interval: the interval of stability (optional, can be pre-computed with interval_of_stability())
	:return:
	"""
	m, n = A.shape
	if A[k, l] == 0:
		raise Exception("You can only perturb non-zero entries: A[%d, %d] is zero." % (k, l))
	# get the region of stability
	if interval is None:
		interval = interval_of_stability(A, Ainv, k, l)
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


def critical_epsilon(Ainv, k, l, i, j):
	"""
	This function finds which epsilon causes the sherman morrison ratio to
	equal 0 (indicating a switch will occur after this value).
	Takes advantage of the monotonicity of the SM ratio.

	:param Ainv: Inverse input matrix
	:param k: index (row being perturbed)
	:param l: index (column being perturbed)
	:param i: index (entry looking for a switch)
	:param j: index (entry looking for a switch)
	:return: scalar (epsilon value)
	"""
	denom = Ainv[i, j]*Ainv[l, k] - Ainv[i, k] * Ainv[l,j]
	if denom == 0:
		res = np.inf
	else:
		res = -Ainv[i, j] / float(denom)
	return res


def num_switch_from_crit_eps(crit_epsilon_array, stab_int_array, k, l):
	"""
	This function will get the full description of the num_switch function
	(using interval notation).

	:param crit_epsilon_array: array of critical epsilon (tensor indexed by (k, l, i, j))
	:param stab_int_array: array of intervals of stability (tensor index by (k, l, start, end))
	:param k: index (row to perturb)
	:param l: index (column to perturb)
	:return: list with entries (num_switch_value, (start, end))
	"""
	crit_eps_in_range = []
	n = crit_epsilon_array.shape[0]  # number of rows/columns
	stab_start = stab_int_array[k, l, 0]  # stability interval start
	stab_end = stab_int_array[k, l, 1]  # stability interval end
	# get the epsilons in the range
	for i in range(n):
		for j in range(n):
			crit_eps = crit_epsilon_array[k, l, i, j]
			if crit_eps < stab_end and crit_eps > stab_start:
				crit_eps_in_range.append(crit_eps)
	# get the epsilon values in the stability range
	# both lists start from zero outward
	# pad with stability endpoints in order to get the endpoints right
	pos_crit_eps_in_range = np.append(np.sort([i for i in crit_eps_in_range if i > 0]), stab_end)
	neg_crit_eps_in_range = np.append(np.sort([i for i in crit_eps_in_range if i < 0])[::-1], stab_start) #less to more negative
	num_switch_func = []
	# do the positives first
	is_first = True
	has_zero = 0 in crit_eps_in_range
	zero_count = np.sum(crit_eps_in_range == 0)
	if pos_crit_eps_in_range.size:  # if the list is not empty
		current = pos_crit_eps_in_range[0]  # need to start somewhere
		counter = zero_count + 1 # if zero is in the list, any perturbation will cause a switch
		# start at 1 otherwise (since there's at least one critical epsilon so one switch
		for index in range(1, len(pos_crit_eps_in_range)):
			next = pos_crit_eps_in_range[index]  # get the next guy in the list
			if next == current:  # if it's the same, just increment the counter
				counter += 1
			else:
				if is_first and has_zero:  # if zero is in the list, need to start the interval at zero
					num_switch_func.append((counter, (0, next)))
					is_first = False
				else:
					# otherwise, add this interval in to the return list
					num_switch_func.append((counter, (current, next)))
				current = next
				counter += 1  # you crossed another crit eps, so increment
	# Then do negatives (same as above, just a different list)
	is_first = True
	if neg_crit_eps_in_range.size:
		current = neg_crit_eps_in_range[0]
		counter = zero_count + 1
		for index in range(1, len(neg_crit_eps_in_range)):
			next = neg_crit_eps_in_range[index]
			if next == current:
				counter += 1
			else:
				if is_first and has_zero:
					num_switch_func.append((counter, (next, 0)))
					is_first = False
				else:
					num_switch_func.append((counter, (next, current)))
				current = next
				counter += 1
	return num_switch_func


def exp_num_switch_from_crit_eps(n, k, l, num_switch_funcs, stab_intervals, dist=None):
	"""
	Compute the expectation using the num_switch function

	:param n: dimension of the matrix
	:param k: index (pert row)
	:param l: index (pert column)
	:param num_switch_funcs: The shape of the num_switch_function from num_switch_from_crit_eps
	:param stab_intervals: all the intervals of stability
	:param dist: probability distribution to use
	:return: float (expected value)
	"""
	# get the appropriate function
	num_switch_func = num_switch_funcs[k, l]
	# get the interval of stability
	interval = stab_intervals[k, l]

	# get the right distribution function
	if not dist:
		my_std = (interval[1] - interval[0])/2.
		my_mean = 0
		a, b = (interval[0] - my_mean) / my_std, (interval[1] - my_mean) / my_std
		dist = st.truncnorm(a, b, 0, my_std)
	# Now do the integration
	exp_value = 0
	for (value, (start, end)) in num_switch_func:
		exp_value += value*(dist.cdf(end) - dist.cdf(start))
	return exp_value/float(n*n)


def num_switch_to_step(num_switch_funcs, intervals, k, l):
	"""
	Helper function to get input arguments for matplotlib step function plot

	:param num_switch_funcs: the full num switch functions in interval notation
	:param intervals: asympt stab intervals
	:param k: row you're looking at
	:param l: column you're looking at
	:return: (x,y) suitable for input to plt.step(x,y)
	"""
	num_switch_func = num_switch_funcs[k, l]
        #num_switch_func = num_switch_funcs["(" + str(k) + ", " + str(l) + ")"]
        # in R dict keys get converted to char
        #strr = "'(" + str(k) + ", " + str(l) + ")'"
        #num_switch_func = num_switch_funcs[k, l]
	start, stop = intervals[k, l]
	if not num_switch_func:
		x = [start, stop]
		y = [0, 0]
	else:
		num_switch_func_sorted = sorted(num_switch_func, key=lambda x: x[1][0])  # sort based off of increasing x-value
		val, (first_int_start, first_int_end) = num_switch_func_sorted[0]
		val, (last_int_start, last_int_end) = num_switch_func_sorted[-1]
		num_switch_func_sorted.insert(0, (-1, (start, 0)))
		num_switch_func_sorted.append((-1, (stop, 0)))
		#return num_switch_func_sorted
		for i in range(len(num_switch_func_sorted)):
			val, (int_start, int_stop) = num_switch_func_sorted[i]
			if int_start > 0:
				zero_start = num_switch_func_sorted[i - 1][1][1]
				zero_stop = num_switch_func_sorted[i][1][0]
				break
		if first_int_start > start:
			zero_start = start
		if last_int_end < stop:
			zero_stop = stop
		try:  # in case they weren't defined, num_switch doesn't cross zero
			zero_start
			zero_stop
		except NameError:
			zero_start = 0
			zero_stop = np.inf
			for temp_val, (temp_start, temp_end) in num_switch_func:
				if np.abs(temp_start) < zero_stop:
					zero_stop = temp_start
				if np.abs(temp_end) < zero_stop:
					zero_stop = temp_end

		full_switch_func = sorted(num_switch_func + [(0, (zero_start, zero_stop))], key=lambda x: x[1][0])
		# get rid of zero length intervals
		full_switch_func_clean = []
		for val, (int_start, int_stop) in full_switch_func:
			if int_stop - int_start != 0:
				full_switch_func_clean.append((val,(int_start, int_stop)))
		full_switch_func = full_switch_func_clean
		x = [i[1][0] for i in full_switch_func] + [full_switch_func[-1][1][1]]
		y = [full_switch_func[0][0]] + [i[0] for i in full_switch_func]
	return (x, y)

def import_matrix(file_path):
	"""
	Import a matrix, look for row/column labels, use pandas if they are, otherwise, use numpy

	:param file_path: input matrix file path
        :return: (matrix:numpy.array, row_labels:Array[string], column_labels:Array[string]
	"""
	if not os.path.exists(file_path):
		raise Exception("The file %s does not appear to exist.")
	try:  # csv plain
		A = np.loadtxt(file_path, delimiter=',')
		row_names = ['%d' % i for i in range(A.shape[0])]
		column_names = ['%d' % i for i in range(A.shape[1])]
		return (A, row_names, column_names)
	except ValueError:
		try:  # tsv plain
			A = np.loadtxt(file_path, delimiter='\t')
			row_names = ['%d' % i for i in range(A.shape[0])]
			column_names = ['%d' % i for i in range(A.shape[1])]
			return (A, row_names, column_names)
		except ValueError:
			try:  # csv with labels
				df = pd.read_csv(file_path, header=0, index_col=0)
				if df.empty:  # might instead be tsv
					df = pd.read_csv(file_path, header=0, index_col=0, sep='\t')
				A = df.values
				column_names = ['%s' % i for i in df.columns.values]
				row_names = ['%s' % i for i in df.index.values]
				return (A, row_names, column_names)
			except:
				raise Exception("Could not read in the file %s. Maybe not tsv or csv or otherwise mangled?" % file_path)


# This has been checked against Mathematica
# a, b = (-.5 - 0) / .75, (1 - 0) / .75
#dist = st.truncnorm(a, b, 0, .75)
#x = np.linspace(-.6, 1.1, 100)
#fig, ax = plt.subplots(1, 1)
#ax.plot(x, dist.pdf(x),'r-',label='norm pdf')
#fig.show()
# dist.pdf(0)

# stick these constants in for easy access later down the road/testing
Atri = np.array([[-0.237, -1, 0, 0], [0.1, -0.015, -1, 0], [0, 0.1, -0.015, -1], [0, 0, 0.1, -0.015]])
Aigp = np.array([[-0.237, -1, 0, 0], [0.1, -0.015, -1, -1], [0, 0.1, -0.015, -1], [0, .045, 0.1, -0.015]])
Atriinv = np.linalg.inv(Atri)
Aigpinv = np.linalg.inv(Aigp)

#######################################################################
# test cases
# test matrix
def fast_tests():
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
	interval = interval_of_stability(Aigp, Aigpinv, 0, 0)
	assert np.abs(interval[0] - -0.08298659074229853) < .001
	assert np.abs(interval[1] - 0.10901820241962716) < .001
	#interval = interval_of_stability_crawl(Aigp, Aigpinv, 0, 0, step_size=.0001)
	#assert np.abs(interval[0] - -0.08298659074229853) < .001
	#assert np.abs(interval[1] - 0.10901820241962716) < .001
	interval = interval_of_stability(Aigp, Aigpinv, 1, 0)
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


def test_num_switch_from_crit_eps():
	A = np.array([[-0.237, -1, 0, 0], [0.1, -0.015, -1, -1], [0, 0.1, -0.015, -1], [0, .045, 0.1, -0.015]])
	Ainv = np.linalg.inv(A)
	n = A.shape[0]
	crit_epsilon_array = np.zeros((n, n, n, n))
	for k in range(n):
		for l in range(n):
			for i in range(n):
				for j in range(n):
					crit_epsilon_array[k, l, i, j] = critical_epsilon(Ainv, k, l, i, j)

	stab_int_array = np.zeros((n, n, 2))
	for k in range(n):
		for l in range(n):
			if True:  # A[k, l] != 0:
				interval = interval_of_stability(A, Ainv, k, l, max_bound=10)
				stab_int_array[k, l, 0] = interval[0]
				stab_int_array[k, l, 1] = interval[1]

	# generated using TestingNewSwitchIntervalTestCases.nb
	known_answers =[[(1, (0.00801341, 0.109018))],
					[(1, (-9.999, -0.03))],
					[(1, (0.0694652, 1.99872))],
					[],
					[(1, (0.0040656, 10.))],
					[(2, (-9.999, -0.312)), (1, (-0.312, -0.013))],
					[(1, (0.0298758, 0.720876)), (2, (0.720876, 0.985876)), (4, (0.985876, 1.))],
					[(1, (-7.11695, -2.93395))],
					[],
					[(1, (0.272959, 0.919474))],
					[(3, (-9.999, -0.984)), (1, (-0.984, -0.692))],
					[(1, (0.720384, 0.967384)), (3, (0.967384, 0.985384)), (5, (0.985384, 1.))],
					[(1, (0.000361142, 0.0103611)), (2, (0.0103611, 0.0123611)), (3, (0.0123611, 0.014284))],
					[(4, (-0.044, -0.043)), (2, (-0.043, -0.031)), (1, (-0.031, -0.001))],
					[(1, (0.00357428, 0.272574)), (2, (0.272574, 0.63448))],
					[(3, (-9.999, -0.434)), (1, (-0.434, -0.312))]]
	iter = 0
	for k in range(n):
		for l in range(n):
			known_answer = known_answers[iter]
			res = num_switch_from_crit_eps(crit_epsilon_array, stab_int_array, k, l)
			# get rid of the super small intervals that mathematica misses
			new_res = []
			for index in range(len(res)):
				val, (start, stop) = res[index]
				if stop - start > 0.000001:
					new_res.append((val, (start, stop)))
			res = new_res
			res = sorted(res, key=lambda x: x[1][0])
			#print("k=%d, l=%d" % (k, l))
			for index in range(len(known_answer)):
				known_val, (known_start, known_end) = known_answer[index]
				test_val, (test_start, test_end) = res[index]
				assert known_val == test_val
				assert np.abs(known_start - test_start) < 0.02
				assert np.abs(known_end - test_end) < 0.02
			iter += 1

	# Tri-diagonal case
	A = np.array([[-0.237, -1, 0, 0], [0.1, -0.015, -1, 0], [0, 0.1, -0.015, -1], [0, 0, 0.1, -0.015]])
	Ainv = np.linalg.inv(A)
	n = A.shape[0]
	crit_epsilon_array = np.zeros((n, n, n, n))
	for k in range(n):
		for l in range(n):
			for i in range(n):
				for j in range(n):
					crit_epsilon_array[k, l, i, j] = critical_epsilon(Ainv, k, l, i, j)

	stab_int_array = np.zeros((n, n, 2))
	for k in range(n):
		for l in range(n):
			if True:  # A[k, l] != 0:
				interval = interval_of_stability(A, Ainv, k, l, max_bound=10)
				stab_int_array[k, l, 0] = interval[0]
				stab_int_array[k, l, 1] = interval[1]

	# generated using TestingNewSwitchIntervalTestCases.nb
	# note that some of these are not empty because I'm letting it perturb a zero entry
	known_answers = [[],
					[],
					[(2, (2.37059, 2.52559)), (3, (2.52559, 4.32057))],
					[(1, (0.356, 9.978)), (2, (9.978, 10.))],
					[],
					[],
					[],
					[(2, (0.150145, 0.301145)), (3, (0.301145, 4.37015)), (4, (4.37015, 4.52924))],
					[(2, (0.0243159, 0.0253159)), (3, (0.0253159, 0.0432057))],
					[],
					[],
					[],
					#[(1, (-0.00973428, 0.000265719)), (1, (0.0102657, 0.0123244))],  # step size was too large in mathematica
					[(4, (-0.0107243, -0.0103543)), (3, (-0.0103543, -0.0100143)), (2, (-0.0100143, -0.00997428)), (1, (-0.00997428, -0.000354281)), (1, (0.0100057, 0.0123244))],
					[(2, (0.00242148, 0.00342148)), (3, (0.00342148, 0.0452924))],
					[],
					[]]
	iter = 0
	for k in range(n):
		for l in range(n):
			#print("k=%d, l=%d" % (k, l))
			known_answer = known_answers[iter]
			res = num_switch_from_crit_eps(crit_epsilon_array, stab_int_array, k, l)
			# get rid of the super small intervals that mathematica misses
			new_res = []
			for index in range(len(res)):
				val, (start, stop) = res[index]
				if stop - start > 0.000001:
					new_res.append((val, (start, stop)))
			res = new_res
			res = sorted(res, key=lambda x: x[1][0])
			# print("k=%d, l=%d" % (k, l))
			for index in range(len(known_answer)):
				known_val, (known_start, known_end) = known_answer[index]
				test_val, (test_start, test_end) = res[index]
				#print('known: %f, test: %f' % (known_val, test_val))
				#print(res)
				assert known_val == test_val
				assert np.abs(known_start - test_start) < 0.02
				assert np.abs(known_end - test_end) < 0.02
			iter += 1


def test_exp_num_switch_from_crit_eps():

	Atri = np.array([[-0.237, -1, 0, 0], [0.1, -0.015, -1, 0], [0, 0.1, -0.015, -1], [0, 0, 0.1, -0.015]])
	Aigp = np.array([[-0.237, -1, 0, 0], [0.1, -0.015, -1, -1], [0, 0.1, -0.015, -1], [0, .045, 0.1, -0.015]])
	Atriinv = np.linalg.inv(Atri)
	Aigpinv = np.linalg.inv(Aigp)

	# IGP case
	A = Aigp
	Ainv = Aigpinv
	n = Aigp.shape[0]
	crit_epsilon_array = np.zeros((n, n, n, n))
	for k in range(n):
		for l in range(n):
			for i in range(n):
				for j in range(n):
					crit_epsilon_array[k, l, i, j] = critical_epsilon(Ainv, k, l, i, j)

	stab_int_array = np.zeros((n, n, 2))
	for k in range(n):
		for l in range(n):
			if A[k, l] != 0:
				interval = interval_of_stability(A, Ainv, k, l, max_bound=10)
				stab_int_array[k, l, 0] = interval[0]
				stab_int_array[k, l, 1] = interval[1]

	num_switch_funcs = dict()
	for k in range(n):
		for l in range(n):
			if True:  # A[k, l] != 0:
				num_switch_func = num_switch_from_crit_eps(crit_epsilon_array, stab_int_array, k, l)
				num_switch_funcs[k, l] = num_switch_func

	exp_value = exp_num_switch_from_crit_eps(n, 3, 2, num_switch_funcs, stab_int_array)
	assert np.abs(exp_value - 0.081) < 0.01
	exp_value = exp_num_switch_from_crit_eps(n, 1, 2, num_switch_funcs, stab_int_array)
	assert np.abs(exp_value - 0.048) < 0.01
	exp_value = exp_num_switch_from_crit_eps(n, 0, 0, num_switch_funcs, stab_int_array)
	assert np.abs(exp_value - 0.032) < 0.01
	exp_value = exp_num_switch_from_crit_eps(n, 1, 1, num_switch_funcs, stab_int_array)
	assert np.abs(exp_value - 0.12) < 0.01
	exp_value = exp_num_switch_from_crit_eps(n, 2, 2, num_switch_funcs, stab_int_array)
	assert np.abs(exp_value - 0.16) < 0.01
	exp_value = exp_num_switch_from_crit_eps(n, 3, 3, num_switch_funcs, stab_int_array)
	assert np.abs(exp_value - 0.17) < 0.01


def run_all_tests():
	"""
	Runs all the tests.

	:return: None
	"""
	test_num_switch_from_crit_eps()
	fast_tests()
	test_exp_num_switch_from_crit_eps()