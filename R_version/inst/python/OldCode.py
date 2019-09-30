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


def num_switched_by(Ainv, eps, i, j):
	"""
	When perturbing each of the k,l entries by eps, how many times does the (i,j)
	entry switch signs?
	TODO: I don't know about this one, since the eps should be different for
	TODO: each of the (k, l) entries...
	:param Ainv: inverse input matrix, numpy array
	:param eps: perturbation size (scalar)
	:param i: index
	:param j: inex
	:return: Natural number
	"""
	n = Ainv.shape[0]
	if any([index >= n for index in [i, j]]):
		raise Exception("Matrix is only %dx%d, invalid choice of subscripts: %d, %d." % (n, n, k, l))
	ns_sum = 0
	for k in range(n):
		for l in range(n):
			ns_sum += ind_switch(Ainv, eps, i, j, k, l)
	return ns_sum