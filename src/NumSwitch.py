# This script will implement equations 3.4, 3.5, and 3.6
import numpy as np


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
