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
		# TODO: put custom pdf in here



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


exp_num_switch(Atri, Atriinv, 0, 0)

A = Atri
Ainv = Atriinv
k = 0
l = 0

class my_pdf(st.rv_continuous):
	c = A[k, l]
	if A[k, l] > 0:
		max_pert = np.inf
		min_pert = -A[k, l]
	else:
		max_pert = -A[k, l]
		min_pert = -np.inf

	def _pdf(self, x):
		a = self.min_pert
		b = self.max_pert
		c = self.c
		if a < x and x <= b:
			# this is the result of Mathematica: PDF[TruncatedDistribution[{a, b}, NormalDistribution[c, 1]], x]
			numerator = np.exp(-1 / 2. * (-c + x) ** 2)
			denom = np.sqrt(2 * np.pi) * (
						-1 / 2. * special.erfc((-a + c) / (np.sqrt(2))) + 1 / 2. * special.erfc((-b + c) / np.sqrt(2)))
			value = numerator / float(denom)
			return value
		else:
			return 0

my_cv = my_pdf(a=my_pdf.min_pert, b=my_pdf.max_pert, name='my_pdf')

print(integrate.quad(my_pdf, -np.inf, np.inf))