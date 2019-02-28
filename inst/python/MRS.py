# this script will implement equation 3.21
import numpy as np


def MRS(A):
	"""
	This function implements equation 3.21 from the manuscript

	:param A: a non-empty matrix
	:return: MRS(A) (a scalar)
	"""
	M_i, M_j = np.nonzero(A)  # location of nonzero
	M = len(M_i)  # number of nonzeros
	Ainv = np.linalg.inv(A)  # matrix inverse of A
	Ainv_T = np.sum(np.abs(Ainv))
	n = A.shape[0]  # square matrix, this is the dimension
	assert n != 0  # this should be nonzero, otherwise things below are undefined
	# now for the main computation
	outer_sum = 0  # numerator is 1
	for k, l in zip(M_i, M_j):
		if Ainv[l, k] != 0:
			outer_sum_const = 1 / np.abs(Ainv[l, k])
			inner_sum1 = 0
			for i in range(n):
				inner_sum1_const = np.abs(Ainv[i, k])
				inner_sum2 = 0
				for j in range(n):
					inner_sum2 += np.abs(Ainv[l, j])
				inner_sum1 += inner_sum1_const * inner_sum2
			outer_sum += outer_sum_const * inner_sum1

	MRS = 1/(M * Ainv_T) * outer_sum
	return MRS

def quant_sens(A, k ,l):
	"""
	Quantitative sensitivity of a matrix (as epsilon goes to infinity). See equation 3.20

	:param A: input numpy matrix
	:param k: row index (int)
	:param l: column index (int)
	:return: float
	"""
	M_i, M_j = np.nonzero(A)  # location of nonzero
	M = len(M_i)  # number of nonzeros
	Ainv = np.linalg.inv(A)  # matrix inverse of A
	Ainv_T = np.sum(np.abs(Ainv))
	n = A.shape[0]  # square matrix, this is the dimension
	assert n != 0  # this should be nonzero, otherwise things below are undefined
	# now for the main computation
	outer_sum = 0  # numerator is 1
	if Ainv[l, k] != 0:  # if the entry is zero, can't do anything with it
		outer_sum_const = 1 / np.abs(Ainv[l, k])
		inner_sum1 = 0
		for i in range(n):
			inner_sum1_const = np.abs(Ainv[i, k])
			inner_sum2 = 0
			for j in range(n):
				inner_sum2 += np.abs(Ainv[l, j])
			inner_sum1 += inner_sum1_const * inner_sum2
		outer_sum += outer_sum_const * inner_sum1

	return outer_sum / float(Ainv_T)


# test cases
def tests():
	"""
	Run all the test cases.

	:return: None
	"""
	Atri = np.array([[-0.237, -1, 0, 0], [0.1, -0.015, -1, 0], [0, 0.1, -0.015, -1], [0, 0, 0.1, -0.015]])
	Aigp = np.array([[-0.237, -1, 0, 0], [0.1, -0.015, -1, -1], [0, 0.1, -0.015, -1], [0, .045, 0.1, -0.015]])
	assert np.abs(MRS(Atri) - 3.52748) < .00001
	assert np.abs(MRS(Aigp) - 1.62801) < .00001
	assert np.abs(quant_sens(Aigp, 0, 0) - 1.08) < 0.01
	assert np.abs(quant_sens(Aigp, 0, 1) - .8) < 0.01
	assert np.abs(quant_sens(Aigp, 2, 3) - 8.53) < 0.01
	assert np.abs(quant_sens(Atri, 0, 0) - 3.43) < 0.01
	assert np.abs(quant_sens(Atri, 1, 2) - 24.5) < 0.1



