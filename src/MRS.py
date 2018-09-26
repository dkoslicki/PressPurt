# this script will implement equation 3.21
import numpy as np
A = np.array([[-0.237,-1,0,0],[0.1,-0.015,-1,0],[0,0.1,-0.015,-1],[0,0,0.1,-0.015]])  # initial input matrix
M_i, M_j = np.nonzero(A)  # location of nonzero
M = len(M_i)  # number of nonzeros
Ainv = np.linalg.inv(A)  # matrix inverse of A
Ainv_T = np.sum(np.abs(Ainv))
n = A.shape[0]  # square matrix, this is the dimension
assert n != 0  # this should be nonzero, otherwise things below are undefined
# now for the main computation
outer_sum = 1  # numerator is 1
for k, l in zip(M_i, M_j):
	outer_sum /= np.abs(Ainv[l, k])
	for i in range(n):
		inner_sum1 = np.abs(Ainv[i, k])
		inner_sum2 = 0
		for j in range(n):
			inner_sum2 += np.abs(Ainv[l, j])
		inner_sum1 *= inner_sum2
	outer_sum *= inner_sum1

MRS = 1/(M * Ainv_T) * outer_sum
