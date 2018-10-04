import numpy as np
import scipy as sp

# create random asymptotically stable matrices
n = 15
eig = np.diag(-np.random.rand(n))
v = sp.linalg.orth(2*(np.random.rand(n, n)-.5))
A = np.transpose(v).dot(eig.dot(v))
np.savetxt('/home/dkoslicki/Desktop/PressPurtCoreAlg/ExampleJacobians/big_test.csv', A, delimiter=',')
Ainv = np.linalg.inv(A)


# Test the num switch stuff
import NumSwitch
reload(NumSwitch)
A = NumSwitch.Aigp
Ainv = NumSwitch.Aigpinv
n = Aigp.shape[0]
crit_epsilon_array = np.zeros((n, n, n, n))
for k in range(n):
	for l in range(n):
		for i in range(n):
			for j in range(n):
				crit_epsilon_array[k, l, i, j] = NumSwitch.critical_epsilon(Ainv, k, l, i, j)

stab_int_array = np.zeros((n, n, 2))
for k in range(n):
	for l in range(n):
		if A[k, l] != 0:
			interval = NumSwitch.interval_of_stability(A, Ainv, k, l, max_bound=10, num_sample=1000)
			stab_int_array[k, l, 0] = interval[0]
			stab_int_array[k, l, 1] = interval[1]

reload(NumSwitch)
print(NumSwitch.num_switch_from_crit_eps(crit_epsilon_array, stab_int_array, 0, 0))

