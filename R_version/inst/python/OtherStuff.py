import numpy as np
import scipy as sp

# create random asymptotically stable matrices
n = 20
eig = np.diag(-np.random.rand(n))
v = sp.linalg.orth(2*(np.random.rand(n, n)-.5))
A = np.transpose(v).dot(eig.dot(v))
np.savetxt('/home/dkoslicki/Desktop/PressPurt/ExampleJacobians/big_test.csv', A, delimiter=',')
#np.savetxt('/home/dkoslicki/Dropbox/Repositories/PressPurt/ExampleJacobians/big_test.csv', A, delimiter=',')
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
		if True:#A[k, l] != 0:
			interval = NumSwitch.interval_of_stability(A, Ainv, k, l, max_bound=10, num_sample=1000)
			stab_int_array[k, l, 0] = interval[0]
			stab_int_array[k, l, 1] = interval[1]

reload(NumSwitch)
res = NumSwitch.num_switch_from_crit_eps(crit_epsilon_array, stab_int_array, 0, 0)
print(sorted(res, key=lambda x: x[1][0]))

num_switch_funcs = dict()
for k in range(n):
	for l in range(n):
		if True:#A[k, l] != 0:
			num_switch_func = NumSwitch.num_switch_from_crit_eps(crit_epsilon_array, stab_int_array, k, l)
			num_switch_funcs[k, l] = num_switch_func




test=[[(1, (0.00801341, 0.109018))],
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




# Let's try a binary search for the interval of stability

def largest_root(A, k, l, eps):
	zero_matrix = np.zeros(A.shape)
	zero_matrix[k, l] = eps
	[s, _] = np.linalg.eig(A + zero_matrix)
	return max([np.real(i) for i in s])

A = NumSwitch.Aigp
Ainv = NumSwitch.Aigpinv
k = 0
l = 0
max_bound = 10
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

# start on the positive side
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