import numpy as np
import scipy as sp
n = 8
eig = np.diag(-np.random.rand(n))
v = sp.linalg.orth(2*(np.random.rand(n, n)-.5))
A = np.transpose(v).dot(eig.dot(v))
