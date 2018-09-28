# this will implement equation 3.40 and 3.42 in the manuscript
# volume of sign sensitivity region when perturbing multiple entries
import numpy as np
import sympy as sp
import scipy.stats as st
from NumSwitch import interval_of_stability
#entries_to_perturb = np.ones((4,4))
entries_to_perturb = np.array([[1,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]])
Aigp = np.array([[-0.237, -1, 0, 0], [0.1, -0.015, -1, -1], [0, 0.1, -0.015, -1], [0, .045, 0.1, -0.015]])
Aigpinv = np.linalg.inv(Aigp)
pert_locations_i, pert_locations_j = np.where(entries_to_perturb)
symbol_string = ""
for i,j in zip(pert_locations_i, pert_locations_j):
	symbol_string += "eps_%d_%d" % (i, j)
symbol_tup = sp.symbols(symbol_string)

# Next, get the intervals of stability for each of these
