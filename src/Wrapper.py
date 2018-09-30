# This script will run through all the other functions, generating the results and storing them
import numpy as np
import scipy as sp
import scipy.stats as st
import MRS
import SS
import NumSwitch
import matplotlib.pyplot as plt

######################
# input matrix
A = np.array([[-0.237, -1, 0, 0], [0.1, -0.015, -1, -1], [0, 0.1, -0.015, -1], [0, .045, 0.1, -0.015]])
Ainv = np.linalg.inv(A)
######################
# Generate figure 2
k = 3
l = 2
padding = .2
interval = NumSwitch.interval_of_stability(A, Ainv, k, l, max_bound=10, num_sample=1000)
x_range = np.linspace(interval[0] - padding, interval[1] + padding, 500)
ns_values = [NumSwitch.NS(Ainv, eps, k, l) for eps in x_range]
# for the fun of it, overlay the distribution too
my_std = (interval[1] - interval[0])/2.
my_mean = 0
a, b = (interval[0] - my_mean) / my_std, (interval[1] - my_mean) / my_std
dist = st.truncnorm(a, b, 0, my_std)
dist_vals = [dist.pdf(eps) for eps in x_range]

fig, ax1 = plt.subplots()
ax1.plot(x_range, ns_values, 'b-')
ax1.set_xlabel('Epsilon value')
# Make the y-axis label, ticks and tick labels match the line color.
ax1.set_ylabel('NumSwitch(eps, %d, %d)' % (k + 1, l + 1), color='b')
ax1.tick_params('y', colors='b')

ax2 = ax1.twinx()
ax2.plot(x_range, dist_vals, 'tab:gray')
ax2.set_ylabel('Probability', color='tab:gray')
ax2.tick_params('y', colors='tab:gray')
ax2.fill(x_range, dist_vals, 'tab:gray', alpha=0.5)

fig.tight_layout()
plt.show()
#plt.figure()
#plt.plot(x_range, ns_values)
#plt.show()
