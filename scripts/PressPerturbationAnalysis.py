import argparse
import numpy as np
import scipy.stats as st
import MRS
import SS
import NumSwitch
import matplotlib.pyplot as plt
import os

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description="This script all indices from the Koslicki & Novak (2018) paper", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('input_file', type=str, help="Input comma separated file for the jacobian matrix.")
	parser.add_argument('--max_bound', type=int, help="some of the matrices are unbounded stable towards one end, this is the limit the user imposes", default=10)
	parser.add_argument('--num_sample', type=int, help="number of points to sample when looking for the region of asymptotic stability of the matrix", default=1000)
	parser.add_argument('--num_points_plot', type=int, help="number of points to plot in the first figure", default=500)
	parser.add_argument('--num_iterates', help="number of Monte-Carlo points to sample for the SS index", default=5000)

	# read in the arguments
	args = parser.parse_args()
	input_file = os.path.abspath(args.input_file)
	if not os.path.exists(input_file):
		raise Exception("It appears that the file %s does not exits." % input_file)
	max_bound = args.max_bound
	assert max_bound > 0
	num_sample = args.num_sample
	assert num_sample > 10
	num_points_plot = args.num_points_plot
	assert num_points_plot > 10
	num_iterates = args.num_iterates
	assert num_iterates > 100

	# read in the input matrix
	A = np.loadtxt(input_file, "rb", delimiter=",")
	Ainv = np.linalg.inv(A)
	m, n = A.shape
	######################
	# Generate figure 2
	k = 3
	l = 2
	padding = .2
	interval = NumSwitch.interval_of_stability(A, Ainv, k, l, max_bound=max_bound, num_sample=num_sample)
	x_range = np.linspace(interval[0] - padding, interval[1] + padding, num_points_plot)
	ns_values = [NumSwitch.NS(Ainv, eps, k, l) for eps in x_range]
	# for the fun of it, overlay the distribution too
	my_std = (interval[1] - interval[0]) / 2.
	my_mean = 0
	a, b = (interval[0] - my_mean) / my_std, (interval[1] - my_mean) / my_std
	dist = st.truncnorm(a, b, 0, my_std)
	dist_vals = [dist.pdf(eps) for eps in x_range]

	# plot both simultaneously on the same graph (two y-axis plot)
	fig, ax1 = plt.subplots()
	ax1.plot(x_range, ns_values, 'b-')
	ax1.set_xlabel('Epsilon value')
	ax1.set_ylabel('NumSwitch(eps, %d, %d)' % (k + 1, l + 1), color='b')
	ax1.tick_params('y', colors='b')

	ax2 = ax1.twinx()
	ax2.plot(x_range, dist_vals, 'tab:gray')
	ax2.set_ylabel('Probability', color='tab:gray')
	ax2.tick_params('y', colors='tab:gray')
	ax2.fill(x_range, dist_vals, 'tab:gray', alpha=0.5)

	fig.tight_layout()
	plt.draw()
	plt.pause(0.01)
	# plt.show(block=False)
	# plt.figure()
	# plt.plot(x_range, ns_values)
	# plt.show()

	#####################
	# Get the expected number of sign switches, in a table
	exp_num_switch_array = np.zeros((m, n))
	for i in range(m):
		for j in range(n):
			if A[i, j] != 0:  # only perturb non-zero entries
				exp_num_switch_array[i, j] = NumSwitch.exp_num_switch(A, Ainv, i, j, num_sample=num_sample, dist=None,
																	  interval=None)
			else:
				exp_num_switch_array[i, j] = 0

	#####################
	# Generate figure 3
	fig, ax = plt.subplots()
	# im = ax.imshow(exp_num_switch_array, cmap=plt.get_cmap('seismic'))
	im = ax.imshow(exp_num_switch_array, cmap=plt.get_cmap('YlOrBr'))
	# We want to show all ticks...
	ax.set_xticks(np.arange(m))
	ax.set_yticks(np.arange(n))
	# ... and label them with the respective list entries
	ax.set_xticklabels(['%d' % (i + 1) for i in range(n)])
	ax.set_yticklabels(['%d' % (j + 1) for j in range(m)])
	ax.set_xlabel('l')
	ax.set_ylabel('k')
	# Rotate the tick labels and set their alignment.
	plt.setp(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")

	# Loop over data dimensions and create text annotations.
	for i in range(m):
		for j in range(n):
			text = ax.text(j, i, '%.2f%%' % (100 * exp_num_switch_array[i, j]), ha="center", va="center", color="k")

	ax.set_title("Expected number of sign switches")
	fig.tight_layout()
	plt.draw()
	plt.pause(0.01)
	# plt.show(block=False)

	######################
	# Compute MRS

	print('Quantitative sensitivity of A: %f' % MRS.MRS(A))

	######################
	# Recreate figure 5
	quant_sens_values = np.zeros((m, n))
	for i in range(m):
		for j in range(n):
			if A[i, j] != 0:  # only perturb non-zero entries
				quant_sens_values[i, j] = MRS.quant_sens(A, i, j)
			else:
				quant_sens_values[i, j] = 0

	#####################
	# Generate figure 3
	fig, ax = plt.subplots()
	# im = ax.imshow(exp_num_switch_array, cmap=plt.get_cmap('seismic'))
	# im = ax.imshow(quant_sens_values, cmap=plt.get_cmap('YlOrBr'))
	im = ax.imshow(quant_sens_values, cmap=plt.get_cmap('Wistia'))
	# We want to show all ticks...
	ax.set_xticks(np.arange(m))
	ax.set_yticks(np.arange(n))
	# ... and label them with the respective list entries
	ax.set_xticklabels(['%d' % (i + 1) for i in range(n)])
	ax.set_yticklabels(['%d' % (j + 1) for j in range(m)])
	ax.set_xlabel('l')
	ax.set_ylabel('k')
	# Rotate the tick labels and set their alignment.
	plt.setp(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")
	# Loop over data dimensions and create text annotations.
	for i in range(m):
		for j in range(n):
			text = ax.text(j, i, '%.2f' % quant_sens_values[i, j], ha="center", va="center", color="k")

	ax.set_title("Quantitative sensitivity")
	fig.tight_layout()
	plt.draw()
	plt.pause(0.01)
	# plt.show(block=False)

	#####################
	# Section 3.4, perturbing multiple entries
	interval_length = 0.01
	ss_val = SS.SS(A, num_iterates=num_iterates, interval_length=interval_length)
	print('The percent of uniform perturbations over an interval of length %.2f is: %f' % (interval_length, ss_val))

	input("Press any key to quit")
