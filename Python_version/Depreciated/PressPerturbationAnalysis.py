import argparse
import numpy as np
import scipy.stats as st
import os
import sys
import matplotlib.pyplot as plt
import timeit
# import stuff in the src folder
try:
	import MRS
except ImportError:
	sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'src'))
	import MRS
try:
	import SS
except ImportError:
	sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'src'))
	import SS
try:
	import NumSwitch
except ImportError:
	sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'src'))
	import NumSwitch

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description="This script all indices from the Koslicki & Novak (2018) paper", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('input_file', type=str, help="Input comma separated file for the jacobian matrix.")
	parser.add_argument('-o', '--output_folder', type=str, help="Output folder")
	parser.add_argument('-p', '--prefix', help="Prefix of output files, if you so choose.")
	parser.add_argument('-a', '--all_numswitch_plots', action='store_true', help="Include this flag if you want all the num switch plots (it coul be large)")
	parser.add_argument('-l', '--list_of_numswitch_to_plot', nargs='+', help="List of entries you want visualized with num switch. Eg. -l 1 1 1 2 to plot the (1,1) and (1,2) entries.")
	parser.add_argument('--max_bound', type=int, help="some of the matrices are unbounded stable towards one end, this is the limit the user imposes", default=10)
	parser.add_argument('--num_sample', type=int, help="number of points to sample when looking for the region of asymptotic stability of the matrix", default=100)
	parser.add_argument('--num_points_plot', type=int, help="number of points to plot in the first figure", default=500)
	parser.add_argument('--num_iterates', help="number of Monte-Carlo points to sample for the SS index", default=5000)
	parser.add_argument('-asf', '--asymp_stability_file', type=str, help="location of where to save/load the intervals of asymptotic stability. File extension .npy must be used.")
	parser.add_argument('-sasf', '--save_asymp_stability_file', action='store_true', help="Flag to include if you want to save the intervals of stability. This will clobber the file in -asf.")
	parser.add_argument('-ss', '--run_global_sign_sensitivity', action='store_true', help="Flag to include the calculation of the global sign sensitivity.")
	parser.add_argument('--ss_interval_length', type=float, help="Interval length over wich global sign sensitivity will be calculated.", default=0.01)

	# read in the arguments
	args = parser.parse_args()
	input_file = os.path.abspath(args.input_file)
	output_folder = args.output_folder
	if output_folder:
		output_folder = os.path.abspath(output_folder)
	prefix = args.prefix
	all_numswitch_plots = args.all_numswitch_plots
	if args.list_of_numswitch_to_plot:
		list_of_numswitch_to_plot = [int(i)-1 for i in list(args.list_of_numswitch_to_plot)]
	else:
		list_of_numswitch_to_plot = []
	max_bound = int(args.max_bound)
	num_sample = int(args.num_sample)
	num_points_plot = int(args.num_points_plot)
	num_iterates = int(args.num_iterates)
	save_asymp_stability_file = args.save_asymp_stability_file
	asymp_stability_file = args.asymp_stability_file
	if asymp_stability_file:
		asymp_stability_file = os.path.abspath(asymp_stability_file)
	run_global_sign_sensitivity = args.run_global_sign_sensitivity
	interval_length = float(args.ss_interval_length)

	# check for read-writability of files
	if asymp_stability_file and not save_asymp_stability_file:
		if not os.access(asymp_stability_file, os.R_OK):
			raise Exception("The provided asymptotic stability file %s does not exist/is not readable." % asymp_stability_file)
	if save_asymp_stability_file:
		if not os.access(os.path.dirname(asymp_stability_file), os.W_OK):
			raise Exception("The provided directory %s is not writable." % os.path.dirname(asymp_stability_file))
	if not os.path.exists(input_file):
		raise Exception("It appears that the input file %s does not exits." % input_file)

	# check for sanity of input parameters
	if not max_bound > 0:
		raise Exception("max_bound must be larger than 0; provided value: %d." % max_bound)
	if not num_sample > 10:
		raise Exception("num_sample must be larger than 10; provided value: %d." % num_sample)
	if not num_points_plot > 10:
		raise Exception("num_points_plot must be larger than 10; provided value: %d." % num_points_plot)
	if not num_iterates > 100:
		raise Exception("num_iterates must be larger than 100; provided value: %d." % num_iterates)
	if not interval_length > 0:
		raise Exception("ss_interval_length must be larget than 0; provided value: %f" % interval_length)
	if all_numswitch_plots and list_of_numswitch_to_plot:
		raise Exception("The flags -a and -l are mutually exclusive.")
	if len(list_of_numswitch_to_plot) % 2 != 0:
		raise Exception("Only an even length list can be passed to -l.")

	# read in the input matrix
	A = np.loadtxt(input_file, delimiter=",")
	Ainv = np.linalg.inv(A)
	m, n = A.shape

	# if appropriate, get the indices to plot
	indices_to_plot = []
	# quick way to take entries two at a time
	it = iter(list_of_numswitch_to_plot)
	for i, j in zip(*[it]*2):
		indices_to_plot.append((i, j))

	# make sure the original matrix is itself asymptotically stable
	if not SS.is_stable(A):
		raise Exception("Sorry, the input matrix is not stable itself (all eigenvalues must have negative real part). Please try again.")

	# compute/load the intervals of stability
	if asymp_stability_file and not save_asymp_stability_file:
		intervals = np.load(asymp_stability_file)
	else:  # you need to actually compute it
		intervals = np.zeros((m, n, 2))
		for k in range(m):
			for l in range(n):
				if A[k, l] != 0:
					intervals[k, l, :] = NumSwitch.interval_of_stability(A, Ainv, k, l, max_bound=max_bound)

	# save these if that was asked for
	if save_asymp_stability_file:
		print("Saving asymptotic stability to: %s" % asymp_stability_file)
		np.save(asymp_stability_file, intervals)

	######################
	# Generate figure 2
	# let's do a mxn grid of these figures
	# TODO: should really make this into its own function
	if all_numswitch_plots:
		padding = .2
		big_fig, axarr = plt.subplots(m, n)
		big_fig.suptitle("Number of mis-predictions versus perturbation value, \n overlaid with distribution over stable perturbation values")
		for k in range(m):
			for l in range(n):
				if A[k, l] != 0:
					interval = intervals[k, l, :]
					x_range = np.linspace(interval[0] - padding, interval[1] + padding, num_points_plot)
					ns_values = [NumSwitch.NS(Ainv, eps, k, l) for eps in x_range]
					# for the fun of it, overlay the distribution too
					my_std = (interval[1] - interval[0]) / 2.
					my_mean = 0
					a, b = (interval[0] - my_mean) / my_std, (interval[1] - my_mean) / my_std
					dist = st.truncnorm(a, b, 0, my_std)
					dist_vals = [dist.pdf(eps) for eps in x_range]

					# plot both simultaneously on the same graph (two y-axis plot)
					ax1 = axarr[k, l]
					ax1.title.set_text('(%d,%d) entry' % (k+1, l+1))
					ax1.plot(x_range, ns_values, 'b-')
					#ax1.set_xlabel('Epsilon value')
					#ax1.set_ylabel('NumSwitch(eps, %d, %d)' % (k + 1, l + 1), color='b')
					ax1.tick_params('y', colors='b')

					ax2 = ax1.twinx()
					ax2.plot(x_range, dist_vals, 'tab:gray')
					#ax2.set_ylabel('PDF value', color='tab:gray')
					ax2.tick_params('y', colors='tab:gray')
					ax2.fill(x_range, dist_vals, 'tab:gray', alpha=0.5)
				else:
					axarr[k, l].axis('off')  # don't show the ones we are not perturbing
		plt.tight_layout(pad=0.1, w_pad=.1, h_pad=.9)
		big_fig.text(0.5, 0.01, 'Epsilon value', ha='center', va='center')
		big_fig.text(0.03, 0.5, 'Number of incorrect predictions', ha='center', va='center', rotation='vertical', color='b')
		big_fig.text(.99, 0.5, 'Probability density', ha='center', va='center', rotation='vertical', color='tab:gray')
		plt.subplots_adjust(top=.9)
		plt.draw()
		plt.pause(0.01)
	elif indices_to_plot:  # you only want to plot individual entries
		padding = .2
		grid_size = int(np.ceil(np.sqrt(len(indices_to_plot))))
		big_fig, axarr = plt.subplots(grid_size, grid_size)
		big_fig.suptitle(
			"Number of mis-predictions versus perturbation value, \n overlaid with distribution over stable perturbation values")
		num_plotted = 0
		for k,l in indices_to_plot:
			num_plotted += 1
			if A[k, l] != 0:
				interval = intervals[k, l, :]
				x_range = np.linspace(interval[0] - padding, interval[1] + padding, num_points_plot)
				ns_values = [NumSwitch.NS(Ainv, eps, k, l) for eps in x_range]
				# for the fun of it, overlay the distribution too
				my_std = (interval[1] - interval[0]) / 2.
				my_mean = 0
				a, b = (interval[0] - my_mean) / my_std, (interval[1] - my_mean) / my_std
				dist = st.truncnorm(a, b, 0, my_std)
				dist_vals = [dist.pdf(eps) for eps in x_range]

				# plot both simultaneously on the same graph (two y-axis plot)
				ax1 = axarr[k, l]
				ax1.title.set_text('(%d,%d) entry' % (k + 1, l + 1))
				ax1.plot(x_range, ns_values, 'b-')
				# ax1.set_xlabel('Epsilon value')
				# ax1.set_ylabel('NumSwitch(eps, %d, %d)' % (k + 1, l + 1), color='b')
				ax1.tick_params('y', colors='b')

				ax2 = ax1.twinx()
				ax2.plot(x_range, dist_vals, 'tab:gray')
				# ax2.set_ylabel('PDF value', color='tab:gray')
				ax2.tick_params('y', colors='tab:gray')
				ax2.fill(x_range, dist_vals, 'tab:gray', alpha=0.5)
			else:
				axarr[k, l].axis('off')  # don't show the ones we are not perturbing
		axes_flat = axarr.flatten()
		for i in range(grid_size**2-num_plotted, grid_size**2):
			axes_flat[i].axis('off')
		plt.tight_layout(pad=0.1, w_pad=.1, h_pad=.9)
		big_fig.text(0.5, 0.01, 'Epsilon value', ha='center', va='center')
		big_fig.text(0.03, 0.5, 'Number of incorrect predictions', ha='center', va='center', rotation='vertical',
					 color='b')
		big_fig.text(.99, 0.5, 'Probability density', ha='center', va='center', rotation='vertical', color='tab:gray')
		plt.subplots_adjust(top=.9)
		plt.draw()
		plt.pause(0.01)
	#####################
	# Get the expected number of sign switches, in a table
	exp_num_switch_array = np.zeros((m, n))
	for i in range(m):
		for j in range(n):
			if A[i, j] != 0:  # only perturb non-zero entries
				exp_num_switch_array[i, j] = NumSwitch.exp_num_switch(A, Ainv, i, j, num_sample=num_sample, dist=None, interval=intervals[i, j, :])
			else:
				exp_num_switch_array[i, j] = 0

	if output_folder:
		np.savetxt(os.path.join(output_folder, prefix + '_Expected_num_switch_array.csv'), exp_num_switch_array, delimiter=',')

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
	if m <= 10:
		for i in range(m):
			for j in range(n):
				text = ax.text(j, i, '%.2f%%' % (100 * exp_num_switch_array[i, j]), ha="center", va="center", color="k")
	else:
		fig.colorbar(im)

	ax.set_title("Expected number of mis-predictions when perturbing the (k,l) entry.")
	fig.tight_layout()
	plt.draw()
	plt.pause(0.01)

	num_non_zero = len(np.where(exp_num_switch_array)[0])
	ave_expected_num_sign_switches = exp_num_switch_array.sum()/float(num_non_zero)
	print("Average expected percentage of mis-predictions (perturbing each edge individually): %f%% (i.e. %f entries)" % (100*ave_expected_num_sign_switches, ave_expected_num_sign_switches*m*n))
	# TODO: check with Mark if I should be multiplying by n^2 or num_non_zero

	######################
	# Compute MRS
	print('Average quantitative sensitivity of A when perturbing each each individually by an arbitrarily large amount : %f' % MRS.MRS(A))

	######################
	# Recreate figure 5
	quant_sens_values = np.zeros((m, n))
	for i in range(m):
		for j in range(n):
			if A[i, j] != 0:  # only perturb non-zero entries
				quant_sens_values[i, j] = MRS.quant_sens(A, i, j)
			else:
				quant_sens_values[i, j] = 0

	if output_folder:
		np.savetxt(os.path.join(output_folder, prefix + '_Quantitative_sensitivity_array.csv'), quant_sens_values, delimiter=',')


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
	if m <= 10:
		for i in range(m):
			for j in range(n):
				text = ax.text(j, i, '%.2f' % quant_sens_values[i, j], ha="center", va="center", color="k")
	else:
		fig.colorbar(im)

	ax.set_title("Quantitative sensitivity\n when perturbing the (k,l) entry by an arbitrarily large value.")
	fig.tight_layout()
	plt.draw()
	plt.pause(0.01)

	#####################
	# Section 3.4, perturbing multiple entries
	if run_global_sign_sensitivity:
		ss_val = SS.SS(A, num_iterates=num_iterates, interval_length=interval_length)
		print('The percent of multi-perturbation space (over an interval of length %.2f) that causes some qualitative mis-prediction is: %f' % (interval_length, ss_val))

	input("Press any key to quit")
