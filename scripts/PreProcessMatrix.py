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
	parser.add_argument('output_folder', type=str, help="Output folder. A number of files will be created in the form 'output_folder/<prefix>_*.npy'")
	parser.add_argument('-p', '--prefix', help="Prefix of output files, if you so choose.", default=None)
	parser.add_argument('-m', '--max_bound', type=int, help="some of the matrices are unbounded stable towards one end, this is the limit the user imposes", default=10)
	parser.add_argument('-n', '--num_sample', type=int, help="number of points to sample when looking for the region of asymptotic stability of the matrix", default=1000)
	parser.add_argument('-z', '--zero_perturb', action='store_true', help="Flag to indicated you want to pertub the zero entries.", default=False)

	# read in the arguments
	args = parser.parse_args()
	input_file = os.path.abspath(args.input_file)
	output_folder = args.output_folder
	if output_folder:
		output_folder = os.path.abspath(output_folder)
	prefix = args.prefix
	max_bound = int(args.max_bound)
	num_sample = int(args.num_sample)
	pert_zero = args.zero_perturb
	if not os.access(output_folder, os.W_OK):
		raise Exception("The provided directory %s is not writable." % output_folder)

	if prefix:
		asymp_stab_file = os.path.join(output_folder, prefix + "_asymptotic_stability.npy")
		num_switch_file = os.path.join(output_folder, prefix + "_num_switch_funcs.npy")
	else:
		asymp_stab_file = os.path.join(output_folder, "asymptotic_stability.npy")
		num_switch_file = os.path.join(output_folder, "num_switch_funcs.npy")

	# check for sanity of input parameters
	if not max_bound > 0:
		raise Exception("max_bound must be larger than 0; provided value: %d." % max_bound)
	if not num_sample > 10:
		raise Exception("num_sample must be larger than 10; provided value: %d." % num_sample)

	# read in the input matrix
	A = np.loadtxt(input_file, delimiter=",")
	Ainv = np.linalg.inv(A)
	m, n = A.shape

	# make sure the original matrix is itself asymptotically stable
	if not SS.is_stable(A):
		raise Exception("Sorry, the input matrix is not stable itself (all eigenvalues must have negative real part). Please try again.")

	# compute the intervals of stability
	intervals = np.zeros((m, n, 2))
	# TODO: parallelize this
	for k in range(m):
		for l in range(n):
			if A[k, l] != 0:
				intervals[k, l, :] = NumSwitch.interval_of_stability(A, Ainv, k, l, max_bound=max_bound, num_sample=num_sample)
			elif pert_zero:
				intervals[k, l, :] = NumSwitch.interval_of_stability(A, Ainv, k, l, max_bound=max_bound, num_sample=num_sample)

	# save these
	print("Saving asymptotic stability to: %s" % asymp_stab_file)
	np.save(asymp_stab_file, intervals)

	# Compute the num switch functions
	# TODO: parallelize this
	crit_epsilon_array = np.zeros((n, n, n, n))
	for k in range(n):
		for l in range(n):
			for i in range(n):
				for j in range(n):
					crit_epsilon_array[k, l, i, j] = NumSwitch.critical_epsilon(Ainv, k, l, i, j)

	num_switch_funcs = dict()
	for k in range(n):
		for l in range(n):
			res = NumSwitch.num_switch_from_crit_eps(crit_epsilon_array, intervals, k, l)
			num_switch_funcs[k,l] = res

	# Save it
	print("Saving shape of num switch functions to: %s" % num_switch_file)
	fid = open(num_switch_file, 'w')
	for k in range(n):
		for l in range(n):
			key = (k,l)
			dict_val = num_switch_funcs[key]
			# TODO: switch based on zero entries
	#for key, dict_val in num_switch_funcs.items():
		#print(key)
		#print(dict_val)
			fid.write("%d\t%d\t" % (key[0], key[1]))
			for i in range(len(dict_val) - 1):
				val, (start, stop) = dict_val[i]
				fid.write("%f\t%f\t%f\t" % (val, start, stop))
			if dict_val:
				val, (start, stop) = dict_val[-1]
				fid.write("%d\t%f\t%f\n" % (val, start, stop))
			else:
				fid.write("\n")
	fid.close()