# This script will take the intervals of stability, the num_switch_funcs, and a distribution and compute the
# expected num switch array. Then save it
import argparse
import numpy as np
import os
import sys
import pickle

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
	parser.add_argument('input_folder', type=str, help="Input folder. The location of the files created by PreProcessMatrix.py. eg 'output_folder/<prefix>_asymptotic_stability.npy'. This is also where the expected num swith array will be saved.")
	parser.add_argument('-n', type=int, help="Dimension of matrix. Eg 4x4 matrix means you use -n 4", required=True)
	parser.add_argument('-p', '--prefix', help="Prefix of output files, if you so choose.", default=None)

	# TODO: add custom distribution options
	# read in the arguments
	args = parser.parse_args()
	input_folder = args.input_folder
	output_folder = os.path.abspath(input_folder)
	prefix = args.prefix
	n = int(args.n)

	if not os.access(output_folder, os.R_OK):
		raise Exception("The provided directory %s is not readable." % output_folder)

	if prefix:
		asymp_stab_file = os.path.join(output_folder, prefix + "_asymptotic_stability.npy")
		num_switch_file = os.path.join(output_folder, prefix + "_num_switch_funcs.pkl")
		exp_num_switch_file = os.path.join(output_folder, prefix + "_expected_num_switch.csv")
	else:
		asymp_stab_file = os.path.join(output_folder, "asymptotic_stability.npy")
		num_switch_file = os.path.join(output_folder, "num_switch_funcs.pkl")
		exp_num_switch_file = os.path.join(output_folder, "expected_num_switch.csv")

	# read in the files
	asymp_stab = np.load(asymp_stab_file)
	num_switch_funcs = pickle.load(open(num_switch_file, 'rb'))

	# Do the expectation calculation
	exp_num_switch_array = np.zeros((n, n))
	for k in range(n):
		for l in range(n):
			try:
				val = NumSwitch.exp_num_switch_from_crit_eps(n, k, l, num_switch_funcs, asymp_stab, dist=None)
				exp_num_switch_array[k, l] = val
			except KeyError:
				pass

	# export the results
	np.savetxt(exp_num_switch_file, exp_num_switch_array, delimiter=',')
