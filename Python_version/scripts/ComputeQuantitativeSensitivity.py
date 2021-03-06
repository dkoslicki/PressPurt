import argparse
import numpy as np
import os
import sys
import pandas as pd

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


def get_parser():
	parser = argparse.ArgumentParser(description="This script Generates the quantitative sensitivity: perturbing each single entry off to infinity.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('input_file', type=str, help="Input comma separated file for the jacobian matrix.")
	parser.add_argument('output_folder', type=str, help="Output folder. A file <prefix>_quantitative_sensitivity.csv' will be put here.")
	parser.add_argument('-p', '--prefix', help="Prefix of output files, if you so choose.", default=None)
	return parser


if __name__ == '__main__':
	parser = get_parser()

	# read in the arguments
	args = parser.parse_args()
	input_file = os.path.abspath(args.input_file)
	output_folder = args.output_folder
	if output_folder:
		output_folder = os.path.abspath(output_folder)
	prefix = args.prefix

	if not os.access(output_folder, os.W_OK):
		raise Exception("The provided directory %s is not writable." % output_folder)

	if prefix:
		quant_sens_file = os.path.join(output_folder, prefix + "_quantitative_sensitivity.csv")
		MRS_file = os.path.join(output_folder, prefix + "_MRS.csv")
	else:
		quant_sens_file = os.path.join(output_folder, "quantitative_sensitivity.csv")
		MRS_file = os.path.join(output_folder, "MRS.csv")

	# read in the input matrix
	#A = np.loadtxt(input_file, delimiter=",")
	A, row_names, column_names = NumSwitch.import_matrix(input_file)
	Ainv = np.linalg.inv(A)
	m, n = A.shape

	# make sure the original matrix is itself asymptotically stable
	if not SS.is_stable(A):
		raise Exception("Sorry, the input matrix is not stable itself (all eigenvalues must have negative real part). Please try again.")

	# do the computation
	quant_sens_values = np.zeros((m, n))
	for i in range(m):
		for j in range(n):
			if A[i, j] != 0:  # only perturb non-zero entries
				quant_sens_values[i, j] = MRS.quant_sens(A, i, j)
			else:
				quant_sens_values[i, j] = 0

	#np.savetxt(quant_sens_file, quant_sens_values, delimiter=',')
	df = pd.DataFrame(quant_sens_values)
	df.index = row_names
	df.columns = column_names
	df.to_csv(quant_sens_file)
	mrs = MRS.MRS(A)
	np.savetxt(MRS_file, [mrs], delimiter=',')
