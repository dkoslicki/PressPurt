import argparse
import numpy as np
import os
import sys
import matplotlib.pyplot as plt
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
	parser = argparse.ArgumentParser(description="This script Generates the quantitative sensitivity figure.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('input_folder', type=str, help="Input folder containing the output of ComputeQuantitativeSensitivity.py")
	parser.add_argument('-p', '--prefix', help="Prefix of output files, if you so choose.", default=None)
	return parser


if __name__ == '__main__':
	parser = get_parser()

	# read in the arguments
	args = parser.parse_args()
	input_folder = os.path.abspath(args.input_folder)
	prefix = args.prefix

	if not os.access(input_folder, os.R_OK):
		raise Exception("The provided directory %s is not writable." % input_folder)

	if prefix:
		quant_sens_file = os.path.join(input_folder, prefix + "_quantitative_sensitivity.csv")
		MRS_file = os.path.join(input_folder, prefix + "_MRS.csv")
	else:
		quant_sens_file = os.path.join(input_folder, "quantitative_sensitivity.csv")
		MRS_file = os.path.join(input_folder, "MRS.csv")

	#quant_sens_values = np.loadtxt(quant_sens_file, delimiter=',')
	df = pd.read_csv(quant_sens_file, header=0, index_col=0)
	quant_sens_values = df.values
	m, n = quant_sens_values.shape

	# print out a statistic
	mrs = np.loadtxt(MRS_file, delimiter=',')
	print('Average quantitative sensitivity of A when perturbing each each individually by an arbitrarily large amount : %f' % mrs)

	# generate the figure
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
	input("Press any key to quit")
