# This script will take the intervals of stability, the num_switch_funcs, and a distribution and compute the
# expected num switch array. Then save it
import argparse
import numpy as np
import os
import sys
import pickle
import timeit
import scipy.stats as st
from scipy.stats import rv_continuous

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

class trunc_lognorm():
	"""
	Truncated log normal
	"""
	def __init__(self, a, b, s, loc, scale):
		self.a = a
		self.b = b
		self.s = s
		self.loc = loc
		self.scale = scale
		self.nrm = st.lognorm.cdf(self.b, self.s, loc=self.loc, scale=self.scale)-st.lognorm.cdf(self.a, self.s, loc=self.loc, scale=self.scale)
	def cdf(self, x):
		return st.lognorm.cdf(x, self.s, loc=self.loc, scale=self.scale)/self.nrm
	def pdf(self, x):
		return st.lognorm.pdf(x, self.s, loc=self.loc, scale=self.scale) / self.nrm

class custom_beta():
	"""
	custom beta function
	"""
	def __init__(self, a, b, loc, scale):
		self.a = a
		self.b = b
		self.loc = loc
		self.scale = scale
	def cdf(self, x):
		return st.beta.cdf(x, self.a, self.b, loc=self.loc, scale=self.scale)
	def pdf(self, x):
		return st.beta.pdf(x, self.a, self.b, loc=self.loc, scale=self.scale)

known_distributions = ['truncnorm', 'uniform', 'trunc_lognorm', 'beta']

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description="This script all indices from the Koslicki & Novak (2018) paper", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('input_folder', type=str, help="Input folder. The location of the files created by PreProcessMatrix.py. eg 'output_folder/<prefix>_asymptotic_stability.npy'. This is also where the expected num swith array will be saved.")
	parser.add_argument('-n', type=int, help="Dimension of matrix. Eg 4x4 matrix means you use -n 4", required=True)
	parser.add_argument('-p', '--prefix', help="Prefix of output files, if you so choose.", default=None)
	parser.add_argument('-d', '--distribution_type', type=str, help="Kind of distribution to use. Valid choices are: %s." % (', '.join(known_distributions)), default='truncnorm')
	parser.add_argument('-a', type=float, help="First parameter to the distribution you choose. For truncnorm, this is the mean.", default=0)
	parser.add_argument('-b', type=float, help="First parameter to the distribution you choose. For truncnorm, this is the variance. Using a negative value indicates you want the standard deviation to be the length of the interval divided by the absolute value of the input parameter.", default=-2)

	# TODO: add custom distribution options
	# read in the arguments
	args = parser.parse_args()
	input_folder = args.input_folder
	output_folder = os.path.abspath(input_folder)
	prefix = args.prefix
	n = int(args.n)
	distribution_type = args.distribution_type
	input_a = float(args.a)
	input_b = float(args.b)

	if not os.access(output_folder, os.R_OK):
		raise Exception("The provided directory %s is not readable." % output_folder)

	if prefix:
		asymp_stab_file = os.path.join(output_folder, prefix + "_asymptotic_stability.npy")
		num_switch_file = os.path.join(output_folder, prefix + "_num_switch_funcs.pkl")
		exp_num_switch_file = os.path.join(output_folder, prefix + "_expected_num_switch.csv")
		distribution_file = os.path.join(output_folder, prefix + "_distributions.pkl")
	else:
		asymp_stab_file = os.path.join(output_folder, "asymptotic_stability.npy")
		num_switch_file = os.path.join(output_folder, "num_switch_funcs.pkl")
		exp_num_switch_file = os.path.join(output_folder, "expected_num_switch.csv")
		distribution_file = os.path.join(output_folder, "distributions.pkl")

	if distribution_type not in known_distributions:
		raise Exception("You can only choose between the following distributions: %s. You provided '%s'." % (", ".join(known_distributions), distribution_type))
	# read in the files
	asymp_stab = np.load(asymp_stab_file)
	num_switch_funcs = pickle.load(open(num_switch_file, 'rb'))

	# create all the distributions
	dists = dict()
	for k in range(n):
		for l in range(n):
			interval = asymp_stab[k, l]
			if interval[1] - interval[0] > 0:
				if distribution_type == 'truncnorm':  # normal distribution truncated to interval of stability
					if input_b < 0:  # negative indicates to scale the variance based on the size of the interval of stability
						my_std = (interval[1] - interval[0]) / float(np.abs(input_b))
					else:
						my_std = input_b  # otherwise use a fixed variance
					my_mean = input_a
					a, b = (interval[0] - my_mean) / my_std, (interval[1] - my_mean) / my_std
					dist = st.truncnorm(a, b, 0, my_std)
					dists[k, l] = dist
				elif distribution_type == 'uniform':  # uniform on the interval of stability
					loc = interval[0]
					scale = interval[1] - interval[0]
					dist = st.uniform(loc=loc, scale=scale)
					dists[k, l] = dist
				elif distribution_type == 'trunc_lognorm':
					a = interval[0]
					b = interval[1]
					loc = 0
					s = input_a
					scale = input_b
					dists[k, l] = trunc_lognorm(a, b, s, loc, scale)
				elif distribution_type == 'beta':
					loc = interval[0]
					scale = interval[1] - interval[0]
					a = input_a
					b = input_b
					dists[k, l] = custom_beta(a, b, loc, scale)

	# Save the distributions
	fid = open(distribution_file, 'wb')
	pickle.dump(dists, fid)
	fid.close()

	# Do the expectation calculation
	exp_num_switch_array = np.zeros((n, n))
	for k in range(n):
		for l in range(n):
			#t0 = timeit.default_timer()
			try:
				val = NumSwitch.exp_num_switch_from_crit_eps(n, k, l, num_switch_funcs, asymp_stab, dist=dists[k, l])
				exp_num_switch_array[k, l] = val
			except KeyError:
				pass
			#t1 = timeit.default_timer()
			#print(t1 - t0)

	# export the results
	np.savetxt(exp_num_switch_file, exp_num_switch_array, delimiter=',')
