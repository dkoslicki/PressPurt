import numpy as np
import scipy.linalg as sl
import argparse
import os

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description="Generate a random symmetric, asymptotically stable matrix.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('output_file', type=str, help="Name of file to create (csv).")
	parser.add_argument('-n', '--num_rows', type=int, help="Number of rows/columns of the matrix.", default=10)
	parser.add_argument('-r', '--range', type=float, help="(sqrt-ish) of the range of the random matrix.", default=2)

	# read in the arguments
	args = parser.parse_args()
	output_file = os.path.abspath(args.output_file)
	if not os.access(os.path.dirname(output_file), os.W_OK):
		raise Exception("The directory %s is not writeable." % os.path.dirname(output_file))

	n = int(args.num_rows)
	r = float(args.range)
	# create random asymptotically stable matrices
	eig = np.diag(-r*np.random.rand(n))  # negative eigenvalues
	v = sl.orth(r*(np.random.rand(n, n)-.5))
	A = np.transpose(v).dot(eig.dot(v))
	np.savetxt(output_file, A, delimiter=',')
