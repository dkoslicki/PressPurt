import argparse
import os
import sys
import numpy as np
import multiprocessing

# import stuff in the src folder
try:
    import NaiveSS
except ImportError:
    sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'src'))
    import NaiveSS
try:
    import NumSwitch
except ImportError:
    sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'src'))
    import NumSwitch


def get_parser():
    parser = argparse.ArgumentParser(description="This script pre-processes a matrix by figuring out what the intervals of asymptotic stability are, as well as finding which perturbation values lead to a sign switch.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('input_file', type=str, help="Input comma separated file for the jacobian matrix.")
    parser.add_argument('-n', '--num_iterates', type=int, help="Number of iterates in the Monte Carlo sampling to perform.", default=10000)
    parser.add_argument('-l', '--interval_length', type=float, help="Interval length over which to make the perturbations.", default=0.01)
    parser.add_argument('-t', '--threads', type=int, help="Number of threads to use.", default=multiprocessing.cpu_count())
    return parser

def run_MultiEntry(input_file, num_iterates, interval_length, threads):
    if __name__ == '__main__':
        #parser = get_parser()

        # read in the arguments
        #args = parser.parse_args()
        #num_threads = int(args.threads)
        num_threads = int(threads)
        #input_file = os.path.abspath(args.input_file)
        input_file = input_file
        #num_iterates = int(args.num_iterates)
        num_iterates = int(num_iterates)
        #interval_length = float(args.interval_length)
        interval_length = float(interval_length)

        # check for sanity of input parameters
        if not interval_length > 0:
            raise Exception("interval_length must be larger than 0; provided value: %d." % interval_length)

        # read in the input matrix
        #A = np.loadtxt(input_file, delimiter=",")
        A, row_names, column_names = NumSwitch.import_matrix(input_file)
        expectation = NaiveSS.naive_SS(A, num_iterates=num_iterates, interval_length=interval_length, num_threads=num_threads)
        return(expectation)
