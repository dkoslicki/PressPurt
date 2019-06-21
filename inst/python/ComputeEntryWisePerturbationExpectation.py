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
from multiprocessing import Pool  # Much faster without dummy (threading)
import multiprocessing
import pandas as pd
import decimal
# import stuff in the src folder
import MRS
import SS
import NumSwitch

# TODO: note, these must be present if you attempt to import the pickling of the distributions and a custom one was chosen
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

class Helper(object):
    def __int__(self):
        pass

    def helper(self, n, k, l, num_switch_funcs, asymp_stab, dists):
        try:
            val = NumSwitch.exp_num_switch_from_crit_eps(n, k, l, num_switch_funcs, asymp_stab, dist=dists[k, l])
        except KeyError:
            val = None
        return (val, k, l)
    def helper_star(self, arg):
        return self.helper(*arg)

def run_EntryWise(input_folder, prefix, distribution_type,
                  input_a, input_b, threads, **kwargs):
    if __name__ == '__main__':
        save = False
        prefix = prefix
        distribution_type = distribution_type
        input_a = float(input_a)
        input_b = float(input_b)
        num_threads = int(threads)
        if distribution_type not in known_distributions:
            raise Exception("You can only choose between the following distributions: %s. You provided '%s'." %(", ".join(known_distributions), distribution_type))
        if input_folder is not None:
            # If input file is specified, get path to files, check if they
            # exist, read in files
            save = True
            input_folder = input_folder
            output_folder = os.path.abspath(input_folder)
            if not os.access(output_folder, os.R_OK):
                raise Exception("The provided directory %s is not readable." % output_folder)
            if prefix:
                asymp_stab_file = os.path.join(output_folder, prefix + "_asymptotic_stability.npy")
                num_switch_file = os.path.join(output_folder, prefix + "_num_switch_funcs.pkl")
                exp_num_switch_file = os.path.join(output_folder, prefix + "_expected_num_switch.csv")
                distribution_file = os.path.join(output_folder, prefix + "_distributions.pkl")
                matrix_size_file = os.path.join(output_folder, prefix + "_size.npy")
                row_names_file = os.path.join(output_folder, prefix + "_row_names.txt")
                column_names_file = os.path.join(output_folder, prefix + "_column_names.txt")
            else:
                asymp_stab_file = os.path.join(output_folder, "asymptotic_stability.npy")
                num_switch_file = os.path.join(output_folder, "num_switch_funcs.pkl")
                exp_num_switch_file = os.path.join(output_folder, "expected_num_switch.csv")
                #exp_num_test = os.path.join(output_folder, "expected_num_switch.pkl")
                distribution_file = os.path.join(output_folder, "distributions.pkl")
                matrix_size_file = os.path.join(output_folder, "size.npy")
                row_names_file = os.path.join(output_folder, "row_names.txt")
                column_names_file = os.path.join(output_folder, "column_names.txt")

            # check if files exist
            required_files = [asymp_stab_file, num_switch_file,
                              row_names_file, column_names_file, matrix_size_file]
            for file in required_files:
                if not os.access(file, os.R_OK):
                    raise Exception("Missing files. Please first run PreProcessMatrix.py and ComputeEntryWisePerturbationExpectation.py. Files should include: asymptotic_stablity.npy, num_switch_funcs.pkl, expected_num_switch.csv, and distributions.pkl")

            # read in the files
            n = int(np.load(matrix_size_file))
            asymp_stab = np.load(asymp_stab_file)
            num_switch_funcs = pickle.load(open(num_switch_file, 'rb'))

            row_names = []
            column_names = []
            with open(row_names_file, 'r') as fid:
                for line in fid.readlines():
                    line_strip = line.strip()
                    row_names.append(line_strip)

            with open(column_names_file, 'r') as fid:
                for line in fid.readlines():
                    line_strip = line.strip()
                    column_names.append(line_strip)
        else:
            n = kwargs["matrix_size"]
            asymp_stab = kwargs["asymp_stab"]
            num_switch_funcs = kwargs["num_switch"]
            #num_switch_file = "test_r/test3/num_switch_funcs.pkl"
            #num_switch_file = "test_r/test3/test_num_switch.pkl"
            #num_switch_funcs = pickle.load(open(num_switch_file, 'rb'))
            row_names = kwargs["row_names"]
            column_names = kwargs["col_names"]

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
                        #loc = 0
                        s = input_a
                        #scale = input_b
                        loc = interval[0]
                        scale = interval[1] - interval[0]
                        dists[k, l] = trunc_lognorm(a, b, s, loc, scale)
                    elif distribution_type == 'beta':
                        loc = interval[0]
                        scale = interval[1] - interval[0]
                        a = input_a
                        b = input_b
                        dists[k, l] = custom_beta(a, b, loc, scale)

        # Save the distributions
        if save is True:
            fid = open(distribution_file, 'wb')
            pickle.dump(dists, fid)
            fid.close()

        # Do the expectation calculation
        exp_num_switch_array = np.zeros((n, n))
        if n < 10:  # in serial
            for k in range(n):
                for l in range(n):
                    try:
                        val = NumSwitch.exp_num_switch_from_crit_eps(n, k, l, num_switch_funcs, asymp_stab, dist=dists[k, l])
                        exp_num_switch_array[k, l] = val
                    except KeyError:
                        pass
        else:  # do it in parallel
            # try the same thing, but parallelized
            #def helper(k, l):
            #    """
            #    helper function to make it easier to parallelize
            #    :param k: index
            #    :param l: index
            #    :return: none or float
            #    """
            #    try:
            #        val = NumSwitch.exp_num_switch_from_crit_eps(n, k, l, num_switch_funcs, asymp_stab, dist=dists[k, l])
            #    except KeyError:
            #        val = None
            #        return (val, k, l)

            #def helper_star(arg):
            #    """
            #    unwrap the tuple
            #    :param arg: tuple of (k,l) indices
            #    :return: none or float
            #    """
            #    return helper(*arg)

            # get all the arguments to compute
            to_compute_args = []
            for k in range(n):
                for l in range(n):
                    to_compute_args.append((n, k, l, num_switch_funcs, asymp_stab, dists))

            # start the pool and do the computations
            helper_obj = Helper()
            pool = Pool(processes=num_threads)
            res = pool.map(helper_obj.helper_star, to_compute_args)
            #return(res)
            # collect the results
            for val, k, l in res:
                if val is not None:
                    exp_num_switch_array[k, l] = val

        # export the results
        #np.savetxt(exp_num_switch_file, exp_num_switch_array, delimiter=',')
        #exp_test = exp_num_switch_array
        df = pd.DataFrame(exp_num_switch_array)
        df.index = row_names
        df.columns = column_names
        #exp_num_test = pickle.load(open("test_r/test3/expected_num_switch.pkl", 'rb'))
        if save is True:
            df.to_csv(exp_num_switch_file)
        if save is False:
            return(dists, df)
            #return(kwargs["num_switch"])


#infold = "test_r/test3"
#tips = run_EntryWise(input_folder = infold, prefix = None, distribution_type = "truncnorm",
#                     input_a = 0, input_b = -2, threads = 1)
#print(tips)
#print(type(tips))

