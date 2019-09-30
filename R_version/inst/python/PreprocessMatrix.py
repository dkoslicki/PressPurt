import argparse
import numpy as np
import os
import sys
import pickle
import timeit
from multiprocessing import Pool  # Much faster without dummy (threading)
import multiprocessing
import itertools


import NumSwitch

class Helper(object):
    def __init__(self):
        pass

    def helper(self,k ,l, A, Ainv, max_bound):
        val = NumSwitch.interval_of_stability(A, Ainv, k, l, max_bound=max_bound)
        return (val, k, l)
    def helper_star(self, arg):
        return self.helper(*arg)


def run_preproc(input_file, output_folder, prefix, max_bound, zero_peturb, threads, verbose):
    if __name__ == '__main__':
        save = False
        verbose = verbose
        num_threads = int(threads)
        input_file = os.path.abspath(input_file)
        output_folder = output_folder
        if output_folder is not None:
            output_folder = os.path.abspath(output_folder)
            save = True
            if not os.access(output_folder, os.W_OK):
                raise Exception("The provided directory %s is not writable." % output_folder)
        prefix = prefix
        max_bound = float(max_bound)
        pert_zero = zero_peturb

        if save is True:
            if prefix:
                asymp_stab_file = os.path.join(output_folder, prefix + "_asymptotic_stability.npy")
                num_switch_file = os.path.join(output_folder, prefix + "_num_switch_funcs.pkl")
                matrix_size_file = os.path.join(output_folder, prefix + "_size.npy")
                row_names_file = os.path.join(output_folder, prefix + "_row_names.txt")
                column_names_file = os.path.join(output_folder, prefix + "_column_names.txt")
                num_nonzero_file = os.path.join(output_folder, prefix + "_num_non_zero.npy")
            else:
                asymp_stab_file = os.path.join(output_folder, "asymptotic_stability.npy")
                num_switch_file = os.path.join(output_folder, "num_switch_funcs.pkl")
                matrix_size_file = os.path.join(output_folder, "size.npy")
                row_names_file = os.path.join(output_folder, "row_names.txt")
                column_names_file = os.path.join(output_folder, "column_names.txt")
                num_nonzero_file = os.path.join(output_folder, "num_non_zero.npy")

        # check for sanity of input parameters
        if not max_bound > 0:
            raise Exception("max_bound must be larger than 0; provided value: %d." % max_bound)

        # read in the input matrix
        #A = np.loadtxt(input_file, delimiter=",")
        A, row_names, column_names = NumSwitch.import_matrix(input_file)

        # save the number of non-zero entries
        num_non_zero = len(np.where(A)[0])

        Ainv = np.linalg.inv(A)
        m, n = A.shape

        if save is True:
            with open(row_names_file, 'w') as fid:
                for item in row_names:
                    fid.write("%s\n" % item)
            with open(column_names_file, 'w') as fid:
                for item in column_names:
                    fid.write("%s\n" % item)
            np.save(num_nonzero_file, num_non_zero)
            np.save(matrix_size_file, n)

        # make sure the original matrix is itself asymptotically stable
        if not NumSwitch.is_stable(A):
            raise Exception("Sorry, the input matrix is not stable itself (all eigenvalues must have negative real part). Please try again.")

        # compute the intervals of stability
        intervals = np.zeros((m, n, 2))
        #def helper(k,l):
        #    val = NumSwitch.interval_of_stability(A, Ainv, k, l, max_bound=max_bound)
        #    return (val, k, l)
        #def helper_star(arg):
        #    return helper(*arg)

        to_compute_args = []
        for k in range(m):
            for l in range(n):
                if A[k, l] != 0:
                    #intervals[k, l, :] = NumSwitch.interval_of_stability(A, Ainv, k, l, max_bound=max_bound)
                    to_compute_args.append((k, l, A, Ainv, max_bound))
                elif pert_zero:
                    #intervals[k, l, :] = NumSwitch.interval_of_stability(A, Ainv, k, l, max_bound=max_bound)
                    to_compute_args.append((k, l, A, Ainv, max_bound))
        # start the pool and do the computations
        helper_obj = Helper()
        pool = Pool(processes=num_threads)
        res = pool.map(helper_obj.helper_star, to_compute_args)
        # collect the results
        for val, k, l in res:
            intervals[k, l, :] = val

        # save these
        if verbose:
            print("Saving asymptotic stability to: %s" % asymp_stab_file)
        if save is True:
            np.save(asymp_stab_file, intervals)

        # Compute the num switch functions
        # Looks like up to matrices of size 50, the parallelization doesn't help. NumSwitch.critical_epsilon is CPU non-intensive
        crit_epsilon_array = np.zeros((n, n, n, n))
        for k in range(n):
            for l in range(n):
                for i in range(n):
                    for j in range(n):
                        crit_epsilon_array[k, l, i, j] = NumSwitch.critical_epsilon(Ainv, k, l, i, j)

        num_switch_funcs = dict()
        # Looks like up to matrices of size 50, the parallelization doesn't help. Here again, the computation is CPU non-intensive
        for k in range(n):
            for l in range(n):
                if A[k, l] != 0:
                    res = NumSwitch.num_switch_from_crit_eps(crit_epsilon_array, intervals, k, l)
                    num_switch_funcs[k, l] = res
                elif pert_zero:
                    res = NumSwitch.num_switch_from_crit_eps(crit_epsilon_array, intervals, k, l)
                    num_switch_funcs[k, l] = res

        # Save it
        if verbose:
            print("Saving shape of num switch functions to: %s" % num_switch_file)
        if save is True:
            fid = open(num_switch_file, 'wb')
            pickle.dump(num_switch_funcs, fid)
            fid.close()
        if save is False:
            return(A, n, column_names, row_names, num_non_zero, num_switch_funcs, intervals)

        # This is a text format
       # for k in range(n):
       #     for l in range(n):
       #         key = (k, l)
       #         dict_val = num_switch_funcs[key]
       #         # TODO: switch based on zero entries
       #         fid.write("%d\t%d\t" % (key[0], key[1]))
       #         for i in range(len(dict_val) - 1):
       #             val, (start, stop) = dict_val[i]
       #             fid.write("%f\t%f\t%f\t" % (val, start, stop))
       #             if dict_val:
       #                 val, (start, stop) = dict_val[-1]
       #                 fid.write("%d\t%f\t%f\n" % (val, start, stop))
       #             else:
       #                 fid.write("\n")
       # fid.close()


#infile = "ExampleJacobians/Modules/IGP.csv"
#out_folder = "test_2/"

#tp = run_preproc(infile, out_folder, prefix=None, max_bound=10, zero_peturb=False, threads=2, verbose=False)
#print(tp)


