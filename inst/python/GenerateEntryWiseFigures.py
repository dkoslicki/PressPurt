import argparse
import numpy as np
import scipy.stats as st
import os
import sys
#import matplotlib.pyplot as plt
import pickle
#import matplotlib.gridspec as gridspec
import pandas as pd


import MRS
import SS
import NumSwitch

# Required for the pickle import with a custom class
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


def make_plot(input_folder, prefix, all_numswitch_plots, list_of_numswitch_to_plot):
    if __name__ == '__main__':
        #parser = get_parser()
        # read in the arguments
        #args = parser.parse_args()
        #input_folder = args.input_folder
        input_folder = input_folder
        output_folder = os.path.abspath(input_folder)
        #prefix = args.prefix
        prefix = prefix
        #all_numswitch_plots = args.all_numswitch_plots
        all_numswitch_plots = all_numswitch_plots
        #if args.list_of_numswitch_to_plot:
        #    list_of_numswitch_to_plot = [int(i) - 1 for i in list(args.list_of_numswitch_to_plot)]
        #else:
        #    list_of_numswitch_to_plot = []
        if list_of_numswitch_to_plot:
            list_of_numswitch_to_plot = [int(i) - 1 for i in list(list_of_numswitch_to_plot)]
        else:
            list_of_numswitch_to_plot = []



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
            num_nonzero_file = os.path.join(output_folder, prefix + "_num_non_zero.npy")
        else:
            asymp_stab_file = os.path.join(output_folder, "asymptotic_stability.npy")
            num_switch_file = os.path.join(output_folder, "num_switch_funcs.pkl")
            exp_num_switch_file = os.path.join(output_folder, "expected_num_switch.csv")
            distribution_file = os.path.join(output_folder, "distributions.pkl")
            matrix_size_file = os.path.join(output_folder, "size.npy")
            row_names_file = os.path.join(output_folder, "row_names.txt")
            column_names_file = os.path.join(output_folder, "column_names.txt")
            num_nonzero_file = os.path.join(output_folder, "num_non_zero.npy")

        n = int(np.load(matrix_size_file))
        m = n
        required_files = [asymp_stab_file, num_switch_file, exp_num_switch_file, distribution_file, row_names_file, column_names_file, matrix_size_file, num_nonzero_file]
        for file in required_files:
            if not os.access(file, os.R_OK):
                raise Exception("Missing files. Please first run PreProcessMatrix.py and ComputeEntryWisePerturbationExpectation.py. Files should include: asymptotic_stablity.npy, num_switch_funcs.pkl, expected_num_switch.csv, and distributions.pkl")

        if all_numswitch_plots and list_of_numswitch_to_plot:
            raise Exception("The flags -a and -l are mutually exclusive.")
        if len(list_of_numswitch_to_plot) % 2 != 0:
            raise Exception("Only an even length list can be passed to -l.")

        # load the row/column names
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

        # if appropriate, get the indices to plot
        indices_to_plot = []
        # quick way to take entries two at a time
        it = iter(list_of_numswitch_to_plot)
        for i, j in zip(*[it]*2):
            indices_to_plot.append((i, j))

        # load the intervals of stability
        intervals = np.load(asymp_stab_file)

        # load the num_switch functions
        fid = open(num_switch_file, 'rb')
        num_switch_funcs = pickle.load(fid)
        fid.close()

        # load the distributions
        fid = open(distribution_file, 'rb')
        dists = pickle.load(fid)
        fid.close()

        ######################
        # Generate figure 2
        # let's do a mxn grid of these figures
        # TODO: should really make this into its own function
        #padding = .1
        if all_numswitch_plots:
            big_fig, axarr = plt.subplots(m, n)
            big_fig.suptitle("Number of mis-predictions versus perturbation value, \n overlaid with distribution over stable perturbation values")
            for k in range(m):
                for l in range(n):
                    if (k, l) in dists:
                        interval = intervals[k, l, :]
                        padding = (interval[1] - interval[0])/float(100)
                        x_range = np.linspace(interval[0] - padding, interval[1] + padding, 250)
                        dist = dists[k, l]
                        dist_vals = [dist.pdf(eps) for eps in x_range]
                        (nsx, nsy) = NumSwitch.num_switch_to_step(num_switch_funcs, intervals, k, l)
                        # plot both simultaneously on the same graph (two y-axis plot)
                        ax1 = axarr[k, l]
                        #ax1.title.set_text('(%d,%d) entry' % (k+1, l+1))
                        ax1.title.set_text('(%s, %s) entry' % (column_names[k], row_names[l]))
                        ax1.step(nsx, nsy, 'b-')
                        ax1.tick_params('y', colors='b')

                        ax2 = ax1.twinx()
                        ax2.plot(x_range, dist_vals, 'tab:gray')
                        ax2.tick_params('y', colors='tab:gray')
                        ax2.fill_between(x_range, dist_vals, color='tab:gray', alpha=0.5)
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
            grid_size = int(np.ceil(np.sqrt(len(indices_to_plot))))
            big_fig, axarr = plt.subplots(grid_size, grid_size)
            gs1 = gridspec.GridSpec(grid_size, grid_size)
            gs1.update(wspace=0.025, hspace=0.05)
            try:
                axarr_flat = axarr.flatten()
            except AttributeError:
                axarr_flat = [axarr]
            big_fig.suptitle("Number of mis-predictions versus perturbation value, \n overlaid with distribution over stable perturbation values")
            num_plotted = 0
            for k, l in indices_to_plot:
                ax1 = axarr_flat[num_plotted]
                if (k, l) in dists:
                    interval = intervals[k, l, :]
                    padding = (interval[1] - interval[0]) / float(100)
                    x_range = np.linspace(interval[0] - padding, interval[1] + padding, 250)
                    dist = dists[k, l]
                    dist_vals = [dist.pdf(eps) for eps in x_range]
                    (nsx, nsy) = NumSwitch.num_switch_to_step(num_switch_funcs, intervals, k, l)
                    # plot both simultaneously on the same graph (two y-axis plot)
                    #ax1 = axarr[k, l]
                    #ax1.title.set_text('(%d,%d) entry' % (k + 1, l + 1))
                    ax1.title.set_text('(%s, %s) entry' % (column_names[k], row_names[l]))
                    ax1.step(nsx, nsy, 'b-')
                    ax1.tick_params('y', colors='b')

                    ax2 = ax1.twinx()
                    ax2.plot(x_range, dist_vals, 'tab:gray')
                    ax2.tick_params('y', colors='tab:gray')
                    ax2.fill_between(x_range, dist_vals, color='tab:gray', alpha=0.5)
                    num_plotted += 1
                else:
                    axarr_flat[num_plotted].axis('off')  # don't show the ones we are not perturbing
            try:
                axes_flat = axarr.flatten()
            except AttributeError:
                axes_flat = [axarr]
            for i in range(num_plotted, grid_size**2):
                axes_flat[i].axis('off')
            #plt.tight_layout(pad=0.1, w_pad=.1, h_pad=.9)
            big_fig.text(0.5, 0.01, 'Epsilon value', ha='center', va='center')
            big_fig.text(0.03, 0.5, 'Number of incorrect predictions', ha='center', va='center', rotation='vertical', color='b')
            big_fig.text(.99, 0.5, 'Probability density', ha='center', va='center', rotation='vertical', color='tab:gray')
            plt.subplots_adjust(top=.9)
            plt.draw()
            plt.pause(0.01)

        #####################
        # Get the expected number of sign switches, in a table
        #exp_num_switch_array = np.loadtxt(exp_num_switch_file, delimiter=',')
        df = pd.read_csv(exp_num_switch_file, header=0, index_col=0)
        exp_num_switch_array = df.values

        #####################
        # Generate figure 3
        fig, ax = plt.subplots()
        # im = ax.imshow(exp_num_switch_array, cmap=plt.get_cmap('seismic'))
        im = ax.imshow(exp_num_switch_array, cmap=plt.get_cmap('YlOrBr'))
        # We want to show all ticks...
        ax.set_xticks(np.arange(m))
        ax.set_yticks(np.arange(n))
        # ... and label them with the respective list entries
        #ax.set_xticklabels(['%d' % (i + 1) for i in range(n)])
        #ax.set_yticklabels(['%d' % (j + 1) for j in range(m)])
        ax.set_xticklabels(column_names)
        ax.set_yticklabels(row_names)
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

        num_non_zero = np.load(num_nonzero_file)
        if num_non_zero != 0:
            ave_expected_num_sign_switches = exp_num_switch_array.sum()/float(num_non_zero)
        else:
            ave_expected_num_sign_switches = 0
        print("Average expected percentage of mis-predictions (perturbing each edge individually): %f%% (i.e. %f entries)" % (100*ave_expected_num_sign_switches, ave_expected_num_sign_switches*m*n))
        # TODO: check with Mark if I should be multiplying by n^2 or num_non_zero
        input("press any key to quit")



#infolder = "test_2/"

#make_plot(input_folder = infolder, prefix = None, all_numswitch_plots = True, list_of_numswitch_to_plot = False)


