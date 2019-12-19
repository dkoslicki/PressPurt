# PressPurt

This repository contains an *R* implementation of all the results contained in the Koslicki &amp; Novak JMB paper [1].  See [Python_version](https://github.com/dkoslicki/PressPurt/tree/master/Python_version) for an alternative implementation.

- [Installation](#installation)
  * [Required python installation](#required-python-installation)
    + [Conda](#conda)
        * [Linux](#linux)
        * [MacOS](#macos)
        * [Windows](#windows)
      - [Install python dependencies](#install-python-dependencies)
        * [From R](#from-r)
        * [Via the command line](#via-the-command-line)
    + [Virtualenv](#virtualenv)
      - [Installing python dependencies](#installing-python-dependencies)
        * [From R](#from-r-1)
- [Usage](#usage)
  * [Quick start](#quick-start)
  * [Detailed usage](#detailed-usage)
    + [Single edge uncertainty](#single-edge-uncertainty)
      - [Qualitative analysis](#qualitative-analysis)
      - [Quantitative analysis](#quantitative-analysis)
    + [Multiple edge uncertainty](#multiple-edge-uncertainty)
      - [Qualitative analysis](#qualitative-analysis-1)
- [Complete documentation](#complete-documentation)
- [Citations](#citations)


# Installation
`PressPurt` depends on R >= 3.1.0. To check your R version run `version` from the R console.

This package also depends on the following R libraries:

* data.table
* ggplot2
* gridExtra
* reticulate >= 1.11

If not already installed, these packages will be installed upon the installation of `PressPurt`.


## Required python installation
This package is dependent on python and uses the R package `reticulate` under the hood to use python. Thus python and the dependant packages must be installed.

There are two suggested approaches:

* Conda (Anaconda/miniconda) which supports Linux/Mac/Windows
* Virtual environments which supports Linux/Mac OS

It is also possible to set a specific version of python with `reticulate::use_python("/path/to/python")` before you run 
any of the PressPurt commands but it is not recommended. You’ll also have to ensure that the dependent python packages 
are installed. For more information on using `reticulate directly`: [see the reticulate docs](https://rstudio.github.io/reticulate/articles/python_packages.html).

### Conda
You will need to install Anaconda or Miniconda. Conda is a package and environment manager that is open source. The main 
difference between Anaconda and Miniconda is that Anaconda comes with a bundle of pre-installed packages so it takes 
longer to download and install. Both Miniconda and Anaconda come with python but you can specify a specific version of 
python later as well.

This document will show you how to install Miniconda via the command-line. For more information on installation and how 
to install it graphically (no command-line) [see the conda docs](https://docs.conda.io/projects/continuumio-conda/en/latest/user-guide/install/index.html).

##### Linux
1. [Download the Miniconda installer for Linux](https://docs.conda.io/en/latest/miniconda.html)
1. In a terminal window run: bash Miniconda3-latest-Linux-x86_64.sh
1. Follow the prompts - accepting defaults should be fine.
1. Close your Terminal window and re-open a new one.
1. Test your installation by running: conda list
    * This should result in a list of installed packages.

##### MacOS
1. [Download the Miniconda installer for macOS](https://docs.conda.io/en/latest/miniconda.html)
1. In a terminal window run: bash Miniconda3-latest-MacOSX-x86_64.sh
1. Follow the prompts - accepting defaults should be fine.
1. Close your Terminal window and re-open a new one.
1. Test your installation by running: conda list
    * This should result in a list of installed packages.
1. Download, install, and activate (by launching at least once after installation) [Xcode](https://developer.apple.com/xcode/) via the Apple App Store.

##### Windows
1. [Download the Miniconda installer for Windows](https://docs.conda.io/en/latest/miniconda.html)
1. Run the .exe file by double clicking.
1. Follow the prompts - accepting defaults should be fine.
1. Test your installation by opeing the Anacond Prompt from the Start menu and running running: conda list
    * This should result in a list of installed packages.

#### Install python dependencies

You will now need to download a number of required python packages.

##### From R

First you’ll need to load `PressPurt`, check which versions of python are found and if any Conda environments exists.

```R
# Install PressPurt
install.packages("devtools")
# you may need to manually install libcurl4-openssl-dev via apt-get if using linux
devtools::install_github("dkoslicki/PressPurt", subdir = "R_version", build_vignettes = FALSE)
# load the library
library(PressPurt)
# Find the versions of python, conda environments and virtualenvs available on your machine:
find_python()
#> Default Python:
#>  /home/gibbond/anaconda3/envs/r-reticulate/bin/python 
#> 
#>  Python versions found:
#>  /home/gibbond/anaconda3/envs/r-reticulate/bin/python /home/gibbond/anaconda3/bin/python /usr/bin/python /usr/bin/python3 /home/gibbond/.virtualenvs/PressPurt/bin/python /home/gibbond/.virtualenvs/Test/bin/python /home/gibbond/anaconda3/envs/lpi-asms/bin/python /home/gibbond/anaconda3/envs/picrust2/bin/python /home/gibbond/anaconda3/envs/test_r_versions/bin/python /home/gibbond/anaconda3/envs/test-env/bin/python 
#> 
#>  List of condaenvironments:
#>  anaconda3: /home/gibbond/anaconda3/bin/python, lpi-asms: /home/gibbond/anaconda3/envs/lpi-asms/bin/python, picrust2: /home/gibbond/anaconda3/envs/picrust2/bin/python, r-reticulate: /home/gibbond/anaconda3/envs/r-reticulate/bin/python, test-env: /home/gibbond/anaconda3/envs/test-env/bin/python, test_press: /home/gibbond/anaconda3/envs/test_press/bin/python, test_r_versions: /home/gibbond/anaconda3/envs/test_r_versions/bin/python, 
#> 
#>  List of virtualenvs:
#>  PressPurt Test
```
This shows you your default Python path, other available python paths, a list of your conda environments and virtualenvs.

The `set_python` function lets one specify a conda environment to use. If you want to make a new conda environment, 
you’ll need to specify which python version to use in addition to the conda environment name.

You may have python installed under `/usr/bin/python3`, `Anaconda2` and/or `Anaconda3` installed so it’s important to 
specify which python to use. In this case I’ll use the anaconda3 python.

```R
set_python(condaenv="test-r", version="/home/gibbond/anaconda3/bin/python",verbose = TRUE)
```
Next you’ll need to install the python dependencies in your conda environment. 
This will install the following into your `test-r` conda env:

* numpy
* scipy
* matplotlib
* sympy
* pathos
* pandas

```R
py_depend(condaenv = "test-r")
```

If the dependencies don’t install correctly, see below for instructions on how to install via the command-line.

If you want to use the `test-r` environment in a new R session you can access it via:
```R
set_python(condaenv="test-r")
```

##### Via the command line
You can also install your Python dependencies via the command line with conda or with the Anaconda Prompt (Windows).

```bash
# Check that conda is running
conda --version
# List conda envs
conda info --envs
# Create conda environment "test-r"
conda create --name test-r
# activate test-r
source activate test-r
# check installed packages
conda list
# install python dependencies
conda install matplotlib numpy pandas scipy sympy
pip install pathos
# Check to make sure they were installed
conda list
```
Now that you have created the conda environments and installed dependencies, you can use the new env in R! Run:
```R
# load the library
library(PressPurt)
# You should see "test-r" listed as a conda env
find_python()
# Set your conda env
set_python(condaenv="test-r")
```

### Virtualenv
Instead of using conda, you can use a virtual python environment instead. If needed, install `virtualenv` with `pip`:
```bash
# install virtualenv with pip
pip install virtualenv
# test your installation
virtualenv --version
# create a virtualenv
virtualenv venv
```

#### Installing python dependencies
As above, you will need to install the python dependencies in this new virtual environment.

##### From R
First, you will need to install the R package. In R, run
```R
# Install PressPurt
install.packages("devtools")
devtools::install_github("dkoslicki/PressPurt", subdir = "R_version", build_vignettes = FALSE)
# load the library
library(PressPurt)
```
The `set_python_virtual` function lets one specify a virtual environment to use. If you want to make a new virtual 
environment, you’ll need to specify which python version to use in addition to the virtual environment name.

You may check which version(s) of python you have in a fresh R session with `find_python()`. You should use the sample 
one you used to install `virtualenv`.
```R
set_python_virtual(version = "/usr/bin/python3", virtualenv = "r-reticulate", verbose = TRUE)
```

Next you’ll need to install the python dependencies in your virtual environment. This will install the following into 
your `PressPurt` virtual env:

* numpy
* scipy
* matplotlib
* sympy
* pathos
* pandas

```R
py_depend(virtualenv = "PressPurt")
```

If you want to use the `PressPurt` environment in a new R session you can access it via:
```R
set_python_virtual(virtualenv = "PressPurt")
```

Now that you have created the virtual environment and installed dependencies, you can use the new env in R! Run:

```R
# load the library
library(PressPurt)
# You should see "test-r" listed as a conda env
find_python()
# Set your virtualenv
set_python_virtual(virtualenv = "PressPurt")
```

Information on how to remove/delete an environment may be found in the following [conda documentation](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#removing-an-environment).

# Usage
**Please note! After following the installation instructions above, you will need to restart your R session to complete the installation**

**Please also note! A detailed guide about how to use the R version (including visualizations and how to interact with the 
data structures) is included [in this vignette](https://htmlpreview.github.io/?https://github.com/dkoslicki/PressPurt/blob/master/R_version/vignettes/basic_tutorial.html).**

## Quick start

* Load the library
* Check available python versions
* Set python version and conda environment. You can set up a new conda env or use a previous one.

```
library(PressPurt)
find_python()
set_python("r-reticulate", verbose = TRUE)  # or whatever version of python you want to install the dependencies on
```

* install python dependencies

```
py_depend("r-reticulate")
```

* Run PreprocessMatrix
* Run ComputeEntryWisePerturbationExpectation


```
# Path to your matrix
infile <- "../Python_version/ExampleJacobians/Modules/IGP.csv"
PreProsMatrix <- PreprocessMatrix(input_file = infile, output_folder = NULL, max_bound = 10, threads = 2)
Entrywise <- ComputeEntryWisePerturbationExpectation(PreProsMatrix = PreProsMatrix, 
                                        distribution_type="truncnorm", 
                                        input_a=0, input_b=-2, threads=1)
```

* Plot
    * specific figures
    * all figures

```
list_of_numswitch_to_plot <- list(c(1, 1), c(1, 2))
GenerateEntryWiseFigures(EntryWise=Entrywise, 
                         all_numswitch_plots = FALSE, 
                         list_of_numswitch_to_plot=list_of_numswitch_to_plot)
GenerateEntryWiseFigures(EntryWise=Entrywise, 
                         all_numswitch_plots = TRUE)
```

## Detailed usage
This section contains details about what the various scripts are actually computing and mirrors the information contained 
in the [Python_version Readme](https://github.com/dkoslicki/PressPurt/blob/master/Python_version/README.md).


In general, there are three kinds of analyses:
1. Single edge uncertainty
    - Qualitative analysis
    - Quantitative analysis
2. Multiple edge uncertainty
    - Qualitative analysis

Each analysis requires a numeric Jacobian matrix evaluated at an equilibrium point. This must be 
provided as a CSV file containing column and row names that are either non-numeric ("A","B","C",...) 
or numeric (but not including zero) entries.  If you created your matrix within R, you will have to
export it to a CSV file before proceeding.

Before any analysis can be performed, you must first pre-process a matrix so that the intervals 
of asymptotic stability can be determined, as well as a calculation that will determine which 
uncertainty values on which edges will result in a sign change in the net effects matrix. This 
can be accomplished via:
```R
# Use the build in example data, otherwise make this point to your Jacobian CSV file
infile <- system.file("extdata", "Modules", "IGP.csv", package = "PressPurt")  
# preprocess the IGP.csv example and store the results in memory
PreProsMatrix <- PreprocessMatrix(input_file = infile, 
                                  output_folder = NULL, 
                                  max_bound = 10, 
                                  threads = 2)   
```
You will notice that your workspace now contains a number of files that are required for downstream analysis. Details about 
the workspace data structures can be found [in this vignette](https://htmlpreview.github.io/?https://github.com/dkoslicki/PressPurt/blob/master/R_version/vignettes/basic_tutorial.html).

Note that you need only preprocess the Jacobian once, as changing parameters in the following analyses does not require 
you to re-run `PreprocessMatrix`.

### Single edge uncertainty
This set of analyses are for the situation where you wish to assess how the presence of uncertainty 
in the strength of the direct effect (edge weight) between a single focal pair of species (nodes) alters 
the signs of the predicted press perturbation net effects between all pairs of species in the network.

#### Qualitative analysis
The function `ComputeEntryWisePerturbationExpectation` will compute the expected number of sign 
switches/mispredictions in the net effects matrix when perturbing each entry of the Jacobian matrix individually.
This can be run with:
```R
# run the script on the preprocessed matrix contained in the Results folder assuming the edge uncertainties are distributed according to a truncated normal distribution using a mean of 0 and a variance equal to (the length of the interval divided by the absolute value of the input parameter)^2
Entrywise <- ComputeEntryWisePerturbationExpectation(PreProsMatrix = PreProsMatrix, 
                                        distribution_type="truncnorm", 
                                        input_a=0, input_b=-2, threads=1)   
```
The results will also be placed your workspace. Details about 
the workspace data structures can be found [in this vignette](https://htmlpreview.github.io/?https://github.com/dkoslicki/PressPurt/blob/master/R_version/vignettes/basic_tutorial.html).
In short, the matrix `Entrywise$expected_num_switch` is a matrix whose (i,j) entry shows the expected fraction of sign switches/mispredictions in the net effects 
matrix when the edge between nodes i and j experiences uncertainty according to the distribution specified. You can then 
make your own heat map of this matrix if you wish. 

Averaging the non-zero entries in `Entrywise$expected_num_switch` calculates the expected percentage of mispredictions
across the entire net effects matrix when each edge individually experiences uncertainty.

Use `?ComputeEntryWisePerturbationExpectation` to view the different kinds of uncertainty distributions available.

You can then visualize the results using the function `GenerateEntryWiseFigures`:
```R
# List of a few entries to plot
list_of_numswitch_to_plot <- list(c(1, 1), c(1, 2))
# plot just those entries
GenerateEntryWiseFigures(EntryWise=Entrywise, 
                         all_numswitch_plots = FALSE, 
                         list_of_numswitch_to_plot=list_of_numswitch_to_plot)
# plot everything
GenerateEntryWiseFigures(EntryWise=Entrywise, 
                         all_numswitch_plots = TRUE)
```
This command generates a plot of the (edge uncertainty value) versus 
(number of sign switches/mispredictions) in the net effects matrix, overlaid with the edge uncertainty distribution (truncated by the 
region of asymptotic stability). The option `list_of_numswitch_to_plot` specifies that only the (1,1) and (1,2) edges should be shown. 
If you use `-all_numswitch_plots = TRUE` instead, this will show similar plots for *all* edges (warning: this may be uninterpretable if you matrix is large). 


#### Quantitative analysis
Instead of asking whether a qualitative misprediction occurs, one might be interested in quantifying by how much you mispredict 
the change in abundance of other nodes when edge weights (individually) are under uncertainty.

As of yet, this analysis can only be performed using the [Python version](https://github.com/dkoslicki/PressPurt/tree/master/Python_version). 


### Multiple edge uncertainty

While the above analysis computes the sensitivity of the net effects matrix when edges individually experience
uncertainty, the following allows multiple edges to experience uncertainty simultaneously.

#### Qualitative analysis
The function `ComputeMultiEntryPerturbationExpectation` performs a Monte Carlo sampling of edge uncertainties according 
to a uniform distribution (within a specified interval length) and reports the average fraction of sampled uncertainties 
that lead to a sign switch/misprediction in the net effects matrix (considering only uncertainties that result in a stable Jacobian).

The command can be called with
```R
MultiPert <- ComputeMultiEntryPerturbationExpectation(input_file, num_iterates = 1000,
  interval_length = 0.01, threads = 4)
```
The above command samples uncertainties 1000 times on all edges (non zero entries of the input Jacobian `IGP.csv`) 
centered around the respective Jacobian value +-0.005, and reports the fraction of samples that led to a 
misprediction. For example, a returned value of 0.306 would indicate that 306 of the 1000 samples led to a
misprediction in the net effects matrix.



# Complete documentation
The complete documentation of the underlying python code can be found at: [http://math.oregonstate.edu/~koslickd/PressPurtCoreAlg/html/](http://math.oregonstate.edu/~koslickd/PressPurt/html/).

# Citations
1. Koslicki, D., & Novak, M. (2018). Exact probabilities for the indeterminacy of complex networks as perceived through press perturbations. Journal of mathematical biology, 76(4), 877-909. DOI: [ https://doi.org/10.1007/s00285-017-1163-0]( https://doi.org/10.1007/s00285-017-1163-0)
