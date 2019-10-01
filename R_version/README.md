# PressPurtCoreAlg

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
This repository contains a *R* implementation of all the results contained in the Koslicki &amp; Novak JMB paper [1].

In short, this is a computational package designed to identify the most sensitive interactions 
within a network which must be estimated most accurately in order to produce qualitatively 
robust predictions to a press perturbation. The package produces data and visualization when 
uncertainty is associated to one or more edges in the network and according to a variety of 
distributions. The software requires the network to be described by a system of differential 
equations but only requires as input a numerical Jacobian matrix.

Use cases include modeling a food web using a system of differential equations and asking the question: 
if I were cause a sustained increase in abundance of one organism, how robust would my predictions 
be about the net change in abundance of other organisms given that the interaction strengths between 
organisms is uncertain? Which edge interaction strengths would I need to accurately quantify in order for 
my predictions to be robust?


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

If the dependencies don’t install correctly see below for instructions on how to install via the command-line.

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
library(PressPurtCoreAlg)
# You should see "test-r" listed as a conda env
find_python()
# Set your virtualenv
set_python_virtual(virtualenv = "PressPurt")
```









## R Usage

### Quick start

* Load the library
* Check available python versions
* Set python version and conda environment. You can set up a new conda env or use a previous one.

```
library(PressPurtCoreAlg)
find_python()
set_python("r-reticulate", verbose = TRUE)
```

* install python dependencies

```
py_depend("r-reticulate")
```

* Run PreprocessMatrix
* Run ComputeEntryWisePerturbationExpectation


```
# Path to your matrix
infile <- "../ExampleJacobians/Modules/IGP.csv"
PreProsMatrix <- PreprocessMatrix(input_file = infile, output_folder = NULL, max_bound = 10, threads = 2)
Entrywise <- ComputeEntryWisePerturbationExpectation(PreProsMatrix = PreProsMatrix, 
                                        distribution_type="truncnorm", 
                                        input_a=0, input_b=-2, threads=1)
```

* Plot
    * specific figures
    * all figures

```
list_of_numswitch_to_plot <- list(c(0, 0), c(0, 1))
GenerateEntryWiseFigures(EntryWise=Entrywise, 
                         all_numswitch_plots = FALSE, 
                         list_of_numswitch_to_plot=list_of_numswitch_to_plot)
GenerateEntryWiseFigures(EntryWise=Entrywise, 
                         all_numswitch_plots = TRUE)
```


# Complete documentation
The complete documentation can be found at: [http://math.oregonstate.edu/~koslickd/PressPurtCoreAlg/html/](http://math.oregonstate.edu/~koslickd/PressPurtCoreAlg/html/).

# Citations
1. Koslicki, D., & Novak, M. (2018). Exact probabilities for the indeterminacy of complex networks as perceived through press perturbations. Journal of mathematical biology, 76(4), 877-909. DOI: [ https://doi.org/10.1007/s00285-017-1163-0]( https://doi.org/10.1007/s00285-017-1163-0)