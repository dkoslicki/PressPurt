# PressPurtCoreAlg
This repository contains an implementation of all the results contained in the Koslicki &amp; Novak JMB paper [1].

# Installation
I suggest doing this in a [virtual environment](https://docs.python.org/3/library/venv.html). This repository is created for python3, so to create a virtual environment, enter these commands at the command line:

```bash
 python3 -m venv <name>
 source <name>/bin/activate
```
where ``<name>`` is the name of the folder you want to put the virtual environment in.

To install this repository, clone it and then install with pip
```bash
git clone https://github.com/dkoslicki/PressPurtCoreAlg.git
cd PressPurtCoreAlg
pip install -U pip
pip install -r requirements.txt
```
The command line programs are located in the ``scripts`` folder.

## R installation

PressPurtCoreAlg is also avaialable as an R package. You'll need to have python and conda installed. 

Install as an R package:

```
install.package("devtools")
devtools::install_github("dkoslicki/PressPurtCoreAlg", 
                         ref = "gibbond", build_vignettes = TRUE)
```

# Usage
In the ``scripts`` folder are a number of command line scripts that can be used.
The basic steps are:
1. Preprocess the matrix 
2. Compute statistics on it 
3. Visualize the results 

When any script is called, it will create a number of auxiliary files that are used by other scripts.
A user will be interested in the produced ``*.csv`` files, as these contain the statistical information about their network.


# Example
All these examples are run in the `scripts` directory.

For single entry perturbation statistics, use the following:
```bash
python PreprocessMatrix.py ../ExampleJacobians/IGP.csv ../ExampleJacobians
python ComputeEntryWisePerturbationExpectation.py ../ExampleJacobians
python GenerateEntryWiseFigures.py ../ExampleJacobians -a
```
For any script, you can run use the ``-h`` flag to see all the ways in which you can use the script, as each allows for custom options to be set. 
For example: ``python PreprocessMatrix.py -h``.

To compute quantitative stability statistics, use the following:
```bash
python ComputeQuantitativeSensitivity.py ../ExampleJacobians/IGP.csv ../ExampleJacobians/
python GenerateQuantitativeSensitivityFigure.py ../ExampleJacobians/
```

For multi-entry perturbation statistics, use:
```bash
python ComputeMultiEntryPerturbationExpectation.py ../ExampleJacobians/IGP.csv 
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