This repository contains a *python* implementation of all the results contained in the Koslicki &amp; Novak JMB paper [1].  See  [R_version](https://github.com/dkoslicki/PressPurt/tree/master/R_version) for an alternative implementation.

# PressPurt
- [Installation](#installation)
- [Usage](#usage)
  * [Quick Python Usage](#quick-python-usage)
  * [Detailed Python Usage](#detailed-python-usage)
    + [Single edge uncertainty](#single-edge-uncertainty)
      - [Qualitative analysis](#qualitative-analysis)
      - [Quantitative analysis](#quantitative-analysis)
    + [Multiple edge uncertainty](#multiple-edge-uncertainty)
      - [Qualitative analysis](#qualitative-analysis-1)
- [Complete documentation](#complete-documentation)
- [Citations](#citations)


# Installation
I suggest doing this in a [virtual environment](https://docs.python.org/3/library/venv.html). This repository is created for python3, so to create a virtual environment, enter these commands at the command line:

```bash
 python3 -m venv <name>
 source <name>/bin/activate
```
where ``<name>`` is the name of the folder you want to put the virtual environment in.

To install this repository, clone it and then install with pip
```bash
git clone https://github.com/dkoslicki/PressPurt.git
cd PressPurt/Python_version
pip install -U pip
pip install -r requirements.txt
```
The command line programs are located in the ``scripts`` folder.

# Usage
In the ``scripts`` folder are a number of command line scripts that can be used.
The basic steps are:
1. Preprocess the matrix 
2. Compute statistics on it 
3. Visualize the results 

When any script is called, it will create a number of auxiliary files that are used by other scripts.
A user will be interested in the produced ``*.csv`` files, as these contain the statistical information about their network.


##  Quick Python Usage
What follows is a quick example using some data included in this repository.

All these examples are run in the `Python_version/scripts` directory.

For single entry perturbation statistics, use the following:
```bash
python PreprocessMatrix.py ../ExampleJacobians/Modules/IGP.csv ../ExampleJacobians
python ComputeEntryWisePerturbationExpectation.py ../ExampleJacobians
python GenerateEntryWiseFigures.py ../ExampleJacobians -a
```
For any script, you can run use the ``-h`` flag to see all the ways in which you can use the script, as each allows for custom options to be set. 
For example: ``python PreprocessMatrix.py -h``.

To compute quantitative stability statistics, use the following:
```bash
python ComputeQuantitativeSensitivity.py ../ExampleJacobians/Modules/IGP.csv ../ExampleJacobians/
python GenerateQuantitativeSensitivityFigure.py ../ExampleJacobians/
```

For multi-entry perturbation statistics, use:
```bash
python ComputeMultiEntryPerturbationExpectation.py ../ExampleJacobians/Modules/IGP.csv 
```

## Detailed Python Usage

In general, there are three kinds of analyses:
1. Single edge uncertainty
    - Qualitative analysis
    - Quantitative analysis
2. Multiple edge uncertainty
    - Qualitative analysis

Each analysis requires a Jacobian matrix evaluated at an equilibrium point. This must be provided as
a CSV file containing column and row names that are either non-numeric ("A","B","C",...) or numeric
(but not including zero) entries.

The following examples will use data provided in this repository.

Before any analysis can be performed, you must first pre-process a matrix so that the intervals 
of asymptotic stability can be determined, as well as a calculation that will determine which 
uncertainty values on which edges will result in a sign change in the net effects matrix. This 
can be accomplished via:
```python
cd Python_version/ExampleJacobians  # go to the examples
mkdir Results  # create an output directory
python ../scripts/PreprocessMatrix.py Modules/IGP.csv Results  # preprocess the IGP.csv example and store the results in the Results folder  
```
You will notice that the `Results` folder now contains a number of files that are required for downstream analysis.

Note that you need only preprocess the Jacobian once, as changing parameters in the following analyses does not require 
you to re-run `PreprocessMatrix.py`.

### Single edge uncertainty
This set of analyses are for the situation where you wish to assess how the presence of uncertainty 
in the strength of the direct effect (edge weight) between a single focal pair of species (nodes) alters 
the signs of the predicted press perturbation net effects between all pairs of species in the network.

#### Qualitative analysis
The script `ComputeEntryWisePerturbationExpectation.py` will compute the expected number of sign 
switches/mispredictions in the net effects matrix when perturbing each entry of the Jacobian matrix individually.
The script can be run with:
```python
python ../scripts/ComputeEntryWisePerturbationExpectation.py --distribution_type 'truncnorm' -a 0 -b -2 Results  # run the script on the preprocessed matrix contained in the Results folder assuming the edge uncertainties are distributed according to a truncated normal distribution using a mean of 0 and a variance equal to (the length of the interval divided by the absolute value of the input parameter)^2
```
The results will also be placed in the `Results` subfolder. The file ``Results/expected_num_switch.csv`` 
is a matrix whose (i,j) entry shows the expected fraction of sign switches/mispredictions in the net effects 
matrix when the edge between nodes i and j experiences uncertainty according to the distribution specified.

Use `ComputeEntryWisePerturbationExpectation.py -h` to view the different kinds of distributions available and a more detailed 
description of the flags `-a` and `-b`.

You can then visualize the results using the script `GenerateEntryWiseFigures.py`:
```python
python ../scripts/GenerateEntryWiseFigures.py Results --list_of_numswitch_to_plot 1 1 1 2
```
This command generates a heatmap of the `expected_num_switch.csv` file as well as a plot of the (edge uncertainty value) versus 
(number of sign switches/mispredictions) in the net effects matrix overlaid with the edge uncertainty distribution (truncated by the 
region of asymptotic stability). The flag `--list_of_numswitch_to_plot 1 1 1 2` specifies that only the (1,1) and (1,2) edges should be shown. 
If you use `--all_numswitch_plot` instead, this will show similar plots for *all* edges (warning: this may be uninterpretable if your matrix is large). 
This script also calculates the expected percentage of mispredictions when each edge individually experiences uncertainty 
(i.e. averaging the matrix `Results/expected_num_switch.csv`)

#### Quantitative analysis
Instead of asking whether a qualitative misprediction occurs, one might be interested in quantifying by how much you mispredict 
the change in abundance of other nodes when edge weights (individually) are under uncertainty. The script `ComputeQuantitativeSensitivity.py` allows you to do this. However, a decision needs to be made 
about what level of uncertainty the edge interactions have experienced (i.e. have they been changed to their mean uncertainty value, 
median value, maximal value, etc). The script `ComputeQuantitativeSensitivity.py` takes a "worst case" analysis and assumes 
that each edge (individually) has experienced a limited value of uncertainty (i.e. the maximal/limiting value). The script 
is called with:
```python
python ../scripts/ComputeQuantitativeSensitivity.py Modules/IGP.csv Results
```  
This runs the computation on the `IGP.csv` Jacobian and stores the results in the `Results` subfolder. 
This folder will now contain a file called `MRS.csv` which gives the mean (averaged over all individually perturbed entries) 
factor by which your net effects matrix predictions will be off when an (averaged) single edge edge experiences maximal uncertainty. 
The file `quantitative_sensitivity.csv` also depicts the mean factor by which the net effects matrix predictions will be off by when a 
single edge (i,j) experiences maximal/limiting uncertainty.

The script `GenerateQuantitativeSensitivityFigure.py` visualizes the quantitative analysis by creating a heat map of the 
`quantitative_sensitivity.csv` file:
 ```python
python ../scripts/GenerateQuantitativeSensitivityFigure.py Results
```
The value returned is the same as the one contained in `MRS.csv`.

### Multiple edge uncertainty

While the above analysis computes the sensitivity of the net effects matrix when edges individually experience
uncertainty, the following allows multiple edges to experience uncertainty simultaneously.

#### Qualitative analysis
The script `ComputeMultiEntryPerturbationExpectation.py` performs a Monte Carlo sampling of edge uncertainties according 
to a uniform distribution (within a specified interval length) and reports the average fraction of sampled uncertainties 
that lead to a sign switch/misprediction in the net effects matrix (considering only uncertainties that result in a stable Jacobian).

The command can be called with
```python
python ../scripts/ComputeMultiEntryPerturbationExpectation.py Modules/IGP.csv --num_iterates 1000 --interval_length 0.01 --threads 4
```
The above command samples uncertainties 1000 times on all edges (non zero entries of the input Jacobian `IGP.csv`) 
centered around the respective Jacobian value +-0.005, and reports the fraction of samples that led to a 
misprediction. For example, a returned value of 0.306 would indicate that 306 of the 1000 samples led 
to a misprediction in the net effects matrix.


# Complete documentation
The complete documentation can be found at: [http://math.oregonstate.edu/~koslickd/PressPurtCoreAlg/html/](http://math.oregonstate.edu/~koslickd/PressPurtCoreAlg/html/).

# Citations
1. Koslicki, D., & Novak, M. (2018). Exact probabilities for the indeterminacy of complex networks as perceived through press perturbations. Journal of mathematical biology, 76(4), 877-909. DOI: [ https://doi.org/10.1007/s00285-017-1163-0]( https://doi.org/10.1007/s00285-017-1163-0)
