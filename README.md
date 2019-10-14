# PressPurt
This repository contains an implementation of all the results contained in the Koslicki &amp; Novak JMB paper [1].


PressPurt is a computational package for identifying the interactions (edges) within a network 
whose strengths (edge weights) must be estimated most accurately in order to produce qualitatively 
robust predictions of the network's response to press perturbations. The package provides methods 
for calculating and visualizing these edge-specific sensitivities (tolerances) when estimation-uncertainty 
is associated to one or more edges according to a variety of different error distributions
(e.g., uniform, truncated normal). The software requires the network to be described by a system of 
differential equations and only requires as input a numerical Jacobian matrix and a specification of the 
presumed error distribution.

Use cases include the study of food web dynamics in the context of fisheries management, asking: 
If the harvest of a focal species were to be increased, how robust would qualitative predictions of net 
change in the abundance of all other species be given that their estimated interaction strength are 
all subject to uncertainty? Which interaction strengths would need to be quantified most accurately for 
predictions to be robust?

# About

There are two flavors of the code (in case you find one easier to work with than the other):
1. [Python](https://github.com/dkoslicki/PressPurt/tree/master/Python_version)
2. [R](https://github.com/dkoslicki/PressPurt/tree/master/R_version)

Please click on the folders [Python_version](https://github.com/dkoslicki/PressPurt/tree/master/Python_version) or [R_version](https://github.com/dkoslicki/PressPurt/tree/master/R_version) to see the installation and usage instructions for each.

# Citations
1. Koslicki, D., & Novak, M. (2018). Exact probabilities for the indeterminacy of complex networks as perceived through press perturbations. Journal of mathematical biology, 76(4), 877-909. DOI: [ https://doi.org/10.1007/s00285-017-1163-0]( https://doi.org/10.1007/s00285-017-1163-0)
