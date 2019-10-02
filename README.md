# PressPurt
This repository contains an implementation of all the results contained in the Koslicki &amp; Novak JMB paper [1].

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

# About

There are two flavors of the code (in case you find one easier to work with than the other):
1. [Python](https://github.com/dkoslicki/PressPurt/tree/master/Python_version)
2. [R](https://github.com/dkoslicki/PressPurt/tree/master/R_version)

Please click on the folders [Python_version](https://github.com/dkoslicki/PressPurt/tree/master/Python_version) or [R_version](https://github.com/dkoslicki/PressPurt/tree/master/R_version) to see the installation and usage instructions for each.

# Citations
1. Koslicki, D., & Novak, M. (2018). Exact probabilities for the indeterminacy of complex networks as perceived through press perturbations. Journal of mathematical biology, 76(4), 877-909. DOI: [ https://doi.org/10.1007/s00285-017-1163-0]( https://doi.org/10.1007/s00285-017-1163-0)