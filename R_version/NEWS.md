# Version 1.0.1

**Date:** 05/07/2020

## Updated per CRAN comments

* DESCRIPTION: 
    * Description: added Koslicki, D., & Novak, M. (2017)<doi:10.1007/s00285-017-1163-0>
    * Version: 1.0.0 -> 1.0.1
* Added `@return None` -> `\value{None}` to:
    * find_python
    * create_conda_env
    * set_python_conda
    * create_virtual_env
    * set_python_virtual
    * py_depend
* Added @return to: 
    * get_distributions_single
    * ns_to_step
    * GenerateEntryWiseFigures
* Changed `cat`/`print` -> `message` at the following lines in `helper_funcs.R`:
    * 25, 32, 35, 38, 40, 73, 80, 84, 88, 93, 94, 122, 125, 130, 
    * 164, 171, 175, 179, 183, 184, 211, 214, 219, 318, 585, 643
* Changed `T -> TRUE` or `F -> FALSE`
    * ComputeMultiEntryPerturbationExpectation.R:
        * Lines: 27, 31, 34
    * ComputeEntryWisePerturbationExpectation.R:
        * Lines: 55, 58, 61, 64
    * helper_funcs.R:
        * Lines: 77, 86, 127, 168, 177, 216, 321, 386, 395, 400, 417, 637
        * Lines: 305, 640
    * PreprocessMatrix.R:
        * Lines: 34, 37
* Update NEWS.md and cran-comments.md

# PressPurt 1.0.0

Submitting version 1.0.0 to Cran on 2/21/2020.
