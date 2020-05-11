# Resubmission

This is the resubmission addressing the comments after the first submission of `PressPurt`


# Addressing CRAN comments:

* DESCRIPTION:
    * Description: added Koslicki, D., & Novak, M.
      (2017)<doi:10.1007/s00285-017-1163-0>
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
* Did not find any variables named `T` or `F`

# First time submission

This is the first submission for `PressPurt`

# Test environments

* No ERRORs, WARNINGs or NOTEs when `devtools::check()` was run
* Tested locally on:
  * R version 3.6.2 (2019-12-12), Platform: x86_64-pc-linux-gnu (64-bit), Running under: CentOS Linux 7 (Core)
  * R version 3.6.1 (2019-07-05), Platform: x86_64-pc-linux-gnu (64-bit), Running under: Ubuntu 18.04.2 LTS
  * R version 3.5.2 (2018-12-20), Platform: x86_64-apple-darwin15.6.0 (64-bit), Running under: macOS High Sierra 10.13.6
* Tested with:
  * `devtools::check_win_release()`: 1 NOTE: New submission
  * `devtools::check_win_devel()`: 1 NOTE: New submission
  * `devtools::check_win_oldrelease()`: 1 NOTE: New submission
* Tested using rhub:
  * `rhub::check_for_cran()`: 1 NOTE: New submission

# Downstream Dependencies

This package is dependent on python, using `reticulate`, which is why the vignettes aren't built in this submission.