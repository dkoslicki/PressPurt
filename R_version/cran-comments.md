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