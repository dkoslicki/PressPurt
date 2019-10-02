# Resources:
## Building R packages
### https://kbroman.org/pkg_primer/
### http://r-pkgs.had.co.nz/
### https://hilaryparker.com/2014/04/29/writing-an-r-package-from-scratch/
## background on RStudio projects
### https://r4ds.had.co.nz/workflow-projects.html

# For Continued package development:
## From RStudio:
### Open the R project: file -> open project -> navigate to project folder -> click on PressPurt.Rproj

## Important files/folders (from project top directory)
### R scripts: R/
### python scripts: inst/python
### example data: inst/extdata/
### directory of R documentation files for each function (.Rd): man/
### vignettes: vignettes/
### NAMESPACE: what to import/export, edited automatically
### DESCRIPTION: Information about the package

## Summary of commands
### load the package (do this first):
devtools::load_all()
### Now you have access to the package's functions and docs
### Update documentation & NAMESPACE
devtools::document()
### Install the package to your local machine from the project directory
devtools::install()
### Build the package to make a .tar.gz file
devtools::build()
### To check for CRAN errors/warnings/notes
devtools::check()

## Documentation
### The documentation is created from the text you have above an exported function.
### You need a title, write a Description,
### specify @param, @export, @import (which packages to import),
### @examples (you can show examples of how to use a function)
#### In this case, you don't want to actually run the examples
#### on cran because they all require python to run.
### There are other options as well.
### See any of the .R scripts for examples



