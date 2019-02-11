# Contributing

_A quick overview of the contents of this repository and how to contribute._


## Repository File Structure

### `R/` Folder

This is where functions are placed in order to be documented (i.e. have documentation pop up when a user uses ?name_of_function) and to be installed as part of the package.

### `man/` Folder

This is an auto-generated folder that contains (also auto-generated) documentation for functions and other objects. Discouraged to manually edit these files; instead, do so via the ROxygen methods detailed below.


### CTMv6

Version 6.0 of the Cohort Theory Model. Currently built in Excel VBA (visual basic for applications)

_work needed_: VBA developers translating into R. Currently [in progress by JH](https://github.com/tilbud/MEMs/tree/master/CTMv6/R)

### CTMv5

Version 5.0 of the Cohort Theory Model. Translated into R, still under development.

_work needed:_ R developers [improving and expanding current code](https://github.com/tilbud/MEMs/tree/master/CTMv5/R), further documentation needed


## Documentation

Descriptions of functions, function parameters, etc. stored in `man/` directory.


Make sure that (if you’re using Studio) under Build > Configure Build Tools, “Generate documentation with ROxygen is checked.

*Resources*:  
[Hadley Wickham, R packages objection documentation](http://r-pkgs.had.co.nz/man.html)


### Rbuildignore

List of files and directories that are NOT to be included in the built package.This means that we can include documentation, folders of other data files, etc. that may be helpful for communication and for developing the package, but that we don’t want (or can’t) have in the package itself.

Standard notation is to include a ^ in front of the file and its extension ([see here](http://r-pkgs.had.co.nz/package.html), go to “Bundled Packages” for more)


### Namespace file

Current template is a placeholder. [See here for more](https://stat.ethz.ch/pipermail/r-package-devel/2016q2/000862.html).


## Building Packages

Ultimately, we will likely want to generate bundled and binary versions of the package. The standard package that most R users are familiar with (those that are installed via `install.packages()`) are binary packages. [See here](http://r-pkgs.had.co.nz/package.html) for a description of each of these package versions.

Build source: can be conducted via Build > Build Source Package
Build binary: can be conducted via Build > Build Binary Package
