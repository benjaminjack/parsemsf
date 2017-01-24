# parsemsf - Parse Thermo MSF files and estimate protein abundances

The main purpose of the ParseMSF package is to parse proprietary Thermo MSF files into a format readable by R. This package makes it easy to view individual peptide information, including peak areas, and to map peptides to locations within the parent protein sequence. This package also estimates protein abundances from peak areas and across multiple technical replicates.

Currently, the ParseMSF package provides functions for parsing Thermo MSF files produced by Proteome Discoverer 1.4.x only.

# Installation

This latest stable release of this package can be installed from [CRAN](https://cran.r-project.org/package=parsemsf) by running the following command in the R console:

```
install.packages("parsemsf")
```

If you want to install the development version this package from Github by running the following command in the R console:

```
devtools::install_github("benjaminjack/parsemsf")
```

# Usage

To get an introduction to the functions in ParseMSF, please see the [Introduction vignette on CRAN](https://cran.r-project.org/package=parsemsf).

# Acknowledgements

Some of the SQL queries in this package come from the now-defunct [paRseMSF package](https://github.com/ashokapol/parsemsf) by Ashoka Polpitiya.
