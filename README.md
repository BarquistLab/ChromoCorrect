
# ChromoCorrect

<!-- badges: start -->
<!-- badges: end -->

The goal of ChromoCorrect is to ...

## Installation

You can install the development version of ChromoCorrect from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("gerisullivan/ChromoCorrect")
```

## Example

# If current working directory is folder with all .csv files:
readcounts <- structure_rc()

# If wanting to specify path:
readcounts <- structure_rc(csvpath = "~/path/to/files")

# If not wanting locus information:
readcounts <- structure_rc(csvpath = "~/path/to/files", getLocusInfo = FALSE)

# If column names for control are "MH_1" and "MH_2":
normalise_CB(x, control = "MH", windowSize = "auto", minrc = 10, \cr writePlots = TRUE, locusInfo = TRUE, path = "~/path/to/files")

``` r
library(ChromoCorrect)
## basic example code
```

