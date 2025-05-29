# Rglimclim: a multisite, multivariate weather generator based on generalised linear models

This package is designed to fit and simulate Generalised Linear Models to daily climate sequences from a network of sites (e.g. weather stations, or model grid nodes). It can be used to analyse historical data and to provide simulations of future climate scenarios, for example to provide input to climate change impact assessment studies.

A key feature of the package is its ability to impute missing historical observations, and hence to quantify (some of) the uncertainty in properties of past climate.

The underpinning methodology is summarised in 

> Chandler, R.E. (2020). Multisite, multivariate weather generation based on generalised linear models. _Environmental Modelling and Software_ **134**, 104867. doi: [10.1016/j.envsoft.2020.104867](https://doi.org/10.1016/j.envsoft.2020.104867). 

## Requirements

The package requires [R version 4.4 or later](https://www.r-project.org/), as well as:

* The `devtools` package in R. If you don't have this already, it can be installed via the Tools menu in [RStudio](https://posit.co/download/rstudio-desktop/) or via `install.packages("devtools", lib=<whatever>)` from an R console).
* Tools for compiling packages from source, notably a `Fortran` compiler. `Windows` users will need the version of [Rtools](https://cran.r-project.org/bin/windows/Rtools/) appropriate to their version of `R`; `Mac` users may also need to ensure that the relevant [compilation tools](https://mac.r-project.org/tools/) are available, depending on their setup. `Linux` users will presumably know what they're doing. 

## Installation

With the tools above in place, `Rglimclim` can be installed using

```
library(devtools)
install_github("Richard-Chandler/Rglimclim")
```

If this succeeds, you are now ready to start. `help("Rglimclim-package")` may be a good entry point. There is also a worked example in Section 5 of the PDF package manual (available from the `R` help page), and a more extensive one in the "Supplementary data" for the [paper cited above](https://doi.org/10.1016/j.envsoft.2020.104867). 
