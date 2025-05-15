\name{Rglimclim-package}
\alias{Rglimclim-package}
\alias{Rglimclim}
\docType{package}
\title{
Generalized Linear Modelling for Daily Climate Time Series
}
\description{
Modelling and simulation of daily multi-site time series, based on generalised 
linear models and designed in particular for climatological downscaling 
applications. The references below give an overview of the theory upon which 
much of this work is based. 

The mathematical specification of the models can be found in Section 3 and in 
the Appendix of the \href{../doc/manual.pdf}{PDF package manual}.
}
\details{This package is designed to allow flexibility in fitting and simulating fairly complicated models for daily climate/weather data, and to avoid some of the tedious manipulation of data files that would be necessary if standard software packages were used to perform this kind of analysis. The package is written to take advantage of the fact that all exercises involving fitting, and simulating, GLMs for daily weather series must inevitably share common features, as follows:
\itemize{
\item{Models are to be fitted over a network of sites, for which various
attributes (e.g location) are known. Flexibility in the placement and
known attributes of sites can be achieved via the use of a database
containing site information - this can be referenced as required by the fitting and simulation programs
}
\item{Typically, possible covariates fall into a small number of
categories, as follows:
\itemize{
\item{A constant term}
\item{Site effects}
\item{`Year' effects (e.g. long-term trends)}
\item{`Month' effects (e.g. seasonality)}
\item{`Day' effects (e.g. day-to-day temporal autocorrelation)}
\item{Interactions}
}
The software exploits this small number of categories, and treats each
separately. In each category a variety of choices can be made
regarding parametrisation by selecting from a `menu' of choices. The
software is, hopefully, written in a sufficiently modular way that
users can customise these menus if they so desire. 

`External' effects (such as ENSO or the North Atlantic Oscillation)  
are dealt with under the appropriate timescale - for example, if 
you wanted to use a monthly ENSO index as a covariate in your model, 
this would be counted as a `monthly' effect; if you wanted to use 
a `winter NAO' series (one value per year), it would count as a `yearly' 
effect.
}
}

Climate datasets often have other unusual features relating to 
measurement methods. For example, in daily raingauge data any non-zero
amount that is less than some small threshold may be recorded as a 
`trace' amount, because it is too small to be measured accurately. Such
features pose potential problems for statistical analysis. The software
aims to provide methods for dealing with them.

The package is essentially an \R port of the original \code{GLIMCLIM} package (Chandler, 2002) which was written in \code{Fortran 77}. The \code{Fortran} code for model fitting and simulation has been retained since it is considerably faster than, for example, the \code{\link{glm}} command in \R itself. However, the \R interface makes it much easier to save and update models, and to transfer and visualise information. In particular, \code{GLIMCLIM} made extensive use of definition files to define sites and models. All of these had to be prepared manually, with a consequent risk of error. \code{Rglimclim} contains routines that will read old \code{GLIMCLIM} definition files; however, it is intended that the primary means of defining sites, models etc. will be via manipulation of \R objects, and hence that definition files will eventually become defunct. 
}

\seealso{See the help for \code{\link{GLCfit}} (the main model fitting routine) to get started. A full description of the package can be found in the \href{../doc/manual.pdf}{PDF package manual}.
}
\author{Richard E. Chandler

Maintainer: Richard E. Chandler <r.chandler@ucl.ac.uk>
}
\references{
Chandler, R.E. (2020). Multisite, multivariate weather generation based on generalised linear models. \emph{Environmental Modelling and Software} \bold{134}: 104867. \href{https://dx.doi.org/10.1016/j.envsoft.2020.104867}{doi: 10.1016/j.envsoft.2020.104867}.

Chandler, R.E. (2005). On the use of generalized linear models for interpreting climate variability. \emph{Environmetrics} \bold{16(7)}: 699-715.

Chandler, R.E. and Bate, S. (2007). Inference for clustered data using the independence log-likelihood. \emph{Biometrika} \bold{94}: 167-183.

Yan, Z., Bate, S., Chandler, R.E., Isham, V. and Wheater, H. (2002). An analysis of daily maximum windspeed in northwestern Europe using generalized linear models. \emph{J. Climate} \bold{15(15)}: 2073-2088.

Yang, C., Chandler, R.E., Isham, V. and Wheater, H.S. (2005). Spatial-temporal rainfall simulation using Generalized Linear Models. \emph{Water Resources Research} \bold{41}, doi:10.1029/2004WR003739. 
}
\keyword{ package }

