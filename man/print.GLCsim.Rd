\name{Rglimclim simulation objects}

\alias{print.GLCsim}
\alias{summary.GLCsim}
\alias{print.summary.GLCsim}
\alias{plot.summary.GLCsim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Methods for Rglimclim simulation objects}
\description{Simulations produced by the \code{\link{GLCsim}} routine are stored in separate files, which may be large especially if daily data are stored. Information about each simulation is stored in a small \R object of class \code{GLCsim}; the methods documented here are intended to provide the user with a quick means of summarising and visualising simulation results. 
}
\usage{
\method{print}{GLCsim}(x, ...)
\method{summary}{GLCsim}(object,which.variables,which.sites,which.regions,
                           which.timescales,thresholds,season.defs,...)
\method{print}{summary.GLCsim}(x, ...)
\method{plot}{summary.GLCsim}(x,imputation,quantiles,
                                which.variables,which.sites,which.regions,
                                which.timescales,which.stats,which.seasons,
                                plot.titles,ylabs,colours.sim="greyscale",
                                colour.obs="black",...)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x, object}{An object of class \code{GLCsim}, resulting from a call to \code{\link{GLCsim}} or \code{summary.GLCsim} (see below).}
  \item{which.variables}{A vector selecting the variables for which to produce plots and summaries. If omitted, plots and summaries will be produced for all available variables. \strong{Note} that the specification of \code{which.variables} differs between the \code{GLCsim} methods and the \code{summary.GLCsim} method (currently just \code{plot.summary.GLCsim}). For the former, \code{which.variables} should be a numeric vector of variable numbers, the numbering corresponding to that in the original data file from which the simulations were initialised. For the latter, \code{which.variables} should be a \emph{character} vector of variable names as stored in the relevant \code{dimnames} attribute and revealed by the \code{print} method for the object.}
  \item{which.sites}{An optional character vector of 4-character site codes. If supplied, plots and summaries will be restricted to the corresponding sites; otherwise they will be produced for all sites.}
  \item{which.regions}{An optional vector selecting regions for which to produce plots or summaries. If omitted, plots and summaries will be produced for all available regions. Like \code{which.variables}, the specification differs between the \code{GLCsim} and \code{summary.GLCsim} methods: for the former, it should be a numeric vector of region numbers (0 corresponding to the whole area), and for the latter it should be a character vector of region names. Regions must be defined in the \code{siteinfo} component of the simulation object; use \code{\link{read.regiondef}} or \code{\link{define.regions}} to achieve this.} 
  \item{which.timescales}{Either \code{"daily"}, \code{"monthly"} or \code{c("daily","monthly")}. This selects the simulation output files to be processed. If omitted, all available files will be processed.}
  \item{thresholds}{A numeric vector the same length as \code{which.variables}. If present and if the element corresponding to a particular variable is non-\code{NA}, the \code{summary} method will calculate the proportion of exceedances of the corresponding threshold for that variable, as well as the mean and standard deviation conditional on threshold exceedance.}
  \item{season.defs}{A list of numeric vectors, defining months to be grouped together to form "seasons" when processing monthly output files. If not specified, separate summaries are produced for each month of the year. 
  }
  \item{...}{Other arguments to generic methods.}
  \item{imputation}{For the \code{plot} method, an optional object of class \code{summary.GLCsim} that is treated as containing summaries from a set of imputations conditioned on all available data. If supplied, an imputation envelope will be overlaid on the plots.}
  \item{quantiles}{For the \code{plot} methods, a vector of quantiles controlling which quantiles of simulated distributions will be plotted. Defaults to \code{c(0,0.01,0.05,0.1,0.25,0.5,0.75,0.9,0.95,0.99,1)}.}
  \item{which.stats}{For the \code{plot} methods, a character vector specifying the summary statistics to be plotted for daily data (if requested via \code{which.timescales}). These should exactly match the names of the statistics as output by the \code{summary} method, and as stored in the \code{dimnames} attribute of the corresponding elements of a \code{summary.GLCsim} object (see "Value" section below).}
  \item{which.seasons}{For the \code{plot} methods, a character vector specifying the names of the seasons for which plots are to be produced. Again, these names should match the \code{dimnames} from the corresponding \code{summary.GLCsim} object.}
  \item{plot.titles}{For the \code{plot} methods, a character vector of titles for the plots. It is the user's responsibility to ensure that the length of this vector matches the number of plots actually produced. If omitted, the routine will construct plot titles automatically.}
  \item{ylabs}{For the \code{plot} methods, a character vector of y-axis labels for the plots. Similar comments apply as for \code{plot.titles}.}
  \item{colours.sim}{Either \code{"greyscale"} (the default), \code{"colour"} or a vector of valid colour specifiers, of length \code{length(quantiles)-1}. These colours will be used to shade the simulated distributions. \code{"greyscale"} produces plots that are suitable for inclusion in printed material; \code{"colour"} uses a default colour scale generated using the \code{\link{rainbow}} command.}
  \item{colour.obs}{The colour to use for plotting an imputation envelope if required. Defaults to \code{"black"}.}
}
\details{The \code{print} method for the \code{GLCsim} class produces a summary of the simulation settings: input file name, which variables were simulated, the simulation period and numbers of simulations, as well as details of which output has been stored and where the output files are located. 

The \code{summary} method is used to produce summary statistics that can be plotted or used for further analysis. The objects that it produces are lists (see "Value" section below) that have their own \code{summary.GLCsim} class.

The \code{print} method for the \code{summary.GLCsim} class produces a concise printout of the summaries that have been calculated and can be useful for (e.g.) finding the relevant names that should be used in a call to the \code{plot} method for this class (i.e. the available choices for the \code{which.variables}, \code{which.regions} arguments etc.). 

There is no \code{plot} method for objects of class \code{GLCsim}; plotting is done on objects of class \code{summary.GLCsim} (see "Value" section below). Usually, the most efficient way to plot results of a simulation is precompute all summary statistics via a single call to the \code{summary} method, and then to use multiple calls to \code{plot.summary.GLCsim} to produce plots for subsets of these statistics without requiring any further computations. It may be useful to know the order in which this \code{plot} method produces its plots: it loops over variables and, for each variable, first produces all plots for daily summary statistics (if requested) and then the plots for monthly or seasonal means. For the daily summary statistics, the routine first loops over all sites for which plots have been requested (and, within this, over all statistics that have been requested) and then over all regions. For monthly or seasonal means, the routine just loops over the seasons that have been defined. 

If an \code{imputation} argument is supplied to the \code{plot} method, a band will be superimposed on the plots of simulated distributions, showing the range of values found within the \code{imputation} object. This is useful to assess whether the simulations can capture the observed behaviour, taking account of uncertainty in the observations due to missing values. It is worth noting that if the \code{imputation} object is derived from a simulation containing 39 imputations, then the resulting band will be a 95\% confidence interval for the actual values of the quantities plotted (because if we could pool the actual value with the 39 imputations, there would be a 1/40 chance that the actual value would be the maximum of the pooled sample and a 1/40 chance that it would be the minimum). 
}

\value{The \code{print} and \code{plot} methods produce \code{NULL} values. 

The \code{summary} method produces a list object of class \code{summary.GLCsim}, with components \code{Daily} and \code{Monthly} (which are \code{NULL} unless the corresponding summaries have been selected via the \code{which.timescales} argument). The \code{Daily} component itself contains components \code{Sites} and \code{Regions}, each of which is a 5-dimensional array of summary statistics: the dimensions are \code{Simulation} (simulation number), \code{Month}, \code{Site} or \code{Region}, \code{Variable} and \code{Statistic}. The \code{Monthly} component is a 5-dimensional array of seasonal means for selected regions, with dimensions \code{Simulation}, \code{Region}, \code{Variable}, \code{Year} and \code{Season (months)}; the individual \code{dimnames} for the \code{Season (months)} dimension are of the form (for example) \code{"XXX, YYY, ZZZ"} where \code{XXX}, \code{YYY} and \code{ZZZ} are elements of \code{\link{month.abb}}. Note that (a) the seasonal means are computed as a straight average of the corresponding monthly values (so they do not account for different numbers of days in each month) (b) if a season is defined via \code{season.defs} as, for example, \code{c(12,1,2)} then, for year \code{YYYY}, the routine will calculate the mean of each variable for month 12 in year \code{YYYY} and months 1 and 2 in year \code{YYYY+1}. 
}
%%\references{
%% ~put references to the literature/web site here ~
%%}
\author{
Richard Chandler (\email{r.chandler@ucl.ac.uk})
}
\note{Because the \code{summary} method produces objects with their own print method, the user cannot see the values of the summaries simply by calling \code{print} on a \code{summary} object. To see the contents of a \code{summary.GLCsim} object in its entirety, use \code{\link{unclass}}. Specific parts of the object can be accessed directly using the component names. 

A potential source of confusion with these routines is that the \code{which.variables} and \code{which.regions} are specified differently in the \code{GLCsim} and \code{summary.GLCsim} methods. In the \code{GLCsim} methods, they should be supplied as numeric vectors; in the \code{summary.GLCsim} method(s) they should be supplied as character vectors, corresponding to the required \code{dimnames} of the \code{GLCsim} object to which they are applied. 
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{GLCsim}.}
}
%%\examples{
%% Need an example involving the extraction of a time series of maxima for an EV analysis. 
%%}

