\name{Rglimclim model objects}

\alias{GLC.modeldef}
\alias{print.GLC.modeldef}
\alias{plot.GLC.modeldef}
\alias{anova.GLC.modeldef}
\alias{logLik.GLC.modeldef}
\alias{summary.GLC.modeldef}
\alias{print.anova.GLC.modeldef}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Methods for Rglimclim model objects}
\description{Many of the objects produced by the \code{Rglimclim} package are lengthy lists of class \code{GLC.modeldef}, needed for interfacing with the underlying \code{Fortran} code. The list structures are not intended to be manipulated by the user; in general, the objects should be visualised and accessed using the functions documented here. To find out more about the individual components of a \code{GLC.modeldef} object, use the \code{\link{names}} command.
}
\usage{
\method{print}{GLC.modeldef}(x, scr.width=NULL, global.warn=TRUE, 
mean.only=FALSE, which.se="robust",...)
\method{plot}{GLC.modeldef}(x, which.plots = 1:2, sd.plots = TRUE, 
site.options = list(add.to.map = FALSE, scale = NULL, axis.labels = NULL, 
                    coord.cols = 1:2, site.labels="all"),
                    titles=TRUE,distance.units,plot.cols=gray(c(0.4,0)), ...)
\method{anova}{GLC.modeldef}(object, ..., include.constraints = FALSE, nonlin.warning = TRUE) 
\method{logLik}{GLC.modeldef}(object, ...)
\method{summary}{GLC.modeldef}(object, tables=c("month","site","year"), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x, object}{An object of (or inheriting from) class \code{GLC.modeldef}, for example resulting from a call to \code{\link{read.modeldef}} or \code{\link{GLCfit}}.}
  \item{scr.width}{Desired width of screen output. Labels are truncated so that the width of each line of output is \code{scr.width}. Defaults to \code{as.numeric(options("width")) - 2} - this should be adequate if \code{options("width")} is set correctly. If lines get wrapped on output, reduce this value accordingly.}
  \item{global.warn}{Controls whether or not warnings are issued if the user fails to define 'global' quantities that might be expected to be present for different choices of model - such as `trace' thresholds for logistic regression models (because these models might relate to daily rainfall for which such thresholds are usually defined.
}
  \item{mean.only}{For the \code{print} method, controls whether or not to restrict output just to the mean component of the model (i.e. to omit information about dispersion, spatial structure etc.}
  \item{tables}{Controls whether tables of Pearson residual summaries are produced by the summary method. The default is to produce summary tables broken down separately by month, site and year (in that order). 
}
  \item{which.se}{For the \code{print} method, controls whether \code{"naive"} or
  \code{"robust"} standard errors are printed. The naive ones are computed under the assumption that the model is correct and all sites are independent. Robust ones use a ``sandwich'' covariance matrix estimator which assumes only that residuals from different days are independent (this can be assured by incorporating an adequate representation of dependence on lagged values within the model). Default is \code{"robust"}.
  
  \strong{Note:} the calculation of robust covariance matrices currently does not take account of differences in dispersion parameter estimates, for models that have them. This is on my "to-do" list, and may remain there for some time. 
  }
  \item{which.plots}{A numeric vector controlling which residual plots are produced in the \code{plot} method. If the vector contains the value 1, plots will be produced by month; if it contains the value 2, plots will be produced by year; if it contains the value 3, plots will be produced by site; if it contains the value 4 then (where appropriate) a quantil-quantile plot of suitably-defined residuals will be prouced to check the distributional assumptions of the model; and if it contains the value 5 then a plot of estimated inter-site (residual or latent) correlations will be produced, with opacity used to indicate the sample size for each estimated correlation (correlations that are estimated very imprecisely are shown more faintly) and with any fitted correlation model superimposed. Defaults to \code{1:2} so that both monthly and annual plots are produced. 
  }
  \item{sd.plots}{A logical scalar controlling whether the \code{plot} method will plot the standard deviations of Pearson residuals, if monthly or annual summaries are requested. Defaults to \code{TRUE}.
  }
  \item{site.options}{A list controlling the behaviour of the \code{plot} method when producing residual plots by site (i.e. when the \code{which.plots} argument contains a 3). Possible entries are:
  \describe{ 
   \item{add.to.map}{A logical scalar. If \code{TRUE}, a bubble map of mean Pearson residuals will be added to an existing plot; if \code{FALSE}, a new plot will be created. The idea is that the user may wish to produce a topographic map of a study area (e.g. using \code{\link{image}}) and then overlay the residuals onto this map. \strong{Note:} there is an implicit assumption that if \code{add.to.map=TRUE}, the plot of residuals by site is the only plot being produced by this particular command: any attempt to add site residuals to another plot is likely to produce nonsense.}
   \item{scale}{A scalar controlling the relative sizes of the circles drawn. The
   default attempts to find a visually appealing scale, such that the circles are
   reasonably well separated (and hence distinguishable) but at the same time 
   clearly visible. The value of \code{scale} scales the radii of the circles 
   relative to this default (so that setting \code{scale=2} gives circles that
   have twice the radii of the default).}
   \item{axis.labels}{Labels for the horizontal and vertical axes of a bubble map showing mean residuals by site. Ignored if \code{add.to.map} is \code{TRUE}.}
   \item{coord.cols}{Indicates which columns of the \code{siteinfo} component of \code{x} should be taken as \code{x}- and \code{y}-axis co-ordinates respectively. Defaults to \code{1:2}. Note that this is appropriate for an `Eastings, Northing' co-ordinate system but not for `Latitude, Longitude'. For the latter, it is necessary to specify \code{coord.cols=c(2,1)} (assuming that the first two site attributes in \code{sitiinfo} represent the \code{x}- and \code{y}-coordinates of spatial location.}
   \item{site.labels}{Controls whether to print site codes at the centre of each bubble. If equal to \code{"all"} (the default), the site code will be printed for all bubbles. If equal to \code{"significant"}, codes will only be printed for bubbles where the mean residual differs significantly from zero at the 5\% level (this can be useful when there are many sites and you just want to identify the ones that are potentially problematic). Otherwise no codes will be printed on the map. For any other value, site labels will not be printed.}
  }
  \strong{Note:} the \code{plot} method currently makes some strong assumptions about the way in which it will be called --- e.g. that if a spatial plot is required and if \code{add.to.map} is \code{TRUE}, the spatial plot is the \emph{only} plot required. These assumptions may be relaxed in future versions. 
  }
  \item{titles}{For the \code{plot} method, controls whether or not titles are printed above each plot. The default is \code{TRUE}, although users may want to set it to \code{FALSE} in order to suppress printing of titles (e.g. to add their own plot titles later).
  }
  \item{distance.units}{For the \code{plot} method, when plotting inter-site correlations (if plot 5 is selected using \code{which.plots}) the default label for the x-axis is "Distance". If \code{distance.units} is set to a character scalar, the label will be changed: for example, if \code{distance.units="km"} then the label will be "Distance (km)".
  }
  \item{plot.cols}{A vector of two colour definitions, used to plot respectively the points and lines representing observations and modelled values if plots 4 or 5 are selected using \code{which.plots}.}
  \item{include.constraints}{For the \code{anova} method, indicates whether to include a summary of the constraints imposed at each stage as one moves through a sequence of nested models. Although the default is \code{FALSE}, it would be nice if users could set this to \code{TRUE} to check that model nestings are identified correctly by the code \ldots{}  
  }
  \item{nonlin.warning}{For the \code{anova} method, if Model A contains a term involving a nonlinear transformation of a covariate with parameters that must be estimated, and that term is missing from Model B, then the software will proceed as though Model B has a zero coefficient attached to that term and nonlinear parameter values equal to those in Model A. For most purposes this makes no difference whatsoever to anything, but it could make a small difference to the robust likelihood ratio test results. Therefore, by default, the software issues a warning whenever this occurs. To turn off the warning, set \code{nonlin.warning=FALSE}.
  
  }
  \item{...}{further arguments to be passed to or from other methods. They are ignored by the \code{print} and \code{logLik} methods here.
}
}
\details{The \code{print} method for objects of class \code{GLC.modeldef} produces a table showing the model structure and coefficient / parameter values, along with standard errors, test statistics and p-values where relevant if these are available (typically, they will be available if \code{model} is the result of a call to \code{\link{GLCfit}}). Standard errors are derived from a robust estimate of the parameter covariance matrix (a sandwich estimate accounting for inter-site dependence), \emph{unless} the model definition uses a spatial ``independence'' structure (which is the default).  %% 

The \code{summary} method extracts relevant diagnostic information and prints it. It returns a \code{NULL} value, invisibly.

The \code{plot} method produces plots of mean Pearson residuals for different subsets of the observations. \strong{Note} that at present this command is useful for model-building and checking, but it may not be suitable for the generation of publication quality graphics where finer control may be needed.

The \code{anova} method produces an analyis of variance / deviance table for comparing nested models, analogous to that from \code{\link{anova.lm}}. Any number of objects may be passed to the method; the code will sort them into order according to their degrees of freedom. Comparisons are based on chi-squared distributions for the log-likelihood ratio statistic, with dispersion parameters handled via their maximum likelihood estimates under an assumption of independence between sites. The ANOVA table reports both a "naive" test which assumes that the log-likelihoods are correctly specified; and a "robust" test which corrects for mis-specification due to inter-site dependence and other failures of model assumptions. The correction is based on the "vertical" adjustment proposed in Chandler and Bate (2007). See notes below for more information. 
}

\value{The \code{print} method returns the value of \code{model}, invisibly. 

The \code{logLik} method returns the `independence' log-likelihood (i.e. the log-likelihood calculated under the assumption that data from different sites are independent), with attributes representing the numbers of observations and parameters estimated (including dispersion parameters - where these have been estimated, maximum likelihood estimates are used for the log-likelihood calculations). 
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...

See comment above regarding the accuracy / fitness for purpose of this routine.
}
\references{
Chandler, R.E. and S. Bate (2007). Inference for clustered data using the independence loglikelihood. \emph{Biometrika}, \strong{94}, 167-183.
%% ~put references to the literature/web site here ~
}
\author{
Chiara Ambrosino and Richard Chandler (\email{r.chandler@ucl.ac.uk})
}
\note{
The degrees of freedom reported by the \code{summary} and \code{anova} methods differ slightly from the usual degrees of freedom associated with linear models: they count all parameters estimated, including constant terms and dispersion
parameters. The reason is that the technique used to adjust likelihood ratio tests for inter-site dependence (and other aspects of model mis-specification) cannot be used to adjust F-tests, and the software therefore reports tests
based on the chi-squared distribution with all parameters set at their
maximum likelihood values. The degrees of freedom therefore must correspond
to the number of parameters estimated via maximum likelihood. This is unlikely to have a major impact on results, because the datasets arising in this kind of application are typically so large that the chi-squared and F-tests will be almost exactly equivalent in any case. 

In the \code{anova} method, the robust likelihood ratio calculations currently do not take account of differences in dispersion parameter estimates, for models that have them. This is on my "to-do" list, and may remain there for some time. 

In the \code{anova} method, when comparing more than two models the adjustments required to obtain the "robust" likelihood ratio tests are all based on the robust covariance matrix from the largest model: the results when comparing models M2 with M3 (say) may be slightly different, therefore, depending on whether or not a larger model M1 is included simultaneously in the \code{anova} command. However, if the difference is \emph{very} large then this in itself is a sign that M1 produces very different results from  either M2 or M3. The implication is that neither M2 nor M3 captures all of the structure in the data, and hence that M1 should be preferred in such settings.  

%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{read.modeldef}}, \code{\link{GLCfit}}
}
\examples{
require(Rglimclim)
##
##       Define site information (see make.siteinfo for details)
##
data(GLCdemo)
##
##       Use the "SimpleModel" object supplied with the demo
##
print(names(SimpleModel))               # Shows names of underlying list
old.warn <- options()$warn
options(warn=1)                         # Show warnings as they occur
print(SimpleModel)                      # Show model in interpretable format
print(SimpleModel,global.warn=FALSE)    # Suppress warning message
options(warn=old.warn)                  # Reset handling of warnings
##
##      Diagnostics for a fitted model object
##
summary(DemoModel,tables=NULL)
opar <- par(no.readonly=TRUE)
par(mfrow=c(2,2))
plot(DemoModel,which.plots=1:2)
par(mfrow=c(1,1))
plot(DemoModel,which.plots=3,site.options=list(site.labels="none"),
     titles=FALSE)
title("My own plot title")
}

