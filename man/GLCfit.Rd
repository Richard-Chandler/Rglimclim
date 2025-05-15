\name{GLCfit}
\alias{GLCfit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Fits a variety of generalized models to daily climate time series}
%
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This routine fits generalised linear models (GLMs) to daily climate time series. It is designed to allow flexibility in fitting fairly complicated models to daily climate/weather data, and to avoid some of the tedious manipulation of data files which would be necessary if standard software packages were used to perform this kind of analysis. It also has some nonlinear fitting capabilities that are not available in standard GLM routines such as \code{\link{glm}}.
}
\usage{
GLCfit(model.def, dispersion.def,response.check=1, max.iter=Inf,  
       data.file="gaugvals.dat", missval=-99.99, which.response=1, 
       nprev.required, nvar.check=TRUE, allow.incomplete.averages=FALSE,
       siteinfo, model.type, 
       external.files = c("yr_preds.dat", "mn_preds.dat", "dy_preds.dat"), 
       diagnostics = 2, verbosity=2, resid.file, cor.file)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{model.def}{
        An object of class \code{GLC.modeldef}, created for example
        via a call to \code{\link{read.modeldef}}.
 }
 \item{dispersion.def}{
        An object of class \code{GLC.modeldef} defining the dispersion
        model structure when \code{model="normal-heteroscedastic"} (see
        documentation for \code{model} argument below).
 }
 \item{response.check}{
       Code to determine whether to check the values of the response
        variable. Options:
        \describe{
          \item{-1}{No checking performed}
          \item{0}{Checking performed, but invalid values will lead to 
          warnings rather than errors; the routine will continue to 
          attempt model fitting}
          \item{1}{Checking performed and response variables 
                calculated as follows:
                \itemize{
                       \item If \code{model="logistic"}, any 
                         strictly positive value is set to 1 and 
                         all other values are set to zero
                       \item If \code{model="gamma"}, only strictly 
                         positive values are retained
                        }
                This is motivated by applications involving daily 
                rainfall modelling, in which typically one wants to 
                use logistic regression to model rainfall occurrence
                and a gamma GLM for wet-day rainfall intensities.
                }
          \item{2}{Checking performed and execution halted if any 
                   invalid response values are found}
        }
  }
  \item{max.iter}{Maximum number of iterations for the iterative weighted least squares fitting algorithm. Negative values, along with the default of \code{Inf}, place no limit on the number of iterations allowed.  
  }
  \item{data.file}{
   Name of file containing daily weather / climate data (defaults to 
   \code{gaugvals.dat}, for compatibility with previous versions of 
   \code{GLIMCLIM}). This must be an ASCII file, and must be 
   structured as follows:
   \itemize{
       \item Each line represents one day's data at one location.
       \item The file is sorted by date, and within that by site code.
       \item Each record is in fixed format, with fields as follows: 
        \describe{
         \item{Positions 1-4}{contain the year (e.g. 1978).}
         \item{Positions 5-6}{contain the month (1-12).}
         \item{Positions 7-8}{contain the day.}
         \item{Positions 9-12}{contain the four-character site identifier, as
          defined in \code{siteinfo}.}
         \item{The remaining positions}{contain the data values, in 
         fields of width 6 and each with two decimal places.}
        }
   The FORTRAN format for reading each record is I4,I2,I2,A4,NF6.2 where
   N is the number of variables in the data file. Thus, in a data file 
   containing information on three variables, the records will have the
   form \code{YYYYMMDD$$$$111.11222.22333.33} where \code{YYYY} is the year,
   \code{MM} the month, \code{DD} the day, \code{$$$$} the site identifier
   and \code{111.11} to \code{333.33} the data values. 
   
   The \code{\link{write.GLCdata}} command can be used to create such a file 
   from within \R, although this could be slow for very large datasets 
   running into millions of records (because \R itself is slow for such 
   large datasets). For very large datasets therefore, it may be preferable 
   to use another software package / environment to create data files outside
   \R.    
  }
  For sites which record non-negative values as `trace' below a 
  certain threshold, trace values should be entered as any strictly
  positive number less than the trace threshold (the trace threshold
  itself is defined to the system within \code{model.def}).
  
  If a site provides no data for any variables on a particular
  day, the corresponding record should be omitted from the file. If a 
  site provides data for some variables but not others on a particular 
  day, the missing values should be coded as \code{missval} (see below).

  The underlying FORTRAN code reserves values less than \code{-1e100} to
  keep track of missing data (since FORTRAN does not have a NA value). 
  Problems may be expected, therefore, if such values are provided in 
  the data file. This is unlikely to pose a serious constraint, however! 
  }
 \item{missval}{A numeric value used to indicate missing data in \code{data.file}. Should be set to a value that does not occur as a genuine value for any of the variables in the file. Default is -99.99.
 }
 \item{which.response}{Indicates which of the variables in \code{data.file}
  should be taken as the response. Default is 1.
 } 
 \item{nprev.required}{Used to select a subset of the cases in \code{data.file} for model fitting. If the value is a scalar, models will be fitted only to cases for which the current and \code{nprev.required} previous days' values of \emph{all} variables in \code{data.file} are present. If the value is a vector with length equal to the number of variables in \code{data.file}, then models will be fitted to cases for which the required number of days is present for each variable. If the model contains covariates representing lagged values at lags larger than those specified in \code{nprev.required}, of course, then these additional previous days' values must also be present for a case to be used in the model fitting. 
 
The rationale is that when comparing models, for example using likelihood ratio tests via the \code{\link[=anova.GLC.modeldef]{anova}} method, the models to be compared must be fitted to the same data set. This can cause problems when comparing models with different numbers of lagged terms: setting this argument is a way of ensuring that the models \emph{are} fitted to the same data. This is such an important point, and so easily overlooked, that the software does not provide a default value: the user must decide on this before fitting any model.

\strong{Note:} to completely ignore a particular variable when deciding whether to include a case in the fitting, set the corresponding element of \code{nprev.required} to a negative number. To fit to all possible observations, set \code{nprev.required} to a negative scalar.
 }
 \item{nvar.check}{Controls whether to check that the number of variables in 
  \code{data.file} matches the number of variable names defined in 
  \code{model.def}. The default is \code{TRUE}, in which case a warning is issued in the event of a mismatch.
 } 
 \item{allow.incomplete.averages}{When fitting models, cases can only be
  included if all required covariates are present. In multivariate settings,
  this can lead to large amounts of data being discarded in situations where
  a covariate is a lagged or contemporaneous value of a variable other 
  than the response, and where different stations record different subsets 
  of the variables (for example: if rainfall and temperature are recorded
  at different stations then temperature cannot be included as 
  a covariate in a model for rainfall). However, if the covariate 
  is a (possibly weighted) average over all sites then it may still be 
  reasonably well-defined for the purposes of model fitting, even if it
  was not recorded at the same site as the response. Setting
  \code{allow.incomplete.averages=TRUE} allows the inclusion of such cases
  when fitting models. The results may be slightly unreliable if data are
  sparse however (for example, if precipitation at a site depends only
  on temperatures in a small neighbourhood of that site, and if there
  are no neighboring stations with temperature data). Therefore this 
  argument defaults to \code{FALSE}. 
  
  \strong{Note:} this argument has no effect on lagged values of the 
  response itself. 
  There are three reasons for this: (i) there is less potential 
  for discarding such large quantities of data here, because by definition 
  the lagged values are measured at the same set of stations as the response
  (ii) the results may be slightly unsatisfactory as indicated above 
  (iii) dependence upon lagged values of the response is often strong and 
  may also be localised in space, so that any problems due to the slightly 
  unsatisfactory nature of the solution are likely to be exacerbated.   
  } 
 \item{siteinfo}{
        An object containing site information, for example the result of a call to \code{\link{read.siteinfo}} or \code{\link{make.siteinfo}}.
 }
  \item{model.type}{
   One of \code{"logistic"}, \code{"normal"}, \code{"gamma"} or
   \code{normal-heteroscedastic}, respectively to fit a logistic 
   regression model, linear regression model and gamma GLM with a 
   log link and a linear regression model with logarithm of error
   variance also represented as a linear combination of covariates. 
   Defaults to the \code{model.type} component of \code{model.def}. 
   If equal to \code{normal-heteroscedastic}, \code{model.def} is taken
   to define the mean part of the model: the dispersion part can be
   defined \emph{either} by supplying an appropriate value for 
   \code{dispersion.def} \emph{or} via the \code{dispersion} component of 
   \code{model.def} (if this is an object of class \code{GLC.modeldef} 
   then it will be taken to define the dispersion part of the model).    
  }
 \item{external.files}{
        3-element character vector giving names of files containing
        "external" covariates at annual, monthly and daily resolution
        (in that order). Each element is used only if the model
        requires it. The format of these files is tightly specified; use
        \code{\link{write.GLCexternal}} to create them.
 }
 \item{diagnostics}{
        Controls the amount of diagnostic information (residual analysis 
        etc. aimed at checking the model assumptions) produced for the fit.
        Options:
        \describe{
          \item{0}{No diagnostics produced}
          \item{1}{Produce some standard diagnostics and store these
                   in the result of the call; the \code{plot} and 
                   \code{summary} methods can be used to access these (
                   see documentation for the \code{\link{GLC.modeldef}}
                   class for more details).
                }
          \item{2}{Produce and store standard diagnostics and,
                   additionally, write detailed information to 
                   file \code{resid.file} (see below). 
        }
        The default value (and value assumed if any value other than 0 
        or 1 is supplied) is 2. 
        }
 }
 \item{verbosity}{Controls the amount of information written to screen as
        the fitting progresses. To work silently, set this to zero; to 
        provide brief progress updates, set it to 1; and to provide 
        detailed information set it to 2 (the default).
 }
 \item{resid.file}{
        Optional (produced only if \code{diagnostics=2} output file, 
        containing residual information for 
        each case in the fitting database so that the user can carry out
        further residual analyses if desired. The file has 7 columns and
        a header row. The columns are \code{SITE}, \code{YEAR}, 
        \code{MONTH}, \code{DAY}, \code{OBSERVED}, \code{PREDICTED} and
        \code{ETA}. Most of these are self-explanatory. If the 
        data file contained trace values, these are replaced in the
        \code{OBSERVED} by their approximate conditional expectation
        under the fitted model. The \code{PREDICTED} column contains the 
        fitted values, and the \code{ETA} column contains the linear
        predictors. Records from these files can be read using the
        FORTRAN format 
        \code{A4,1X,I4,T14,I2,T17,I3,T21,F8.4,T31,F8.4,T41,F8.4}.
 }
 \item{cor.file}{
        This output file is produced when the user chooses to 
        incorporate inter-site dependence via correlations (see 
        manual for details). The files can be used directly as 
        input to the \code{\link{GLCsim}} simulation routine 
        if required. They contain a single header row, and then a 
        record for each pair of sites. Each record has six entries:
        4-character site codes (as defined in \code{siteinfo}) for each
        site, then the observed correlation, the number of pairs of
        residuals used in the calculation and finally the site separation 
        in terms of the first two attributes from \code{siteinfo} 
        (the assumption is that these attributes represent 
        geographical coordinates; if they are not defined, the site
        separations are simply output as zero). The FORTRAN format 
        for reading each record is \code{A4,4X,A4,4X,F7.4,1X,I5}.
        Correlations that could not be computed due to a lack of 
        available observations are coded as \code{-9.999}.
 }
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{The function returns an object of class \code{\link{GLC.modeldef}}
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Richard Chandler (\email{r.chandler@ucl.ac.uk})
}
\note{When fitting models involving gamma distributions (currently those with \code{model.type} equal to \code{gamma} or \code{normal-heteroscedastic}), the routine sometimes issues error and warning messages such as 
\preformatted{
****ERROR**** In GAMMLN, you've asked me to calculate ln[GammFn(x)] for x <= 0.
***WARNING*** Can't compute digamma(####) to required accuracy
}
where \code{####} is some small number. These messages are associated with poor starting values for the model coefficients. The routine can usually recover from these, in which case the messages will disappear after the first few iterations and should not be a cause for concern. If the fitting runs for a very large number of iterations and these messages persist, it may be worth trying to change the starting values (in this situation, sometimes it can be helpful simply to start with all of the regression coefficients set to zero). If different starting values still do not help, this may be an indication that the model structure is not appropriate for your dataset. 
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
 \code{\link{GLCdemo}}, \code{\link{write.modeldef}}, 
 \code{\link{GLC.modeldef}}, \code{\link{write.GLCdata}} ,
 \code{\link{write.GLCexternal}}
}
\examples{
require(Rglimclim)
##
##       Load objects for artificial example - help(GLCdemo) for details
##      
data(GLCdemo)
##
##       Collate site information
##
Ashdown.sites <- make.siteinfo(Ashdown.sites,site.names=1,region.col=2)
##
##      Generate data file in current working directory, and fit the
##      simple logistic regression model defined in SimpleModel
##
write.GLCdata(Ashdown.data,file="Ashdown.dat")
print(SimpleModel,global.warn=FALSE)
Model1 <- GLCfit(model.type="logistic",model.def=SimpleModel,
                  siteinfo=Ashdown.sites,diagnostics=1,
                  data.file="Ashdown.dat",nprev.required=0)
print(Model1,global.warn=FALSE)
summary(Model1)
##
##      Summary without residual tables
##
summary(Model1,tables=NULL)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
