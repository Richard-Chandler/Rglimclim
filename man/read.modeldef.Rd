\name{read.modeldef}
\alias{read.modeldef}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{To read model definition from a file.}
\description{This function is intended as a convenient way to read model definitions from old \code{GLIMCLIM} (rather than \code{Rglimclim})  definition files, and to store the result in a format that can be used directly as an argument to other routines such as \code{\link{GLCfit}}.
}
\usage{
read.modeldef(model.file, nhead = 46, model.type, which.part="mean", siteinfo,
var.names, which.response=1, 
external.files = c("yr_preds.dat", "mn_preds.dat","dy_preds.dat"), 
sim = FALSE, oldGlimClim.warning=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model.file}{The name of a model definition file, which should be in the format specified in the Details section below. The file should contain three sections: a header containing \code{nhead} rows (which is ignored); a model title (up to 70 characters) used for labelling output to help with subsequent identification of printouts; and a section defining the model structure. 
}
  \item{nhead}{Number of header rows in \code{model.file}. The default value of 46 is the number of header rows in \code{logistic.def} and \code{gammamdl.def} files for previous versions of \code{GLIMCLIM}. 
}
  \item{model.type}{One of \code{normal}, \code{gamma}, \code{logistic} or \code{normal-heteroscedastic}. 
}
  \item{which.part}{If \code{model.type="normal-heteroscedastic"}, indicates whether we're defining the mean (\code{which.part="mean"}) or dispersion (\code{which.part="dispersion"}) components.}
  \item{siteinfo}{An object of class \code{siteinfo}, resulting from a call to \code{\link{read.siteinfo}} or \code{\link{make.siteinfo}}. This is used for checking that any required site attributes have been defined.}
  \item{var.names}{Vector of variable names (used to construct labels, 
  in multivariate models in particular). If not supplied, this is set to 
  \code{c("Y1","Y2",...,"Ym")} where \code{m} is the maximum variable 
  index found in the definition file.}
  \item{which.response}{Index number of the response variable (used to 
  identify which element of \code{var.names} corresponds to the response).
  Default is 1.
  }
  \item{external.files}{
        3-element character vector giving names of files containing
        "external" covariates at annual, monthly and daily resolution
        (in that order). Each element is used only if such covariates
        are part of the model specification: thus, for example, if the model 
	requires monthly external covariates but not daily or annual ones, then
        only the second element will be used (but it is still necessary
	to supply a 3-element vector, for example by setting the first and
	third elements to \code{NA}). See manual for details.  
 }
  \item{sim}{Logical value, indicating whether we're defining a model for simulation or fitting purposes. For simulating some model types, a dispersion parameter should be specified in the model definition file (this is not necessary for model fitting - dispersion parameters will be estimated automatically as necessary). Conversely, if we're fitting models then we'll need to calculate derivatives of log-likelihoods with respect to any parameters that enter nonlinearly intothe "linear predictor". Defaults to \code{FALSE}.}
  \item{oldGlimClim.warning}{Logical value indicating whether to issue a warning if the model definition file contains a value for \code{CODE3} (see below) in a row relating to lagged daily values. The reason is that the meaning of this code is different in \code{Rglimclim} to that in the original \code{GlimClim} package: users of old definition files will therefore need to make the appropriate changes.  
  }
}
\details{
Model definition files are used to define covariates for
a GLM (along with their associated coefficients), and to 
select spatial dependence structures. For the model fitting programs,
the parameter values in these files are treated as initial estimates 
for the numerical estimation algorithm.

The final section of the file (see documentation for argument \code{model.file} above) defines the model structure proper. In this section, each row corresponds to a single parameter in the model, and contains
five entries in fixed-width format:
                                                                          
\preformatted{COMPONENT       VALUE    CODE1   CODE2   CODE3   TEXT}
                                                                            
The field widths are 5,10,5,5,5,40. An explanation of the entries is as follows:
\describe{
 \item{\code{COMPONENT}}{occupies the first 5 positions in the record,
 and is used to identify the type of quantity being defined. Valid  values are:
\describe{
\item{\code{0}}{if the record relates to the constant term in the
  linear predictor (i.e. the regression part of the model).}
\item{\code{1}}{if the record relates to a site effect in the linear predictor.}
\item{\code{2}}{if the record relates to a `year' effect in the
  linear predictor.}
\item{\code{3}}{if the record relates to a `month' effect in the
  linear predictor.}
\item{\code{4}}{if the record relates to a `day' effect, in the
  linear predictor. This includes previous days' values,
  as well as any other covariate that varies on a daily timescale (except those that vary with both day and site).}
\item{\code{5}}{if the record relates to a 2-way interaction in the
  linear predictor.}
\item{\code{6}}{if the record relates to a 3-way interaction in the
  linear predictor.}
\item{\code{7}}{if the record relates to the nonlinear transformation
  of one of the covariates in the linear predictor.}
\item{\code{8}}{if the record relates to a global quantity such as a
`trace' threshold for raingauges. Such quantities are not strictly
part of the model, but must be defined to the system somehow.}
\item{\code{9}}{if the record relates to a dispersion parameter for the
  model (not all models have a dispersion parameter, and the parameter does not need to be specified when fitting models - see documentation for argument \code{sim} above).}
\item{\code{10}}{if the record relates to the spatial structure of the model.}
}
The rows of a model definition file must be ordered according to the
value of \code{ COMPONENT}; if the rows are out of order, an error message will result.}
\item{\code{VALUE}}{is the value of the parameter / coefficient being defined, occupying positions 6-15 of the record.}
\item{\code{CODE1}, \code{CODE2} and \code{CODE3}}{are used to define 
the precise details of the model to the system. Their interpretation 
depends on the value of \code{COMPONENT}. In general, it is not
necesary to define all three codes. If they are all defined,
\code{CODE1} occupies positions 16-20 of the record, \code{CODE2} occupies positions 21--25 and \code{CODE3} occupies positions 26-30. See Tables 2 to 8 of the \href{../doc/manual.pdf}{PDF package manual} for full details on coding.}
\item{\code{TEXT}}{contains descriptive text for this record, and
appears after position 31. It is intended to help make the definition file
readable by the user, and is not required or used by the program.}
}
}
\value{An object of class \code{\link{GLC.modeldef}}. This is not intended to be manipulated by the user: it is a list containing a large number of components, all of which are passed to the underlying \code{Fortran} code when fitting and simulating models.
}
\author{
Richard Chandler (\email{r.chandler@ucl.ac.uk})
}


%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
require(Rglimclim)
##
##       Define site information (see make.siteinfo for details)
##
data(GLCdemo)
siteinfo <- make.siteinfo(Ashdown.sites,site.names=1,region.col=2,
                          attr.names=c(
               "Eastings (inches from left of 11\" wide map)",
               "Northings (inches from bottom of 8\" high map)"))
##
##       Valid logistic.def file required in current working directory
##
#logistic.modeldef <- read.modeldef("logistic.def",
#                                   model.type="logistic",
#                                   siteinfo=siteinfo)
}

