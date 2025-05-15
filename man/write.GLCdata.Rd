\name{write.GLCdata}
\alias{write.GLCdata}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Write an Rglimclim data file}
\description{To convert a data frame into a data file of potential response variables,
 for input into the \code{Rglimclim} model fitting and simulation routines.
}
\usage{
write.GLCdata(x, date.cols = 1:3, site.col = 4, data.cols, missval = -99.99, 
              file, check.file = TRUE)
}
\arguments{
  \item{x}{A data frame containing a row for each day:site combination for which data are available. The data frame should contain columns for year (4-digit integer), month (2 digits), day (2 digits), site (4 characters) and data values (one column for each variable required). The 4-character site codes should match those defined in \code{\link{make.siteinfo}} or \code{\link{read.siteinfo}} - note in particular that alignment is important in case the codes contain leading or trailing spaces. 
}
  \item{date.cols}{A numeric vector of length 3, giving the numbers of the columns in \code{x} corresponding respectively to  year, month and day of each observation. The default is \code{1:3}.
}
  \item{site.col}{A numeric scalar, giving the number of the column in \code{x} containing 4-character site codes (this will be converted to a character vector if it is not already one). The default is \code{4}.
}
  \item{data.cols}{A numeric vector giving the number(s) of the column(s) in \code{x} containing the data values to be written. The default is all columns except those defined in \code{date.cols} and \code{site.col}. 
}
  \item{missval}{The numeric value that will be used to replace any \code{NA} values in \code{x} (default \code{-99.99}). This is needed for compatibility with the underlying \code{Fortran} code which does not recognise a separate \code{NA} value. 
}
  \item{file}{Name of the output file.
}
  \item{check.file}{If \code{TRUE} (the default), a check will be made to see if \code{file} already exists and, if so, the user will be asked whether they wish to overwrite the existing file. If \code{FALSE}, any existing file will be overwritten without warning.
}
}
\details{The \code{\link{GLCfit}} and \code{\link{GLCsim}} routines require data files in a tightly specified format, as described in the help for the \code{data.file} argument of the \code{\link{GLCfit}} routine). The present routine provides an easy way to generate them from an existing data frame in \R.
}
\value{The routine returns a \code{NULL} value; its sole effect is to write the required data to file.
}
\author{
Richard Chandler (\email{r.chandler@ucl.ac.uk})
}
\note{This routine has not been tested extensively. Moreover, it may be slow and memory-intensive when writing very large datasets containing hundreds of thousands of cases for several variables. Finally: it is the user's responsibility to check that the rows of  \code{x} are in chronological order - the routine does not check this (the model fitting and simulation routines will pick up any problems later, however).
}

\seealso{\code{\link{GLCfit}}, \code{\link{read.GLCdata}}.}
%\examples{
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
