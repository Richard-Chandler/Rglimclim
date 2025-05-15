\name{write.GLCexternal}
\alias{write.GLCexternal}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Write files containing 'external' covariate data
}
\description{To convert a data frame of external covariate data into a file suitable for input into the \code{Rglimclim} model fitting or simulation routines. 
}
\usage{
write.GLCexternal(x, var.names, file, check.file = TRUE)
}
\arguments{
  \item{x}{A data frame containing predictor information at a daily, monthly or 
  yearly time scale. If the time scale is daily then \code{x} must have columns named
  \code{Year}, \code{Month} and \code{Day} (case-sensitive) in addition to the 
  predictor data. If monthly, it must have columns \code{Year} and \code{Month}; 
  and if yearly it must have a column \code{Year}. All other columns will be taken 
  to contain data values.
}
  \item{var.names}{An optional vector of names for the non-date variables in \code{x}. 
  The default is to use the corresponding column names from \code{x} itself.
}
  \item{file}{Name of the output file.
}
  \item{check.file}{If \code{TRUE} (the default), a check will be made to see if \code{file} already exists and, if so, the user will be asked whether they wish to overwrite the existing file. If \code{FALSE}, any existing file will be overwritten without warning.
}
}
\details{If the software is to be used for downscaling then files of external predictor (covariate) data are required (see the documentation for \code{\link{GLCfit}}, in particular the \code{external.files} argument to that function). The format of these files is tightly specified; the present routine provides an easy way to generate them.

External predictor data for the fitting and simulation routines are stored in different files depending on whether they vary on a daily, monthly or annual timescale. The appropriate timescale for the data in \code{x} is inferred from the variable names (so if \code{x} contains a column called \code{Year} and column called \code{Month} but no variable called \code{Day}, the inference is that the data are monthly). 
}
\value{The routine returns a \code{NULL} value; its sole effect is to write the required data to \code{file}.
}
\note{This routine has not been tested extensively. Note also that it is the user's responsibility to check that the dates in \code{x} are contiguous and in the right order - the routine does not check this (the model fitting and simulation routines will pick up any problems later, however).
}
\author{
Richard Chandler (\email{r.chandler@ucl.ac.uk})
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{GLCfit}}.
}
%%\examples{
%%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
