\name{write.modeldef}
\alias{write.modeldef}
\title{Write a model definition file}
\description{To write a definition file for an object of class \code{\link{GLC.modeldef}}, that can be read subsequently using \code{\link{read.modeldef}}.
}
\usage{
write.modeldef(x, file, check.file=TRUE, mean.only=FALSE)
}
\arguments{
 \item{x}{An object of class \code{\link{GLC.modeldef}}, for example resulting from a call to \code{\link{GLCfit}}. }
 \item{file}{Name(s) of output file(s). If \code{x} is a joint mean-variance model then this should be a vector of length two (one file for the mean component and one for the dispersion); otherwise a single character string.
}
 \item{check.file}{If \code{TRUE} (the default), the routine checks that the output file(s) in \code{file} does / do not already exist; if it does / they do, user will be prompted before proceeding.}
 \item{mean.only}{If \code{TRUE}, the routine will only write information relating to the mean component of the model (no dispersion info, spatial structure etc.). Defaults to \code{FALSE}.
}
}
\details{This function is provided for back-compatibility with \code{GLIMCLIM} (rather than \code{Rglimclim}). See the documentation for \code{\link{read.modeldef}} for more details. 
}
\value{The function returns a \code{NULL} value: its effect is solely to write a model specification to file.}

\author{
Richard Chandler (\email{r.chandler@ucl.ac.uk})
}

\seealso{\code{\link{read.modeldef}}, and \code{\link{GLCfit}} for examples.}

