\name{read.GLCdata}
\alias{read.GLCdata}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Read an Rglimclim data file}
\description{To read a daily data file of response variables in 
\code{Rglimclim} format (see the help for the \code{data.file} 
argument of the \code{\link{GLCfit}} routine). This can also be used 
to read daily output files from simulations produced by 
\code{\link{GLCsim}}. 
}
\usage{
read.GLCdata(file, missval = -99.99)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file}{Character string giving the name of the data file}
  \item{missval}{The number used to denote missing values in the data file}
}
\details{
Although \code{Rglimclim} data files are in fixed-width format,
the use of \code{\link{read.fwf}} to read them can be extremely 
slow. The present command is much faster.  
}
\value{The function returns a data frame, with named columns
\code{Year}, \code{Month}, \code{Day}, \code{Site} and 
\code{Var1} to \code{VarN}, where \code{N} is the number of
variables present in \code{file}. Values of \code{missval}
will be replaced with \code{NA}.
}
\author{Richard Chandler (\email{r.chandler@ucl.ac.uk})}
\note{
Daily weather datasets can be extremely large. Indiscriminate use of
this routine may cause slow running or memory problems, therefore.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{write.GLCdata}}, \code{\link{GLCfit}}}
%%\examples{
%%##---- Should be DIRECTLY executable !! ----
%%##-- ==>  Define data, use random,
%%##--	or do  help(data=index)  for the standard data sets.
%%
%%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
