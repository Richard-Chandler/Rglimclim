\name{Rglimclim example datasets}
\alias{GLCdemo}
\alias{Ashdown.sites}
\alias{Ashdown.regions}
\alias{Ashdown.data}
\alias{Ashdown.IntensityModel.Initial}
\alias{ConstantModel}
\alias{DemoModel}
\alias{SimpleModel}
\docType{data}
\title{Rglimclim example datasets}
\description{These data sets relate to an artificial example involving daily 
rainfall from six sites in Ashdown Forest in Sussex, England, from 1970 to 1989. 
A map of the area can be found in Milne (1958). The data were actually generated 
by simulating a GLM fitted elsewhere, with appropriate modifications.
These simulated data exhibit many typical features of rainfall sequences in 
northwestern Europe. To make things more realistic, rainfall amounts less than 
0.1mm have been set to `trace' values and appear in the files as values of 
0.05mm. Moreover, approximately 20\% of the values are missing (at random).}
\usage{data(GLCdemo)}
\format{See Details section below.}
\details{This data library contains several objects used in the examples throughout the \code{Rglimclim} package. The objects are as follows:
\describe{
 \item{\code{Ashdown.sites}}{A data frame with six rows and three columns, containing details of six locations (one per row) in Ashdown Forest. The columns are as follows:
 \tabular{ll}{
  \code{Name} \tab The site name \cr
  \code{Region} \tab The region (1 or 2) in which the site is located - see "Ashdown.regions" below \cr
  \code{Eastings} \tab Eastings co-ordinate (as measured, in inches, from the inside jacket cover of Milne 1958) \cr
  \code{Northings} \tab Northings co-ordinate, measured similarly
 }
 The \code{rownames} of the data frame contain 4-character site identifiers. \bold{NB} the first two characters of each identifier are spaces! 
 }
 \item{\code{Ashdown.regions}}{A data frame with three rows and two columns, defining the 
 names of the region and subregions of Ashdown Forest. The columns are:
 \tabular{ll}{
  \code{Region} \tab 0, 1 or 2 (o corresponding to the whole area) \cr
  \code{Name} \tab The region name, to be used in labelling simulation output \cr
 }
 }
 \item{\code{Ashdown.data}}{A data frame containing artificial daily rainfall 
 time series from six sites in Ashdown Forest, from 1970 to 1989. The columns 
 are \code{Year}, \code{Month}, \code{Day}, \code{Site} (a four-character site 
 identifier) and \code{Rain}.
 }
 \item{\code{Ashdown.IntensityModel.Initial}}{An object of class \code{GLC.modeldef}
 defining the structure for the rainfall intensity model used in Section 5.3.7 of the 
 \href{../doc/manual.pdf}{PDF package manual}.
 }
 \item{\code{ConstantModel}}{An object of class \code{GLC.modeldef}
  defining the structure for a model consisting of a constant term and nothing else.
  This can be used to generate an initial model definition file at the start of
  an analysis, using \code{\link{write.modeldef}}.
  }
 \item{\code{SimpleModel}}{An object of class \code{GLC.modeldef}
  defining the structure for a simple logistic regression model consisting 
  of a constant term, a Fourier representation of seasonality and indicators
  for non-zero amounts on three previous days. Again, this can perhaps be used to 
  generate an initial model definition file.
  }
 \item{\code{DemoModel}}{An object of class \code{GLC.modeldef}
  defining the structure for the final logistic regression model considered in
  Chapter 5 of the \href{../doc/manual.pdf}{PDF package manual}. This object
  also contains residual information for the model, as fitted to the data in
  \code{Ashdown.data}; it is used in the examples for \code{\link{GLCfit}}.
  }
}
}
\seealso{the \href{../doc/manual.pdf}{PDF package manual}, in particular the example in Chapter 5.
}
%%\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
%%}
\references{
Milne, A. A. (1958). \emph{The World of Pooh (the complete Winnie-the-Pooh and The House at Pooh Corner)}. Methuen, London.
}
\keyword{datasets}
