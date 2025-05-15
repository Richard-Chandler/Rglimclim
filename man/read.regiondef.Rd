\name{read.regiondef}
\alias{read.regiondef}
\alias{define.regions}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Defining (sub)region names}
\description{\code{read.regiondef} is intended as a convenient way to read region definitions from old \code{GLIMCLIM} (rather than \code{Rglimclim}) definition files.

\code{define.regions} achieves the same thing, but taking its input directly from \R objects.
}
\usage{
read.regiondef(region.file, nhead = 19)
define.regions(names, codes = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{region.file}{
A file containing region codes and names, in the format specified in the Details section below. The file should contain two sections: a header containing \code{nhead} rows (which is ignored); and a section defining the individual regions. 
}  
  \item{nhead}{Number of header rows in \code{region.file}. The default value of 19 is the number of header rows in \code{regions.def} files for previous versions of \code{GLIMCLIM}. 
}
  \item{names}{A character vector of region names, including an entry giving the name of the entire study area}
  \item{codes}{A numeric vector of region codes corresponding to the entries in \code{names}: region zero corresponds to the entire study area, and other (sub)regions should be numbered from 1 to \code{nreg}, where \code{nreg} is the number of subregions being defined. The default assumes that the first element of \code{names} gives the name for the entire study area, with subsequent elements giving names for regions 1 through \code{nreg}.
}
}
\details{The definition of regions is required solely for the simulation routines in \code{Rglimclim}, which offer the option to output summary information for separate regions. In particular, region information is ignored by the model fitting routines. 

For the \code{read.regiondef} function, the structure of the \code{region.file} file is tightly specified (a \code{Fortran 77} relic). After the header, each line of the file defines one region, and looks something like

\code{NUMBER TEXT}

where \code{NUMBER} is a 2-digit region code and \code{TEXT} is a name for the region. \code{NUMBER} occupies positions 1-5 of the line, and \code{TEXT} occupies the remainder. The first region - number 0 - always represents the entire study area and \emph{must} be present. Subsequent (sub)regions must be numbered consecutively from 1 upwards.
}
\value{A data frame with two columns: \code{Region} (the region number, starting with zero) and \code{Name} (the region name). 
}
\author{
Richard Chandler (\email{richard.chandler@ucl.ac.uk})
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%%\seealso{\code{\link{}}
%%}
\examples{
require(Rglimclim)
##
##       Names defined via a character vector
##
region.names <- c("Ashdown Forest",
                  "Pooh and Piglet's side of the forest",
                  "Christopher Robin's side of the forest")
regions <- define.regions(region.names)
print(regions)
##
##       Alternative (if there is a valid regions.def file 
##       in the current working directory - uncomment to run
##       this version):
##
# regions <- read.regiondef("regions.def")
}
