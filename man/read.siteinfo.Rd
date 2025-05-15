\name{read.siteinfo}
\alias{read.siteinfo}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{To read site information from a file.}
\description{This function is intended as a convenient way to read site information from old \code{GLIMCLIM} (rather than \code{Rglimclim}) site definition files, and to store the result in a format that can be used directly as an argument to other routines such as \code{\link{GLCfit}}.
}
\usage{
read.siteinfo(site.file, nhead = 38, which.coords=NULL, regions=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{site.file}{A file containing site information, in the format specified in the Details section below. The file should contain three sections: a header containing \code{nhead} rows (which is ignored); a section defining labels for the attributes of each sites; and a section defining the individual sites and their attributes. 
}
  \item{nhead}{Number of header rows in \code{site.file}. The default value of 38 is the number of header rows in \code{siteinfo.def} files for previous versions of \code{GLIMCLIM}. 
}
  \item{which.coords}{For many purposes (e.g. fitting spatial correlation structures, calculating distance-weighted averages of lagged values) it is necessary to calculate distances between sites. This argument is used to define which attributes will be used to calculate (Euclidean) distances, if more than two attributes are defined (if one or two attributes are defined, these are automatically used for distance calculations). The default value is \code{1:2}, the rationale being that users will usually define the spatial co-ordinates as the first two site attributes. No more than two attributes may be used to define inter-site distances.
} 
  \item{regions}{An optional data frame giving names of regions or subareas, in the format returned by \code{\link{read.regiondef}} or \code{\link{define.regions}}.
}
}
\details{The structure of the \code{site.file} file is tightly specified (a \code{Fortran 77} relic). The first line after the header should contain a single integer (\code{n.attr}, say) which is the number of site attributes to be defined. The subsequent \code{n.attr} lines contain text describing each attribute: this will be used in labelling model output. 
After the attribute definition section of the file, there is a line which is ignored by this routine (in previous versions of \code{GLIMCLIM} this contained text "****END OF ATTRIBUTE DEFINITION****"); this is followed by the final section of the file which is used to specify the details of individual sites. Each site is defined on two lines: the first contains a text string giving the site name and the second looks something like:

\code{  CODE  REGION  ATTRIB[1]  ...  ATTRIB[2]  ...  ATTRIB[n.attr]}

where \code{CODE} is a 4-character code used to identify the site, \code{REGION} is the number of the region in which the site lies (as defined using \code{\link{read.regiondef}} or \code{\link{define.regions}} - this field is optional, with a value of 0 or blank associating the site with no region), and \code{ATTRIB[]} are the values of the various attributes for the site (in the order specified in the `attribute definition' section). \code{CODE} occupies the first 4 positions of the line; \code{REGION} occupies positions 5-9; and each successive attribute occupies the first 10 free positions of the line starting with position 10.

There is no limit to the number of sites that can be defined, or to the number of attributes that can be defined for each site.}
\value{See help for \code{\link{make.siteinfo}}.}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Richard Chandler (\email{r.chandler@ucl.ac.uk})
}
%\note{
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{make.siteinfo}} to define sites based on information held in a data frame, \code{\link{GLCfit}}
}
\examples{
require(Rglimclim)
##
##       Valid siteinfo.def file required in current working directory
##
# siteinfo <- read.siteinfo("siteinfo.def")
}

