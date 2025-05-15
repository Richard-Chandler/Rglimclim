\name{make.siteinfo}
\alias{make.siteinfo}
\alias{print.siteinfo}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{To convert site information from a data frame}
\description{To extract site information from a data frame and produce a \code{siteinfo} object suitable for passing to \code{\link{GLCfit}}.
}
\usage{
make.siteinfo(site.data, coord.cols = NULL, site.codes = NULL, 
              site.names = NULL, region.col=NULL, attr.cols = NULL, 
              attr.names = NULL, regions = NULL)
\method{print}{siteinfo}(x,...)
}
\arguments{
  \item{site.data}{
A data frame containing information on the attributes of each site to be defined and, optionally, site names and short 4-character identification codes. There should be one row for each site. See \code{\link{GLCdemo}} for an example). 
}
  \item{coord.cols}{
A numeric vector containing at most two elements, identifying the number(s) of the column(s) in \code{site.data} that contain spatial co-ordinates of the sites. If unspecified, the first two columns that do not contain names or site codes are taken as the spatial co-ordinates (if there is just one candidate column, this is taken as a single spatial co-ordinate).
}
  \item{site.codes}{
Either a character vector of short site identifiers, or a single number indicating which column of \code{site.data} contains these identifiers. Each identifier should be exactly 4 characters long (including spaces if necessary). The default is to take \code{rownames(site.data)}: in this case the the row names should themselves be 4 characters in length. The codes are used to identify observations when reading data files.}
  \item{site.names}{
Either a character vector of site names, or a single number indicating which column of \code{site.data} contains the site names. This must be supplied and is used to label software output. 
}
  \item{region.col}{
A single number indicating which (if any) column of \code{site.data} contains the numbers of the regions in which each site is located (codes should correspond to regions defined via \code{\link{read.regiondef}} or \code{\link{define.regions}}). Default is not to associate sites with any region except the entire area. 
}
  \item{attr.cols}{
Numeric vector indicating which columns of \code{site.data} contain site attribute information. Default is all columns that do not contain site codes, site names or region codes.
}
  \item{attr.names}{
Character vector of descriptive text for each of the site attributes being defined. The default is to take the corresponding column names of \code{site.data}.
}
  \item{regions}{An optional data frame giving names of regions or subareas, in the format returned by \code{\link{read.regiondef}} or \code{\link{define.regions}}.
}
  \item{x}{For the \code{print} method, an object of class \code{siteinfo}.}
  \item{...}{Additional arguments to the \code{print} method (not used here).}
}
\details{There is no limit to the number of sites that can be defined, or to the number of attributes that can be defined for each site.}
\value{The function returns a \code{list} result of class \code{siteinfo}, with components
\item{Nsites}{The number of sites defined in file \code{site.data}}
\item{Site.names}{A vector of character strings giving the names of the sites (i.e. the values of \code{site.names})}
\item{Site.codes}{A vector of 4-character identifiers for each site (this is used to identify observations when reading data files)}
\item{Regions}{A numeric vector indicating the region in which each site is located (codes correspond to regions defined via \code{\link{define.regions}}) or \code{\link{read.regiondef}}. Contains zeroes for sites where no region information is given. Regions are defined solely for the purpose of producing regional summary statistics from simulated weather sequences.}
\item{Region.names}{A character vector giving the names of the regions. Again, this is used only for labelling simulation output.}
\item{Nattr}{The number of attributes defined for each site.} 
\item{Attribute.names}{A vector of character strings describing each attribute (this is used for labelling model output). The attributes corresponding to \code{coord.cols} are listed first.}
\item{Attribute.values}{A matrix containing the attribute values for each site: rows correspond to sites (labelled using \code{Site.codes}) and columns to attributes (labelled using truncated \code{Attribute.names}). The attributes corresponding to \code{coord.cols} appear first.}

The \code{print} method shows what attributes have been defined and in what order. 
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Richard Chandler (\email{r.chandler@ucl.ac.uk})
}
%\note{
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{read.siteinfo}} to define sites based on an old \code{GLIMCLIM} site definition file, \code{\link{GLCdemo}},
\code{\link{GLCfit}}
}
\examples{
require(Rglimclim)
#
#       Load example data and view the Ashdown.sites dataframe
#
data(GLCdemo)
print(Ashdown.sites)
#
#       Convert the data frame to the required format, taking the 
#       site names from the first column, the region definitions
#       from the second (and, by default, taking the 
#       remaining columns as site attributes) and defining
#       some extended descriptive text for each of the 
#       attributes.
#       
siteinfo <- make.siteinfo(Ashdown.sites,site.names=1,region.col=2,
                          attr.names=c(
               "Eastings (inches from left of 11\" wide map)",
               "Northings (inches from bottom of 8\" high map)"))
print(siteinfo)                 # Brief summary on screen
print(unclass(siteinfo))        # To see the gory details
}
