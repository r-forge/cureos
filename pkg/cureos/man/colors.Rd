\name{colors}
\alias{get.colors}
\alias{colorbar}
\alias{show.color.matrix}
\alias{col.chromosome.bands}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Convert numeric values into colors. }
\description{
  Several functions to help in displaying numeric or categorical data as
  colors.
}
\usage{
get.colors(x, col="rainbow", ncol=NULL, decreasing=FALSE, lim=NULL)
colorbar(col, lim, ..., ncol=NULL, decreasing=FALSE, border=FALSE, xlab="")
show.color.matrix(colors)
col.chromosome.bands(cytobands)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ Vector of values to be represented by colors. }
  \item{col}{ Either vector of colorscodes, function or one of:
  \code{heat},\code{terrain},\code{blue},\code{red},\code{green},
  \code{gray},\code{grey},\code{bluered},\code{redgreen},
  \code{bluegreen},\code{greenred},\code{rainbow}.  
  \item{ncol}{ Number of colors to generate }
  \item{decreasing}{ Values should be mapped to colors in decreasing order. }
  \item{lim}{ Limits for color scale }
  \item{border}{ Display borders. }
  \item{xlab}{ Label of X axis. }
  \item{cytobands}{ Table with cytoband positions and annotations. }
  \item{\dots}{ Further arguments passed through. }

}
\value{
  Vector of colors in as hex-coded character vectors.
}
\author{ Tim Beissbarth }
\seealso{ code{\link{read.ucsc.track}}, code{\link{hcl}}}
\examples{
get.colors(-30:100, "greenred", 100, lim=c(-100,100))
colorbar("greenred", lim=c(-100,100))
colorbar(heat.colors, ncol=100, lim=c(0,100))
#zz<-gzcon(url("http://hgdownload.cse.ucsc.edu/goldenPath/hg18/database/cytoBand.txt.gz"))
#f=tempfile()
#write(readLines(zz), f)
#close(zz)
#cytobands=read.ucsc.track(f, "http://hgdownload.cse.ucsc.edu/goldenPath/hg18/database/cytoBand.sql") 
#cytobands=col.chromosome.bands(cytobands)
}
\keyword{ color }

