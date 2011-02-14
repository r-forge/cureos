\name{chrompos}
\alias{split.chrompos}
\alias{paste.chrompos}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Change formats of chromosome position descriptions. }
\description{
  Convert between chromosome positions as individual values for
  \code{chrom}, \code{start}, \code{end} to character strings in the
  format chr12:012345678-123456789.
  
}
\usage{
split.chrompos(chrompos)
paste.chrompos(chrom, start, end)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{chrompos}{ Character vector with chromosome positions. }
  \item{chrom}{ Character vector with chromosomes. }
  \item{start}{ Integer vector with start positions. }
  \item{end}{ Integer vector with end positions. }
}
\value{
  Returns either character vector with chromosome positions
  or data.frame with colums \code{chrom}, \code{start}, \code{end}.
}
\author{ Tim Beissbarth }
\seealso{\code{\link{read.ucsc.track}}}}
\examples{
chrompos = c("chr16:011282193-011282693", "chr1:151278824-151280218",
             "chr1:000058953-000059871",  "chr1:000357521-000358458",
             "chr1:000357521-000358458")
table = split.chrompos(chrompos)
table
paste.chrompos(table$chrom, table$start, table$end)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ misc }

