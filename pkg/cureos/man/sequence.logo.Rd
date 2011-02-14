\name{sequence.logo}
\alias{simple.sequence.logo.plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Generate Sequence Logo from Count matrices. }
\description{
  No documentation yet.
}
\usage{
simple.sequence.logo.plot(matrix, e = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{matrix}{ Count matrix, numeric matrix with one row for each
    letter and columns for positions.}
  \item{e}{ error term. }
}
\details{
  No documentation yet.
}
\value{
  ~Describe the value returned
  If it is a LIST, use
  \item{comp1 }{Description of 'comp1'}
  \item{comp2 }{Description of 'comp2'}
  ...
}
\references{ Schneider et al, NAR, 1990 }
\author{ Tim Beissbarth }
\seealso{ \code{\link{read.matrices}}}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function(matrix, e=0) {
  entropy <- -colSums(matrix *log2(matrix), na.rm=TRUE) # Shannon entropy
  R <- log2(nrow(matrix)) - entropy + e
  hights <- t(t(matrix) * R)
  barplot(hights, col=rainbow(nrow(matrix)))
  }
}
\keyword{misc}% at least one, from doc/KEYWORDS
