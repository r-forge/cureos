\name{moving.average}
\alias{window.average}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Computes moving averages. }
\description{
  Computes moving averages over a sequence of numbers, either of windows
  of a certain width (like n consecutive numbers) or using an additional
  position vector (like positions on a chromosome), where the windows
  are defined by the distances in the position vector.  
}
\usage{
moving.average(numbers, window.size=10)
window.average(numbers, positions, window.size=1000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{numbers}{ Numeric vector. }
  \item{positions}{ Numeric vector of the same length as numbers. }
  \item{window.size}{ Number of elements to take averages or distances
    in position vector. }
}
\details{
  The function \code{moving.average} computes the moving averages of
  windows of \code{n=window.size} consequtive numbers in a numerical
  vector. The result is a new vector which is \code{n-1} elements shorter than
  the old one.

  The function \code{window.average} uses an additional position vector
  (like e.g. chromosomal positions). The positions vector must be
  ordered, otherwise averaging will start new at each break of the
  ordering (like e.g. several different chromosomes). A numeric vector
  with window averages of the same length as the numbers and positions
  vector is returned. The windows are computed based on the positions
  vector and the numbers, where the window is centered at each of the
  positions and the numbers to the left and the right, which fall into
  the region, are averaged.
  
}
\value{
  \item{numbers}{ Vector with window averages. }
}
\references{ ~put references to the literature/web site here ~ }
\author{ Tim Beissbarth }

 ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ \code{\link{pwm.scan}} }
\examples{
library(cureos)
moving.average(1:30, window.size=10)

window.average(rep(5, 30), c(1:10/2, 1:10, 1:10*2), window.size=10)

}
\keyword{ misc }% at least one, from doc/KEYWORDS
