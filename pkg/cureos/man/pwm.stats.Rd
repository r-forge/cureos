\name{PWM statistics}
\alias{pwm.scoredist}
\alias{conserved.pwm.scoredist}
\alias{pwm.exact.scoredist}
\alias{pwm.quality}
\alias{scores2pvalues}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Compute score distributions for Position Weight Matrices (PWMs) }
\description{
  No documentation yet.
}
\usage{
pwm.scoredist(pwm, sequence = NULL, npermutations = 1000)
conserved.pwm.scoredist(pwm, alignment, nsample = 1000)
pwm.exact.scoredist(pwm, background.dist = NULL, eps = 0.001)
pwm.quality(pwm, profile.matrix, background.dist = NULL, eps = 0.001,
            window = 500, instances = 1, plot=TRUE)
scores2pvalues(scores, scoredist, method = "scoredist")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pwm}{ A numeric matrix, where rows represent letters in Sequence
    and colums positions. }
  \item{profile.matrix}{ A numeric matrix, where rows represent letters in Sequence
    and colums positions. }
  \item{background.dist}{ A vector for the background distribution of
    letters or a profile matrix. }
  \item{eps}{ Rounding granularity. }
  \item{sequence}{ A DNA or protein sequence. }
  \item{alignment}{ A DNA or protein alignment. }
  \item{npermutations,nsample}{ Integer number of random tests used in
    permutations to compute maximum scores in permuted sequences. }
  \item{scores}{ Scores resulting from PWM search. }
  \item{scoredist}{ Distribution of scores a PWM reaches in random sequences. }
  \item{method}{ Either equal to "scoredist" or "precomputed.pval". }
  \item{window}{ Length of sequence window to compute significance in. }
  \item{instances}{ Number of instances of a PWMs in window to be
    computed to be significant. }
  \item{plot}{ Produce a plot summarizing pwm.quality. }

}
\details{
  No documentation yet.
}
\value{
  \item{pwm.scoredist}{ An ordered numeric vector with maximum scores for PWM computed in
    permutations of sequence. }
  \item{pwm.exact.scoredist}{List with minimim, maximum and expected
    scores for a PWM, and vectors with the probability and density
    distributions of the of the scores at a computed rounding
    granularity. }
  \item{pwm.quality}{ Data Frame.}
  \item{scores2pvalues}{ Numeric vector with p-values for PWM scores
    based on score distribution. }
}
\references{ Rahmann et all, On the Power of Profiles for Transcription
  Factor Binding Site Detection, 2003 }
\author{ Tim Beissbarth}
\seealso{\code{\link{pwm.scan}} \code{\link{conserved.pwm.scan}} \code{\link{read.matrices}} }
\examples{
library(cureos)

# Read a couple of PWMs
transfac.file <- system.file("data/transfac40.dat", package="cureos")
count.matrices <- read.matrices(transfac.file)
reg.count.matrices <- regularize.count.matrices(count.matrices)
profile.matrices<-count2profile(reg.count.matrices)
pwms<-profile2score(profile.matrices)

# Get matrix for AML1
aml<-grep("AML",names(pwms))
amlprofile <- profile.matrices[aml]
amlcount <- reg.count.matrices[aml]
amlpwm <- pwms[aml]

# Some statitics for matrix
pwmscores <- pwm.scoredist(amlpwm[[1]])
hist(pwmscores,"FD")
qual<-pwm.quality(amlpwm[[1]],amlprofile[[1]])

# Test in searching
sequences <- random.strings(c(1000, 500, 2000))
testmatches1 <- pwm.scan(sequences,amlpwm[1],cut=0.8, pmethod=list(pwmscores))
testmatches2 <- pwm.scan(sequences,amlpwm[1],cut=0.8, pmethod=list(qual$type1))
testscores <- pwm.scores(sequences[1],amlpwm[[1]])
testpvalues1 <- scores2pvalues(testscores, pwm.scores)
testpvalues2 <- scores2pvalues(testscores, qual$type1, "precomputed.pval")

# Test different function. Compare pvalue distributions.
# Here all should give fairly similar results.
qual<-pwm.quality(amlpwm[[1]],amlprofile[[1]], plot=FALSE)
randtest1 <- pwm.scoredist(amlpwm[[1]],random.strings(500),1000)
randtest2 <- conserved.pwm.scoredist(amlpwm[[1]], random.strings(500), 1000)
plot(as.numeric(names(qual$type1)), qual$type1, type="l", col="green")
lines(randtest1, scores2pvalues(randtest1,randtest1), col="red")
lines(randtest2[,1], scores2pvalues(randtest2[,1],randtest2[,1]),col="blue")
lines(randtest2[,2], scores2pvalues(randtest2[,2],randtest2[,2]),col="brown")



}
\keyword{misc}% at least one, from doc/KEYWORDS
