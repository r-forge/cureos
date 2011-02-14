\name{PWM Searching}
\alias{pwm.scan}
\alias{pwm.scores}
\alias{pattern.scan}
\title{ Search in DNA or protein sequences using scoring matrices (PWMs). }
\description{
  Some functions used to search with a Position Weight Matrix (PWM) within a
  DNA or Protein sequence and to compute significance of hits.
}
\usage{
pwm.scan(sequences, pwms, rev.strand = TRUE, cutoff = 0.1, pmethod = 1000)
pwm.scores(sequence, pwm)
pattern.scan(sequences, patterns, rev.strand = TRUE, cutoff = 0.1, pmethod = 1000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sequences}{ A character vector with DNA or protein sequences. }
  \item{pwms}{ A list of pwm matrices. }
  \item{patterns}{ Sequence motive given by strings or regular expressions. }
  \item{pwm}{ A numeric matrix, where rows represent letters in Sequence
    and colums positions.}
  \item{rev.strand}{ TRUE or FALSE whether to search with reverse
    complement PWMs (only useful for DNA). }
  \item{cutoff}{ Define cutoff for scores or pvalues. }
  \item{pmethod}{ Most commonly integer number of random test used to
    compute pvalues. Can also be a list of numeric vectors with pvalue or
    score distribributions for each pwm. If NULL no pvalues are
    computed. In pattern scan only permutation tests are feasable.}
}
\details{
  No documentation yet.
}
\value{
  \item{pwm.scan}{Data frame with match info.}
  \item{pwm.scores}{ A numeric vector with length
    length(sequence)-ncol(pwm)+1 returning all the scores of a PWM on a sequence.}
}
\references{ http://cureos.wehi.edu.au }
\author{ Tim Beissbarth }
\seealso{\code{\link{read.matrices}}, \code{\link{pwm.scoredist}}}
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

# Generate a couple of random Sequences
sequences <- random.strings(c(1000, 500, 2000))

# Search by foot
testscores <- pwm.scores(sequences[1],amlpwm[[1]])
pwmscores <- pwm.scoredist(amlpwm[[1]])
testpvalues <- scores2pvalues(testscores, pwmscores)

# Easy search interface
testmatches <- pwm.scan(sequences,amlpwm[1],cut=0.8, pmethod=1000)
s <- substr(sequences[testmatches$Sequence.Name],
testmatches$Match.StartPos, testmatches$Match.EndPos)
testmatches$seqs <- s

}
\keyword{misc}% at least one, from doc/KEYWORDS
