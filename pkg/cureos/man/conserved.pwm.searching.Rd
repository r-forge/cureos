\name{Conserved PWM Searching}
\alias{conserved.pwm.scores}
\alias{conserved.pwm.scan}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Search in DNA or protein alignment for conserved high scoring PWMs. }
\description{
  Some functions used to search with a Position Weight Matrix (PWM)
  within a  DNA or Protein sequence alignment and to compute significance of hits.
}
\usage{
pwm.scan(alignment, pwms, rev.strand = TRUE, cutoff = 0.1, pmethod = 1000)
conserved.pwm.scores(alignment, pwm)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{alignment}{ A character vector with DNA or protein sequences in
    alignment. All sequences have same length. Gap character "-".}
  \item{pwms}{ A list of pwm matrices. }
  \item{pwm}{ A numeric matrix, where rows represent letters in Sequence
    and colums positions.}
  \item{rev.strand}{ TRUE or FALSE whether to search with reverse
    complement PWMs (only useful for DNA). }
  \item{cutoff}{ Define cutoff for scores or pvalues. }
  \item{pmethod}{ Integer for number of random test used to
    compute pvalues.}
}
\details{
  No documentation yet.
}
\value{
  \item{conserved.pwm.scan}{Data frame with match info.}
  \item{conserved.pwm.scores}{ A numeric matrix with nchar(alignment[1])
    rows and length(alignment)+1 columns. The columns contain the
    individual scores of the searched pwm at each position. Column
    length(alignment)+1 contains the conserved scores of the pwm,
    i.e. the sum of the individual scores (min 0) weighted by the percentege of
    conserved positions in the alignment.}
}
\references{ http://cureos.wehi.edu.au }
\author{ Tim Beissbarth }
\seealso{ \code{\link{read.matrices}}, \code{\link{read.fasta}}, \code{\link{conserved.pwm.scoredist}}, \code{\link{pwm.scan}} }
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

# Read example alignment
alignment <- read.fasta(system.file("data/testali.fa", package="cureos"))

# Search by foot
testscores <- conserved.pwm.scores(alignment,amlpwm[[1]])
pwmscores <- conserved.pwm.scoredist(amlpwm[[1]],alignment)
testpvalues <- scores2pvalues(testscores[,1], pwmscores[,1])
conspvalues <- scores2pvalues(testscores[,length(alignment)+1],
                              pwmscores[,length(alignment)+1])

# Compare score distribution of best scores of random alignment with
# scores in test alignment
plot(density(pwmscores[,4]),  xlab="Scores", 
     xlim=range(c(testscores[testscores[,4]>0,4], pwmscores[,4])),
     main="Best scores in random alignments vs. real scores.")
lines(density(testscores[testscores[,4]>0,4]), col="red")

# Easy search interface
testmatches <- conserved.pwm.scan(alignment,amlpwm[1],cut=0.8,pmethod=1000)

}
\keyword{misc}% at least one, from doc/KEYWORDS
