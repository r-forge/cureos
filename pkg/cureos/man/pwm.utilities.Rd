\name{PWM Utilities}
\alias{read.matrices}
\alias{regularize.count.matrix}
\alias{regularize.count.matrices}
\alias{count2profile}
\alias{profile2score}
\alias{revcomp.pwms}
\alias{pwm2consensusSequence}
\alias{pwm2pattern}
\alias{pattern2countmatrix}
\alias{sequence2pattern}
\alias{pattern2sequence}
\alias{random.profile.sequence}
\alias{submatrices}
\alias{core.matrices}
\alias{unwhich}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Utilities to read and create PWMs. }
\description{
  No documentation yet.
}
\usage{
read.matrices(file)
regularize.count.matrix(matrix, method = "Hdiff", background.dist = NULL, weight = 1.5)
regularize.count.matrices(matrices, method = "Hdiff", background.dist = NULL, weight = 1.5)
count2profile(matrices, add.pseudocounts = 0)
profile2score(matrices, background.dist = NULL)
revcomp.pwms(pwms)
pwm2consensusSequence(matrix)
pwm2pattern(pwm)
pattern2countmatrix(pattern, letters="dna", nobs=1)
sequence2pattern(sequence, letters="dna")
pattern2sequence(pattern, letters="dna")
random.profile.sequence(profile.matrix, nseq=1)
submatrices(matrices, starts, ends=NULL, lengths=NULL)
core.matrices(matrices, core.length=5)
unwhich(indices, length)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file}{ File name }
  \item{matrix,profile.matrix,pwm}{ A numeric matrix, where rows
    represent letters in sequence and colums positions. }
  \item{matrices,pwms}{ List of matrices. }
  \item{method}{ Used in regularize either "Hdiff" or "Hfactor". Entropy of
    matrix will change by diff or factor and value in weight. }
  \item{background.dist}{ Distribution of letters in bacground
    model. Numeric vector. }
  \item{weight}{ Used in regularize. Specifies how much should Entropy of
    resulting matrix will change. }
  \item{add.pseudocounts}{ For count2profile. Adds a pseudo.cound before
    computing profile. Alternatively use regularized.matrices. }
  \item{nobs}{ Used in pattern2countmatrix. Total number of observations
    generated in countmatrix. }
  \item{sequence}{ DNA or peptide sequence using IUPAC symbols. }
  \item{pattern}{ Pattern (Regular expression) describing sequence motive. }
  \item{letters}{ Can be either "dna" or "peptide" or a string
    specifying the desired letters. }
  \item{nseq}{ Number of sequences to create. }
  \item{starts,ends,lengths,core.length}{ Start and end positions or lengths for
    submatrices. }
}
\details{
  No documentation yet.
}
\value{
  \item{comp1 }{Description of 'comp1'}
  \item{comp2 }{Description of 'comp2'}
  ...
}
\references{ Rahmann et all, On the Power of Profiles for Transcription
  Factor Binding Site Detection, 2003 }
\author{ Tim Beissbarth }
\seealso{ \code{\link{pwm.scan}} }
\examples{
library(cureos)
# Read a file of matrices
transfac.file <- system.file("data/transfac40.dat", package="cureos")
jaspar.file <- system.file("data/jaspar-08-Jan-2004.dat", package="cureos")
count.matrices <- read.matrices(transfac.file)
# Convert between different types of matrices
reg.count.matrices <- regularize.count.matrices(count.matrices)
profile.matrices<-count2profile(reg.count.matrices)
pwms<-profile2score(profile.matrices)
# Get sequence pattern approximating PWM
pwm2consensusSequence(pwms[[1]])
pwm2consensusSequence(revcomp.pwms(pwms)[[1]])
pwm2pattern(pwms[[1]])
random.profile.sequence(profile.matrices[[1]], nseq=10)
# Generate a PWM from a sequence pattern
pattern<-sequence2pattern("TGWGGY")
pattern2sequence(pattern)
m<-list(from.pattern=pattern2countmatrix(pattern, nobs=10))
profile2score(count2profile(regularize.count.matrices(m)))

}
\keyword{misc}% at least one, from doc/KEYWORDS
