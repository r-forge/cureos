\name{sequence.utilities}
\alias{read.fasta}
\alias{read.subsequences}
\alias{random.strings}
\alias{permute.sequences}
\alias{background.dist.from.seq}
\alias{rev.strings}
\alias{revcomp}
\alias{randomize.alignment}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Utilities to work with DNA or protein sequences. }
\description{
  No documentation yet.
}
\usage{
read.fasta(file)
read.clustalw(file)
read.subsequences(files, starts=1, ends=NULL, lengths=NULL)
random.strings(lengths, background.dist=NULL, letters="dna")
permute.sequences(sequences)
background.dist.from.seq(sequence, letters="dna")
rev.strings(strings)
revcomp(sequences)
randomize.alignment(alignment, letters="dna")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{files,file}{ Character vector with file names. }
  \item{starts}{ Start positions of sequence to be read
                [index starts with 1]}
  \item{ends}{ End positions of sequence to be read. [alternative: length] }
  \item{lengths}{ Length of sequences. }
  \item{background.dist}{ Background distribution of nucleotides,
    aminoacids or letters. }
  \item{letters}{ Letters can be "dna" or "peptide" or character
    string. }
  \item(strings,sequences,sequence){ Character vector. }
  \item{alignment}{ Character vector with aligned sequences to be randomized. }
}
\details{
  No documentation yet.
}
\value{
  \item{sequence}{ Character vector, with sequences or alignment. }
  \item{background.dist}{ Integer vector describing distribution of
    nucleotides of aminoacids. }
}
\references{ ~put references to the literature/web site here ~ }
\author{ Tim Beissbarth }

 ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ \code{\link{pwm.scan}} }
\examples{
library(cureos)
# Tools for DNA or peptide sequences
seqs<-random.strings(c(20,20), letters="dna")
revcomp(seqs)
rev.strings(seqs)
permute.sequences(seqs)
background.dist.from.seq(random.strings(10000,background.dist=c(0.3,0.4,0.2,0.1)))
# Randomize alignment sequences, keeping conservation and gaps
alignment<-c("aa-cgtagc", "aaa-gtagc", "aa-cgcagc")
as.matrix(alignment)
as.matrix(randomize.alignment(alignment))
# Reading sequences from FASTA files.
chrM <- system.file("data/humanMitochondrium.fa", package="cureos")
read.subsequence(rep(chrM,4), starts=c(1,20,90,16000), ends=NULL, lengths=100)
ali <- system.file("data/testali.fa", package="cureos")
alignment <- read.fasta(ali)
substr(alignment, 1, 100)
}
\keyword{ misc }% at least one, from doc/KEYWORDS
