\name{genome.annotation.tracks}
\alias{read.ucsc.track}
\alias{regions.to.track}
\alias{track.region}
\alias{plot.tracks}
\alias{lines.track}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Utilities to read and plot genome annotations from UCSC genome
  browser tracks.}
\description{
  Several functions to help in the display of information and positions
  from the UCSC genome browser.
}
\usage{
read.ucsc.track(track, sql = NULL)
regions.to.track(table, colnames = c(name="name", chrompos="chrompos"))
track.region(track, chrom, start, end)
plot.tracks(tracks, chrom, start, end, track.offset=1, gene.offsets=0.6,
            col="black", width=0.4, chrom.pos=FALSE, chrom.scale=1000, type="b", ...)
lines.track(track, at = 0, offset = 1, width = 0.6, col = "black",
            track.name = NA, names = TRUE, lty = 1, line = 0, las = 2,
            adj = 1, padj = 0.5, type = "b", chrom.scale = 1000,
            border = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{track}{ Data.frame in UCSC genome browser format with columns:
    \code{name}, \code{chrom}, \code{start}, \code{end}
    Possible further columns (comma separated):
    \code{exonStarts}, \code{exonEnds}, \code{strand}, \code{color}
  }
  \item{sql}{ Description of database tabe of UCSC track. }
  \item{table}{ Table in format:
    \code{name}, \code{chrom}, \code{start}, \code{end}
    Multiple exons and respective information are given in
    seperate rows.
  }
  \item{colnames}{ Description of column names in table to translate
    formats. Column names must have:
    \code{name}, \code{chrom}, \code{start}, \code{end}, or \code{chrompos}
  }
  \item{chrom}{ Character vector describing chromosomes,
    e.g. \code{chrX}, \code{chr11} }
  \item{start}{ Integer vector indication chromosomal start positions. }
  \item{end}{ Integer vector indication chromosomal start positions. }
  \item{tracks}{ List of tracks }
  \item{track.offset}{ Gap between tracks. }
  \item{gene.offsets}{ Gap between different genes (indicated by name)
    in a track. }
  \item{col}{ Colors }
  \item{width}{ Display width of boxes in individual track lane. }
  \item{chrom.pos}{ Logical, whether to display xaxis with chromosome positions. }
  \item{chrom.scale}{ Rescale chromosome position }
  \item{type}{ Can be either: r=rectangles, s=segments, b=both }
  \item{at}{ y position of track start. }
  \item{offset}{ Offset between individual lines in track. }
  \item{track.name}{ Name of track to display. }
  \item{names}{ Logical, should gene names (from column name) be displayed }
  \item{lty}{ Line type of segments. }
  \item{line}{ line to display gene names. }
  \item{las}{ label orientation }
  \item{adj}{ label adjustment }
  \item{padj}{ label adjustment }
  \item{border}{ Logical, whether to display border of boxes. }
  \item{\dots}{ Arguments passed to other functions }
}
\details{
  A genome annotation track is stored as a data.frame.
  Annotation tracks can be obtained for example via
  http://hgdownload.cse.ucsc.edu/downloads.html. Lists of
  tracks can be displayed as box and line plots similar to a
  genome browser for different chromosomal regions using the
  \code{plot.tracks} function. Essential columns in the track format
  are: \code{name}, \code{chrom}, \code{start}, \code{end}. Further
  columns are useful, for example to display genes with several
  exons and introns: \code{exonStarts}, \code{exonEnds}, \code{color}.
  These columns are comma separated lists of either integer or colorcode
  vectors of equal length.
}
\value{
  Usually a data.frame or list of data.frames with genome.annotation
  track information.
}
\author{ Tim Beissbarth }
\examples{
#tracknames = c("refGene.txt", "cytoBand.txt",, "cpgIslandExt.txt",
#               "phastConsElements28way.txt")

#tracks=list()
#for(i in tracknames)
#  tracks[[i]]=read.ucsc.track(i)

#tracks$refGene.txt$color = c("darkgreen", "darkred")[factor(tracks$refGene.txt$strand)]
#tracks$cpgIslandExt$color=get.colors(tracks$cpgIslandExt$obsExp,
                                     col="gray", ncol=20, decreasing=TRUE, lim=c(0,2))
#tracks$phastConsElements28way.txt$color=get.colors(tracks$phastConsElements28way.txt$score,
#  col="gray", ncol=20, decreasing=TRUE, lim=c(0,1000))

#plot.tracks(tracks, "chr10", 35303095, 35763633, gene.offsets=c(0.6,0,0,0))

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ misc }

