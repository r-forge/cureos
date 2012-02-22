.First.lib <- function(lib, pkgname, where) {
  library.dynam(pkgname, pkgname, lib)
}

# Read from Transfact matrix.dat style format in a count, score, or profile matrix
read.matrices <- function(file) {
  input <- readLines(file, -1)
  input <- grep("(^ID)|(^[0-9][0-9]*)", input, value=TRUE)
  idlines <- logical(length(input))
  idlines[grep("^ID", input)] <- TRUE
  if(sum(idlines) == 0) {
    stop("Matrix file in wrong format. No rows starting with ID found.")
  }
  ids <- sub("^ID  *([^ ]*).*", "\\1", input[idlines])
  counts <- as.matrix(as.data.frame(strsplit(input[!idlines], "  *")))
  if(nrow(counts) == 4) {
    counts <- matrix(as.numeric(counts),nrow=4)
    rownames(counts) <- c("A", "C", "G", "T")
  } else if(nrow(counts)==5 || nrow(counts)==6) {
    counts <- matrix(as.numeric(counts[2:5,]),nrow=4)
    rownames(counts) <- c("A", "C", "G", "T")
  } else if(nrow(counts)==20) {
    counts <- matrix(as.numeric(counts),nrow=20)
    rownames(counts) <- c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")
  } else if(nrow(counts)==21 || nrow(counts)==22) {
    counts <- matrix(as.numeric(counts[2:21,]),nrow=20)
    rownames(counts) <- c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")
  } else {
    stop("Matrix file in wrong format. Number of columns has to be 4-6 for DNA or 20-22 for Peptide matrices")
  }
  counts2 <- matrix(0, nrow(counts), length(input))
  counts2[,!idlines] <- counts
  rownames(counts2) <- rownames(counts)
  idlines <- c(which(idlines), length(input)+1)
  matrices <- list()
  for(i in 1:(length(idlines)-1)) {
    matrices[[i]] <- counts2[,(idlines[i]+1):(idlines[i+1]-1)]
  }
  names(matrices) <- ids
  matrices
}

# Read from Transfac matrix.dat style format in a count, score, or profile matrix
write.matrices <- function(x, file) {
  write(paste("VV Exported matrices from CUREOS R package", date(), "\nXX\n//"), file)
  if(!is.list(x)) stop("Must be a list of matrices")
  for(i in 1:length(x)) {
    write(paste("XX\nID  ", names(x[i]), "\nXX"), file, append=TRUE)
    tm<-rbind(c("PO", formatC(rownames(x[[i]]), width=12), ""),
              t(rbind(formatC(1:ncol(x[[i]]), width=2, format="d", flag="0"),
                      formatC(x[[i]], width=12, digits=7, format="f"),
                      unlist(strsplit(pwm2consensusSequence(x[[i]]),"")))))
    write.table(tm, file, append=TRUE, sep=" ", col.names=FALSE, row.names=FALSE, quote=FALSE)
    write(paste("XX\n//"), file, append=TRUE)

  }
}

# get regions from matrices
submatrices <- function(matrices, starts, ends=NULL, lengths=NULL) {
  starts<-pmax(1,rep(starts,length.out=length(matrices)),na.rm=TRUE)
  if(is.null(ends)&&!is.null(lengths)) ends<-starts+lengths-1
  if(is.null(ends)) stop("No end positions specified\n")
  ends<-pmin(sapply(matrices, ncol), rep(ends,length.out=length(matrices)))
  mat.out<-matrices
  for(i in 1:length(matrices)) {
    mat.out[[i]] <- (matrices[[i]])[,starts[i]:ends[i]]
  }
  mat.out
}

# find most conserved bases in pwm matrices
get.core.matrices <- function(matrices, core.length=5) {
  cs<-numeric(length(matrices))+NA;
  mat.out<-matrices
  for(i in 1:length(matrices)) {
    matrix <- matrices[[i]]
    cmax<-matrix[cbind(max.col(t(matrix)),1:ncol(matrix))]
    if(any(cmax==-Inf)&&!all(cmax==-Inf)) cmax[cmax==-Inf] <- min(cmax[cmax>-Inf])-100000
    if(ncol(matrix) > core.length) {
      vec <- numeric(ncol(matrix)-core.length+1)
      for(j in 1:(ncol(matrix)-core.length+1)) {
        vec[j]<-sum(cmax[j:(j+core.length-1)])
      }
      cs[i]<-match(TRUE, max(vec,na.rm=TRUE)-vec<.Machine$double.eps^0.5)
      mat.out[[i]] <- (matrices[[i]])[,cs[i]:(cs[i]+core.length-1)]
    }
  }
  names(mat.out) <- paste(names(matrices), "_CORE_", cs, "-", cs+core.length-1, sep="")
  mat.out
}

# Read fasta file 
read.fasta <- function(file) {
  lines <- readLines(file, -1)
  names <- grep("^>", lines)
  sequences <- character(length(names))
  borders <- c(names, length(lines)+1)
  for(i in 1:length(names)) {
    sequences[i] <- paste(lines[(borders[i]+1):(borders[i+1]-1)], collapse="");
  }
  names<-lines[names]
  names <- sub("^> *([^ ]*).*", "\\1", names)
  names(sequences) <- names
  sequences
}

# Read clustalw file 
read.clustalw <- function(file){
  lines=readLines(file)
  if(!is.na(pmatch("CLUSTAL", lines[1]))) lines[1]=""
  lines=sub('^[ :.]*$', '', lines, perl=TRUE)
  lines=sub('^#.*', '', lines, perl=TRUE)
  breaks=lines==""
  blocks=cumsum(breaks)
  lines=lines[!breaks]
  blocks=blocks[!breaks]  
  lens=table(blocks)
  if(!all(lens==lens[1])) error("Not Clustal format. Blocks with different numbers of sequences")
  annos=sub("^(.*[^ ])  *[^ ]*$", "\\1", lines[blocks==names(lens)])
  lines=sub("^.* ([^ ]*)$", "\\1", lines)
  s=rep(1:lens[1], length(lens))
  lines=split(lines, s)
  lines=unlist(lapply(lines, paste, collapse=""), use.names=FALSE)
  names(lines)=annos
  lines
}

# Extract sequences from multiple chromosome style fasta files
read.subsequences <- function(files, starts=1, ends=NULL, lengths=NULL) {
  starts <- rep(starts, length.out=length(files))
  if(is.null(ends) && !is.null(lengths)) {
    ends <- rep(starts+lengths-1, length.out=length(files))
  } else if(!is.null(ends)) {
    ends <- rep(ends, length.out=length(files))
  } else {
    ends <- rep(-1, length.out=length(files))
  }
  ends[is.na(ends)] <- -1
  files <- as.character(files)
  sequences <- character(length(files))
  names <- character(length(files))
  for(f in unique(files)) {
    infile <- which(f == files)
    con <- file(f,"r")
    lines <- readLines(con, n = 2, ok = FALSE)
    seqstart <- nchar(lines[1])+1
    linelength <- nchar(lines[2])
    names[infile] <- rep(sub("^>([^ ]*)", "\\1", lines[1]),length.out=length(infile))
    s <- starts[infile] + seqstart + floor((starts[infile]-1)/linelength) - 1
    l <- ends[infile] + seqstart + floor((starts[infile]-1)/linelength) - s
    l[ends[infile]==-1] <- -1
    sequences[infile] <- apply(cbind(s,l), 1, function(p) {
      p[1]<-seek(con, where=p[1], origin = "start")
      if(p[2]==-1) {
        return(paste(readLines(con, n = -1, ok = TRUE),collapse=""))
      } else {
        return(gsub("[^a-z_-]","",readChar(con, p[2]),ignore.case=TRUE))
      }
    })
    close(con)
  }
  ends[ends==-1] <- nchar(sequences[ends==-1]) - starts[ends==-1] + 1
  names(sequences) <- paste(names, ":", starts, "-", ends, sep="")
  return(sequences)
}

# Regularize Matrices with the methods of Rahmann et al. 2003
regularize.count.matrices <- function(matrices, method="Hdiff", background.dist=NULL, weight=1.5) {
  if(is.null(background.dist))
    background.dist <- rep(1,nrow(matrices[[1]]))/nrow(matrices[[1]])
  lapply(matrices, regularize.count.matrix, method=method, background.dist=background.dist, weight=weight)
}

# Method can be constant, Hfactor or Hdiff
regularize.count.matrix <- function(matrix, method="Hdiff", background.dist=NULL, weight=1.5) {
  osums <- colSums(matrix)
  nrow<-nrow(matrix)
  if(is.null(background.dist))
    background.dist <- rep(1,nrow(matrix))/nrow
  if(method=="constant") {
    m<-matrix+(background.dist*weight)
  } else {
    vdist <- matrix/rep(osums,each=nrow)
    oent <- 2*colSums(matrix)*colSums(vdist*log(vdist/background.dist), na.rm=T) # Rahmann entropy
    tent <- oent - weight
    if(method=="Hfactor") tent <- oent * weight 
    wtH<-rep(1, length(tent))
    wtL<-rep(0, length(tent))
    wtg<-rep(0.5, length(tent))
    for(i in 1:50) {
      ndist<- (1-rep(wtg,each=nrow))*vdist + rep(wtg,each=nrow)*background.dist
      nent <- 2*colSums(matrix)*colSums(ndist*log(ndist/background.dist), na.rm=T)
      inc <- nent > tent
      wtL[inc] <- wtg[inc]
      wtH[!inc] <- wtg[!inc]
      wtg = (wtL+wtH)/2
    }
    m<-(1-rep(wtg,each=nrow))*vdist + rep(wtg,each=nrow)*background.dist
  }
  nsums <- colSums(m)
  m <- m * rep(osums/nsums, each=nrow)
  m
}

count2profile <- function(matrices, add.pseudocounts=0) {
  lapply(matrices, function(m) {
    (m+add.pseudocounts)/rep(colSums(m+add.pseudocounts), each=nrow(m))
  })
}

profile2score <- function(matrices, background.dist=NULL) {
  if(is.null(background.dist))
    background.dist <- rep(1,nrow(matrices[[1]]))/nrow(matrices[[1]])
  lapply(matrices, function(m) {
    log2(m/background.dist)
  })
}

pwm2consensusSequence <- function(matrix) {
  if(is.null(rownames(matrix))) stop("Wrong format for matrix. Needs rownames.\n")
  m<-matrix==matrix(matrix[cbind(max.col(t(matrix)),1:ncol(matrix))],nrow(matrix),ncol(matrix),byrow=TRUE)
  m<-paste(rep(c("[",rownames(matrix),"]"),ncol(matrix))[rbind(TRUE,m,TRUE)],collapse="")
  if(nrow(matrix)==4) {
    m<-pattern2sequence(m,"dna")
  } else {
    m<-pattern2sequence(m,"peptide")
  }
  m
}

pwm2pattern <- function(pwm) {
  if(is.null(rownames(pwm))) stop("Wrong format for matrix. Needs rownames.\n")
  gsub("\\[(.)\\]","\\1",paste(rep(c("[",rownames(pwm),"]"),ncol(pwm))[rbind(TRUE,pwm>0,TRUE)],collapse=""))
}

background.dist.from.seq <- function(sequence, letters="dna") {
  letters<-toupper(letters)
  if(letters == "DNA") letters="ACGT"
  if(letters == "PEPTIDE") letters="ACDEFGHIKLMNPQRSTVWY"
  t<-table(strsplit(gsub(paste("[^",letters,"]",collapse=""),"",toupper(sequence)),""))
  t/sum(t)
}

random.profile.sequence <- function(profile.matrix, nseq=1) {
  if(is.null(rownames(profile.matrix))) {
    if(nrow(profile.matrix)==4) letters="acgt"
    if(nrow(profile.matrix)==20) letters="ACDEFGHIKLMNPQRSTVWY"
  } else {
    letters<-as.character(paste(rownames(profile.matrix),collapse=""))[1]
  }
  nletters=as.integer(nchar(letters))
  if(nletters != nrow(profile.matrix) || !identical(all.equal(colSums(profile.matrix), rep(1,ncol(profile.matrix))), TRUE))
    stop("Profile matrix  needs to be have column sums equal to 1 and rownames indicating letters.")
  randseq<-matrix("n",nseq,ncol(profile.matrix))
  if(!nseq>0) stop("At least one sequence needs to be generated.")
  for(i in 1:ncol(profile.matrix)) {
    randseq[,i] <- .C("random_strings", character(nseq), as.integer(rep.int(1,nseq)), as.integer(nseq),
                      letters, as.numeric(profile.matrix[,i]), nletters)[[1]]
  }
  apply(randseq,1,paste,collapse="")
}

pattern2countmatrix <- function(pattern, letters="dna", nobs=1) {
  letters<-toupper(letters)
  if(letters == "DNA") {
    letters="ACGT"
    sequence2pattern(pattern, "dna")
  }
  if(letters == "PEPTIDE") {
    letters="ACDEFGHIKLMNPQRSTVWY"
    sequence2pattern(pattern, "peptide")
  }
  letters<-unique(unlist(strsplit(letters,"")))
  pattern <- gsub("\\[([^\\]]*)\\]","[\\#\\1]", pattern)
  subpattern1<-unlist(strsplit(pattern, "[\\]\\[]"))
  subpattern2<-strsplit(subpattern1,"")
  subpattern2[grep("\\#",subpattern1)] <- gsub("\\#","",grep("\\#",subpattern1, value=TRUE))
  pattern<-strsplit(unlist(subpattern2),"")
  matrix<-matrix(0,length(letters),length(pattern))
  rownames(matrix)<-letters
  for(i in 1:length(pattern)) {
    matrix[match(pattern[[i]],letters),i]<-nobs/length(pattern[[i]])
  }
  matrix
}

pattern2sequence <- function(pattern, letters="dna") {
  letters<-toupper(letters)
  sequence<-toupper(pattern)
  sequence <- gsub("\\[(.)\\]","\\1",sequence)
  if(letters == "DNA" || letters == "ACGT") {
    sequence <- gsub("\\[AC\\]","M",sequence)
    sequence <- gsub("\\[CA\\]","M",sequence)
    sequence <- gsub("\\[AG\\]","R",sequence)
    sequence <- gsub("\\[GA\\]","R",sequence)    
    sequence <- gsub("\\[AT\\]","W",sequence)
    sequence <- gsub("\\[TA\\]","W",sequence)
    sequence <- gsub("\\[CG\\]","S",sequence)
    sequence <- gsub("\\[GC\\]","S",sequence)
    sequence <- gsub("\\[CT\\]","Y",sequence)
    sequence <- gsub("\\[TC\\]","Y",sequence)
    sequence <- gsub("\\[GT\\]","K",sequence)
    sequence <- gsub("\\[TG\\]","K",sequence)
    sequence <- gsub("\\[ACG\\]","V",sequence)
    sequence <- gsub("\\[ACT\\]","H",sequence)
    sequence <- gsub("\\[AGT\\]","D",sequence)
    sequence <- gsub("\\[CGT\\]","B",sequence)
    sequence <- gsub("\\.","N",sequence)
    sequence <- gsub("\\[[^\\]]*\\]","N",sequence)
  } else if(letters == "PEPTIDE" || letters == "ACDEFGHIKLMNPQRSTVWY") {
    sequence <- gsub("\\[DN\\]","B",sequence)
    sequence <- gsub("\\[ND\\]","B",sequence)
    sequence <- gsub("\\[EQ\\]","Z",sequence)
    sequence <- gsub("\\[QE\\]","Z",sequence)
    sequence <- gsub("\\.","X",sequence)
    sequence <- gsub("\\[[^\\]]\\]","X",sequence)
  } else {
    warning("Error in input: letters must indicate \"dna\" or \"peptide\".")
  }
  sequence
}

sequence2pattern <- function(sequence, letters="dna") {
  letters<-toupper(letters)
  sequence<-toupper(sequence)
  sequence<-gsub("[^\\]\\[\\)\\(\\}\\{A-Z]","",sequence)
  if(letters == "DNA" || letters == "ACGT") {
    sequence <- chartr("U(){}","T[][]",sequence)
    sequence <- gsub("M","[AC]",sequence)
    sequence <- gsub("R","[AG]",sequence)
    sequence <- gsub("W","[AT]",sequence)
    sequence <- gsub("S","[CG]",sequence)
    sequence <- gsub("Y","[CT]",sequence)
    sequence <- gsub("K","[GT]",sequence)
    sequence <- gsub("V","[ACG]",sequence)
    sequence <- gsub("H","[ACT]",sequence)
    sequence <- gsub("D","[AGT]",sequence)
    sequence <- gsub("B","[CGT]",sequence)
    sequence <- gsub("N",".",sequence)
    sequence <- gsub("X",".",sequence)    
  } else if(letters == "PEPTIDE" || letters == "ACDEFGHIKLMNPQRSTVWY") {
    sequence <- chartr("(){}","[][]",sequence)
    sequence <- gsub("B","[DN]",sequence)
    sequence <- gsub("Z","[EQ]",sequence)
    sequence <- gsub("X",".",sequence)
  } else {
    warning("Error in input: letters must indicate \"dna\" or \"peptide\".")
  }
  sequence
}

# Sequence Logos as defined by Schneider et al, NAR, 1990
simple.sequence.logo.plot <- function(matrix, e=0) {
  entropy <- -colSums(matrix *log2(matrix), na.rm=TRUE) # Shannon entropy
  R <- log2(nrow(matrix)) - entropy + e
  hights <- t(t(matrix) * R)
  barplot(hights, col=rainbow(nrow(matrix)))
}

random.strings <- function(lengths, background.dist=NULL, letters="dna") {
  if(letters == "dna") letters="acgt"
  if(letters == "peptide") letters="ACDEFGHIKLMNPQRSTVWY"
  nletters=as.integer(nchar(letters))
  if(is.null(background.dist)) background.dist <- rep(1,nletters)/nletters
  if(nletters != length(background.dist) || !identical(all.equal(sum(background.dist), 1),TRUE) )
    stop("dist needs to be a parameter with the same length as number of letters and sum equal to 1, or NULL for uniform distributions.")
  .C("random_strings", character(length(lengths)), as.integer(lengths), length(lengths),
     as.character(letters)[1], as.numeric(background.dist), nletters)[[1]]
}

# Compute score distribution of a PWM according to Rahmann et al 2003
pwm.exact.scoredist <- function(pwm, background.dist=NULL, eps=0.001) {
  if(is.null(background.dist)) background.dist <- rep(1, length(pwm))/nrow(pwm)
  background.dist <- as.numeric(rep(background.dist, length.out=length(pwm)))
  maxscore <- sum(pwm[cbind(max.col(t(pwm)), 1:ncol(pwm))])
  minscore <- sum(pwm[cbind(max.col(t(-pwm)), 1:ncol(pwm))])
  if(eps>1) eps <- (maxscore-minscore)/eps
  len <- as.integer(1 + round((maxscore-minscore)/eps, 0));
  if(len>50000) eps <- (maxscore-minscore)/50000
  if(len<100) eps <- (maxscore-minscore)/100
  dens <- .C("pwmexactscoredist", as.numeric(pwm), background.dist, nrow(pwm), ncol(pwm), as.numeric(eps)[1],  
               numeric(len), len)[[6]]
  scores <- round(minscore,0) + (1:len)*eps
  return(list(minscore=minscore, maxscore=maxscore, expect=sum(dens*scores), scores=scores, dens=dens, prob=cumsum(dens),
              pwm=pwm, background.dist=matrix(background.dist, nrow(pwm)), eps=eps))
}

pwm.quality <- function(pwm, profile.matrix, background.dist=NULL, eps=0.001, window=500, instances=1, plot=TRUE) {
  negative <- pwm.exact.scoredist(pwm, background.dist, eps)
  positive <- pwm.exact.scoredist(pwm, profile.matrix, eps)
  scores<-negative$scores
  e1 <- 1-exp(window*log(pmax(0,negative$prob)))
  e2 <- 1-exp(instances*log(instances * (pmax(0,1-positive$prob))))
  names(e1) <- scores
  names(e2) <- scores
  pbal <- min(c(which(e1<e2), length(e1)), na.rm=T)
  Qbal = 1-e2[pbal]
  tbal = scores[pbal]
  if(plot) {
    par(mfrow=c(1,2))
    plot(e1,e2, col="blue", type="l", main=paste("ROC curve (Qbal=", round(Qbal,3), ")", sep=""),
         xlab="Type I error", ylab="Type II error", xlim=c(0,1), ylim=c(0,1))
    abline(v=0.05, lty="dotted")
    abline(0, 1, lty="dotted")
    abline(h=0.05, lty="dotted")
    text(e1[pbal], e2[pbal], paste("tbal=", round(tbal,2), sep=""))
    plot(scores,e2, col="green", type="l", main=paste("Scoredist (m=", window, ", n=", instances, ")", sep=""),
         xlab="Scores", ylab="Errors", ylim=c(0,1))
    lines(scores,e1, col="red")
    legend(scores[1], 1, c("Type I errors", "Type II errors"), pch=20, col=c("red", "green"), bg="white")
  }
  return(list(minscore=negative$minscore, maxscore=negative$maxscore, Qbal=Qbal, tbal=tbal,
              expect1=negative$expect, expect2=positive$expect, type1=e1, type2=e2))
}
         
pwm.scoredist <- function(pwm, sequence=NULL, npermutations=1000) {
  if(is.character(sequence)) {
  } else if(is.numeric(sequence) && length(sequence)==5) {
    sequence <- random.strings(sequence[1], nseq=1, background.dist=sequence[2:5], letters="dna")
  } else if(is.numeric(sequence) && length(sequence)==21) {
    sequence <- random.strings(sequence[1], nseq=1, background.dist=sequence[2:5], letters="peptide")
  } else {
    len <- 500
    if(is.numeric(sequence) && length(sequence)==1) len <- as.integer(sequence)
    if(nrow(pwm)==4) {
      sequence <- paste(rep(c("acgt"), floor(len/4)), collapse="")
    } else {
      sequence <- paste(rep("ACDEFGHIKLMNPQRSTVWY", floor(len/20)), collapse="")
    }    
  }
  sequence <- as.character(sequence)[1]
  npermutations <- as.integer(npermutations)[1]
  .C("pwmscoredist", sequence, nchar(sequence), pwm, nrow(pwm), ncol(pwm),
     numeric(npermutations), npermutations)[[6]]

}

pwm.scores <- function(sequence, pwm) {
  sequence <- as.character(sequence)[1];
  if(nchar(sequence)<ncol(pwm)) stop("Sequence shorter than matrix\n");
  scores<-.C("pwm_scores", sequence, nchar(sequence), as.numeric(pwm), nrow(pwm), ncol(pwm),
                  numeric(nchar(sequence)))[[6]]
  scores[1:(nchar(sequence)-ncol(pwm)+1)]
}

conserved.pwm.scores <- function(alignment, pwm) {
  if(length(unique(nchar(alignment)))!=1) stop("Sequences in alignment must all have the same lengths.");
  if(nchar(alignment[1])<ncol(pwm)) stop("Sequence shorter than matrix\n");
  scores <- .C("conserved_pwm_scores", as.character(alignment), nchar(alignment[1]), length(alignment),
               as.numeric(pwm), nrow(pwm), ncol(pwm), numeric(nchar(alignment[1])*(length(alignment)+1)))[[7]]
  rn <- names(alignment)
  if(is.null(rn)) rn <- 1:length(alignment)
  dim(scores)<-c(nchar(alignment[1]),(length(alignment)+1));
  colnames(scores) <- c(rn, "Cons.Scores")
  scores;
}

conserved.pwm.scoredist <- function(pwm, alignment, nsample=1000) {
  if(length(unique(nchar(alignment)))!=1) stop("Sequences in alignment must all have the same lengths.");
  scores <- .C("conserved_pwm_scoredist", as.character(alignment), nchar(alignment[1]), length(alignment),
               as.numeric(pwm), nrow(pwm), ncol(pwm), as.integer(nsample)[1],
               numeric(as.integer(nsample)[1]*(length(alignment)+1)))[[8]]
  rn <- names(alignment)
  if(is.null(rn)) rn <- 1:length(alignment)
  dim(scores)<-c(as.integer(nsample)[1], (length(alignment)+1));
  colnames(scores) <- c(rn, "Cons.Scores")
  scores;
}

scores2pvalues <-function(scores, scoredist, method="scoredist") {
  if(method=="precomputed.pval") {
    pvals <- as.numeric(scoredist)
    scoredist <- as.numeric(names(scoredist))
    use.pvals=as.integer(1)
  } else {
    pvals <- numeric(1)
    scoredist <- as.numeric(scoredist)
    use.pvals=as.integer(0)
  }
  scores[is.na(scores)|scores==-Inf] <- scoredist[1]-1
  .C("scores2pvalues", as.numeric(scores), length(scores), scoredist, pvals, length(scoredist),
     use.pvals, numeric(length(scores)))[[7]]
}

revcomp.pwms <- function(pwms) {
  lapply(pwms, function(m) {
    matrix(rev(m),nrow=nrow(m),dimnames=dimnames(m))
  })
}

rev.strings <- function(strings) {
  lstrings=as.character(strings) ## make copy here - something wrong with duplicate function in .Call
  lstrings[is.na(lstrings)]=NA
  .Call("revstrings", strings)
}

revcomp <- function(sequences) {
  chartr("ACGTUMRYKVHDBacgtumrykvhdb", "TGCAAKYRMBDHVtgcaakyrmbdhv", rev.strings(sequences))
}

filter.multimatch.pwms <- function(matches, keep.best=50) {
  keep <-logical(nrow(matches))|TRUE
  pwmscores<-aggregate(matches$PWM.Score, list(matches$PWM.Name, matches$Match.StartPos), max)    
  names(pwmscores) <- c("mat", "pos", "score")
  o <- order(pwmscores$score, pwmscores$mat, decreasing=TRUE)
  pwmscores <- pwmscores[o,]
  frequentfactors<-names(which(table(pwmscores$mat)>keep.best))
  if(!is.null(matches$Cons.Score)) {
    consscores<-aggregate(matches$Cons.Score, list(matches$PWM.Name, matches$Match.StartPos), max)    
    names(consscores) <- c("mat", "pos", "score")
    o <- order(consscores$score, consscores$mat, decreasing=TRUE)
    consscores <- consscores[o,]
    frequentfactors<-unique(c(frequentfactors, names(which(table(consscores$mat)>keep.best))))
  }
  for(factor in frequentfactors) {
    bestpos<-as.integer(as.character(pwmscores$pos[pwmscores$mat == factor]))
    if(length(bestpos)>keep.best) bestpos<-bestpos[1:keep.best]
    tmp <- matches$PWM.Name == factor
    keep[tmp] <- !is.na(match(matches$Match.StartPos[tmp], bestpos))
    if(!is.null(matches$Cons.Score)) {
      bestpos<-as.integer(as.character(consscores$pos[consscores$mat == factor]))
      if(length(bestpos)>keep.best) bestpos<-bestpos[1:keep.best]
      tmp <- matches$PWM.Name == factor
      keep[tmp] <- keep[tmp] | !is.na(match(matches$Match.StartPos[tmp], bestpos))
    }
  }
  matches <- matches[keep,]
  return(matches)
}

pattern.scan <- function(sequences, patterns, rev.strand = TRUE, cutoff = 0.1, pmethod = 1000,
                         perl=FALSE, fixed=FALSE) {
  res <- list(Sequence.Name=character(), Pattern=character(), Match.StartPos=integer(),
              Match.EndPos=integer(), Pvalue=numeric(), Match.Ori=character(), seqs=character())
  t2 <- 0
  ns = as.character(1:length(sequences))
  if(!is.null(names(sequences))) ns = names(sequences)
  np=patterns
  if(!is.null(names(patterns))) np = names(patterns)
  if(rev.strand) {
    r <- revcomp(sequences)
    nr <- nchar(r)
  }
  for(s in 1:length(sequences)) {
    if(pmethod>0) rand = permute.sequences(rep(sequences[s], pmethod))
    if(rev.strand && pmethod>0) randr = revcomp(rand)
    for(p in 1:length(patterns)) {
      ts = integer()
      te = integer()
      pval = c(1,1)
      ori = character()
      lts=c(0,0)
      t1 = gregexpr(patterns[p], sequences[s], perl=perl, fixed=fixed)[[1]]
      if(t1[1]>-1) {
        ts = as.integer(t1)
        te = ts + as.integer(attr(t1,"match.length")) - 1
        lts[1] = length(ts)
        ori = rep("+", lts[1])
      }
      if(rev.strand) {
        t1 = gregexpr(patterns[p], r[s], perl=perl, fixed=fixed)[[1]]
        if(t1[1]>-1) {
          ts = c(ts, nr[s] - as.integer(t1) - as.integer(attr(t1,"match.length")) + 2)
          te = c(te, nr[s] - as.integer(t1) +1)
          lts[2] = length(t1)
          ori <- c(ori, rep("-", lts[2]))
        }
      }
      if(sum(lts)<1) next
      if(pmethod > 0) {
        if(lts[1] > 0)
          pval[1] = max(c(0,length(grep(patterns[p], rand, fixed=fixed, perl=perl))), na.rm=TRUE) / pmethod
        if(rev.strand)
          pval[2] = max(c(0,length(grep(patterns[p], randr, fixed=fixed, perl=perl))), na.rm=TRUE) / pmethod
      }
      if(!is.na(cutoff) && cutoff>0 && sum(pval<=cutoff)<2) {
        t1=rep.int(pval<=cutoff,times=lts)
        ts=ts[t1]
        te=te[t1]
        ori=ori[t1]
        pval=rep.int(pval, times=lts)[t1]
      } else {
        pval=rep.int(pval, times=lts)
      }
      if(length(ts)>0) {
        if(lts[2] > 0) {
          o=order(ts)
          ts=ts[o]
          te=te[o]
          ori=ori[o]
          pval=pval[o]
        }
        lts=length(ts)
        res$Sequence.Name <- c(res$Sequence.Name, rep(ns[s], lts))
        res$Pattern <- c(res$Pattern, rep(np[p], lts))
        res$Match.StartPos <- c(res$Match.StartPos, ts)
        res$Match.EndPos <- c(res$Match.EndPos, te)
        res$Pvalue <- c(res$Pvalue, pval)
        res$Match.Ori <- c(res$Match.Ori, ori)
        res$seqs <- c(res$seqs, substring(sequences[s], ts, te))
      }
    }
  }
  return(as.data.frame(res))
}

pwm.scan <- function(sequences, pwms, rev.strand=TRUE, cutoff=0.1, pmethod=1000) {
  if(rev.strand) {
    rev.pwms <- revcomp.pwms(pwms)
  } else {
    rev.pwms <- NULL
  }
  as.data.frame(.Call("pwmscan", sequences, pwms, rev.pwms, as.numeric(cutoff), pmethod))
}

conserved.pwm.scan <- function(alignment, pwms, rev.strand=TRUE, cutoff=0.1, pmethod=1000) {
  if(rev.strand) {
    rev.pwms <- revcomp.pwms(pwms)
  } else {
    rev.pwms <- NULL
  }
  as.data.frame(.Call("conserved_pwmscan", alignment, pwms, rev.pwms, as.numeric(cutoff), as.numeric(pmethod)))
}

permute.sequences <- function(sequences) {
  sequences <- as.character(sequences)
  sequences.new <- sequences
  .C("permute_sequences", sequences, sequences.new, nchar(sequences), length(sequences))[[2]]
}

moving.average <- function(numbers, window.size=10) {
  numbers <- as.numeric(numbers)
  window.size=as.integer(window.size)
  if(window.size<=0) stop("window.size needs to be positive integer.")
  if(length(numbers)<1) stop("numbers have length 0.")
  numbers <- .C("moving_average", numbers, length(numbers), window.size)[[1]]
  numbers <- numbers[1:max(1,length(numbers)-window.size+1)]
  numbers
}

window.average <- function(numbers, positions, window.size=1000) {
  numbers <- as.numeric(numbers)
  positions <- as.numeric(positions)
  window.size=as.numeric(window.size)
  if(window.size<=0) stop("window.size needs to be positive number.")
  if(length(numbers)<1) stop("numbers have length 0.")
  if(any(is.na(positions))) stop("No NAs in positions vector allowed.")
  if(length(numbers)!=length(positions)) stop("numbers and positions needs to have the same length.")
  numbers <- .C("window_average", numbers, numbers, positions, length(numbers), window.size)[[2]]
  names(numbers)=positions
  numbers
}

randomize.alignment <- function(alignment, letters="dna") {
  if(is.numeric(letters)) {
    nletters<-as.integer(letters)
  } else {
    letters<-toupper(letters)
    if(letters == "DNA") {
      nletters=4
    } else if(letters == "PEPTIDE") {
      nletters=20
    } else {
      nletters=nchar(letters)
    }
    nletters<-as.integer(nletters)
  }
  alignment <- as.character(alignment)
  ali.new <- alignment
  .C("randomize_alignment", alignment, ali.new, nchar(alignment[1]), length(alignment), nletters)[[2]]
}


# read tables in UCSC Track format
read.ucsc.track <- function(track, sql=NULL) {
  if(is.null(sql)) sql=sub("\\.txt$", ".sql", track)
  sql=readLines(sql)
  start = grep("^CREATE", sql)
  end=min(c(grep("^ *KEY", sql), grep("^ *PRIMARY KEY", sql), grep("^)", sql)))
  if(start+1<end) {
    sql = sql[(start+1):(end-1)]
    sql = sub("^ *`?([^` ]*)`?.*$", "\\1", sql)
    sql = sub("^txStart$", "start", sql)
    sql = sub("^txEnd$", "end", sql)
    sql = sub("^chromStart$", "start", sql)
    sql = sub("^chromEnd$", "end", sql)
  } else {
    sql=NULL
  }
  track<-read.table(track, header=FALSE, sep="\t", quote="",
                    stringsAsFactors=FALSE, comment.char="")
  if(!is.null(sql)) colnames(track) = sql
  track
}

# extract subregion from annotation track
track.region <- function(track, chrom, start, end) {
  track[track[,"chrom"] %in% chrom &
  track[,"end"]>start & track[,"start"]<end,]
}

# plot from annotation track
plot.tracks <- function(tracks, chrom, start, end, track.offset=1, gene.offsets=0.6,
                        col="black", width=0.4, chrom.pos=FALSE, chrom.scale=1000, type="b",
                        track.names=TRUE, ...) {
  if(is.data.frame(tracks)) {
    tracks = track.region(tracks, chrom, start, end)
    n = (nrow(tracks)-1)*gene.offsets[1]
  } else if(is.list(tracks)) {
    track.lengths = numeric(length(tracks))
    for(i in 1:length(tracks)) {
      tracks[[i]] = track.region(tracks[[i]], chrom, start, end)
      track.lengths[i]=nrow(tracks[[i]])
    }
    gene.offsets = rep(gene.offsets, length.out=length(tracks))
    col = rep(col, length.out=length(tracks))
    type = rep(type, length.out=length(tracks))
    ypos = (pmax(0,track.lengths-1)*gene.offsets)+(track.offset*(track.lengths>0))
    n = sum(ypos, na.rm=TRUE)-track.offset
    ypos = cumsum(c(width, ypos))
    if(!track.names) names(tracks)=NULL
  } else {
    warning("tracks must be list or data.frame.\n")
    n = -1
  }
  if(n<0) {
    warning("No genes in region\n")
    plot(1:10, 1:10, type="n", axes=FALSE, xlab="", ylab="")
    return()
  }
  plot(c(start, end)/chrom.scale, type="n", xlab="", ylab="", xlim=c(start, end)/chrom.scale,
       ylim=c(n+2*width, 0), axes=FALSE)
  if(chrom.pos) { axis(1); mtext(paste("Chr. position in", chrom.scale, "bases", sep=" "), 1, line=2) }
  if(is.data.frame(tracks)) {
    lines.track(tracks, at=width, offset=gene.offsets[1], width=width, col=col,
                chrom.scale=chrom.scale, type=type, ...)
  } else {
    for(i in 1:length(tracks)) {
      lines.track(tracks[[i]], at=ypos[i], offset=gene.offsets[i], width=width, col=col[i],
                  track.name=names(tracks)[i], chrom.scale=chrom.scale, type=type[i], ...)
    }
  }
  tracks
}

# plot single track
lines.track <- function(track, at=0, offset=1, width=0.6, col="black", track.name=NA, names=TRUE,
                        lty=1, line=0, las=2, adj=1, padj=0.5, type="b", chrom.scale=1000,
                        border=FALSE, ...) {
  if(nrow(track)<1) return()
  ypos = cumsum(c(at, rep(offset, nrow(track)-1)))
  coords=par("usr")
  if(!("start" %in% colnames(track) && "end" %in% colnames(track)))
    stop("Track needs to have 'start' and 'end' column.", track.name)
  if(names&&offset>0&&"name"%in%colnames(track))
    mtext(as.character(track[,"name"]), side=2, line=line, at=ypos, adj=adj, padj=padj, col=col, las=las)
  if("exonStarts" %in% colnames(track) && "exonEnds" %in% colnames(track)) {
    if(type=="b"||type=="l") {      
      segments(pmax(coords[1], track[,"start"]/chrom.scale), ypos,
               pmin(coords[2], track[,"end"]/chrom.scale), ypos, col=col, lty=lty)
    }
    starts=strsplit(track[, "exonStarts"], split=",")
    l=sapply(starts, length)
    starts=as.numeric(unlist(starts, use.names=FALSE))
    ends=strsplit(track[, "exonEnds"], split=",")
    ends=as.numeric(unlist(ends, use.names=FALSE))
    s = starts>=rep.int(track[,"start"], times=l) & ends<=rep.int(track[,"end"], times=l)
    starts=starts[s]
    ends=ends[s]
    ypos=rep.int(ypos,l)
  } else {
    starts = track[,"start"]
    ends = track[,"end"]
    s=!is.na(starts)
    l= rep(1, nrow(track))
  }
  if("color" %in% colnames(track)) {
    cols=unlist(strsplit(as.character(track[, "color"]), split=","), use.names=FALSE)
    if(length(cols)==nrow(track)) cols=rep.int(cols, times=l)
    cols=cols[s]
  } else {
    cols=col
  }
  if(border==FALSE) border=cols
  if(type=="b"||type=="r") rect(pmax(coords[1], starts/chrom.scale), ypos-width/2,
                  pmin(coords[2], ends/chrom.scale), ypos+width/2, col=cols, border=border)
  if(is.character(track.name)) text(mean(coords[1:2]), at-width, track.name, adj=c(0.5, 0))
}

colorbar <- function(col, lim, ..., ncol=NULL, decreasing=FALSE, border=FALSE, xlab="") {
  if(is.vector(col) && is.null(ncol) && length(col)>1) ncol = length(col)
  if(is.numeric(lim) && length(lim)==2) {
    if(is.null(ncol)) ncol=100
    col = get.colors(seq(lim[1],lim[2],length=ncol), col, ncol, decreasing, lim)
    labels=TRUE
    at=NULL
  } else if(is.vector(lim)) {
    ncol=length(lim)
    col = get.colors(1:ncol, col, ncol, decreasing, NULL)
    labels=lim
    lim=c(0,length(lim))
    at=1:length(labels)-0.5
  } else {
    stop("lim must be range or color labels.")
  }
  plot(0, 0, type="n", xlim=lim, ylim=c(0, 1), axes=FALSE, xlab=xlab, ylab="")
  u=par("usr")
  d=abs(diff(lim))
  l=min(lim, na.rm=TRUE)
  s=((0:(ncol-1))*d/ncol)+l
  e=((1:ncol)*d/ncol)+l
  s[1]=u[1]
  e[length(e)]=u[2]
  rect(s, u[3], e, u[4], col = col, border = border)
  axis(1, at=at, labels=labels, ...)
  box()
}

# convert numeric values to colors
get.colors <- function(x, col="rainbow", ncol=NULL, decreasing=FALSE, lim=NULL, na.col="#CBD5E8") {
  if(is.character(x)) x=as.factor(x)
  if(is.vector(col) && is.null(ncol) && length(col)>1) ncol = length(col)
  if(is.factor(x)) {
    if(is.null(ncol)||!ncol<nlevels(x)) ncol=nlevels(x)
    if(is.null(lim)) lim=c(1, nlevels(x))
    x=as.integer(x)
  }
  if(is.null(ncol)) stop("Need ncol argument!")
  if(is.character(col) && length(col)==1) {
    r=seq(0,1,length=ncol)
    r2=seq(-1,1,length=ncol)
    col=switch(col,
      heat=hcl(h=90-90*r, c=30+50*r^0.2, l=90-60*r^1.5),
      blackheat=hcl(h=60+20*r, c=230+50*r^0.2, l=30+60*r^1.5),
      terrain=hcl(h=0+130*r, c=0+65*r^0.33, l=95-50*r^1.5),
      blackterrain=hcl(h=0+130*r, c=30+65*r^0.33, l=20+40*r^1.5),
      blue=hcl(h=250, c=0+80*r^1.5, l=92-72*r^1.5),
      blackblue=hcl(h=250, c=20+80*r^0.75, l=10+20*r^1.5),
      red=hcl(h=10, c=0+80*r^1.5, l=92-72*r^1.5),
      blackred=hcl(h=90-80*r, c=80+100*r^0.5, l=10+35*r^1.5),
      green=hcl(h=110, c=0+80*r^1.5, l=92-72*r^1.5),
      blackgreen=hcl(h=110, c=80+100*r^0.5, l=10+35*r^1.5),
      gray=hcl(h=0, c=0, l=92-72*r^1.5),
      grey=hcl(h=0, c=0, l=92-72*r^1.5),
      bluered=hcl(h=ifelse(r2>0, 10, 250), c=0+80*abs(r2)^1.5, l=92-72*abs(r2)^1.5),
      redgreen=hcl(h=ifelse(r2>0, 110, 10), c=0+80*abs(r2)^1.5, l=92-72*abs(r2)^1.5),
      greenblue=hcl(h=ifelse(r2>0, 250, 110), c=0+80*abs(r2)^1.5, l=92-72*abs(r2)^1.5),
      greenred=hcl(h=ifelse(r2>0, 90-90*abs(r2), 90+90*abs(r2)), c=30+50*abs(r2)^0.2, l=92-72*abs(r2)^1.5),
      rainbow=hcl(seq(0, 360*(ncol-1)/ncol, length=ncol), 80, 60))
    if(is.null(col)) stop("Not recognized col argument.")
  } else if(is.function(col)) {
    col=col(ncol)
  } else if(is.vector(col)) {
  } else {
    stop("col in wrong format.")
  }
  if(length(col)<ncol) ncol=length(col)
  if(!is.numeric(x)) x=as.numeric(x)
  if(!is.numeric(lim) || length(lim)!=2) lim=range(x)
  x=findInterval(x, seq(lim[1],lim[2],diff(lim)/ncol)[-c(1,ncol+1)])+1
  if(decreasing) col=rev(col)
  res=col[x]
  if(!is.null(na.col)) res[is.na(res)]=na.col
  if(!is.null(dim(x))) dim(res)=dim(x)
  res
}

show.color.matrix <- function(colors) {
  ncol=ncol(colors)
  nrow=nrow(colors)
  plot(0, 0, type = "n", xlim = c(0, ncol), ylim = c(nrow, 0), axes = FALSE,
           xlab = "", ylab = "")
  s = rep(0:(ncol - 1), each=nrow)
  e = s+1
  u = rep(0:(nrow - 1), ncol)
  l=u+1
  rect(s, l, e, u, col = as.vector(t(colors)), border = NA)   
}

col.chromosome.bands <- function(cytobands) {
  if(is.null(cytobands[,"color"])) {
    cytobands[,"color"] = cytobands[,"gieStain"]
    cytobands[,"color"] = sub("^acen$", "blue", cytobands[,"color"])
    cytobands[,"color"] = sub("^stalk$", "red", cytobands[,"color"])
    cytobands[,"color"] = sub("^gvar$", "purple", cytobands[,"color"])
    cytobands[,"color"] = sub("^gneg$", "gray90", cytobands[,"color"])
    g=grep("^gpos", cytobands[,"color"])
    cytobands[g,"color"] = paste("gray",
               round(100-as.numeric(sub("^gpos", "", cytobands[g,"color"])), 0), sep="")
  }
  cytobands
}

# split chromosomal position from format to data frame
split.chrompos <- function(chrompos) {
  if(!is.character(chrompos)) {
    names=names(chrompos)
    chrompos=as.character(chrompos)
    names(chrompos)=names
  }
  test=rep(TRUE, length(chrompos))
  test[grep("^[^:]*:[0-9]*-[0-9]*$", chrompos)]=FALSE
  chrompos[test]=""
  chrom = sub(":.*$", "", chrompos)
  start = sub("^.*:", "", chrompos)
  end =  as.numeric(sub("^.*-", "", start))
  start = as.numeric(sub("-.*$", "", start))
  data.frame(chrom, start, end, row.names=names(chrompos), check.names=FALSE)
}

paste.chrompos <- function(chrom, start, end) {
  paste(chrom, ":", formatC(start, 0,9,"d",flag=0), "-", formatC(end, 0,9,"d",flag=0), sep="")
}

# convert long format with each chromosomal region in one line to track format
regions.to.track <- function(table, colnames=c(name="name", chrompos="chrompos")) {
  names = table[,colnames["name"]]
  colnames = colnames[-match("name", names(colnames))]  
  if("chrompos" %in% names(colnames)) {
    chrompos=split.chrompos(table[,colnames["chrompos"]])
    colnames = colnames[-match("chrompos", names(colnames))]
    chrom = as.character(chrompos$chrom)
    start = chrompos$start
    end = chrompos$end
  } else {
    chrom = table[,colnames["chrom"]]
    colnames = colnames[-match("chrom", names(colnames))]
    start = as.numeric(table[,colnames["start"]])
    colnames = colnames[-match("start", names(colnames))]
    end = as.numeric(table[,colnames["end"]])
    colnames = colnames[-match("end", names(colnames))]  
  }
  exonStarts = unlist(tapply(start, names, paste, collapse=","), use.names=FALSE)
  exonEnds = unlist(tapply(end, names, paste, collapse=","), use.names=FALSE)
  start = unlist(tapply(start, names, min, na.rm=TRUE), use.names=FALSE)
  end = unlist(tapply(end, names, max, na.rm=TRUE), use.names=FALSE)
  chrom = unlist(tapply(chrom, names, head, n=1), use.names=FALSE)
  df=data.frame(name=sort(unique(names)), chrom, start, end, exonStarts, exonEnds)
  for(i in names(colnames)) {
    if(length(unlist(tapply(table[,colnames[i]], names, unique), use.names=FALSE))==nrow(df)) {
      df[[i]] = unlist(tapply(table[,colnames[i]], names, head, n=1), use.names=FALSE)
    } else {
      df[[i]] = unlist(tapply(table[,colnames[i]], names, paste, collapse=","), use.names=FALSE)
    }
  }
  df
}

# convert index vector to logical
unwhich <- function(indices, length) {
  if(length(length)>1) {
    length=length(length)
  } else if(!is.numeric(length) || length(length)!=1) {
    length=max(indices, na.rm=TRUE)
  }
  sel=logical(length)
  sel[indices] = TRUE
  sel
}


# Find position in a given number of numeric intervals -
# in case of several matches only the first is returned
match.intervals <- function(vector, starts, ends) {
  vector <- as.numeric(vector)
  starts <- as.numeric(starts)
  ends <- as.numeric(ends)
  if(length(starts)!=length(ends))
    stop("Interval starts and ends need to have the same length.")
  ttt=is.na(starts)|is.na(ends)
  if(any(ttt)) {
    starts=starts[!ttt]
    ends=ends[!ttt]    
  }
  ttt=is.na(vector)
  if(any(ttt)) {
    vector=vector[!ttt]
  }
  result=integer(length(vector))
  if(length(starts)>0&length(vector)>0) {
    result <- .C("match_intervals", vector, length(vector),
                 starts, ends, length(starts), result)[[6]]
  }
  if(any(ttt)) {
    res2=integer(length(ttt))
    res2[!ttt]=result
    res2[ttt]=NA
    result=res2
  }
  result
}

# match intervals stratified by chromosomes
match.chr.intervals <- function(chr.vector, pos.vector, chrs, starts, ends) {
  if(length(chrs)!=length(starts) || length(chrs)!=length(ends))
    stop("Lengths of chrs, starts, ends do not match")
  int=data.frame(starts, ends)
  int=split(int, chrs)
  result=character(length(pos.vector))
  result[]=NA
  for(i in 1:length(int)) {
    p=chr.vector==names(int)[i]
    res=match.intervals(pos.vector[p], int[[i]][,1], int[[i]][,2])
    # result[p]=ifelse(res>0,paste(names(int)[i], ":", int[[i]][res,1], "-", int[[i]][res,2], sep=""),NA)
    # update (22 Feb 2012, 05:39pm)
    # acces via res is wrong, since res contains entries with 0
    # due to the 0 entries the start and end points of the matched
    # intervals are shifted and lead to wrong entries in the result vector
    res[res==0] <- NA
    result[p] <- ifelse(!is.na(res), paste(names(int)[i], ":", int[[i]][res,1], "-", int[[i]][res,2], sep=""), NA)
  }
  result
}
