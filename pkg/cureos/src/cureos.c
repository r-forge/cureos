#include <R.h>
#include <Rdefines.h>
#include <Rmath.h>

#define SEQ char*
#define BASES "ACGT"
#define AA "ACDEFGHIKLMNPQRSTVWY"

typedef struct {
  int start_pos;
  char strand;
  double score;
  double pvalue;
} PWMsearchResult;

typedef struct {
  int seq;
  int start_pos;
  char strand;
  double score;
  double pvalue;
  double consscore;
  double conspvalue;
} consPWMsearchResult;


/* convert letters in sequence (DNA or Peptide) to numbers 1-N, for letters not specified
   in "letters" put 127, for gap characters '-', '_', put 126 ('~'). If rmgaps not 0
   remove gaps and return modified lengths 
*/
static void seq2num(SEQ seq, int *len, const char *letters, const char nletters, const char rmgaps) {
  int i, lo=0;
  char j, b, c;
  if(nletters>26) error("Error in seq2num. Only letters A-Z supported.\n");
  for(i=0; i<*len; i++) {
    b=seq[i]; c=127;
    if(b=='~' || b=='-' || b=='_') {
      if(rmgaps) continue;
      c='~';
    }else if(b<nletters) {
      continue;
    }else {
      if(b>91) b-=32; /* convert to upper case */;
      for(j=0;j<nletters;j++)
	if(b==letters[j]) {
	  c=j;
	  break;
	}
    }
    seq[lo++]=c;
  }
  *len=lo;
}

/* convert number string into base codes with letters */
static void num2seq(SEQ seq, const int len, const char *letters, const char nletters, const char defaultletter) {
  int i;
  for(i=0; i<len; i++) {
    if(seq[i]=='~') {
      seq[i]='-';
    } else if(seq[i]<nletters) {
      seq[i]=letters[seq[i]];
    } else if(seq[i]<126) {
    } else {
      seq[i]=defaultletter;
    }
  }
}

/* recompute positions in DNA from position vector without gaps putting in gaps */
static inline void adjust_seq_pos(const char *seq_orig, const int length_orig, double *res_mod, int length_mod, const double defval) {
  int i, j=0;
  if(length_mod>length_orig) return;
  for(i=length_orig-1; i>=0,length_mod>0; i--) {
    if(seq_orig[i] == '~' || seq_orig[i] == '-' || seq_orig[i] == '_' || seq_orig[i] == ' ') {
      res_mod[i] = defval;
    } else {
      res_mod[i] = res_mod[--length_mod];
    }
  }
}

/* random number from multinomial distribution */
static inline char random_letter(const char lprob, const double *prob) {
  double r=unif_rand();
  char l=0;
  while(r>prob[l]&&++l<lprob);
  return l;
}

/* distribution of numbers from sequence (converted to numbers before) */
static void get_dist_from_seq(const SEQ seq, const int seq_length, double *dist, const int nletters) {
  int i, len=0;
  for(i=0;i<nletters;i++) dist[i]=0.0;
  for(i=0; i<seq_length; i++) {
    if(seq[i]<nletters) {
      dist[seq[i]]+=1.0;
      len++;
    }
  }
  if(len==0) {
    for(i=0;i<nletters;i++) dist[i]=1/nletters;
  } else {
    for(i=0;i<nletters;i++) dist[i]/=len;
  }
}

/* Create a random permutation of the sequence */
static void permute_seq(const char *seq_in, char *seq_out, const int seq_length, int *tmpvec) {
  int i, j, n=seq_length;
  for (i = 0; i < n; i++)
    tmpvec[i] = i;
  for (i = 0; i < seq_length; i++) {
    j = n * unif_rand();
    seq_out[i] = seq_in[tmpvec[j]];
    tmpvec[j] = tmpvec[--n];
  }
}

/* R interface to permute_seq */
void permute_sequences(char **seq_in, char **seq_out, int *seq_length, int *nseq) {
  int i, *tmpvec;
  GetRNGstate();
  for(i=0; i<nseq[0]; i++) {
    tmpvec = Calloc(seq_length[i], int);
    permute_seq(seq_in[i], seq_out[i], seq_length[i], tmpvec);
    Free(tmpvec);
  }
  PutRNGstate();
}


/* R interface to moving_average */
void moving_average(double *numbers, int *length, int *window_size) {
  double current=0, temp;
  int i, j=0, ws=window_size[0], n=length[0];
  if(ws<=0) return;
  for(i=0; i<ws&&i<n; i++) current += numbers[i];
  temp=current/ws;
  current -= numbers[0];
  numbers[0]=temp;
  for(i=ws; i<n; i++) {
    j++;
    current += numbers[i];
    temp=current/ws;
    current -= numbers[j];
    numbers[j]=temp;
  }
  for(j=j+1; j<n; j++) {
    numbers[j]=0;
  }
}

/* R interface to window_average */
void window_average(double *numbers, double *newnumbers, double *positions, int *length, double *window_size) {
  double current=0, ws=window_size[0];
  int sp=0, ep=0, cp=0, n=length[0], count=0;
  if(ws<=0) return;
  if(n<=0) return;
  for(ep=0; ep<=n; ep++) {
    if(ep==n || (ep>1&&positions[ep]<positions[ep-1])) {
      while(cp<ep) {
	if(count>0) {
	  newnumbers[cp]=current/count;
	} else {
	  newnumbers[cp]=0;
	}
	cp++;
	while(sp<cp&&positions[cp]-positions[sp]>=ws/2) {
	  if(!ISNAN(numbers[sp])) {
	    current -= numbers[sp];
	    count--;
	  }
	  sp++;
	}
      }
      if(ep<n && !ISNAN(numbers[sp])) {
	current=numbers[ep];
	count=1;
      } else {
	current=0;
	count=0;
      }
      sp=ep;
    } else {
      while(cp<ep&&positions[ep]-positions[cp]>ws/2) {
	if(count>0) {
	  newnumbers[cp]=current/count;
	} else {
	  newnumbers[cp]=0;
	}
	cp++;
	while(sp<cp&&positions[cp]-positions[sp]>=ws/2) {
	  if(!ISNAN(numbers[sp])) {
	    current -= numbers[sp];
	    count--;
	  }
	  sp++;
	}
      }
      if(!ISNAN(numbers[ep])) {
	current += numbers[ep];
	count++;
      }
    }
  }
}

/* Create a random sequence from distribution */
static inline void random_seq(SEQ seq_out, const int seq_length, const double *dist, const char nletters) {
  int i;
  for(i=0; i<seq_length; i++)
    seq_out[i] = random_letter(nletters,dist);
}

/* R interface for generating random sequences of distribution dist with letters as in "letters" */
void random_strings(char **strings, int *string_length, int *nseq, char **letters, double *dist, int *nletters) {
  int i, j, b;
  double r;
  char l= (char) nletters[0];
  GetRNGstate();
  for(b=1;b<nletters[0];b++) dist[b]+=dist[b-1];
  for(i=0; i<nseq[0]; i++) {
    strings[i] = R_alloc(string_length[i]+1, sizeof(char));
    for(j=0; j<string_length[i]; j++)
      strings[i][j] = letters[0][random_letter(l,dist)];
    strings[i][j]='\0';
  }
  PutRNGstate();
}

/* compute fractions of conserved letters in sequences excluding gaps (gap character "~" recognized) */
static void alignment_edit_distances(SEQ *seqs, const int len, const int nseq, double *distancematrix) {
  int i,j,s,dist,distlen;
  char *s1, *s2;
  for(i=0;i<nseq;i++) {
    distancematrix[i*nseq+i]=0.0;
    s1=seqs[i];
    for(j=i+1;j<nseq;j++) {
      s2=seqs[j];
      dist=0;
      distlen=0;
      for(s=0;s<len;s++) {
	if(s1[s]=='~' && s2[s]=='~') {
	} else if(s1[s]==s2[s]) {
	  distlen++;
	} else {
	  dist++; distlen++;
	}
      }
      distancematrix[i*nseq+j]=(double)dist/distlen;
      distancematrix[j*nseq+i]=(double)dist/distlen;
    }
  }
}

/* Create a random alignment keeping the structure of conservation and gaps and conserved positions
use base composition from dists[0] (distributions of letters in sequence 1)
*/
static void random_alignment(SEQ *ali_in, SEQ *ali_out, const int seq_length, const int nseq, 
			     double **dists, const char nletters) {
  int s, i, j;
  char b,c,t,n;
  double *d, *dt, *d0, distsum;
  d0 = Calloc(nletters, double);
  dt = Calloc(nletters, double);
  for(s=0; s<nseq; s++) {
    d=dists[s];
    for(n=0;n<nletters;n++) d0[n]=d[n];
    for(n=1;n<nletters;n++) d0[n]+=d0[n-1];
    /* Generate random bases */
    for(i=0; i<seq_length; i++) {
      b=ali_in[s][i];
      if(b=='~') {
	ali_out[s][i] = b;
      } else {
	c=-2;
	for(n=0;n<nletters;n++) dt[n]=d[n];
	for(j=0;j<s;j++) {
	  t=ali_in[j][i];
	  if(t == b) {
	    c=ali_out[j][i];
	    break;
	  } else if(t!='~') {
	    c=-1;
	    dt[ali_out[j][i]] = 0.0;
	  }
	}
	if(c<0) {
	  if(s==0 || c==-2) {
	    c = random_letter(nletters, d0);
	  } else {
	    distsum=0.0;
	    for(n=0;n<nletters;n++) distsum+=dt[n];
	    if(distsum>0) {
	      for(n=0;n<nletters;n++) dt[n]/=distsum;
	    } else {
	      for(n=0;n<nletters;n++) dt[n]=1/nletters;
	    }
	    for(n=1;n<nletters;n++) dt[n]+=dt[n-1];
	    c = random_letter(nletters, dt);
	  }
	}
	ali_out[s][i] = c;
      }
    }
  }
  Free(d0);
  Free(dt);
}

/* R interface for random_alignment */ 
void randomize_alignment(char **ali_in, char **ali_out, int *seq_length, int *nseq, int *nletters) {
  int i, len=seq_length[0];
  double **dists;
  dists=Calloc(nseq[0], double*);
  for(i=0; i<nseq[0]; i++) {
    if(nletters[0]==4) {
      seq2num(ali_in[i],&len,BASES,4,0);
    } else if(nletters[0]==20) {
      seq2num(ali_in[i],&len,AA,20,0);
    } else {
      error("Bad input: Only DNA sequences [ACGT] or Peptide Sequences [ACDEFGHIKLMNPQRSTVWY] are supported.");
    }
    dists[i]=Calloc(nletters[0], double);
    get_dist_from_seq(ali_in[i], len, dists[i], nletters[0]);
  }  
  GetRNGstate();
  random_alignment(ali_in, ali_out, len, nseq[0], dists, (char)nletters[0]);
  PutRNGstate();
  for(i=0; i<nseq[0]; i++) {
    if(nletters[0]==4) {
      num2seq(ali_out[i],len,BASES,4,'N');
    } else {
      num2seq(ali_out[i],len,AA,20,'X');
    }
    Free(dists[i]);
  }
  Free(dists);
}

/* Convert scores to pvalues using either a score and a pvalue vector or a score vector by itself (get fraction of 
   vector larger then score. Uses binary search in score vector.
*/
void scores2pvalues(double *scores, int *len, double *distvec, double *pvec, int *dist_len, int *method, double *result) {
  int i, j, l, m;
  double s;
  if(distvec[dist_len[0]-1]<=distvec[0]) {
    for(i=0; i<*len; i++) result[i]=1.0;
    return;
  }
  for(i=0; i<*len; i++) {
    s = scores[i];
    j = 0; l=0; m=dist_len[0]-1;
    if(distvec[m]<=distvec[l]) {
      error("Bad input: score vector not sorted.");
    } else if(s > distvec[m]) {
      j=m;
    } else if(s<=distvec[l]) {
      j=l;
    } else {
      while(m-l>1) {
	if(distvec[j] >= s) {
	  m=j;
	} else {
	  l=j;
	}
	j = (int)((m-l)/2) + l;
      }
    }
    result[i] = *method==1 ? pvec[j] : 1.0 - ((double) j / dist_len[0]);
  }
}

/* search with a pwm in a sequence and return vector of scores at all positions */
static void pwm_search(SEQ seq, int seq_length, double *pwm, int nrow, int ncol, double *scores) {
  int i,j,k;
  double score;
  for(i=0; i<=seq_length-ncol; i++) {
    score=0.0;
    for(j=i,k=0;j<ncol+i;j++,k+=nrow) {
      if(seq[j] < (char)nrow)
	score += pwm[k+seq[j]];
    }
    scores[i] = score;
  }
  for(; i<seq_length; i++) scores[i] = R_NegInf;
}

/* R interface to pwm_search */
void pwm_scores(char **seq, int *seq_length, double *pwm, int *nrow, int *ncol, double *scores) {
  int i, len=seq_length[0];
  char *sequence=Calloc(seq_length[0], char);
  for(i=0;i<seq_length[0];i++) sequence[i]=seq[0][i];
  if(nrow[0]==4) {
    seq2num(sequence,&len,BASES,4,1);
  } else if(nrow[0]==20) {
    seq2num(sequence,&len,AA,20,1);
  } else {
    error("Bad input: DNA sequences and PWMs with 4 rows or Peptide Sequences and PWMs with 20 rows supported.");
  }
  pwm_search(sequence, len, pwm, nrow[0], ncol[0], scores);
  adjust_seq_pos(seq[0], seq_length[0], scores, len, R_NegInf);
  Free(sequence);
}

/* same as pwm_search but only return maximum score */
static double pwm_search_max(SEQ seq, int seq_length, double *pwm, int nrow, int ncol) {
  int i,j,k;
  double score, max_score=0.0;
  for(i=0; i<=seq_length-ncol; i++) {
    score=0.0;
    for(j=i,k=0; j<ncol+i; j++,k+=nrow) {
      if(seq[j] < (char)nrow)
	score += pwm[k+seq[j]];
    }
    if(score > max_score) 
      max_score=score;
  }
  return max_score;
}

/* Empirical score distribution for PWM with background sequence by permutation */
void pwmscoredist(char **seq, int *seq_length, double *pwm, int *nrow, int *ncol, double *scores, int *nperm) {
  int i, len=seq_length[0], *tmpvec;
  char *sr;
  GetRNGstate();
  if(nrow[0]==4) {
    seq2num(seq[0],&len,BASES,4,1);
  } else if(nrow[0]==20) {
    seq2num(seq[0],&len,AA,20,1);
  } else {
    error("Bad input: DNA sequences and PWMs with 4 rows or Peptide Sequences and PWMs with 20 rows supported.");
  }
  sr = Calloc(len, char);
  tmpvec = Calloc(len, int);
  for(i=0;i<nperm[0];i++) {
    permute_seq(seq[0], sr, len, tmpvec);
    scores[i] = pwm_search_max(sr, len, pwm, nrow[0], ncol[0]);
  }
  R_rsort(scores, nperm[0]);
  PutRNGstate();
  Free(tmpvec);
  Free(sr);
}

/* search with a pwm in alignment and return matrix of scores at all positions and conseved scores*/
static void conserved_pwm_search(SEQ *seq, char *tmpseq, int seq_length, int nseq, double *pwm, int nrow, int ncol, 
				 double *distancematrix, double *scores, double *maxscores) {
  int s,i,j,k,nlen,maxpos,offset;
  double score, maxscore, *tmpscores;
  for(s=0,offset=0; s<nseq; s++,offset+=seq_length) {
    tmpscores = &scores[offset];
    /* remove gaps */
    nlen=0;
    for(i=0; i<seq_length; i++)
      if(seq[s][i]!='~') tmpseq[nlen++] = seq[s][i];
    /* compute pwm scores within sequences and put resulting scores into nseq*seqlength matrix s */
    maxscore=R_NegInf;
    for(i=0; i<=nlen-ncol; i++) {
      score=0.0;
      for(j=i,k=0;j<i+ncol;j++,k+=nrow) {
	if(tmpseq[j] < (char)nrow)
	  score += pwm[k+tmpseq[j]];
      }
      tmpscores[i] = score;
      if(score > maxscore) maxscore=score;
    }
    for(; i<nlen; i++) tmpscores[i] = R_NegInf;
    adjust_seq_pos(seq[s], seq_length, tmpscores, nlen, R_NegInf);
    maxscores[s] = maxscore;
  }
  /* compute conserved scores and put into scores into column nseq+1 in matrix s*/
  maxscore=0.0;
  for(i=0; i<seq_length; i++) {
    score=0.0;
    maxpos=0;
    for(s=i, j=0; j<nseq; s+=seq_length, j++)
      if(scores[s]>score) {
	score=scores[s];
	maxpos=j;
      }
    maxpos *= nseq;
    for(s=i, j=0; j<nseq; s+=seq_length, j++)
      if(scores[s]>0) score+=distancematrix[maxpos+j]*scores[s];
    scores[offset++] = score;
    if(score > maxscore) maxscore=score;
  }
  maxscores[nseq]=maxscore;
}

/* R interface to conserved_pwm_search */
void conserved_pwm_scores(char **seq, int *seq_length, int *nseq, double *pwm, int *nrow, int *ncol, double *scores) {
  int i,len=seq_length[0];
  double *maxscores, *distancematrix;
  char *tmpsequence=Calloc(seq_length[0], char);  
  maxscores=Calloc(nseq[0]+1, double);
  distancematrix=Calloc(nseq[0]*nseq[0], double);
  for(i=0;i<nseq[0];i++) {
    if(nrow[0]==4) {
      seq2num(seq[i],&len,BASES,4,0);
    } else if(nrow[0]==20) {
      seq2num(seq[i],&len,AA,20,0);
    } else {
      error("Bad input: DNA sequences and PWMs with 4 rows or Peptide Sequences and PWMs with 20 rows supported.");
    }
  }
  alignment_edit_distances(seq, len, nseq[0], distancematrix);
  conserved_pwm_search(seq, tmpsequence, len, nseq[0], pwm, nrow[0], ncol[0], distancematrix, scores, maxscores);
  Free(maxscores);
  Free(tmpsequence);
  Free(distancematrix);
}

/* R interface to conserved_pwm_search */
void conserved_pwm_scoredist(char **ali_in, int *seq_length, int *nseq, double *pwm, int *nrow, int *ncol, 
			     int *nsample, double *scoredist) {
  int i,j,k=0,len=seq_length[0];
  double *scores, *maxscores, *distancematrix, **dists;
  char *tmpsequence=Calloc(seq_length[0], char);
  char **tmpalignment=Calloc(nseq[0], char*);
  scores=Calloc((nseq[0]+1)*seq_length[0], double);
  maxscores=Calloc(nseq[0]+1, double);
  distancematrix=Calloc(nseq[0]*nseq[0], double);
  dists=Calloc(nseq[0], double*);
  for(i=0;i<nseq[0];i++) {
    if(nrow[0]==4) {
      seq2num(ali_in[i],&len,BASES,4,0);
    } else if(nrow[0]==20) {
      seq2num(ali_in[i],&len,AA,20,0);
    } else {
      error("Bad input: DNA sequences and PWMs with 4 rows or Peptide Sequences and PWMs with 20 rows supported.");
    }
    tmpalignment[i]=Calloc(seq_length[0], char);
    dists[i]=Calloc(nrow[0], double);
    get_dist_from_seq(ali_in[i], len, dists[i], nrow[0]);
  }
  distancematrix=Calloc(nseq[0]*nseq[0], double);
  alignment_edit_distances(ali_in, len, nseq[0], distancematrix);
  GetRNGstate();
  alignment_edit_distances(ali_in, len, nseq[0], distancematrix);
  for(i=0;i<nsample[0];i++) {
    random_alignment(ali_in, tmpalignment, len, nseq[0], dists, (char)nrow[0]);
    conserved_pwm_search(tmpalignment, tmpsequence, len, nseq[0], pwm, nrow[0], ncol[0], distancematrix, scores, maxscores);
    for(j=0,k=i; j<nseq[0]+1; j++,k+=nsample[0]) scoredist[k]=maxscores[j];
  }
  for(i=0;i<(nseq[0]+1)*nsample[0];i+=nsample[0])
    R_rsort(&scoredist[i], nsample[0]);
  PutRNGstate();
  for(i=0;i<nseq[0];i++) {
    Free(tmpalignment[i]);
    Free(dists[i]);
  }
  Free(tmpalignment);
  Free(dists);
  Free(maxscores);
  Free(scores);
  Free(tmpsequence);
  Free(distancematrix);
}

SEXP pwmscan(SEXP sequences, SEXP pwms, SEXP revpwms, SEXP cutoff, SEXP pmethod) {
  /* Initialize Variables */
  int s, m, i, j, k, p, ncol, nrow, seqlen, modseqlen, pvalmethod, *pvalmethods=NULL, distlen=0, *distlens, nperm=-1, 
    nseq=LENGTH(sequences), npwm=LENGTH(pwms), totalresultlength=0, curresultlen=0, curresultpos=0;
  int *resultlengths, *tmpvec;
  char userevpwms=0, strand='+', ncolout=4, useseqnames=0, usepwmnames=0, npwmrows=0;
  const char *origsequence, *csequence;
  char *letters;
  unsigned char check=9, cutp=0, pvalpresent;
  SEQ sequence;
  SEQ *randomsequences;
  double cut, tmpscore, maxscore;
  double *pwm, *revpwm, *distvec, *pvec=NULL, **distvecs, **pvecs, *scores, *revscores, *pvals;
  SEXP result=R_NilValue, pwmnames=GET_NAMES(pwms), seqnames=GET_NAMES(sequences);
  PWMsearchResult **allresults, *seqresults, tmpresult;
  allresults = Calloc(nseq*npwm, PWMsearchResult*);
  resultlengths = Calloc(nseq*npwm, int);
  GetRNGstate();
  /* Check input */
  if(!IS_CHARACTER(sequences)) error("Wrong input: Sequences needs to be of type character!\n");
  cut = REAL(AS_NUMERIC(cutoff))[0];
  if(!IS_LIST(pwms)) error("Wrong input: PWMS needs to be a list!\n");
  if(isNull(revpwms)) {
    userevpwms=0;
  } else {
    userevpwms=1;
    if(!IS_LIST(pwms)) error("Wrong input: REVCOMP.PWMS needs to be a list or NULL!\n");
    if(npwm != LENGTH(revpwms)) error("Wrong input: %d PWMS and %d REVCOMP.PWMS, numbers do not match.\n", npwm,  LENGTH(revpwms));
  }
  npwmrows=(char)(INTEGER(GET_DIM(VECTOR_ELT(pwms, 0)))[0]);
  if(npwmrows!=4 && npwmrows!=20) error("Wrong input: PWMs have %d rows. Only PWMs with 4 rows for DNA or 20 rows for Peptides are supported!", npwmrows);
  letters = npwmrows==4 ? BASES : AA;
  for(m=0;m<npwm;m++) {
    if(!IS_NUMERIC(VECTOR_ELT(pwms, m)) || LENGTH(GET_DIM(VECTOR_ELT(pwms, m)))!=2)
      error("Wrong input: PWMS in list need to be of type numeric matrix with 2 dimensions!\n");
    if(npwmrows!=INTEGER(GET_DIM(VECTOR_ELT(pwms, m)))[0]) error("Wrong input: All PWM need to have the same number of rows!\n");
    if(userevpwms && (npwmrows!=INTEGER(GET_DIM(VECTOR_ELT(revpwms, m)))[0] ||
		      INTEGER(GET_DIM(VECTOR_ELT(revpwms, m)))[1] != INTEGER(GET_DIM(VECTOR_ELT(pwms, m)))[1]))
      error("Wrong input: All PWMs and REVCOMP.PWMS must have the same size!\n");
  }
  if(!isNull(pmethod)) {
    if(IS_NUMERIC(pmethod) || IS_INTEGER(pmethod) && LENGTH(pmethod)==1) {
      nperm = INTEGER(AS_INTEGER(pmethod))[0];
      distlen = nperm;
      if(nperm<=0) error("Wrong input: Number of permutations if %d, must be a positive integer.\n", nperm);
      distvec = Calloc(nperm, double);
      randomsequences = Calloc(nperm, SEQ);
      pvalmethod=0;
    } else if(IS_LIST(pmethod) && LENGTH(pmethod)==npwm) {
      for(m=0;m<npwm;m++)
	if(!IS_NUMERIC(VECTOR_ELT(pmethod, m))) error("Wrong input: Element %d in list pmethod is not numeric.\n", m);
      pvalmethods=Calloc(npwm,int);
      distvecs=Calloc(npwm, double*);
      distlens=Calloc(npwm, int);
      pvecs=Calloc(npwm, double*);
      for(m=0;m<npwm;m++) {
	if(isNull(GET_NAMES(VECTOR_ELT(pmethod, m)))) {
	  distlens[m] = LENGTH(VECTOR_ELT(pmethod, m));
	  distvecs[m] = REAL(AS_NUMERIC(VECTOR_ELT(pmethod, m)));
	  pvalmethods[m]=0;
	  pvecs[m] = NULL;
	} else {
	  distlens[m] = LENGTH(VECTOR_ELT(pmethod, m));
	  distvecs[m] = REAL(AS_NUMERIC(GET_NAMES(VECTOR_ELT(pmethod, m))));
	  pvecs[m] = REAL(AS_NUMERIC(VECTOR_ELT(pmethod, m)));
	  pvalmethods[m]=1;
	}
      }
    } else {
      error("Wrong input: pmethod must be either single Integer or list with score distributions for each matrix. \n");
    }
  }
  /* Do PWM searches for all sequences */
  for(s=0;s<nseq;s++) { /* start seq loop */
    /* Copy and allocate memory for PWM search in a sequence */
    seqlen=LENGTH(STRING_ELT(sequences, s));
    modseqlen=seqlen;
    origsequence = CHAR(STRING_ELT(sequences, s));
    sequence=Calloc(seqlen, char);
    for(i=0;i<seqlen;i++) sequence[i]=origsequence[i];
    seq2num(sequence,&modseqlen,letters,npwmrows,1);
    if(nperm > 0) {
      tmpvec = Calloc(modseqlen, int);
      for(i=0;i<nperm;i++) {
	randomsequences[i] = Calloc(modseqlen, char);
	permute_seq(sequence, randomsequences[i], modseqlen, tmpvec);
      }
      Free(tmpvec);
    }
    scores = Calloc(seqlen, double);
    if(userevpwms) revscores = Calloc(seqlen, double);
    if(nperm>0||pvalmethods!=NULL) pvals = Calloc(seqlen, double);
    /* Search with all PWMs */
    for(m=0;m<npwm;m++) { /* start pwm loop */
      nrow=INTEGER(GET_DIM(VECTOR_ELT(pwms, m)))[0];
      ncol=INTEGER(GET_DIM(VECTOR_ELT(pwms, m)))[1];
      pwm=REAL(VECTOR_ELT(pwms, m));
      if(modseqlen<ncol) {
	resultlengths[curresultpos]=0;
	allresults[curresultpos++] = NULL;
	continue;
      }
      if(pvalmethods!=NULL) {
	pvalmethod=pvalmethods[m];
	distvec=distvecs[m];
	distlen=distlens[m];
	pvec=pvecs[m];
      }
      for(i=modseqlen-ncol;i<seqlen;i++) scores[i]=R_NegInf;
      if(userevpwms) {
	for(i=modseqlen-ncol;i<seqlen;i++) revscores[i]=R_NegInf;
	revpwm = REAL(VECTOR_ELT(revpwms, m));
      }
      /* Do actual search with PWM */
      maxscore=R_NegInf;
      pwm_search(sequence, modseqlen, pwm, nrow, ncol, scores);
      if(userevpwms) { 
	pwm_search(sequence, modseqlen, revpwm, nrow, ncol, revscores);
	for(i=0; i<modseqlen-ncol+1; i++) {
	  if(scores[i]<revscores[i]) {
	    scores[i]=revscores[i];
	    revscores[i]=1.0;
	  } else {
	    revscores[i]=0.0;
	  }
	}
	adjust_seq_pos(origsequence, seqlen, revscores, modseqlen, 0.0);
      }
      for(i=0; i<modseqlen-ncol+1; i++) if(scores[i]>maxscore) maxscore=scores[i];
      adjust_seq_pos(origsequence, seqlen, scores, modseqlen, R_NegInf);
      /* If no vector for pvalues given but pvalues wanted do permutation test */
      if(nperm > 0) {
	distlen=nperm;
	for(i=0;i<nperm;i++) {
	  distvec[i] = pwm_search_max(randomsequences[i], modseqlen, pwm, nrow, ncol);
	  if(userevpwms) {
	    tmpscore = pwm_search_max(randomsequences[i], modseqlen, revpwm, nrow, ncol);
	    if(tmpscore>distvec[i]) distvec[i] = tmpscore;
	  }
	  /* At some breakpoints check whether random testing can be cut short because p-values are high anyhow */
	  if(i==check) { /* checkpoint */
	    if(i+1==nperm) break;
	    if(check==9) {
	      cutp=0; check=49;
	    } else if(check==49) {
	      cutp=4; check=99;
	    } else if(check==99) {
	      cutp=49; check=199;
	    } else if(check==199) {
	      cutp=179;
	      check=0;
	    }
	    R_rsort(distvec, i+1);
	    if(distvec[cutp]>maxscore) {
	      distlen=i+1;
	      break;
	    }
	  }
	}
	R_rsort(distvec, distlen);
      }
      if(distlen>0)
	scores2pvalues(scores, &seqlen, distvec, pvec, &distlen, &pvalmethod, pvals);
      /* store significant results for output */
      seqresults = Calloc(seqlen, PWMsearchResult);
      curresultlen = 0;
      pvalpresent=distlen>0&&cut<=1&&cut>=0;
      for(i=0;i<seqlen;i++) {
	/* Store results if score is significant  */
	if((pvalpresent&&pvals[i]<=cut)||(!pvalpresent&&scores[i]>=cut)) {
	  seqresults[curresultlen].start_pos=i;
	  seqresults[curresultlen].strand=userevpwms&&revscores[i]>0 ? '-' : '+';
	  seqresults[curresultlen].score=scores[i];
	  seqresults[curresultlen++].pvalue = distlen>0 && pvals[i]>=0 && pvals[i]<=1 ? pvals[i] : 1.0;
	}
      }
      if(curresultlen>0) {
	seqresults = Realloc(seqresults, curresultlen, PWMsearchResult);	
      } else {
	Free(seqresults);
	seqresults=NULL;
      }
      resultlengths[curresultpos]=curresultlen;
      totalresultlength += curresultlen;
      allresults[curresultpos++] = seqresults;
    }
    if(nperm>0) {
      for(i=0;i<nperm;i++)
	Free(randomsequences[i]);
    }
    if(userevpwms) Free(revscores);
    if(nperm>0||pvalmethods!=NULL) Free(pvals);
    Free(sequence);
    Free(scores);
  }
  /* Put Output into R format */
  ncolout=5;
  if(distlen > 0) ncolout++;
  if(userevpwms) ncolout++;
  PROTECT(result = NEW_LIST(ncolout));
  SET_NAMES(result, NEW_CHARACTER(ncolout));
  if(isNull(seqnames)) {
    SET_ELEMENT(result, 0, NEW_INTEGER(totalresultlength));
  } else {
    SET_ELEMENT(result, 0, NEW_CHARACTER(totalresultlength));
    useseqnames=1;
  }
  SET_STRING_ELT(GET_NAMES(result),0,mkChar("Sequence.Name"));
  if(isNull(pwmnames)) {
    SET_ELEMENT(result, 1, NEW_INTEGER(totalresultlength));
  } else {
    SET_ELEMENT(result, 1, NEW_CHARACTER(totalresultlength));
    usepwmnames=1;
  }
  SET_STRING_ELT(GET_NAMES(result),1,mkChar("PWM.Name"));
  SET_ELEMENT(result, 2, NEW_INTEGER(totalresultlength));
  SET_STRING_ELT(GET_NAMES(result),2,mkChar("Match.StartPos"));
  SET_ELEMENT(result, 3, NEW_INTEGER(totalresultlength));
  SET_STRING_ELT(GET_NAMES(result),3,mkChar("Match.EndPos"));
  SET_ELEMENT(result, 4, NEW_NUMERIC(totalresultlength));
  SET_STRING_ELT(GET_NAMES(result),4,mkChar("PWM.Score"));
  ncolout=5;
  if(distlen>0) {
    SET_ELEMENT(result, ncolout, NEW_NUMERIC(totalresultlength));
    SET_STRING_ELT(GET_NAMES(result),ncolout++,mkChar("PWM.Pvalue"));
  }
  if(userevpwms) {
    SET_ELEMENT(result, ncolout, NEW_CHARACTER(totalresultlength));
    SET_STRING_ELT(GET_NAMES(result),ncolout++,mkChar("Match.Ori"));
  }
  i=0;
  k=0;
  for(s=0;s<nseq;s++) { /* start seq loop */ 
    csequence = CHAR(STRING_ELT(sequences, s));
    seqlen=LENGTH(STRING_ELT(sequences, s));
    for(m=0;m<npwm;m++) { /* start pwm loop */
      ncol=INTEGER(GET_DIM(VECTOR_ELT(pwms, m)))[1];
      for(j=0;j<resultlengths[k];j++) {
	tmpresult = allresults[k][j];
	if(useseqnames) {
	  SET_STRING_ELT(VECTOR_ELT(result,0),i,duplicate(STRING_ELT(seqnames, s)));
	} else {
	  INTEGER(VECTOR_ELT(result,0))[i] = s+1;
	}
	if(usepwmnames) {
	  SET_STRING_ELT(VECTOR_ELT(result,1),i,duplicate(STRING_ELT(pwmnames, m)));
	} else {
	  INTEGER(VECTOR_ELT(result,1))[i] = m+1;
	}
	INTEGER(VECTOR_ELT(result,2))[i] = tmpresult.start_pos+1;
	for(p=tmpresult.start_pos,modseqlen=0;p<seqlen&&modseqlen<ncol;p++)
	  if(csequence[p]!='~'&&csequence[p]!='-'&&csequence[p]!='_'&&csequence[p]!=' ') modseqlen++;
	INTEGER(VECTOR_ELT(result,3))[i] = p;
	REAL(VECTOR_ELT(result,4))[i] = tmpresult.score;
	ncolout=5;
	REAL(VECTOR_ELT(result,ncolout++))[i] = tmpresult.pvalue;
	if(userevpwms) {
	  if(tmpresult.strand == '+') {
	    SET_STRING_ELT(VECTOR_ELT(result,ncolout++),i,mkChar("+")); 
	  } else if(tmpresult.strand == '-') {
	    SET_STRING_ELT(VECTOR_ELT(result,ncolout++),i,mkChar("-"));
	  } else {
	    SET_STRING_ELT(VECTOR_ELT(result,ncolout++),i,mkChar(""));
	  }
	}
	i++;
      }
      Free(allresults[k]);
      k++;
    }
  }
  PutRNGstate();
  if(nperm>0) Free(distvec);
  if(nperm>0) Free(randomsequences);
  if(pvalmethods!=NULL) {
    Free(pvalmethods);
    Free(distvecs);
    Free(distlens);
    Free(pvecs);
  }
  Free(allresults);
  Free(resultlengths);
  UNPROTECT(1);
  return(result);
}

SEXP conserved_pwmscan(SEXP alignment, SEXP pwms, SEXP revpwms, SEXP cutoff, SEXP randomtests) {
  int s, r, m, i, j, k, ncolout, nseq=LENGTH(alignment), npwm=LENGTH(pwms), seqlen, distlen=0,
    ntests, allsmall, curresultlen, allresultlen=0, pvalmethod=0, nrow, ncol, *resultlengths;
  char userevpwms=0, npwmrows=0, useseqnames=0, usepwmnames=0;
  unsigned char check=9, cutp=0;
  const char *csequence;
  char *letters;
  double *pwm, *revpwm, *distancematrix, *maxscores, *scores, *revscores, *revmaxscores, 
    *randscores, *randmaxscores, *randrevmaxscores, *randscoredists, **letterdists, *pvals, cut, 
    consscore, conspval;
  SEQ tmpsequence;
  SEQ *originalali;
  SEQ **randomalis;
  SEXP result=R_NilValue, pwmnames=GET_NAMES(pwms), seqnames=GET_NAMES(alignment);
  consPWMsearchResult **allresults, *pwmresults, tmpresult;
  
  /* Check input */
  if(!IS_CHARACTER(alignment)) error("Wrong input: Alignment needs to be of type character!\n");
  seqlen=LENGTH(STRING_ELT(alignment, 0));
  for(s=1;s<nseq;s++)
    if(LENGTH(STRING_ELT(alignment, 0))!=seqlen) error("Wrong input: All sequences in alignment need to be of same length!\n");
  cut = REAL(AS_NUMERIC(cutoff))[0];
  ntests = INTEGER(AS_INTEGER(randomtests))[0];
  if(ntests<=0) error("Wrong input: Number of permutations if %d, must be a positive integer.\n", ntests);
  if(!IS_LIST(pwms)) error("Wrong input: PWMS needs to be a list!\n");
  if(!isNull(revpwms)) {
    userevpwms=1;
    if(!IS_LIST(pwms)) error("Wrong input: REVCOMP.PWMS needs to be a list or NULL!\n");
    if(npwm != LENGTH(revpwms)) error("Wrong input: %d PWMS and %d REVCOMP.PWMS, numbers do not match.\n", npwm,  LENGTH(revpwms));
  }
  npwmrows=(char)(INTEGER(GET_DIM(VECTOR_ELT(pwms, 0)))[0]);
  if(npwmrows!=4 && npwmrows!=20) 
    error("Wrong input: PWMs have %d rows. Only PWMs with 4 rows for DNA or 20 rows for Peptides are supported!", npwmrows);
  letters = npwmrows==4 ? BASES : AA;
  for(m=0;m<npwm;m++) {
    if(!IS_NUMERIC(VECTOR_ELT(pwms, m)) || LENGTH(GET_DIM(VECTOR_ELT(pwms, m)))!=2)
      error("Wrong input: PWMS in list need to be of type numeric matrix with 2 dimensions!\n");
    if(npwmrows!=INTEGER(GET_DIM(VECTOR_ELT(pwms, m)))[0]) error("Wrong input: All PWM need to have the same number of rows!\n");
    if(userevpwms && (npwmrows!=INTEGER(GET_DIM(VECTOR_ELT(revpwms, m)))[0] ||
		      INTEGER(GET_DIM(VECTOR_ELT(revpwms, m)))[1] != INTEGER(GET_DIM(VECTOR_ELT(pwms, m)))[1]))
      error("Wrong input: All PWMs and REVCOMP.PWMS must have the same size!\n");
  }
  /* Allocate/Initialize Memory */
  originalali=Calloc(nseq, SEQ);
  tmpsequence=Calloc(seqlen, char);
  maxscores=Calloc(nseq+1, double);
  scores=Calloc((nseq+1)*seqlen, double);
  resultlengths=Calloc(npwm, int);
  if(userevpwms) { 
    revscores=Calloc((nseq+1)*seqlen, double);
    revmaxscores=Calloc(nseq+1, double);
  }
  distancematrix=Calloc(nseq*nseq, double);
  letterdists=Calloc(nseq, double*);
  allresults=Calloc(npwm, consPWMsearchResult*);
  for(s=0; s<nseq; s++) {
    originalali[s]=Calloc(seqlen, char);
    for(i=0;i<seqlen;i++) originalali[s][i]=CHAR(STRING_ELT(alignment, s))[i];
    seq2num(originalali[s],&seqlen,letters,npwmrows,0);
    letterdists[s]=Calloc(npwmrows, double);
    get_dist_from_seq(originalali[s], seqlen, letterdists[s], npwmrows);
  }
  if(ntests>0) {
    pvals=Calloc((nseq+1)*seqlen, double);
    randscores=Calloc((nseq+1)*seqlen, double);
    randmaxscores=Calloc(nseq+1, double);
    if(userevpwms) randrevmaxscores=Calloc(nseq+1, double);
    randscoredists = Calloc((nseq+1)*ntests, double);
    randomalis = Calloc(ntests, SEQ*);
    GetRNGstate();
    for(r=0;r<ntests;r++) {
      randomalis[r]=Calloc(nseq, SEQ);
      for(s=0;s<nseq;s++) randomalis[r][s]=Calloc(seqlen, char);
      random_alignment(originalali, randomalis[r], seqlen, nseq, letterdists, npwmrows);
    }
    PutRNGstate();  
  }
  /* Cycle through PWMs */
  for(m=0; m<npwm; m++) {
    nrow=INTEGER(GET_DIM(VECTOR_ELT(pwms, m)))[0];
    ncol=INTEGER(GET_DIM(VECTOR_ELT(pwms, m)))[1];
    pwm=REAL(VECTOR_ELT(pwms, m));
    if(seqlen<ncol) {
      allresults[m] = NULL;
      resultlengths[m]=0;
      continue;
    }
    if(userevpwms) revpwm = REAL(VECTOR_ELT(revpwms, m));
    /* Do actual search in alignment */
    alignment_edit_distances(originalali, seqlen, nseq, distancematrix);
    conserved_pwm_search(originalali, tmpsequence, seqlen, nseq, pwm, nrow, ncol, distancematrix, scores, maxscores);
    if(userevpwms) {
      conserved_pwm_search(originalali, tmpsequence, seqlen, nseq, revpwm, nrow, ncol, distancematrix, revscores, revmaxscores);
      for(i=0; i<seqlen*(nseq+1); i++) {
	if(scores[i]<revscores[i]) {
	  scores[i]=revscores[i];
	  revscores[i]=1.0;
	} else {
	  revscores[i]=0.0;
	}
      }
      for(i=0; i<nseq+1; i++)
	if(maxscores[i]<revmaxscores[i]) maxscores[i]=revmaxscores[i];
    }
    /* Test scoredists with random alignments */
    distlen=ntests;
    for(r=0;r<ntests;r++) {
      conserved_pwm_search(randomalis[r], tmpsequence, seqlen, nseq, pwm, nrow, ncol, distancematrix, randscores, randmaxscores);
      if(userevpwms) {
	conserved_pwm_search(randomalis[r], tmpsequence, seqlen, nseq, revpwm, nrow, ncol, distancematrix, randscores, randrevmaxscores);
	for(i=0; i<nseq+1; i++)
	  if(randmaxscores[i]<randrevmaxscores[i]) randmaxscores[i]=randrevmaxscores[i];
      }
      for(j=0,k=r; j<nseq+1; j++,k+=ntests) randscoredists[k]=randmaxscores[j];
      /* At some breakpoints check whether random testing can be cut short because p-values are high anyhow */
      if(r==check) { /* checkpoint */
	if(r+1==ntests) break;
	if(check==9) {
	  cutp=0; check=49;
	} else if(check==49) {
	  cutp=4; check=99;
	} else if(check==99) {
	  cutp=49; check=199;
	} else if(check==199) {
	  cutp=179;
	  check=0;
	}
	allsmall=0;
	for(j=0,k=0;j<nseq+1;j++,k+=ntests) {
	  R_rsort(&randscoredists[k], r+1);
	  if(randscoredists[k+cutp]>maxscores[j]) allsmall++;
	}
	if(allsmall==nseq+1) {
	  distlen=r+1;
	  break;
	}
      }
    }
    /* store significant results for output */
    pwmresults=Calloc(seqlen*nseq, consPWMsearchResult);
    if(ntests>0) {
      pvalmethod=0;
      for(j=0,k=0;j<nseq+1;j++,k+=ntests)
	R_rsort(&randscoredists[k], distlen);
      for(i=0,j=0,k=0;i<nseq+1;i++,j+=seqlen,k+=ntests)
	scores2pvalues(&scores[j], &seqlen, &randscoredists[k], NULL, &distlen, &pvalmethod, &pvals[j]);      
    }
    pvalmethod=distlen>0&&cut<=1&&cut>=0;
    curresultlen=0;
    for(i=0,r=nseq*seqlen;i<seqlen;i++,r++) {
      allsmall = !((pvalmethod&&pvals[r]<=cut)||(!pvalmethod&&scores[r]>=cut));
      consscore = scores[r];
      conspval = distlen<=0 || allsmall || pvals[r]<0 || pvals[r]>1 ? 1.0 : pvals[r];
      for(j=0,k=i;j<nseq;j++,k+=seqlen) {
	/* Store results if consscore is significant and score greater 0, or if score itself is significant */
	if((!allsmall&&scores[k]>0)||(pvalmethod&&pvals[k]<=cut)||(!pvalmethod&&scores[k]>=cut)) {
	  tmpresult.seq = j;
	  tmpresult.start_pos = i;
	  tmpresult.strand = userevpwms&&revscores[k]>0 ? '-' : '+';
	  tmpresult.score = scores[k];
	  tmpresult.pvalue = distlen > 0 && pvals[k]>=0 && pvals[k]<=1 ? pvals[k] : 1.0;
	  tmpresult.consscore = consscore;
	  tmpresult.conspvalue = conspval;
	  pwmresults[curresultlen++] = tmpresult;
	}
      }
    }
    if(curresultlen>0) {
      allresults[m] = Realloc(pwmresults, curresultlen, consPWMsearchResult);
      allresultlen+=curresultlen;
    } else {
      allresults[m] = NULL;
      Free(pwmresults);
    }
    resultlengths[m]=curresultlen;
  }
  /* Free some stuff */
  for(s=0;s<nseq;s++) {
    Free(originalali[s]);
    Free(letterdists[s]);
  }
  Free(originalali);
  Free(tmpsequence);
  Free(maxscores);
  Free(scores);
  Free(distancematrix);
  Free(letterdists);
  if(userevpwms) {
    Free(revscores);
    Free(revmaxscores);
  }
  if(ntests>0) {
    Free(pvals);
    Free(randscores);
    Free(randmaxscores);
    if(userevpwms) { Free(randrevmaxscores); }
    Free(randscoredists);
    for(r=0;r<ntests;r++) {
      for(s=0;s<nseq;s++) Free(randomalis[r][s]);
      Free(randomalis[r]);
    }
    Free(randomalis);
  }
  /* Put Output into R format */
  ncolout=6;
  if(ntests > 0) ncolout+=2;
  if(userevpwms) ncolout++;
  PROTECT(result = NEW_LIST(ncolout));
  SET_NAMES(result, NEW_CHARACTER(ncolout));
  if(isNull(seqnames)) {
    SET_ELEMENT(result, 0, NEW_INTEGER(allresultlen));
  } else {
    SET_ELEMENT(result, 0, NEW_CHARACTER(allresultlen));
    useseqnames=1;
  }
  SET_STRING_ELT(GET_NAMES(result),0,mkChar("Sequence.Name"));
  if(isNull(pwmnames)) {
    SET_ELEMENT(result, 1, NEW_INTEGER(allresultlen));
  } else {
    SET_ELEMENT(result, 1, NEW_CHARACTER(allresultlen));
    usepwmnames=1;
  }
  SET_STRING_ELT(GET_NAMES(result),1,mkChar("PWM.Name"));
  SET_ELEMENT(result, 2, NEW_INTEGER(allresultlen));
  SET_STRING_ELT(GET_NAMES(result),2,mkChar("Match.StartPos"));
  SET_ELEMENT(result, 3, NEW_INTEGER(allresultlen));
  SET_STRING_ELT(GET_NAMES(result),3,mkChar("Match.EndPos"));
  SET_ELEMENT(result, 4, NEW_NUMERIC(allresultlen));
  SET_STRING_ELT(GET_NAMES(result),4,mkChar("PWM.Score"));
  SET_ELEMENT(result, 5, NEW_NUMERIC(allresultlen));
  SET_STRING_ELT(GET_NAMES(result),5,mkChar("Cons.Score"));
  ncolout=6;
  if(ntests>0) {
    SET_ELEMENT(result, ncolout, NEW_NUMERIC(allresultlen));
    SET_STRING_ELT(GET_NAMES(result),ncolout++,mkChar("PWM.Pvalue"));
    SET_ELEMENT(result, ncolout, NEW_NUMERIC(allresultlen));
    SET_STRING_ELT(GET_NAMES(result),ncolout++,mkChar("Cons.Pvalue"));
  }
  if(userevpwms) {
    SET_ELEMENT(result, ncolout, NEW_CHARACTER(allresultlen));
    SET_STRING_ELT(GET_NAMES(result),ncolout++,mkChar("Match.Ori"));
  }
  i=0;
  k=0;
  for(m=0;m<npwm;m++) { /* start pwm loop */
    ncol=INTEGER(GET_DIM(VECTOR_ELT(pwms, m)))[1];
    for(j=0;j<resultlengths[k];j++) {
      tmpresult = allresults[k][j];
      s=tmpresult.seq;
      csequence = CHAR(STRING_ELT(alignment, s));
      seqlen=LENGTH(STRING_ELT(alignment, s));
      if(useseqnames) {
	SET_STRING_ELT(VECTOR_ELT(result,0),i,duplicate(STRING_ELT(seqnames, s)));
      } else {
	INTEGER(VECTOR_ELT(result,0))[i] = s+1;
      }
      if(usepwmnames) {
	SET_STRING_ELT(VECTOR_ELT(result,1),i,duplicate(STRING_ELT(pwmnames, m)));
      } else {
	INTEGER(VECTOR_ELT(result,1))[i] = m+1;
      }
      INTEGER(VECTOR_ELT(result,2))[i] = tmpresult.start_pos+1;
      for(r=tmpresult.start_pos,s=0;r<seqlen&&s<ncol;r++)
	if(csequence[r]!='~'&&csequence[r]!='-'&&csequence[r]!='_'&&csequence[r]!=' ') s++;
      INTEGER(VECTOR_ELT(result,3))[i] = r;
      REAL(VECTOR_ELT(result,4))[i] = tmpresult.score;
      REAL(VECTOR_ELT(result,5))[i] = tmpresult.consscore;
      ncolout=6;
      REAL(VECTOR_ELT(result,ncolout++))[i] = tmpresult.pvalue;
      REAL(VECTOR_ELT(result,ncolout++))[i] = tmpresult.conspvalue;
      if(userevpwms) {
	if(tmpresult.strand == '+') {
	  SET_STRING_ELT(VECTOR_ELT(result,ncolout++),i,mkChar("+")); 
	} else if(tmpresult.strand == '-') {
	  SET_STRING_ELT(VECTOR_ELT(result,ncolout++),i,mkChar("-"));
	} else {
	  SET_STRING_ELT(VECTOR_ELT(result,ncolout++),i,mkChar(""));
	}
      }
      i++;
    }
    Free(allresults[k]);
    k++;
  }
  Free(allresults);
  Free(resultlengths);
  UNPROTECT(1);
  return(result);
}

void pwmexactscoredist(double *S, double *P, int *nrows, int *ncols, double *eps, double *Rdist, int *length) {
  /* Compute distribution of summed scores from @S when
     position specific score distribution is @P.
     Requires that scores in @S are rounded to granularity eps.
     Returns (distribution vector)
  */
  double colp, minscore=0, maxscore=0, *Rodist, mini, maxi, dz;
  int len=1, nlen, j, z, colstart=0, *offset;
  Rodist = Calloc(length[0], double);
  offset = Calloc(nrows[0], int);
  for(j=0;j<length[0];j++)
    Rdist[j]=0;
  Rdist[0]=1;
  for(colstart=0; colstart<nrows[0]*ncols[0]; colstart+=nrows[0]) { 
    /* Convolute colums into dist via reference Rdist */
    mini=S[colstart]; 
    maxi=S[colstart];
    for(j=0; j<nrows[0]; j++) { 
      if(S[colstart+j]>maxi) { 
	maxi=S[colstart+j];
      } else if(S[colstart+j]<mini) {
	mini=S[colstart+j];
      }
    }
    nlen = 1 + (int)((((maxscore+maxi)-(minscore+mini))/eps[0]) + 0.5);
    if(nlen>length[0]) error("Vector for output not big enough. Require>%d\n", nlen);
    for(j=0;j<len;j++)
      Rodist[j]=Rdist[j]; 
    for(j=0;j<nlen;j++)
      Rdist[j]=0;
    for(j=0; j<nrows[0]; j++)
      offset[j] = (int)(((S[colstart+j]-mini)/eps[0])+0.5);
    for(z=0; z<len; z++) {
      dz=Rodist[z];
      if(dz!=0) {
	for(j=0; j<nrows[0]; j++) {
	  if(z+offset[j]>length[0]) error("Vector for output not big enough. Require >%d\n", z+offset[j]);
	  Rdist[z+offset[j]] += dz*P[colstart+j]; 
	}
      }
    }
    len=nlen; 
    maxscore+=maxi;
    minscore+=mini;
  }
  Free(Rodist);
  Free(offset);
}

/* reverse a string */
SEXP revstrings(SEXP x) {
  SEXP r;
  char *rev;
  int i, j, k, c, n;

  if( !isString(x) )
    error("argument must be a character vector");

  n = length(x);

  PROTECT(r = allocVector(STRSXP, n) );
  for(k=0; k<n; k++ ) {
    rev = strdup(CHAR(STRING_ELT(x, k)));
    for(i = 0, j=strlen(rev)-1; i<j; i++, j--) {
        c = rev[i];
        rev[i] = rev[j];
        rev[j] = c;
    }
    SET_STRING_ELT(r, k, mkChar(rev));
  }
  UNPROTECT(1);
  return(r);
}

/* R interface to match.intervals - finds first interval*/
void match_intervals(double *vec, int *vec_length, double *starts, double *ends, int *nint, int *icol) {
  int i,j;
  for(i=0;i<vec_length[0];i++) {
    icol[i]=0;
    for(j=0;j<nint[0];j++) {
      if(vec[i]>=starts[j]&&vec[i]<=ends[j]) {
	icol[i]=j+1;
	break;
      }
    }
  }
}
