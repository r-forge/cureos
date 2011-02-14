#!/usr/bin/perl

print "VV JASPAR Core, Release 01-March-2006, Download site http://jaspar.cgb.ki.se/DOWNLOAD/MatrixDir\n";
print "XX\n//\n";

@bases = ("A", "C", "G", "T");
open(F, "matrix_list.txt");
while(<F>) {
    @F=split("\t");
    $sysgroupletter = /sysgroup \"(.)/ ? uc($1) : "";
    print "AC  $F[0]\nXX\nID  J$sysgroupletter\$$F[2]\nXX\nNA  $F[2]\nXX\nDE  $F[2]\nXX\n";
    ($seqdb) = /seqdb \"(\w*)\"/;
    ($acc) = /acc \"(\w*)\"/;
    ($species) = /species \"([^\"]*)\"/;
    ($type) = /type \"(\S+)\"/;
    ($medline) = /medline \"(\S+)\"/;
    $ic = $F[1];
    $group = $F[3];
    $acc = "$seqdb:$acc";
    $acc = "" if $acc =~ /:$/;
    $acc = "$acc; $F[2];";
    $acc =~ s/^[;: ]*//;
    print "BF  $acc Species: $species.\nXX\nP0      A      C      G      T\n";
    open(M, "$F[0].pfm");
    $i=0;
    $maxcol=0;
    @matrix=();
    while($l=<M>) {
	$l =~ s/^\s+//;
	$l =~ s/\s+$//;
	@col=split(/\s+/,$l);
	$maxcol=@col if(@col>$maxcol);
	for($j=0;$j<@col;$j++) {
	    $matrix[$j][$i]=$col[$j];
	}
	$i++;
	die("Error in PFM format. $i rows.\n") if $i>4;
    }
    close(M);
    for($i=0;$i<$maxcol;$i++) {
	$max=0;
	for($j=0;$j<4;$j++) {
	    $max = $matrix[$i][$j] if  $matrix[$i][$j] > $max;
	}
	$consens="";
	for($j=0;$j<4;$j++) {
	    $consens .= $bases[$j] if  $matrix[$i][$j] >= $max/2;
	}
	$consens="N" if($consens eq "ACGT");
	$consens="V" if($consens eq "ACG");
	$consens="H" if($consens eq "ACT");
	$consens="D" if($consens eq "AGT");
	$consens="B" if($consens eq "CGT");
	$consens="M" if($consens eq "AC");
	$consens="R" if($consens eq "AG");
	$consens="W" if($consens eq "AT");
	$consens="S" if($consens eq "CG");
	$consens="Y" if($consens eq "CT");
	$consens="K" if($consens eq "GT");
	printf("%02d %6d %6d %6d %6d %6s\n", $i,
	       $matrix[$i][0], $matrix[$i][1], $matrix[$i][2], $matrix[$i][3], $consens);
    }
    print "XX\n";
    if(-r "SITES/$F[0].sites") {
	open(M, "SITES/$F[0].sites");
	$i=-1;
	while($l=<M>) {
	    if($l=~/^>\s*(.*?)\s*$/) {
		$i++;
		$sites[$i] = " $1";
		$sites[$i] =~ s/\s+/; /g;
	    } elsif($l=~/^\s*(..*?)\s*$/) {
		$sites[$i] = "$1$sites[$i]";
	    } else {
		last;
	    }
	}
	close(M);
	$nsites=@sites;
	print "BA  $nsites selected by method \"$type\"\nXX\n";
	foreach $site (@sites) {
	    print "BS  $site.\n";
	}
	print "XX\n";
    }
    print "CC  Total Information Content: $ic.\n";
    print "CC  Factor Class: $group.\nXX\n";
    print "RX  PUBMED: $medline.\nXX\n" if($medline);
    print "//\n";
    
}
close(F);
