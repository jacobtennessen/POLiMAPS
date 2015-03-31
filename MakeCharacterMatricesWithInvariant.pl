#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Std;

##### ##### ##### ##### #####

use vars qw( $opt_t $opt_v $opt_m $opt_o $opt_e $opt_r $opt_n $opt_x);

my $defaultoutputname = "LG_Char_Matrix_With_Invariant";

my $defaultnonmissingthresh = 300;

# Usage
my $usage = "
MakeCharacterMatricesWithInvariant.pl - finds phylogenetic markers sharing sequence reads with linkage map markers
Copyright (C) 2015 by Jacob A Tennessen

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

Usage: perl MakeCharacterMatricesWithInvariant.pl options
 required:
  -t	a tab-delimited table listing all ingroup input files. Each line must have the following format:
	    ParentName	LinkageTable	PileupFile	SamFile
	    where:
	    ParentName is the name of a parent in the cross
	    LinkageTable is a tab-delimited table of linkage map markers with the following format:
		Chromosome	Site	LinkageGroupName
	    PileupFile is a file in pileup format containing the linkage map marker genotypes
	    SamFile is a file in sam format of the parental sequencing reads
  -v	a vcf file containing outgroup sequences
  -m	a comma-delimited list of reference chromosomes, each of which will underlie a distinct matrix
 optional:
  -o	output file name [default is $defaultoutputname]
  -e	outgroup taxa to exclude (comma delimited list)
  -r	required taxa (comma delimited list)
  -n	require non-missing in reference genome [default = T (true)]
  -x	minimum number of nonmissing sites for a linkage group to be included in the matrix [default = $defaultnonmissingthresh]

";

#############

getopts('t:v:m:o:e:r:n:x:');
die $usage unless ($opt_t);
die $usage unless ($opt_v);
die $usage unless ($opt_m);

my ($filetable, $outgroupvcf, $matrixchroms, $output, $excludelist, $reqlist, $reqref, $nonmissingthresh);

$filetable = $opt_t if $opt_t;
$outgroupvcf = $opt_v if $opt_v;
$matrixchroms = $opt_m if $opt_m;
$output = $opt_o ? $opt_o : $defaultoutputname;
$excludelist = $opt_e ? $opt_e : "NA";
$reqlist = $opt_r ? $opt_r : "NA";
$reqref = $opt_n ? $opt_n : "T";
$nonmissingthresh = $opt_x ? $opt_x : $defaultnonmissingthresh;

my %excludelist;

unless ($excludelist =~ /^NA$/) {
    my @excludelist = split ",", $excludelist;
    foreach my $e (@excludelist) {
	$excludelist{$e} = 1;
    }
}

my %reqlist;

unless ($reqlist =~ /^NA$/) {
    my @reqlist = split ",", $reqlist;
    foreach my $r (@reqlist) {
	$reqlist{$r} = 1;
    }
}

my %linkagegroups;

my %pileups;

my %sams;

open(FILETABLE, $filetable) || die "can't open $filetable\n";

while (<FILETABLE>) {
    my $line = $_;
    $line =~ s/\r|\n//g;
    my @data = split "\t", $line;
    $linkagegroups{$data[0]} = $data[1];
    $pileups{$data[0]} = $data[2];
    $sams{$data[0]} = $data[3];
}

close (FILETABLE);

my @matrixchroms = split ",", $matrixchroms;

foreach my $mc (@matrixchroms) {

    my %HoHinformativesites;
    
    my %reference; #record reference base
    
    my %fulllgs;
    
    foreach my $parent (keys %linkagegroups) {
	
	open(LINK, "$linkagegroups{$parent}") || die "can't open $linkagegroups{$parent}\n";
	
	my %markers;
	
	my %HoHchromsperlg;
	
	while (<LINK>) {
	    my $line = $_;
	    $line =~ s/\r|\n//g;
	    my @data = split "\t", $line;
	    #data must be Chrom\tSite\tLG
	    if ($data[0] =~ /^$mc$/) {
		$markers{$data[1]} = $data[2];
		if (defined $HoHchromsperlg{$data[2]}{$data[0]}) {
		    $HoHchromsperlg{$data[2]}{$data[0]} += 1;
		} else {
		   $HoHchromsperlg{$data[2]}{$data[0]} = 1; 
		}
	    }
	}
	
	close (LINK);
	    
	open(PILE, "$pileups{$parent}") || die "can't open $pileups{$parent}\n";
	
	while (<PILE>) {
	    my $line = $_;
	    $line =~ s/\r|\n//g;
	    my @data = split "\t", $line;
	    unless ($data[0] =~ /^$mc$/) {
		next;
	    }
	    my @seqdata = split "", $data[4];
	    my @seqdata2 = split "", $data[7];
	    push @seqdata, @seqdata2;
	    my %bases;
	    my $refbasecount = 0;
	    my $alts = 0;
	    foreach my $s (@seqdata) {
		if (($s =~ /,/)||($s =~ /\./)) {
		    $refbasecount +=1;
		} elsif (($s =~ /A/gi)||($s =~ /C/gi)||($s =~ /G/gi)||($s =~ /T/gi)) {
		    my $uc = uc $s;
		    if (defined $bases{$uc}) {
			$bases{$uc} +=1;
		    } else {
			$bases{$uc} = 1;
			$alts = 1;
		    }
		}
	    }
	    if (defined $markers{$data[1]}) {
		my $mainbase = $data[2];
		my $altbase;
		my $altbasecount = 0;
		foreach my $b (keys %bases) {
		    if ($bases{$b} > $refbasecount) {
			$mainbase = $b;
			$refbasecount = $bases{$b};
		    }
		    if ($bases{$b} > $altbasecount) {
			$altbase = $b;
			$altbasecount = $bases{$b};
		    }
		}
		$markers{$data[1]} = "$markers{$data[1]}\t$data[2]\t$altbase\t$mainbase";
	    }
	    $reference{$data[1]} = $data[2];
	}
	
	close (PILE);
    
	open(SAM, "$sams{$parent}") || die "can't open $sams{$parent}\n";
	
	while (<SAM>) {
	    my $line = $_;
	    $line =~ s/\r|\n//g;
	    if ($line =~ /^@/) {
		next;
	    }
	    my @data = split "\t", $line;
	    unless ($data[2] =~ /^$mc$/) {
		next;
	    }
	    if (($data[5] =~ /I/)||($data[5] =~ /D/)) {
		next;
	    }
	    if ($data[3] =~ /\d/) {
		my @seq = split "", $data[9];
		my %lgs;
		my %localinformativesites;
		for (my $x = 0; $x < (scalar(@seq)); $x++) {
		    my $pos = $data[3] + $x;
		    if (defined $markers{$pos}) {
			my @markerdata = split "\t", $markers{$pos};
			if ($markerdata[1] =~ /$markerdata[3]/) { #if ref base is main base
			    if ($markerdata[2] =~ /$seq[$x]/) {
				$lgs{$markerdata[0]} = 1;
			    }
			} elsif ($markerdata[2] =~ /$markerdata[3]/) { #if alt is main base
			    if ($markerdata[1] =~ /$seq[$x]/) {
				$lgs{$markerdata[0]} = 1;
			    }
			}
		    } elsif (defined $reference{$pos}) {
			unless ($seq[$x] =~ /N/gi) {
			    $localinformativesites{$pos} = $seq[$x];
			}
		    }
		}
		foreach my $lgk (keys %lgs) {
		    my $full_lg = "$parent"."_$lgk";
		    $fulllgs{$full_lg} = 1;
		    for (my $x = 0; $x < (scalar(@seq)); $x++) {
			my $pos = $data[3] + $x;
			my $chromsite = "$data[2]\t$pos";
			if (defined $localinformativesites{$pos}) {
			    if (defined $HoHinformativesites{$chromsite}{$full_lg}) {
				unless ($HoHinformativesites{$chromsite}{$full_lg} =~ /$localinformativesites{$pos}/) {
				    $HoHinformativesites{$chromsite}{$full_lg} = "$HoHinformativesites{$chromsite}{$full_lg}$localinformativesites{$pos}";
				}
			    } else {
				$HoHinformativesites{$chromsite}{$full_lg} = $localinformativesites{$pos};
			    }
			}
		    }
		}
	    }
	}
	
	close (SAM);
    
    }
    
    my %HoHoutgroup;
    
    my @outgroups;
    
    my @ognumbers;
    
    open(OUTGROUP, $outgroupvcf) || die "can't open $outgroupvcf\n";
    
    while (<OUTGROUP>) {
	my $line = $_;
	$line =~ s/\r|\n//g;
	my @data = split "\t", $line;
	my $limit = scalar(@data) - 1;
	if ($data[0] =~ /CHROM/gi) {
	    my @fullnames = @data[9..$limit];
	    my $namecount = 9;
	    foreach my $fn (@fullnames) {
		my @namedata1 = split /\//, $fn;
		my @namedata2 = split /\./, $namedata1[-1];
		unless (defined $excludelist{$namedata2[0]}) {
		    push @outgroups, $namedata2[0];
		    push @ognumbers, $namecount;
		}
		$namecount +=1;
	    }
	    next;
	} elsif ($line=~ /^#/) {
	    next;
	} elsif ($data[7] =~ /INDEL/) {
	    next;
	}
	unless ($data[0] =~ /^$mc$/) {
	    next;
	}
	my $chromsite = "$data[0]\t$data[1]";
	unless (defined $HoHinformativesites{$chromsite}) {
	    next;
	}
	my @inds = @data[@ognumbers];
	my $indcount = 0;
	foreach my $i (@inds) {
	    my @genodata = split ":", $i;
	    if ($genodata[0] =~ /1\/1/) {
		if (length($data[4]) == 1) {
		    $HoHoutgroup{$outgroups[$indcount]}{$data[1]} = $data[4];
		}
	    } elsif (($genodata[0] =~ /^0$/)&&(defined $genodata[1])&&($genodata[1] > 0)) {
		$HoHoutgroup{$outgroups[$indcount]}{$data[1]} = $data[3];
	    } elsif (($genodata[0] =~ /0\/0/)&&(defined $genodata[2])&&($genodata[2] > 0)) {
		$HoHoutgroup{$outgroups[$indcount]}{$data[1]} = $data[3];
	    }
	    $indcount +=1;
	}
    }
    
    close (OUTGROUP);
    
    my @phylotable;
    
    my @lgs = sort (keys %fulllgs);
    
    my $lglist = join "\t", @lgs;
    
    my $outgrouplist = join "\t", @outgroups;
    
    push @phylotable, "Site\tReference\t$outgrouplist\t$lglist";
    
    my %HoAphylip;
    
    my %iupac = (
	"AC" => "M",
	"AG" => "R",
	"AT" => "W",
	"CA" => "M",
	"CG" => "S",
	"CT" => "Y",
	"GA" => "R",
	"GC" => "S",
	"GT" => "K",
	"TA" => "W",
	"TC" => "Y",
	"TG" => "K",
    );
    
    my @infosites = sort by_site (keys %HoHinformativesites);
    
    my %nonmissing;
    
    foreach my $infosite (@infosites) {
	my $nonmissing = 0;
	my @chromdata = split "\t", $infosite;
	my $flag = 0;
	if ($reference{$chromdata[1]} !~ /N/) {
	    $nonmissing = 1;
	} elsif ($reqref =~ /T/) {
	    $flag = 1;
	}
	foreach my $o (@outgroups) {
	    if (defined $HoHoutgroup{$o}{$chromdata[1]}) {
		if ($HoHoutgroup{$o}{$chromdata[1]} !~ /N/) {
		    $nonmissing +=1;
		}
	    } else {
		$HoHoutgroup{$o}{$chromdata[1]} = "N";
	    }
	    if ((defined $reqlist{$o})&&($HoHoutgroup{$o}{$chromdata[1]} =~ /N/)) {
		$flag = 1;
	    }
	}
	if ($flag == 1) {
	    next;
	}
	foreach my $lg ( keys %{ $HoHinformativesites{$infosite} } ) {
	    my @lgbases = split "", $HoHinformativesites{$infosite}{$lg};
	    foreach my $lb (@lgbases) {
		if ($lb !~ /N/) {
		    $nonmissing +=1;
		}
	    }
	}
	if ($nonmissing > 1) {
	    my @tableline;
	    push @tableline, "$chromdata[0]_$chromdata[1]";
	    push @tableline, $reference{$chromdata[1]};
	    push @{$HoAphylip{"Reference"}}, $reference{$chromdata[1]};
	    $nonmissing{"Reference"} = $nonmissingthresh;
	    foreach my $o (@outgroups) {
		push @tableline, $HoHoutgroup{$o}{$chromdata[1]};
		push @{$HoAphylip{$o}}, $HoHoutgroup{$o}{$chromdata[1]};
		$nonmissing{$o} = $nonmissingthresh;
	    }
	    foreach my $fl (@lgs) {
		my $nuc = "N";
		if (defined $HoHinformativesites{$infosite}{$fl}) {
		    $nuc = $HoHinformativesites{$infosite}{$fl};
		}
		if (length($nuc) > 1) {
		    if (defined $iupac{$nuc}) {
		       $nuc = $iupac{$nuc}; 
		    } else {
		       $nuc = "N"; 
		    }
		}
		unless ($nuc =~ /N/) {
		    if (defined $nonmissing{$fl}) {
			$nonmissing{$fl} +=1;
		    } else {
			$nonmissing{$fl} = 1;
		    }
		}
		push @tableline, $nuc;
		push @{$HoAphylip{$fl}}, $nuc;
	    }
	    my $tl = join "\t", @tableline;
	    push @phylotable, $tl;
	}
    }
    
    my @relaxedphylip;
    
    my $seqlength;
    
    my $taxacount;
    
    foreach my $linkgroup (keys %HoAphylip) {
	if ((defined $nonmissing{$linkgroup})&&($nonmissing{$linkgroup} >= $nonmissingthresh)) {
	    my $seq = join "", @{$HoAphylip{$linkgroup}};
	    unless (defined $seqlength) {
		$seqlength = length($seq);
	    }
	    $taxacount +=1;
	    push @relaxedphylip, "$linkgroup.$nonmissing{$linkgroup} $seq";
	}
    }
    
    unshift @relaxedphylip, "  $taxacount $seqlength";
    
    my $relaxedphylipresult = join "\n", @relaxedphylip;
    
    my $relaxedphylipoutputfile = "$output"."_$mc.phyx";
    
    unless ( open(RPHY, ">$relaxedphylipoutputfile") ) {
	print "Cannot open file \"$relaxedphylipoutputfile\" to write to!!\n\n";
	exit;
    }
    print RPHY $relaxedphylipresult;
    
    close (RPHY);
    
    my $phyloresult = join "\n", @phylotable;
    
    my $phylooutputfile = "$output"."_$mc.txt";
    
    unless ( open(PHY, ">$phylooutputfile") ) {
	print "Cannot open file \"$phylooutputfile\" to write to!!\n\n";
	exit;
    }
    print PHY $phyloresult;
    close (PHY);

}

######################

sub by_site {
    my @aarray = split "\t", $a;
    my @barray = split "\t", $b;
    if ($aarray[0] lt $barray[0]) {-1} elsif ($aarray[0] gt $barray[0]) {1} elsif ($aarray[1] < $barray[1]) {-1} elsif ($aarray[1] > $barray[1]) {1} else {0}
}
