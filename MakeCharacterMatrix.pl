#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Std;

##### ##### ##### ##### #####

use vars qw( $opt_t $opt_v $opt_o);

# Usage
my $usage = "
MakeCharacterMatrix.pl - finds phylogenetic markers sharing sequence reads with linkage map markers
Copyright (C) 2014 by Jacob A Tennessen

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

Usage: perl MakeCharacterMatrix.pl options
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
 optional:
  -o	output file name [default is LG_Char_Matrix]

";

#############

getopts('t:v:o:');
die $usage unless ($opt_t);
die $usage unless ($opt_v);

my ($filetable, $outgroupvcf, $output);

$filetable = $opt_t if $opt_t;
$outgroupvcf = $opt_v if $opt_v;
$output = $opt_o ? $opt_o : "LG_Char_Matrix";

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

my %HoHinformativesites;

my %HoHreference; #record reference base

my %fulllgs;

foreach my $parent (keys %linkagegroups) {
    
    open(LINK, "$linkagegroups{$parent}") || die "can't open $linkagegroups{$parent}\n";
    
    my %HoHmarkers;
    
    my %HoHchromsperlg;
    
    while (<LINK>) {
	my $line = $_;
	$line =~ s/\r|\n//g;
	my @data = split "\t", $line;
	#data must be Chrom\tSite\tLG
	$HoHmarkers{$data[0]}{$data[1]} = $data[2];
	if (defined $HoHchromsperlg{$data[2]}{$data[0]}) {
	    $HoHchromsperlg{$data[2]}{$data[0]} += 1;
	} else {
	   $HoHchromsperlg{$data[2]}{$data[0]} = 1; 
	}
    }
    
    close (LINK);

    my %mainchroms;
    
    foreach my $lg (keys %HoHchromsperlg) {
	my $biggestchrom;
	my $biggestchromcount = 0;
	foreach my $chrom ( keys %{ $HoHchromsperlg{$lg} } ) {
	    if ($HoHchromsperlg{$lg}{$chrom} > $biggestchromcount) {
		$biggestchrom = $chrom;
		$biggestchromcount = $HoHchromsperlg{$lg}{$chrom};
	    }
	}
	$mainchroms{$lg} = $biggestchrom;
    }
        
    open(PILE, "$pileups{$parent}") || die "can't open $pileups{$parent}\n";
    
    while (<PILE>) {
	my $line = $_;
	$line =~ s/\r|\n//g;
	my @data = split "\t", $line;
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
	if (defined $HoHmarkers{$data[0]}{$data[1]}) {
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
	    $HoHmarkers{$data[0]}{$data[1]} = "$HoHmarkers{$data[0]}{$data[1]}\t$data[2]\t$altbase\t$mainbase";
	}
	$HoHreference{$data[0]}{$data[1]} = $data[2];
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
	if (($data[5] =~ /I/)||($data[5] =~ /D/)) {
	    next;
	}
	if ($data[3] =~ /\d/) {
	    my @seq = split "", $data[9];
	    my %lgs;
	    my %localinformativesites;
	    for (my $x = 0; $x < (scalar(@seq)); $x++) {
		my $pos = $data[3] + $x;
		if (defined $HoHmarkers{$data[2]}{$pos}) {
		    my @markerdata = split "\t", $HoHmarkers{$data[2]}{$pos};
		    if ($markerdata[1] =~ /$markerdata[3]/) { #if ref base is main base
			if ($markerdata[2] =~ /$seq[$x]/) {
			    $lgs{$markerdata[0]} = 1;
			}
		    } elsif ($markerdata[2] =~ /$markerdata[3]/) { #if alt is main base
			if ($markerdata[1] =~ /$seq[$x]/) {
			    $lgs{$markerdata[0]} = 1;
			}
		    }
		} elsif (defined $HoHreference{$data[2]}{$pos}) {
		    unless ($seq[$x] =~ /N/gi) {
			$localinformativesites{$pos} = $seq[$x];
		    }
		}
	    }
	    foreach my $lgk (keys %lgs) {
		my $full_lg = "$parent"."_$lgk"."_$mainchroms{$lgk}";
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

my %HoHoHoutgroup;

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
	    push @outgroups, $namedata2[0];
	    push @ognumbers, $namecount;
	    $namecount +=1;
	}
	next;
    } elsif ($line=~ /^#/) {
	next;
    } elsif ($data[4] =~ /\./) {
	next;
    } elsif ($data[7] =~ /INDEL/) {
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
		$HoHoHoutgroup{$outgroups[$indcount]}{$data[0]}{$data[1]} = $data[4];
	    }
	} elsif (($genodata[0] =~ /^0$/)&&(defined $genodata[1])&&($genodata[1] > 0)) {
	    $HoHoHoutgroup{$outgroups[$indcount]}{$data[0]}{$data[1]} = $data[3];
	} elsif (($genodata[0] =~ /0\/0/)&&(defined $genodata[2])&&($genodata[2] > 0)) {
	    $HoHoHoutgroup{$outgroups[$indcount]}{$data[0]}{$data[1]} = $data[3];
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

foreach my $infosite (@infosites) {
    my %bases;
    my %doublebases;
    my @chromdata = split "\t", $infosite;
    if ($HoHreference{$chromdata[0]}{$chromdata[1]} !~ /N/) {
	$bases{$HoHreference{$chromdata[0]}{$chromdata[1]}} = 1;
    }
    foreach my $o (@outgroups) {
	if (defined $HoHoHoutgroup{$o}{$chromdata[0]}{$chromdata[1]}) {
	    if ($HoHoHoutgroup{$o}{$chromdata[0]}{$chromdata[1]} !~ /N/) {
		if (defined $bases{$HoHoHoutgroup{$o}{$chromdata[0]}{$chromdata[1]}}) {
		    $doublebases{$HoHoHoutgroup{$o}{$chromdata[0]}{$chromdata[1]}} = 1;
		} else {
		    $bases{$HoHoHoutgroup{$o}{$chromdata[0]}{$chromdata[1]}} = 1;
		}
	    }
	} else {
	    $HoHoHoutgroup{$o}{$chromdata[0]}{$chromdata[1]} = "N";
	}
    }
    foreach my $lg ( keys %{ $HoHinformativesites{$infosite} } ) {
	my @lgbases = split "", $HoHinformativesites{$infosite}{$lg};
	foreach my $lb (@lgbases) {
	    if ($lb !~ /N/) {
		if (defined $bases{$lb}) {
		    $doublebases{$lb} = 1;
		} else {
		    $bases{$lb} = 1;
		}
	    }
	}
    }
    if (((scalar (keys %doublebases)) >= 1)&&((scalar (keys %bases)) > 1)) {
	my @tableline;
	push @tableline, "$chromdata[0]_$chromdata[1]";
	push @tableline, $HoHreference{$chromdata[0]}{$chromdata[1]};
	push @{$HoAphylip{"Reference"}}, $HoHreference{$chromdata[0]}{$chromdata[1]};
	foreach my $o (@outgroups) {
	    push @tableline, $HoHoHoutgroup{$o}{$chromdata[0]}{$chromdata[1]};
	    push @{$HoAphylip{$o}}, $HoHoHoutgroup{$o}{$chromdata[0]}{$chromdata[1]};
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
    my $seq = join "", @{$HoAphylip{$linkgroup}};
    unless (defined $seqlength) {
	$seqlength = length($seq);
    }
    $taxacount +=1;
    push @relaxedphylip, "$linkgroup $seq";
}

unshift @relaxedphylip, "  $taxacount $seqlength";

my $relaxedphylipresult = join "\n", @relaxedphylip;

my $relaxedphylipoutputfile = "$output.phyx";

unless ( open(RPHY, ">$relaxedphylipoutputfile") ) {
    print "Cannot open file \"$relaxedphylipoutputfile\" to write to!!\n\n";
    exit;
}
print RPHY $relaxedphylipresult;

close (RPHY);

my $phyloresult = join "\n", @phylotable;

my $phylooutputfile = "$output.txt";

unless ( open(PHY, ">$phylooutputfile") ) {
    print "Cannot open file \"$phylooutputfile\" to write to!!\n\n";
    exit;
}
print PHY $phyloresult;
close (PHY);

######################

sub by_site {
    my @aarray = split "\t", $a;
    my @barray = split "\t", $b;
    if ($aarray[0] lt $barray[0]) {-1} elsif ($aarray[0] gt $barray[0]) {1} elsif ($aarray[1] < $barray[1]) {-1} elsif ($aarray[1] > $barray[1]) {1} else {0}
}
