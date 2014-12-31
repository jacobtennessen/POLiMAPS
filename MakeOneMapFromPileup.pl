#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Std;

##### ##### ##### ##### #####

use vars qw( $opt_g $opt_f $opt_c $opt_d $opt_m $opt_o);

# Usage
my $usage = "
MakeOneMapFromPileup.pl - reads a pileup file and converts it to OneMap format
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

Usage: perl MakeOneMapFromPileup.pl options
 required:
  -g	a pileup file. First two individuals are the parents (maternal then paternal), and the rest are offspring
 optional:
  -f	a variant must occur with a frequency at least the reciprocal of this value [default = 40]
  -c	a variant must occur at least this many times [default = 2]
  -d	minimum depth to call genotypes [default = 32]
  -m	maximum number of offspring with missing data [default = 1]
  -o	minimum number of heterozygotes and homozygotes observed in offspring [default = 8]

";

#############

# command line processing.
getopts('g:f:c:d:m:o:');
die $usage unless ($opt_g);

my ($pileup, $minallelefractionconstant, $minminorallelecount, $mindepth, $missingthresh, $mingenos);

$pileup	= $opt_g if $opt_g;
$minallelefractionconstant = $opt_f ? $opt_f : 40;
$minminorallelecount = $opt_c ? $opt_c : 2;
$mindepth = $opt_d ? $opt_d : 32;
$missingthresh = $opt_m ? $opt_m : 1;
$mingenos = $opt_o ? $opt_o : 8;

my @pileupfiledata = split /\//, $pileup;

my $genotypes = pop @pileupfiledata;

my $outputfile = "OneMap_$genotypes";

if (defined $pileupfiledata[0]) {
    my $dir = join "/", @pileupfiledata;
    $genotypes = "$dir/$genotypes";
    $outputfile = "$dir/$outputfile";
}

my %parents;

$parents{3} = 1;
$parents{6} = 1;

my $highestallelenumber = 2; #maximum number of alleles allowed per SNP.
    
my $total; #number of columns in pileup file

my @indnumbers; #columns with genotype info for each individual

my $goodline = 0; #count of distinct variable sites

my $indcount; #number of individuals

my @out; #data in OneMap format

open(IN, "$genotypes") || die "can't open $genotypes\n";

while (<IN>) {
    my $line = $_;
    $line =~ s/\r|\n//g;
    my @data = split "\t", $line;
    unless (defined $total) {
	$total = scalar (@data);
	$indcount = ($total - 3)/3;
	for (my $x = 3; $x < $total; $x +=3) {
	    push @indnumbers, $x;
	}
    }
    my %allseenbases;
    $allseenbases{$data[2]} = 1;
    my @genos;
    my %parentgenos;
    my $hetparent = 0;
    my $missingcount = 0;
    my $homocount = 0;
    my $hetcount = 0;
    my $chromsite = "$data[0]_$data[1]";
    my %seengenos;
    foreach my $id (@indnumbers) {
	my $localgeno = "-";
	if ($data[$id] >= $mindepth) {
	    my @seqdata = split "", $data[$id+1];
	    my $extra = 0;
	    my $indel = 0;
	    my $postindel = 0;
	    my %bases;
	    my $allgood = 0;
	    foreach my $s (@seqdata) {
		if ($indel == 1) {
		    if ($s =~ /\d/) {
			$extra = $s;
		    } else {
			next;
		    }
		    $indel = 0;
		    $postindel = 1;
		} elsif ($extra > 0) {
		    if ($postindel == 1) {
			if ($s =~ /\d/) { #indel multiple digits
			    $extra = ($extra*10)+$s;
			} else {
			    $postindel = 0;
			    $extra -=1;
			}
		    } else {
			$extra -=1;
		    }
		} elsif (($s =~ /,/)||($s =~ /\./)) {
		    if (defined $bases{$data[2]}) {
			$bases{$data[2]} +=1;
		    } else {
			$bases{$data[2]} = 1;
		    }
		    $allgood +=1;
		} elsif (($s =~ /A/gi)||($s =~ /C/gi)||($s =~ /G/gi)||($s =~ /T/gi)) {
		    my $uc = uc $s;
		    if (defined $bases{$uc}) {
			$bases{$uc} +=1;
		    } else {
			$bases{$uc} = 1;
		    }
		    $allgood +=1;
		} elsif ($s =~ /\^/) {
		    $extra = 1;
		} elsif (($s =~ /-/)||($s =~ /\+/)) {
		    $indel = 1;
		} elsif (($s =~ /\$/)||($s =~ /\*/gi)||($s =~ /</gi)||($s =~ />/gi)||($s =~ /N/gi)) {
		    next;
		} else {
		    next;
		}
	    }
	    if ($allgood >= $mindepth) {
		$localgeno = "a";
		my $goodbases = 0;
		my $toolowtocall = 0;
		foreach my $b (keys %bases) {
		    if (($bases{$b} >= $allgood/$minallelefractionconstant)&&($bases{$b} >= $minminorallelecount)) {
			unless (defined $allseenbases{$b}) {
			    my $bcount = (scalar(keys %allseenbases))+1;
			    $allseenbases{$b} = $bcount;
			}
			$goodbases +=1;
		    } elsif (($bases{$b} >= $allgood/(2*$minallelefractionconstant))||($bases{$b} > 1)) {
			$toolowtocall = 1;
		    }
		}
		if ($goodbases > 1) {
		    $localgeno = "ab";
		}
		if ($toolowtocall == 0) {
		    $seengenos{$localgeno} = 1;
		} else {
		    $localgeno = "-";
		}
	    }
	}
	if (defined $parents{$id}) {
	    $parentgenos{$id} = $localgeno;
	    if ($localgeno =~ /ab/) {
		$hetparent += 1;
	    } elsif ($localgeno =~ /-/) {
		$missingcount += (10*$missingthresh);
	    }
	} else {
	    push @genos, $localgeno;
	    if ($localgeno =~ /-/) {
		$missingcount +=1;
	    } elsif ($localgeno =~ /ab/) {
		$hetcount +=1;
	    } else {
		$homocount +=1;
	    }
	}
    }
    if ((scalar(keys %allseenbases) <= 1)||(scalar(keys %allseenbases) > $highestallelenumber)) {
	next;
    }
    if ($hetparent == 0) {
	next;
    }
    if ($missingcount > $missingthresh) {
	next;
    }
    if (($hetcount < $mingenos)||($homocount < $mingenos)) {
	next;
    }
    if (($hetparent == 1)&&((scalar(keys %seengenos)) == 3)) {
	next;
    }
    my $model;
    if ($hetparent == 2) {
	$model = "C.8";
	foreach my $g (@genos) {
	    if ($g =~ /ab/) {
		$g = "a";
	    } elsif ($g =~ /a/) {
		$g = "o";
	    }
	}
    } elsif ($parentgenos{3} =~ /ab/) {
	$model = "D1.10";
    } elsif ($parentgenos{6} =~ /ab/) {
	$model = "D2.15";
    }
    my $genoline = join ",", @genos;
    push @out, "*$chromsite\t$model\t$genoline";
    $goodline +=1;
}

close (IN);

my $result = join "\n", @out;

unless ( open(OUT, ">$outputfile") ) {
    print "Cannot open file \"$outputfile\" to write to!!\n\n";
    exit;
}

my $offspringcount = $indcount - 2;

print OUT "$offspringcount $goodline\n$result";

close (OUT);
