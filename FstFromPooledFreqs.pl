#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Std;

##### ##### ##### ##### #####

use vars qw( $opt_f $opt_d $opt_t  );

# Usage
my $usage = "
FstFromPooledFreqs.pl - reads an allele frequency table and calculates pairwise Fst
Copyright (C) 2017 by Jacob A Tennessn 
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
Fst subroutine modeled after Perlymorphism package:
Stajich JE, Hahn MW (2005) Disentangling the effects of demography and selection in human history. Mol Biol Evol. 22:63-73
Usage: perl FstFromPooledFreqs.pl options
 required:
  -f	an allele frequency table (as generated by AlleleFreqsFromPileup.pl)
 optional:
  -d    mimimum depth per sample [default = 8]
  -t    minimum Fst to report [default = 0.2]
";

#############

# command line processing.
getopts('f:d:t:');
die $usage unless ($opt_f);

my ($freqfile, $mindepth, $fstthresh);

$freqfile = $opt_f if $opt_f;
if (defined $opt_d) {
    $mindepth = $opt_d;
} else {
    $mindepth = 8; 
}
if (defined $opt_t) {
    $fstthresh = $opt_t;
} else {
    $fstthresh = 0.2;
}

my @freqfiledata = split /\//, $freqfile;

my $genotypes = pop @freqfiledata;

my $outputfile = "High_Fsts_$genotypes";

if (defined $freqfiledata[0]) {
    my $dir = join "/", @freqfiledata;
    $genotypes = "$dir/$genotypes";
    $outputfile = "$dir/$outputfile";
}

my $total;

my @out;

open(IN, "$genotypes") || die "can't open $genotypes\n";

while (<IN>) {
    my $line = $_;
    $line =~ s/\r|\n//g;
    my @data = split "\t", $line;
    if ($data[0] =~ /^Allele$/) {
        next;
    }
    unless (defined $total) {
        $total = ((scalar(@data)) - 1)/2;
        my @names;
        for (my $s1 = 1; $s1 <= ($total-1); $s1+=1) {
            for (my $s2 = ($s1+1); $s2 <= $total; $s2+=1) {
                push @names, "$s1,$s2";
            } 
        }
        my $namelist = join "\t", @names;
        push @out, "Allele\t$namelist";
    }
    my @fsts;
    my $seengood = 0;
    for (my $s1 = 1; $s1 <= ($total-1); $s1+=1) {
        my $freqpos1 = $s1*2-1;
        my $depthpos1 = $s1*2;
        for (my $s2 = ($s1+1); $s2 <= $total; $s2+=1) {
            my $fstat = "NA";
            if (($data[$depthpos1] =~ /\d/)&&($data[$freqpos1] =~ /\d/)&&($data[$depthpos1] >= $mindepth)) {
                my $freqpos2 = $s2*2-1;
                my $depthpos2 = $s2*2;
                if (($data[$depthpos2] =~ /\d/)&&($data[$freqpos2] =~ /\d/)&&($data[$depthpos2] >= $mindepth)) {
                    my $candfstat = Fst ($data[$depthpos1],$data[$depthpos2],$data[$freqpos1],$data[$freqpos2]);
                    if ($candfstat >= $fstthresh) {
                        $fstat = sprintf "%.4f", $candfstat;
                        $seengood = 1;
                    }
                }
            }
            push @fsts, $fstat;
        }
    }
    if ($seengood == 1) {
        my $fstlist = join "\t", @fsts;
        push @out, "$data[0]\t$fstlist";
    }
}

close (IN);

if (defined $out[0]) {
    my $result = join "\n", @out;
    unless ( open(OUT, ">$outputfile") ) {
        print "Cannot open file \"$outputfile\" to write to!!\n\n";
        exit;
    }
    print OUT $result;
    close (OUT);
}


#########################

sub Fst {
    my ($size1, $size2, $freq1, $freq2) = @_;
    
    my @freqs = ($freq1, $freq2);
    my @samplesizes = ($size1, $size2);
    
    my $num_sub_pops = 2;

    my $Fst;
    my ($TS_sub1,$TS_sub2);
    my @alleles = (0,1);

    my $avg_samp_size         = 0;
    my $avg_allele_freq       = 0;
    my $total_samples_squared = 0;

    foreach (my $c = 0; $c < 2; $c ++) {
	my $s = $samplesizes[$c];
	$avg_samp_size += $s;
	$total_samples_squared += $s**2;
	my $all_freq = $freqs[$c];
	$avg_allele_freq += $s * $all_freq;
    }
    
    my $total_samples =  $avg_samp_size;
    $avg_samp_size /= $num_sub_pops;
    $avg_allele_freq /= $total_samples;

    my $adj_samp_size = ( 1/ ($num_sub_pops - 1)) * ( $total_samples - ( $total_samples_squared/$total_samples));

    my $variance              = 0;
    my $sum_variance          = 0;
    my $i = 0;
    for (my $d = 0; $d <2; $d ++) {
	my $s = $samplesizes[$d];
	$sum_variance += $s * (($freqs[$d] - $avg_allele_freq)**2);
    }
    $variance = ( 1 / (( $num_sub_pops-1)*$avg_samp_size))*$sum_variance;

	       $TS_sub1 = $variance - 
		   ( ( 1/($avg_samp_size-1))*
		     ( ($avg_allele_freq*(1-$avg_allele_freq))-
		       ( (($num_sub_pops-1)/$num_sub_pops)*$variance)));
	       $TS_sub2 = ( (($adj_samp_size-1)/($avg_samp_size-1))*
			      $avg_allele_freq*(1-$avg_allele_freq) ) +
			      ( 1 + ( (($num_sub_pops-1)*
				       ($avg_samp_size-$adj_samp_size))/ 
				      ($avg_samp_size - 1))) * 
				      ($variance/$num_sub_pops);
    
    unless (($freq1 == $freq2)&&(($freq1 == 0)||($freq1 == 1))) {
        $Fst = $TS_sub1 / $TS_sub2;
    } else {
        $Fst = 0;
    }
   
    return $Fst;
}
