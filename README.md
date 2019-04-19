POLiMAPS
========

Perl scripts for analyzing polyploid DNA sequence data (Phylogenetics Of Linkage-Map-Anchored Polyploid Subgenomes)

MakeOneMapFromPileup.pl - reads a pileup file and converts it to OneMap format

MakeOneMapFromPileupNoParents.pl - similar to MakeOneMapFromPileup.pl, but assumes all samples are offspring, parent genotypes unknown.

MakeCharacterMatricesNew.pl - finds phylogenetic markers sharing sequence reads with linkage map markers

AlleleFreqsFromPileup.pl - reads a pileup file and returns within-sample frequencies of putative SNPs

FstFromPooledFreqs.pl - reads an allele frequency table and calculates pairwise Fst

[MakeCharacterMatricesNew.pl includes the functionality previously represented by MakeCharacterMatrix.pl and MakeCharacterMatricesWithInvariant.pl in POLiMAPS v1.0, and also fixes several bugs that were caused by certain input file formats]

[AlleleFreqsFromPileup.pl and FstFromPooledFreqs.pl include the functionality previously represented by GOPOPS, https://github.com/jacobtennessen/GOPOPS]

References:
Tennessen JA, Govindarajulu R, Ashman TL, Liston A. 2014. Evolutionary origins and dynamics of octoploid strawberry subgenomes revealed by dense targeted capture linkage maps. Genome Biol Evol. 6(12):3295-313. doi: 10.1093/gbe/evu261.

Tennessen JA, Wei N, Straub SCK, Govindarajulu R, Liston A, Ashman TL. 2018. Repeated translocation of a gene cassette drives sex-chromosome turnover in strawberries. PLoS Biol. 16(8):e2006062. doi: 10.1371/journal.pbio.2006062

Willoughby JR, Harder AM, Tennessen JA, Scribner KT, Christie MR. 2018. Rapid genetic adaptation to a novel environment despite a genome-wide reduction in genetic diversity. Mol Ecol. 27(20):4041-4051. doi: 10.1111/mec.14726
