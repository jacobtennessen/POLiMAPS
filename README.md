POLiMAPS
========

Perl scripts for analyzing polyploid DNA sequence data (Phylogenetics Of Linkage-Map-Anchored Polyploid Subgenomes)

MakeOneMapFromPileup.pl - reads a pileup file and converts it to OneMap format

MakeOneMapFromPileupNoParents.pl - similar to MakeOneMapFromPileup.pl, but assumes all samples are offspring, parent genotypes unknown.

MakeCharacterMatricesNew.pl - finds phylogenetic markers sharing sequence reads with linkage map markers
[this script includes the functionality previously represented by MakeCharacterMatrix.pl and MakeCharacterMatricesWithInvariant.pl in POLiMAPS v1.1, and also fixes several bugs that were caused by certain input file formats]
