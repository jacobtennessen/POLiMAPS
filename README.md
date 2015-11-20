POLiMAPS
========

Perl scripts for analyzing polyploid DNA sequence data (Phylogenetics Of Linkage-Map-Anchored Polyploid Subgenomes)

MakeOneMapFromPileup.pl - reads a pileup file and converts it to OneMap format

MakeOneMapFromPileupNoParents.pl - similar to MakeOneMapFromPileup.pl, but assumes all samples are offspring, parent genotypes unknown.

MakeCharacterMatrix.pl - finds phylogenetic markers sharing sequence reads with linkage map markers

MakeCharacterMatricesWithInvariant.pl - finds phylogenetic markers sharing sequence reads with linkage map markers; includes all sites even invariants sites and creates a separate matrix for each reference genome chromosome
