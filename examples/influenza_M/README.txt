N1_Genbank.fasta is the set of all N1 NA sequences of influenza A viruses from any hosts, using only full length sequences, removing identical sequences, and excluding lab strains, Downloaded from the Influenza Virus Resource (http://www.ncbi.nlm.nih.gov/genomes/FLU/FLU.html) on 6-21-10.

N1_GISAID.fasta is the set of all N1 NA sequences of influenza A viruses of length at least 1389 as downloaded from GISAID on 6-21-10.

N1.fasta is the merge of all unique sequences in N1_Genbank.fasta and N1_GISAID.fasta, created by running merge_Genbank_and_GISAID.py

N4.fasta contains an N4 neuraminidase to be used as an outgroup.

NA_CA09.fasta is the protein sequence of the neuraminidase from the A/California/04/2009 (H1N1) strain.

PIPS_BUILD_TREE_AND_ALIGNMENT.IN is the input for the pips_build_tree_and_alignment.py script, and PIPS_BUILD_TREE_AND_ALIGNMENT.LOG, seqnames.txt, tree.newick, and renamed_alignment.fasta are the output files.

PIPS_ANALYSIS.IN is the input for pips_analysis.py, and PIPS_ANALYSIS.LOG and PIPS_PREDICTIONS_PRIORS-REGULARIZING.TXT are the output files.

SORTED_PIPS_PREDICTIONS_PRIORS-REGULARIZING.TXT is the result of running analyzing_most_stabilizing.py.

TOP_TWELVE_PREDICTIONS.TXT is the manually annotated set of the top 12 mutations, to be constructed and tested experimentally.
