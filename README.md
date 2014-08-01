vcf_fa_extractor
================

Extract reference and all variant sequences from a vcf file into a multi-fasta file


VCF format
----------

VCF format is a tab-separated text file having the following columns:

0.    Chromosome name
1.    Position
2.    Variant ID
3.    Reference genome
4.    Alternative (i.e. variant)
5.    Quality score
6.    Filter (whether or not the variant passed quality filters)
7.    INFO : Generic information about this variant. SnpEff adds annotation information in this column.

VCF EFF field
-------------
Effects information is added to the INFO field using an 'EFF' tag. There can be multiple effects separated by comma. The format for each effect is:

EFF= Effect ( Effect_Impact | Functional_Class | Codon_Change | Amino_Acid_Change| Amino_Acid_Length | Gene_Name | Transcript_BioType | Gene_Coding | Transcript_ID | Exon_Rank  | Genotype_Number [ | ERRORS | WARNINGS ] )

EFF Sub-field   Meaning
Effect  Effect of this variant. See details here.
Effect impact   Effect impact {High, Moderate, Low, Modifier}. See details here.
Functional Class    Functional class {NONE, SILENT, MISSENSE, NONSENSE}.
Codon_Change / Distance     Codon change: old_codon/new_codon OR distance to transcript (in case of upstream / downstream)
Amino_Acid_Change   Amino acid change: old_AA AA_position/new_AA (e.g. 'E30K')
Amino_Acid_Length   Length of protein in amino acids (actually, transcription length divided by 3).
Gene_Name   Gene name
Transcript_BioType  Transcript bioType, if available.
Gene_Coding     [CODING | NON_CODING]. This field is 'CODING' if any transcript of the gene is marked as protein coding.
Transcript_ID   Transcript ID (usually ENSEMBL IDs)
Exon/Intron Rank    Exon rank or Intron rank (e.g. '1' for the first exon, '2' for the second exon, etc.)
Genotype_Number     Genotype number corresponding to this effect (e.g. '2' if the effect corresponds to the second ALT)
Warnings / Errors   Any warnings or errors (not shown if empty). 
