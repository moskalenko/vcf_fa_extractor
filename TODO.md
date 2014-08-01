vcf_fa_extractor
================

2014-07-08: enhancement

We have talked about this for a while, but the time has come where I need the
codon alignment script to continue the analysis for our Cholera paper.  This
should actually just be an extension of the script you have already wrote (I
included the most recent version).  I attached a copy of an annotated VCF file,
which includes the reference codon and the alternative codon in the INFO field.
I included the annotation format and an example below.  For the script, I just
need it to use the alternative codon designated as the codon after the / in the
Condon_Change field to create a multiple alignment.  As a second part, if
possible, I would like to output a table with the following headers:

Chromosome | Position | Effect | Codon_Change | Amino_acid_change | Gene_Name | Sample 1 (yes/no) | Sample 2 (yes,no) 


Examples:

INFO=<ID=EFF,Number=.,Type=String,Description="Predicted effects for this variant.Format: 'Effect ( Effect_Impact | Functional_Class | Codon_Change | Amino_Acid_change| Amino_Acid_length | Gene_Name | Gene_BioType | Coding | Transcript | Exon  | GenotypeNum [ | ERRORS | WARNINGS ] )' ">

EFF=NON_SYNONYMOUS_CODING(MODERATE|MISSENSE|aCc/aTc|T34I|197|Vch1786_I0077||CODING|Vch1786_I0077|1|1 

ccA/ccT
    ---

Large letter is the base to use.

EFFECT = multiple columns from parsing the EFF tag
(http://snpeff.sourceforge.net/SnpEff_manual.html#output). Parse them out into
separate columns.

Sample = Genotype

1. Use the ALT from the SNPeff instead of the FreeBayes call for the Fasta
sequence.

2. Write an annotation output file with the following header.

Chromosome | Position | EFF | Sample 1 (yes/no) | Sample 2 (yes,no) 

Where EFF = Effect_Impact | Functional_Class | Codon_Change | Amino_Acid_change| Amino_Acid_length | Gene_Name | Gene_BioType | Coding | Transcript | Exon

GenotypeNum results int the last 2 columns (Sample 1 or 2 'yes')
