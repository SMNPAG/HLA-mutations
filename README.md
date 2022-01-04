# HLA-mutations-IAA-PNH

export JAVA_HOME=/usr/java/jre1.8.0_058
export PATH=$PATH: /samtools-1.1
export PATH=$PATH: /picard-tools-1.119
export PATH=$PATH: /novocraft

HLA/IMGT database files: HLA.dat, hla_gen.fa, hla_nuc.fa, HLA alignments 




Current  pipeline for the mutational study of HLA region from target sequencing data.

This workflow is designed to identify both class I and II mutations, in One Sample mode or Tumor/Normal sample mode.
 
Most of the softwares used are publicly available, except for Novoalign/Novocraft: http://www.novocraft.com/products/novoalign/
This is a short read alignment tool, chosen for its versatility especially in highly polymorphic regions. Here the characteristics:

Packed with unsurpassed features:
Mapping with base quality values.
Alignment quality scores using posterior alignment probability.
Paired end alignment.
Mismatches and gaps of up to 50% of read length.
Use of ambiguous codes in reference sequences can be used to reduce allelic bias.
Automatic base quality calibration.
Handles single end and paired end reads up to 950bp/read.
In built adapter trimming and base quality trimming.
Option for amplicon primer trimming.


A java based program has been built for point 2, 8 and 9

Steps for the HLA mutational pipeline:

1) 8-digit Allelic Inference

2) Reference construction from HLA_gen/HLA_nuc Fasta files of IMGT/HLA database based on patient high resolution typing*

3) Index of reference (Picard/Samtool)

4) Alignment to the reference with following criteria (Novoalign): Intermediate clipping, Multialignment, High Seed lenght, High Mismatch and gap penalty, Trimming
   Alternatively other alignment tools can be used, if the combination of parameters listed above is respected, producing a reasonable depth. 

5) GATK pipeline for post-alignment preprocessing (Samtool and Picard: sorting, Markduplicate, Add/Replace read group, Read calibration, Bam validation, Indexing)

6) Samtool Pileup

7) Variant calling with Varscan (This variant caller can be used in Somatic or Tumor only mode). Identification of SNP or Indel from Pileup files and filtering process based on stringent criteria for coverage, Avg quality, pval, strands). Other variant callers can be used (Mutect, deepvariant, but varscan is more flexible in terms of parameter inputs)


8) Filtering for polymorphisms and artefacts (i.e. misaligned reads from the other allele in the same locus).*

9) "Topographic" annotation based on the genomic sequences provided by the IMGT/HLA database *

10) Functional prediction of the mutations.

Evaluation of missense, in frame and frameshift mutations in exons. 
If mutation in UTR regions: computational prediction of transcription factor or miRNA binding sites
If mutation in introns: computational prediction of splicing sites proximity.



A java based program has been compiled for the automatization of points 2, 8 and 9. The filter of polymorphisms is applied to the already filtered variant calling output in order to identify only mutational events non reported in multialignment files provided by the IMGT/HLA database.

Only mutations potentially desruptive are retained (i.e point mutations in introns, if not affecting splicing sites, are discarded)

Genomic positions are calculated basing on the fasta reference corresponding to the mutated allele.


