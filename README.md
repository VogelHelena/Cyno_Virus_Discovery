# Identification of disease associated polyomaviruses from cynomolgus macaques undergoing HSCT

## Background
This repository contains scripts used for de novo assembly and annotation of viral sequences from Illumina short read sequencing of RCA urine DNA from cynomolgus macaques undergoing hematopoietic stem cell transplantation (HSCT).

## Methods Overview

### Sequence alignment
Illumina sequencing reads were trimmed using FASTP (v0.23.2) with default settings. Trimmed reads were aligned with Bowtie2 (v2.5.3) to the cynomolgus macaque host genome ([GCF_012559485.2, MFA1912RKSv2](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_012559485.2/)) and host reads were removed.

### De novo assembly and annotation
All reads not mapping to the host genome were de novo assembled using MEGAHIT (v1.2.9) with default parameters. Assembled contigs were annotated using BLASTn and DIAMOND against the NCBI database for closely related species, and Cenote-Taker2 with default parameters was used to identify more divergent species. Complete polyomavirus genomes assembled through the viral discovery pipeline were used to construct a combined host and viral reference genome. Reads from all samples were realigned to this reference using Bowtie2 to quantify reads matching each polyomavirus genome. A sample was considered positive for a polyomavirus if it had at least 8 reads aligned to the viral reference genome with Bowtie2.

### DNA Virome analysis
CenoteTaker2 contig annotations were joined to coverage data and viral organisms were collapsed into taxonomic groups. Read abundance was normalized to reads per million (RPM) total filtered reads and visualized in R using ggplot2.

### Point mutation variant calling
Genomic variants of cynomolgus macaque polyomaviruses were called using LoFreq (v2) and bcftools (SAMtools v1.21). Variants with amino acid changes were identified with SnpEff (v5.2). AlphaFold (v2.3.2) was used to predict the structure of the VP1 pentamer for MafaPyV2. Inter- and intra-host variants were mapped to the predicted structure.

### Viral motif enrichment, density, and mutation signature analysis
Dinucleotide and trinucleotide motif enrichments for each aligned viral genome were calculated using Markov modeling. Dinucleotide density was calculated across the viral genomes using 100-bp non-overlapping windows. Smoothed fitted lines and 95% confidence intervals of these densities were calculated. Intra-sample variants were characterized by their trinucleotide contexts and single nucleotide substitution type. All mutations detected per virus were normalized by their respective reference motif abundance and visualized using R.

---

## Repository Structure

| Folder | Description |
|--------|-------------|
| [`viral_discovery/`](./viral_discovery/) | Read trimming, host read removal, de novo assembly, and contig annotation |
| [`variant_calling/`](./variant_calling/) | Variant calling, amino acid change identification, and AlphaFold structure prediction |
| [`motif_mutation_analysis/`](./motif_mutation_analysis/) | Motif enrichment, dinucleotide density, and mutation signature analysis |

---

## Dependencies

| Tool | Version | Usage |
|------|---------|-------|
| [FASTP](https://github.com/OpenGene/fastp) | 0.23.2 | Read quality filtering |
| [Bowtie2](https://bowtie-bio.sourceforge.net/bowtie2) | 2.5.3 | Read alignment |
| [SAMtools](http://www.htslib.org/) | 1.21 | BAM processing |
| [BBTools](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/) | 39.06 | Paired-end read repair |
| [bioawk](https://github.com/lh3/bioawk) | 1.0 | SAM/FASTQ parsing |
| [MEGAHIT](https://github.com/voutcn/megahit) | 1.2.9 | De novo assembly |
| [BLAST+](https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html) | 2.15.0 | Nucleotide sequence alignment/annotation |
| [DIAMOND](https://github.com/bbuchfink/diamond) | 2.1.12 | Protein sequence alignment/annotation |
| [CenoteTaker2](https://github.com/mtisza1/Cenote-Taker2) | 2.1.5 | Viral sequence prediction |
| [LoFreq](https://csb5.github.io/lofreq/) | 2 | Variant calling |
| [bcftools](https://samtools.github.io/bcftools/) | 1.21 | Variant calling |
| [SnpEff](https://pcingola.github.io/SnpEff/) | 5.2 | Variant annotation |
| [AlphaFold](https://github.com/google-deepmind/alphafold) | 2.3.2 | Protein structure prediction |
| [R](https://www.r-project.org/) | 4.5.1 | Motif and mutation visualization |

---

## Notes
- The host reference genome (`MFAV2`) is the *Macaca fascicularis* assembly GCF_012559485.2 (MFA1912RKSv2), available from NCBI as `GCF_012559485.2_MFA1912RKSv2_genomic.fna`
- All scripts are designed to run in a SLURM HPC environment and use `$SLURM_CPUS_PER_TASK` and `$SLURM_JOB_ID` environment variables
- Update file paths to match your environment before running
- Scripts should be submitted via `sbatch` with appropriate resource allocations for your cluster
