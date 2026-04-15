# Viral Motif Enrichment, Density, and Mutation Signature Analysis

## Overview

This folder contains scripts for analyzing dinucleotide and trinucleotide motif enrichment and density across cynomolgus macaque polyomavirus genomes (MafaPyV2, MafaPyV3, and SV40 type IIB), and for characterizing intra-host variants by their trinucleotide mutational contexts. Mutations were normalized to genomic trinucleotide frequencies and visualized in R using ggplot2.

---

## Workflow

### 1. Dinucleotide density across the genome — `motif_density_v2.txt`

Dinucleotide density was calculated across each viral genome using 100-bp non-overlapping windows. All complete genomes for a given virus should be combined into a single FASTA file before running.

```bash
module load perl

perl motif_density_v2.txt \
    --kmer 2 \
    --interval 100 \
    --window 100 \
    -f MafaPyV2_completenucleotide.fasta > Mafapyv2motifdensity.txt

perl motif_density_v2.txt \
    --kmer 2 \
    --interval 100 \
    --window 100 \
    -f Mafa_SV40.fasta > MafaSV40motifdensity.txt

perl motif_density_v2.txt \
    --kmer 2 \
    --interval 100 \
    --window 100 \
    -f Mafapyv3.fasta > Mafapyv3motifdensity.txt
```

Output columns: sample ID, window position, motif, count.

---

### 2. Dinucleotide and trinucleotide motif enrichment — `motif_enrichment_v4.txt`

Dinucleotide and trinucleotide motif enrichment was calculated for each aligned genome using Markov modeling. The script computes observed probability, expected probability, and enrichment ratio (observed/expected) for each motif. All complete genomes for a given virus should be combined into a single FASTA file before running.

```bash
module load perl

perl motif_enrichment_v4.txt fasta MafaPyV2_completenucleotide.fasta
perl motif_enrichment_v4.txt fasta Mafa_SV40.fasta
perl motif_enrichment_v4.txt fasta Mafapyv3.fasta
```

Output is written to `<input_filename>.full.txt` in the same directory. Output columns: index, motif, observed probability, expected probability, enrichment ratio, sample ID.

> **Note:** BioPerl must be installed before running these scripts. On Biowulf, load the perl module before running: `module load perl`.

---

### 3. Trinucleotide contexts of intra-host variants — `count_trinuc_muts.py`

Trinucleotide contexts were extracted for each intra-host variant (from the inverted AF VCF files produced in the `variant_calling` pipeline) using the viral reference genome. This was run per sample and then combined by virus.

```bash
python3 count_trinuc_muts.py Mafa_viruses.fasta \
    <sample>_Mafa_invertAF.vars.vcf \
    -o <sample>_Mafa_trinuc_muts.txt
```

Results were then combined by virus across all samples:

```bash
grep MafaPyV2 *_Mafa_trinuc_muts.txt > Combined.Mafapyv2.trinuc.muts.txt
grep MafaPyV3 *_Mafa_trinuc_muts.txt > Combined.Mafapyv3.trinuc.muts.txt
grep Mafa_SV40 *_Mafa_trinuc_muts.txt > Combined.MafaSV40.trinuc.muts.txt
```

Output columns: chromosome, position, 5' tetranucleotide context, 3' tetranucleotide context, trinucleotide context, mutation type, trinucleotide mutation string, strand, 41bp flanking sequence, and counts of C, TC, TCA, TCT, YTCA, and RTCA motifs in the flanking region.

---

### 4. Genomic trinucleotide frequencies — `count_genomic_trints.py`

Genomic trinucleotide and dinucleotide counts were calculated for each viral reference genome to use as the denominator for mutation normalization.

```bash
python3 count_genomic_trints.py Mafapyv2.fasta > Mafapyv2_genomic_trinuc.txt
python3 count_genomic_trints.py Mafapyv3.fasta > Mafapyv3_genomic_trinuc.txt
python3 count_genomic_trints.py Mafa_SV40.fasta > MafaSV40_genomic_trinuc.txt
```

Output columns: motif, count. Reverse complement counts for trinucleotides with C or T at the middle position should be added before calculating scaling factors. Mutation counts per trinucleotide context were then normalized to genomic trinucleotide frequencies and visualized in R using ggplot2.

> **Note:** Ensure you are using the correct Python environment with Biopython installed before running. On Biowulf, activate your conda environment first: `conda activate base`.

---

## Dependencies

| Tool | Version | Usage |
|------|---------|-------|
| Perl | 5.36 | `motif_density_v2.txt`, `motif_enrichment_v4.txt` |
| BioPerl | — | Required by `motif_density_v2.txt` and `motif_enrichment_v4.txt` |
| Python | 3.10 | `count_trinuc_muts.py`, `count_genomic_trints.py` |
| [Biopython](https://biopython.org/) | 1.81 | Required by `count_trinuc_muts.py` and `count_genomic_trints.py` |
| [R](https://www.r-project.org/) | 4.5.1 | Visualization |
| [ggplot2](https://ggplot2.tidyverse.org/) | 4.0.0 | Mutation signature and motif enrichment plots |

---

## Scripts in this folder

| Script | Description |
|--------|-------------|
| `motif_density_v2.txt` | Calculates dinucleotide density across viral genomes in sliding windows |
| `motif_enrichment_v4.txt` | Calculates dinucleotide and trinucleotide motif enrichment using Markov modeling |
| `count_trinuc_muts.py` | Extracts trinucleotide contexts of intra-host variants from VCF files |
| `count_genomic_trints.py` | Counts genomic trinucleotide and dinucleotide frequencies from a reference FASTA |

---

## Notes
- All complete genomes for a given virus should be combined into a single FASTA file before running `motif_density_v2.txt` and `motif_enrichment_v4.txt`
- The `Mafa_viruses.fasta` reference used for `count_trinuc_muts.py` is the same combined viral reference used in the `variant_calling` pipeline
- Update all hardcoded file paths to match your environment before running
