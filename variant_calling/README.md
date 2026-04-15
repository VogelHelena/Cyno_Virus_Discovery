# Point Mutation Variant Calling and Structural Analysis

## Overview

This folder contains scripts for calling and annotating genomic variants of cynomolgus macaque polyomaviruses (MafaPyV2, MafaPyV3, and SV40 type IIB) from Illumina short read sequencing data, and for predicting and visualizing viral protein structures using AlphaFold.

Variants are analyzed in two ways:
- **Inter-host variants**: variants called directly from each sample's BAM file compared against the genotype I consensus reference genome, representing differences between the virus in each animal and the reference
- **Intra-host (subclonal) variants**: identified by inverting allele frequencies (AF) greater than 0.5 in each VCF file, then removing variants with AF=0 after inversion; these represent minority variants circulating within a single host

---

## Workflow

### 1. Build combined reference and align reads

A combined reference FASTA was constructed by concatenating the cynomolgus macaque host genome with the three viral reference genomes assembled with the viral discovery pipeline. A Bowtie2 index was then built from this combined reference.

```bash
cat GCF_012559485.2_MFA1912RKSv2_genomic.fna Mafa_viruses.fasta > MFAV2_Mafa_viruses.fasta
module load bowtie
bowtie2-build --threads 2 MFAV2_Mafa_viruses.fasta MFAV2_Mafa_viruses
```

Reads from two sequencing runs were aligned simultaneously to the combined reference using `run_bowtie2.sh`. FASTP-filtered reads from both runs are provided as comma-separated inputs to Bowtie2.

```bash
#!/bin/bash
module load bowtie
module load samtools

echo "Processing $1"

R1trim=${1}_R1_001.flt.fq.gz,/data/vogelh2/231218_Mafa_urine/Data/fastq/${1}_R1_001.flt.fq.gz
R2trim=${1}_R2_001.flt.fq.gz,/data/vogelh2/231218_Mafa_urine/Data/fastq/${1}_R2_001.flt.fq.gz
R3trim=/data/vogelh2/231218_Mafa_urine/Data/fastq/${1}_R2_unmerged.flt.fq.gz
BAM=${1}.combined.bam
REF=/data/vogelh2/Ref/MFAV2_Mafa_viruses

bowtie2 --very-sensitive --no-unal \
    --rg-id $1 \
    --rg SM:$1 \
    --rg LB:GRC \
    --rg PL:ILLUMINA \
    --rg DS:NextSeq2000 \
    -p $SLURM_CPUS_PER_TASK \
    -x $REF \
    -1 $R1trim \
    -2 $R2trim \
    -U $R3trim | \
    samtools sort -@ 16 -m 4G -o $BAM -

samtools index -@ $SLURM_CPUS_PER_TASK $BAM
```

> **Note:** `--no-unal` suppresses SAM records for reads that failed to align, keeping only aligned reads in the BAM file.

---

### 2. Variant calling — `lofreq.sh`

Variants were called for each sample using LoFreq against the combined reference.

```bash
#!/bin/bash
module load lofreq

lofreq call \
    -f /data/vogelh2/Ref/MFAV2_Mafa_viruses.fasta \
    -o ${1}.vars.vcf \
    ${1}.combined.bam
```

---

### 3. Subset and filter VCF files

VCF files were subset to retain only Mafa virus calls and then split into per-sample files. These steps were performed interactively on the command line without a dedicated script.

```bash
# Subset VCF for Mafa virus calls only
grep 'Mafa' <sample>.vars.vcf > <sample>_Mafa.vars.vcf

# Count variants per virus per sample
grep Mafa *.vars.vcf | grep -v "#" | cut -f1 | sort | uniq -c
```

---

### 4. Invert allele frequencies — `vcfAFinverter.py`

To identify intra-host subclonal variants, allele frequencies above 0.5 were inverted in each per-sample VCF. This swaps the REF and ALT alleles and adjusts the AF and DP4 fields accordingly. After inversion, variants with AF=0 (i.e., those that were fixed at AF=1 before inversion and therefore represent inter-host rather than intra-host variants) were removed. These filtering steps were performed interactively on the command line.

```bash
python3 vcfAFinverter.py <sample>_Mafa.vars.vcf > <sample>_Mafa_invertAF.vars.vcf

# Remove variants with AF=0 after inversion
for f in *_Mafa_invertAF.vars.vcf; do
    out="${f/.vars.vcf/.clean.vars.vcf}"
    awk -F'\t' 'BEGIN{OFS="\t"}
        /^#/ { print; next }
        $8 !~ /(^|;)AF=0(\.0+)?($|;)/ { print }
    ' "$f" > "$out"
done
```

---

### 5. Variant annotation — `snpEff.sh` and `snpEff2.sh`

Variants were annotated using SnpEff to identify those causing amino acid changes. Custom SnpEff databases were built for each virus (MafaPyV2, MafaPyV3, and Mafa_SV40) using complete viral genomes assembled from the viral discovery pipeline. GFF annotation files (`genes.gff`) were created manually for each virus. See the [SnpEff documentation](https://pcingola.github.io/SnpEff/snpeff/build_db/) for instructions on building custom databases.

SnpEff must be run separately for each virus database and variant type (inter-host and intra-host). Representative scripts are provided below.

**Intra-host variants (`snpEff.sh`)** — annotates inverted AF VCF files against the MafaPyV2 database:

```bash
#!/bin/bash
module load snpEff

java -jar snpEff.jar ann Mafapyv2 \
    /data/vogelh2/240124_Mafa_urine2/Data/fastq/${1}Mafa_invertAF.vars.vcf > \
    /data/vogelh2/240124_Mafa_urine2/Data/fastq/${1}Mafa_invertAF.ann.vars.vcf
```

**Inter-host variants (`snpEff2.sh`)** — annotates original VCF files against the MafaPyV2 database:

```bash
#!/bin/bash
module load snpEff

java -jar snpEff.jar ann Mafapyv2 \
    /data/vogelh2/240124_Mafa_urine2/Data/fastq/${1}Mafa.vars.vcf > \
    /data/vogelh2/240124_Mafa_urine2/Data/fastq/${1}Mafa_inter.ann.vars.vcf
```

After annotation, missense and high-impact variants were extracted interactively:

```bash
# Extract moderate and high impact variants across all samples
grep -v "^#" *Mafa_invertAF.ann.vars.vcf | grep "MODERATE\|HIGH" > Mafapyv2.mod.high.ann.vars.vcf
```

> **Note:** When running SnpEff with the MafaPyV2 database, calls for MafaPyV3 and Mafa_SV40 will return `ERROR_CHROMOSOME_NOT_FOUND`. For those viruses, re-run using their respective databases (MafaPyV3 or Mafa_SV40). If a variant had AF inverted (>0.5 before inversion), the annotation will include `WARNING_REF_DOES_NOT_MATCH_GENOME` — this is expected and not an error.

---

### 6. AlphaFold structure prediction

Protein structures were predicted using AlphaFold v2.3.2. The VP1 pentamer of MafaPyV2 was modeled as a multimer using five identical chain sequences. Each prediction was run in two steps: first generating multiple sequence alignments (MSAs), then running the structure prediction using the precomputed MSAs.

Representative scripts for the MafaPyV2 VP1 pentamer are provided. The same approach was applied to VP2 and LT for MafaPyV2, MafaPyV3, and SV40 type IIB using the `monomer` preset for single chains and the `multimer` preset for the LT hexamer.

**MSA generation (`msa_script.multimer.sh`):**

```bash
#!/bin/bash
module load alphafold2/2.3.2

alphafold \
    --model_preset=multimer \
    --fasta_paths=Mafapyv2_VP1.pentamer.fasta \
    --max_template_date=2023-12-31 \
    --msas_only \
    --output_dir=$PWD \
    --num_multimer_predictions_per_model=2
```

**Structure prediction (`model_script.multimer.sh`):**

```bash
#!/bin/bash
module load alphafold2/2.3.2

alphafold \
    --model_preset=multimer \
    --fasta_paths=Mafapyv2_VP1.pentamer.fasta \
    --max_template_date=2023-12-31 \
    --use_precomputed_msas \
    --output_dir=$PWD \
    --num_multimer_predictions_per_model=2
```

Predicted structures were visualized in ChimeraX (v1.7.1) and inter- and intra-host missense variants were mapped onto the predicted VP1 pentamer structure.

---

## Dependencies

| Tool | Version | Usage |
|------|---------|-------|
| [Bowtie2](https://bowtie-bio.sourceforge.net/bowtie2) | 2.5.3 | Read alignment to combined host+viral reference |
| [SAMtools](http://www.htslib.org/) | 1.21 | BAM sorting and indexing |
| [LoFreq](https://csb5.github.io/lofreq/) | 2 | Variant calling |
| [SnpEff](https://pcingola.github.io/SnpEff/) | 5.2 | Variant annotation |
| [AlphaFold](https://github.com/google-deepmind/alphafold) | 2.3.2 | Protein structure prediction |
| [ChimeraX](https://www.cgl.ucsf.edu/chimerax/) | 1.7.1 | Structure visualization and variant mapping |
| Python | 3.10 | `vcfAFinverter.py`, `count_trinuc_muts.py`, `count_genomic_trints.py` |
| [Biopython](https://biopython.org/) | 1.81 | Required by `count_trinuc_muts.py` and `count_genomic_trints.py` |

---

## Scripts in this folder

| Script | Description |
|--------|-------------|
| `run_bowtie2.sh` | Aligns reads from two sequencing runs to the combined host and viral reference |
| `lofreq.sh` | Calls variants per sample using LoFreq |
| `vcfAFinverter.py` | Inverts allele frequencies >0.5 to identify intra-host subclonal variants |
| `snpEff.sh` | Annotates intra-host variants (inverted AF VCFs) using SnpEff |
| `snpEff2.sh` | Annotates inter-host variants (original VCFs) using SnpEff |
| `msa_script.multimer.sh` | Generates MSAs for AlphaFold multimer prediction (MafaPyV2 VP1 pentamer) |
| `model_script.multimer.sh` | Runs AlphaFold structure prediction using precomputed MSAs |

---

## Notes
- Update all hardcoded file paths (e.g., `/data/vogelh2/`) to match your environment before running
- All SLURM scripts are designed to run in an HPC environment using `$SLURM_CPUS_PER_TASK`
- The viral reference FASTA (`Mafa_viruses.fasta`) contains the MafaPyV2 genotype I consensus genome, MafaPyV3 consensus, and SV40 type IIB consensus, all starting at the NCCR going into the late region
- Custom SnpEff databases must be built before running annotation. GFF annotation files were created manually for each virus based on complete genomes assembled in the `viral_discovery` pipeline
