# Viral Discovery

## Overview

This folder contains scripts for de novo assembly and annotation of viral sequences from Illumina short read sequencing of RCA urine DNA from cynomolgus macaques, and for quantifying viral DNA abundance across samples using a DNA virome analysis pipeline.

---

## Workflow

### 1. Quality filtering — `run_fastp.sh`

Filters raw paired-end FASTQ files to remove low-quality bases and adapter sequences.

```bash
#!/bin/bash
module load fastp

fastp -w $SLURM_CPUS_PER_TASK \
    --in1 ${1}_R1_001.fastq.gz \
    --in2 ${1}_R2_001.fastq.gz \
    --out1 ${1}_R1_001.flt.fq.gz \
    --out2 ${1}_R2_001.flt.fq.gz
```

---

### 2. Reference index — Bowtie2 index build

Build a Bowtie2 index from the reference FASTA before running alignment. Only needs to be run once per reference.

> **Reference genome:** `MFAV2.fasta` is derived from the *Macaca fascicularis* reference assembly [GCF_012559485.2 (MFA1912RKSv2)](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_012559485.2/), downloaded as `GCF_012559485.2_MFA1912RKSv2_genomic.fna` from NCBI.

```bash
module load bowtie
bowtie2-build --threads 2 MFAV2.fasta MFAV2
```

---

### 3. Alignment — `run_bowtie2.sh`

Aligns filtered reads to the host reference genome and outputs a sorted, indexed BAM file.

```bash
#!/bin/bash
module load bowtie
module load samtools

echo "Processing $1"

R1trim=${1}_R1_001.flt.fq.gz
R2trim=${1}_R2_001.flt.fq.gz
BAM=${1}.bam
REF=/data/vogelh2/Ref/MFAV2

bowtie2 --very-sensitive \
    --rg-id $1 \
    --rg SM:$1 \
    --rg LB:GRC \
    --rg PL:ILLUMINA \
    --rg DS:NextSeq2000 \
    -p $SLURM_CPUS_PER_TASK \
    -x $REF \
    -1 $R1trim \
    -2 $R2trim | \
    samtools sort -@ 16 -m 4G -o $BAM -

samtools index -@ $SLURM_CPUS_PER_TASK $BAM
```

---

### 4. Host read filtering — `extra_nonhost_pairs.sh`

Extracts non-host read pairs from the BAM file. Retains reads that are unmapped, have an unmapped mate, or do not align to a named reference contig. Repairs broken paired-end reads using BBTools.

```bash
#!/bin/bash
module load bbtools
module load bioawk

samtools view -f 64 ${1}.bam | \
    bioawk -c sam '{if(and($flag,4) || and($flag,8) || $rname !~ /N/) \
    print "@"$qname"\n"$seq"\n+\n"$qual}' | \
    gzip -c > ${1}_nonhost_R1.fastq.gz &

samtools view -f 128 ${1}.bam | \
    bioawk -c sam '{if(and($flag,4) || and($flag,8) || $rname !~ /N/) \
    print "@"$qname"\n"$seq"\n+\n"$qual}' | \
    gzip -c > ${1}_nonhost_R2.fastq.gz
wait

bbtools repair \
    in=${1}_nonhost_R1.fastq.gz \
    in2=${1}_nonhost_R2.fastq.gz \
    out=${1}_nonhost_R1.repair.fastq.gz \
    out2=${1}_nonhost_R2.repair.fastq.gz \
    outs=${1}_nonhost_S.repair.fastq.gz
```

---

### 5. De novo assembly — `run_megahit_nonhost.sh`

Assembles non-host reads into contigs using MEGAHIT.

```bash
#!/bin/bash
module load megahit

megahit \
    -t $SLURM_CPUS_PER_TASK \
    -1 ${1}_nonhost_R1.repair.fastq.gz \
    -2 ${1}_nonhost_R2.repair.fastq.gz \
    -o ${1}_megahit_nonhost \
    --tmp-dir /lscratch/$SLURM_JOB_ID
```

---

### 6. Nucleotide annotation — `run_blastn_nonhost.sh`

Aligns assembled contigs against a nucleotide reference database using BLASTN. DIAMOND was also run in parallel for protein-level annotation and comparison of results.

```bash
#!/bin/bash
module load blast

blastn \
    -max_target_seqs 3 \
    -num_threads $SLURM_CPUS_PER_TASK \
    -query ${1} \
    -out ${1}.blastn.txt \
    -evalue "1e-10" \
    -db ${2} \
    -outfmt "6 qseqid sseqid stitle qstart qend sstart send length mismatch gapopen sstrand evalue bitscore"
```

---

### 7. Protein annotation — `run_blastx_nonhost.sh`

Translates and aligns assembled contigs against a protein database using BLASTX.

```bash
#!/bin/bash
module load blast

blastx \
    -max_target_seqs 10 \
    -num_threads $SLURM_CPUS_PER_TASK \
    -db ${2} \
    -outfmt "6 qseqid sseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq qframe" \
    -query ${1} \
    -out ${1}.blastx.txt
```

---

### 8. Viral sequence prediction — `run_cenote.sh`

Runs CenoteTaker2 to predict and annotate viral sequences from assembled contigs.

```bash
#!/bin/bash
module load prodigal samtools blast bioawk edirect hmmer bowtie trnascan-se bedtools

python /data/vogelh2/Software/Cenote-Taker2/run_cenote-taker2.py \
    -p True \
    -db standard \
    --minimum_length_linear 600 \
    --lin_minimum_hallmark_genes 1 \
    -m 32 \
    -t 2 \
    -c ${1}_megahit/final.contigs.fa \
    -r ${1}_cenote
```

> **Polyomavirus positivity threshold:** Complete polyomavirus genomes assembled through this pipeline were used to construct a combined reference containing the host genome and all viral reference genomes (see `variant_calling/README.md`). Sequencing reads were then realigned to this combined reference using Bowtie2 to enumerate the true number of reads mapping to each polyomavirus genome. A sample was considered positive for a polyomavirus if at least 8 reads aligned to the viral reference genome.

---

## DNA Virome Analysis

Following viral sequence prediction with CenoteTaker2, a DNA virome analysis was performed to quantify viral read abundance per sample and summarize the composition of the viral community.

### 9. Contig coverage — `bowtie2assembly.sh` and `pileup3.sh`

To calculate per-contig read coverage, a Bowtie2 index was built from each sample's MEGAHIT assembly, and reads from both sequencing runs were re-aligned to that per-sample assembly. Coverage statistics were then calculated using BBTools pileup.

```bash
# Build per-sample Bowtie2 index — bowtie2assembly.sh
#!/bin/bash
module load bowtie

bowtie2-build --threads 2 \
    /data/vogelh2/240124_Mafa_urine2/Data/fastq/${1}_megahit_nonhost_combined/final.contigs.fasta \
    ${1}_megahit_combined
```

```bash
# Align reads to per-sample assembly — run_bowtie2.assembly.sh
#!/bin/bash
module load bowtie
module load samtools

echo "Processing $1"

R1trim=${1}_R1_001.flt.fq.gz,/data/vogelh2/231218_Mafa_urine/Data/fastq/${1}_R1_001.flt.fq.gz
R2trim=${1}_R2_001.flt.fq.gz,/data/vogelh2/231218_Mafa_urine/Data/fastq/${1}_R2_001.flt.fq.gz
R3trim=/data/vogelh2/231218_Mafa_urine/Data/fastq/${1}_R2_unmerged.flt.fq.gz
BAM=${1}.combined.assembly.bam
REF=/data/vogelh2/Ref/${1}_megahit_combined

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

```bash
# Calculate contig coverage — pileup3.sh
#!/bin/bash
module load bbtools

pileup.sh in=${1}.combined.assembly.bam out=${1}_pileup3_cov.txt
```

---

### 10. Join coverage to CenoteTaker2 annotations — `cenote2_join.py`

Per-contig coverage from pileup output (`*_pileup3_cov.txt`) was joined to the CenoteTaker2 contig summary files (`*_c_nh2/*_CONTIG_SUMMARY.tsv`) on a per-sample basis, matching by contig name.

```bash
python3 cenote2_join.py
```

Outputs per sample:
- `<sample>_CENOTE2_contig_abundance.tsv` — all contigs with coverage and annotation
- `<sample>_CENOTE2_matched_contigs.tsv` — contigs matched to CenoteTaker2 predictions only

Combined matched contigs across all samples:

```bash
awk 'FNR==1 && NR!=1 { next } { print }' *_CENOTE2_matched_contigs.tsv \
    > SUMMARY_CENOTE2_matched_contigs.tsv
```

---

### 11. Group viral annotations — `group_cenote.py`

Viral organisms identified by CenoteTaker2 were collapsed into taxonomic groups (e.g., Polyomaviridae, Bacteriophage, Cressdnaviricota) and read counts were summed per group per sample.

```bash
python3 group_cenote.py
```

Output per sample: `<sample>_CENOTE2_GROUP.abundance.tsv`

Results were then combined across all samples and sorted:

```bash
awk 'FNR==1 && NR!=1 { next } { print }' *_CENOTE2_GROUP.abundance.tsv \
    > SUMMARY_CENOTE2_GROUP.abundance.tsv
```

The summary file was then split into HSCT recipient and control groups interactively using `awk`, producing `SUMMARY_CENOTE2_GROUP.abundance_HSCT.tsv` and `SUMMARY_CENOTE2_GROUP.abundance_Control.tsv`. Viral read abundance was normalized to reads per million (RPM) total filtered paired-end reads and visualized as stacked bar plots in R using ggplot2.

---

## Dependencies

| Tool | Version | Usage |
|------|---------|-------|
| [FASTP](https://github.com/OpenGene/fastp) | 0.23.2 | Read quality filtering |
| [Bowtie2](https://bowtie-bio.sourceforge.net/bowtie2) | 2.5.3 | Read alignment |
| [SAMtools](http://www.htslib.org/) | 1.21 | BAM processing |
| [BBTools](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/) | 39.06 | Paired-end read repair and contig coverage |
| [bioawk](https://github.com/lh3/bioawk) | 1.0 | SAM/FASTQ parsing |
| [MEGAHIT](https://github.com/voutcn/megahit) | 1.2.9 | De novo assembly |
| [BLAST+](https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html) | 2.15.0 | Nucleotide and protein sequence annotation |
| [DIAMOND](https://github.com/bbuchfink/diamond) | 2.1.12 | Protein sequence annotation |
| [CenoteTaker2](https://github.com/mtisza1/Cenote-Taker2) | 2.1.5 | Viral sequence prediction |
| Python | 3.10 | `cenote2_join.py`, `group_cenote.py` |
| [pandas](https://pandas.pydata.org/) | — | Required by `cenote2_join.py` |
| [R](https://www.r-project.org/) | 4.5.1 | Visualization |
| [ggplot2](https://ggplot2.tidyverse.org/) | 4.0.0 | Stacked bar plots |

---

## Scripts in this folder

| Script | Description |
|--------|-------------|
| `run_fastp.sh` | Quality filters raw paired-end FASTQ files |
| `run_bowtie2.sh` | Aligns filtered reads to host reference genome |
| `extra_nonhost_pairs.sh` | Extracts and repairs non-host read pairs from BAM file |
| `run_megahit_nonhost.sh` | De novo assembles non-host reads into contigs |
| `run_blastn_nonhost.sh` | Annotates contigs against nucleotide database using BLASTN |
| `run_blastx_nonhost.sh` | Annotates contigs against protein database using BLASTX |
| `run_cenote.sh` | Predicts viral sequences from assembled contigs using CenoteTaker2 |
| `bowtie2assembly.sh` | Builds a per-sample Bowtie2 index from MEGAHIT assembly |
| `run_bowtie2.assembly.sh` | Aligns reads from both sequencing runs to per-sample assembly |
| `pileup3.sh` | Calculates per-contig read coverage from BAM file |
| `cenote2_join.py` | Joins contig coverage to CenoteTaker2 annotations per sample |
| `group_cenote.py` | Collapses CenoteTaker2 organisms into taxonomic groups and sums read counts |

---

## Notes
- The host reference genome (`MFAV2`) is the *Macaca fascicularis* assembly GCF_012559485.2 (MFA1912RKSv2), available from NCBI as `GCF_012559485.2_MFA1912RKSv2_genomic.fna`
- All SLURM scripts are designed to run in an HPC environment and use `$SLURM_CPUS_PER_TASK` and `$SLURM_JOB_ID` environment variables
- Update all hardcoded file paths (e.g., `/data/vogelh2/`) to match your environment before running
- The pandas version used by `cenote2_join.py` can be checked with `python3 -c "import pandas; print(pandas.__version__)"`
