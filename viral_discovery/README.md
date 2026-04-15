# Identification of disease associated polyomaviruses from cynomolgus macaques undergoing HSCT

## Background
This repository contains scripts used for de novo assembly and annotation of viral sequences from Illumina short read sequencing of RCA urine DNA from cynomolgus macaques.

## Workflow

1. **Quality filter raw reads** using FASTP to remove low-quality bases and adapter sequences
2. **Align reads to a reference genome** containing host and/or viral sequences using Bowtie2, which generates a BAM file
3. **Filter out host reads** from the BAM file using SAMtools and bioawk to retain non-host read pairs, then repair paired-end reads using BBTools
4. **Assemble non-host reads into contigs** using MEGAHIT de novo assembler
5. **Annotate contigs** by aligning against nucleotide and protein databases using BLASTN and BLASTX, and by predicting viral sequences using CenoteTaker2

---

## Dependencies

| Tool | Version | Usage |
|------|---------|-------|
| [FASTP](https://github.com/OpenGene/fastp) | 0.23.2 | Read quality filtering |
| [Bowtie2](https://bowtie-bio.sourceforge.net/bowtie2) | 2.5.3 | Read alignment |
| [SAMtools](http://www.htslib.org/) | 1.23 | BAM processing |
| [BBTools](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/) | 39.06 | Paired-end read repair |
| [bioawk](https://github.com/lh3/bioawk) | 1.0 | SAM/FASTQ parsing |
| [MEGAHIT](https://github.com/voutcn/megahit) | 1.2.9 | De novo assembly |
| [BLAST+](https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html) | 2.15.0+ | Sequence alignment/annotation |
| [CenoteTaker2](https://github.com/mtisza1/Cenote-Taker2) | 2.1.5 | Viral sequence prediction |

---

## Scripts

### 1. Quality filtering — `run_fastp.sh`

Filters raw paired-end FASTQ files. Takes a sample name as input (without `_R1_001.fastq.gz` suffix).

```bash
#!/bin/bash
module load fastp

fastp -w $SLURM_CPUS_PER_TASK \
    --in1 ${1}_R1_001.fastq.gz \
    --in2 ${1}_R2_001.fastq.gz \
    --out1 ${1}_R1_001.flt.fq.gz \
    --out2 ${1}_R2_001.flt.fq.gz
```

**Usage:**
```bash
sbatch run_fastp.sh <sample_name>
```

---

### 2. Reference index — Bowtie2 index build

Build a Bowtie2 index from the reference FASTA before running alignment. Only needs to be run once per reference.

> **Reference genome:** `MFAV2_Mafapyv2.fasta` is derived from the *Macaca fascicularis* reference assembly [GCF_012559485.2 (MFA1912RKSv2)](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_012559485.2/), downloaded as `GCF_012559485.2_MFA1912RKSv2_genomic.fna` from NCBI.

```bash
module load bowtie
bowtie2-build --threads 2 MFAV2.fasta MFAV2
```

---

### 3. Alignment — `run_bowtie2.sh`

Aligns filtered reads to the reference genome and outputs a sorted, indexed BAM file.

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

**Usage:**
```bash
sbatch run_bowtie2.sh <sample_name>
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

**Usage:**
```bash
sbatch extra_nonhost_pairs.sh <sample_name>
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

**Usage:**
```bash
sbatch run_megahit_nonhost.sh <sample_name>
```

---

### 6. Nucleotide annotation — `run_blastn_nonhost.sh`

Aligns assembled contigs against a nucleotide reference database using BLASTN.

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

**Usage:**
```bash
sbatch run_blastn_nonhost.sh <query.fasta> <blast_db>
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

**Usage:**
```bash
sbatch run_blastx_nonhost.sh <query.fasta> <blast_db>
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

**Usage:**
```bash
sbatch run_cenote.sh <sample_name>
```

---

## Notes
- The host reference genome (`MFAV2`) is the *Macaca fascicularis* assembly GCF_012559485.2 (MFA1912RKSv2), available from NCBI as `GCF_012559485.2_MFA1912RKSv2_genomic.fna`
- All scripts are designed to run in a SLURM HPC environment and use `$SLURM_CPUS_PER_TASK` and `$SLURM_JOB_ID` environment variables
- Update file paths (e.g., `/data/vogelh2/Ref/`, `/data/vogelh2/Software/`) to match your environment before running
- Scripts should be submitted via `sbatch` with appropriate resource allocations for your cluster


