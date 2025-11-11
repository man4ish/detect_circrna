## README.md

# Jakobi2016 RNA-Seq Processing Pipeline

## 1. Overview

This repository provides a reproducible workflow for processing paired-end RNA-seq datasets from **Jakobi et al., 2016**, using the **Mus musculus (GRCm38/mm10)** genome.
The pipeline automates all major preprocessing steps — from downloading SRA reads to generating high-quality STAR-aligned BAM files — designed for HPC environments with Grid Engine (`qsub`).

Key features:

* Automated download and trimming of SRA reads
* Quality filtering using **Flexbar**
* rRNA depletion using **Bowtie2**
* Spliced alignment with **STAR** (paired + unpaired reads)
* Conversion to sorted, indexed BAMs
* Soft-clip removal for cleaner alignments

## 2. Workflow Summary

| Step | Tool               | Description                                          |
| ---- | ------------------ | ---------------------------------------------------- |
| 1    | `fasterq-dump`     | Download SRA paired-end FASTQs                       |
| 2    | `flexbar`          | Adapter trimming and quality filtering               |
| 3    | `bowtie2`          | Filter out rRNA contaminants                         |
| 4    | `STAR`             | Align reads to mm10 genome (paired + unmapped mates) |
| 5    | `samtools` + `awk` | Convert, sort, index, and clean SAM/BAMs             |
| 6    | `qsub`             | Job scheduling across HPC cluster nodes              |

---

## 3. Directory Structure

```bash
in/
  └── reads/
      └── jakobi2016_sra_list.txt
out/
  ├── trimmed/
  ├── bowtie2_filtered_all_samples/
  ├── star_mm10_index/
  ├── STAR_output_all_samples/
  └── rRNA_cluster/
logs/
```

---

## 4. Setup

### Conda Environments

The following environments were used:

| Tool     | Conda Environment | Key Packages |
| -------- | ----------------- | ------------ |
| Flexbar  | base              | flexbar      |
| Bowtie2  | bowtie2_env       | bowtie2      |
| STAR     | star_env          | STAR         |
| Samtools | star_env          | samtools     |

Activate environments before running relevant sections.

### Required Modules (HPC)

```bash
module load samtools
```

---

## 5. Input Data

### SRA Accession List

`in/reads/jakobi2016_sra_list.txt`:

```
SRR7881338
SRR7881276
SRR7881275
SRR7881335
SRR7881337
SRR7881336
SRR7881333
SRR7881334
```

---

## 6. Step-by-Step Workflow

### **Step 1 — Download Reads**

```bash
bash wonderdump.sh --split-files --gzip SRRxxxxxxx
```

or batch run:

```bash
while read sra; do
    ./wonderdump.sh --split-files --gzip "$sra" &
    sleep 60
done < in/reads/jakobi2016_sra_list.txt
```

### **Step 2 — Quality Trimming (Flexbar)**

Uses `flexbar` to trim adapters and low-quality bases.

```bash
flexbar -r SRR*_1.fastq -p SRR*_2.fastq -t out/trimmed/SRRXXXX \
        -n 15 -z GZ -m 30 -q TAIL -qt 28 -as "AGATCGGAAGAG"
```

Each sample is submitted as a `qsub` job for parallel execution.

---

### **Step 3 — rRNA Filtering (Bowtie2)**

Removes reads mapping to rRNA reference index.

```bash
bowtie2 -x out/rRNA_cluster/mus-musculus.rRNA \
        -1 <R1> -2 <R2> \
        --no-unal --un-conc-gz out/bowtie2_filtered_all_samples/SRRXXXX_filtered.fastq.gz
```

---

### **Step 4 — Reference Preparation (STAR Index)**

```bash
wget ftp://ftp.ensembl.org/pub/release-90/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.annotation.gtf.gz
```

Then:

```bash
STAR --runThreadN 8 --runMode genomeGenerate \
     --genomeDir out/star_mm10_index \
     --genomeFastaFiles out/Mus_musculus.GRCm38.dna.primary_assembly.fa \
     --sjdbGTFfile out/gencode.vM25.annotation.nodot.gtf \
     --sjdbOverhang 100
```

---

### **Step 5 — STAR Alignment (Paired Reads)**

Aligns filtered paired-end reads to genome:

```bash
STAR --runThreadN 10 \
     --genomeDir out/star_mm10_index \
     --readFilesIn <READ1> <READ2> \
     --readFilesCommand zcat \
     --quantMode GeneCounts \
     --twopassMode Basic
```

---

### **Step 6 — STAR Alignment (Unmapped Mates)**

Runs STAR separately on unmapped mates (`Unmapped.out.mate1.gz` / `mate2.gz`):

```bash
STAR --readFilesIn <mate1> --quantMode GeneCounts --twopassMode Basic
STAR --readFilesIn <mate2> --quantMode GeneCounts --twopassMode Basic
```

---

### **Step 7 — Post-processing (Samtools + awk)**

Removes soft clipping and converts SAM → sorted BAM:

```bash
samtools view -bS Aligned.noS.sam | samtools sort -@ 10 -o Aligned.noS.bam
samtools index Aligned.noS.bam
```

Custom `awk` script removes soft-clipped bases for cleaner junction analysis.

---

## 7. Expected Outputs

Each sample directory (`out/STAR_output_all_samples/SRRxxxxxxx/`) contains:

```
Aligned.noS.bam
Aligned.noS.bam.bai
Chimeric.noS.bam
Chimeric.noS.bam.bai
Unmapped.out.mate1.gz
Unmapped.out.mate2.gz
Log.final.out
ReadsPerGene.out.tab
```

---

## 8. Parallel Execution

All computationally intensive tasks (Flexbar, Bowtie2, STAR) are run via `qsub`:

```bash
qsub -N <job_name> -pe smp <threads> -cwd -V -j y -o logs/
```

This ensures scalability and fault tolerance across HPC clusters.

