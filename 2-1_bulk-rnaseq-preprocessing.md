# Pre-processing of bulk RNA-seq data

Processing steps from raw FASTQ files to count matrices. 

- Nextflow-based *'nf-core/rnaseq (v3.5)' pipeline is used. 

- Chosen strategy: STAR alignment + Salmon quantification.

- Human genome reference: GRCh38 + GENCODE V39

## Environment preparation

```bash

# Create and activate conda env
conda create --name nf python=3.7 nf-core nextflow
conda activate nf
# Download offline files for the pipeline
nf-core download -c singularity -r 3.5 \
-o /faststorage/project/THOR/tools/nf-core-rnaseq-3.5 \
nf-core/rnaseq

```

## Sequencing batch 1 (234 samples)

```bash

# Check that all files exist
while IFS="," read -r col1 col2 col3 col4
do
    echo "$col1"
        if test -f "$col2"; then
            echo "EXISTS: $col2"
        fi
        if test -f "$col3"; then
            echo "EXISTS: $col3"
        fi
done < <(tail -n +2 samplesheet_kd_r1.csv)
# Run nf-core/rnaseq pipeline
nextflow run /faststorage/project/THOR/rna-seq/nf-core-rnaseq-3.5/workflow \
--input samplesheet_kd_r1.csv \
--aligner star_salmon \
--gencode \
--fasta /faststorage/project/THOR/rna-seq/reference/gencode_r39_GRCh38.p13/GRCh38.p13.genome.fa \
--gtf /faststorage/project/THOR/rna-seq/reference/gencode_r39_GRCh38.p13/gencode.v39.annotation.gtf \
--star_index /faststorage/project/THOR/rna-seq/reference/gencode_r39_GRCh38.p13/index/star \
--salmon_index /faststorage/project/THOR/rna-seq/reference/gencode_r39_GRCh38.p13/index/salmon \
--rseqc_modules 'bam_stat,inner_distance,infer_experiment,junction_annotation,junction_saturation,read_distribution,read_duplication' \
-profile singularity \
-resume

```

## Sequencing batch 2 (99 samples)

```bash

# Check that all files exist
while IFS="," read -r col1 col2 col3 col4
do
    echo "$col1"
        if test -f "$col2"; then
            echo "EXISTS: $col2"
        fi
        if test -f "$col3"; then
            echo "EXISTS: $col3"
        fi
done < <(tail -n +2 samplesheet_kd_r2.csv)
# Run nf-core/rnaseq pipeline
nextflow run /faststorage/project/THOR/rna-seq/nf-core-rnaseq-3.5/workflow \
--input samplesheet_kd_r2.csv \
--aligner star_salmon \
--gencode \
--fasta /faststorage/project/THOR/rna-seq/reference/gencode_r39_GRCh38.p13/GRCh38.p13.genome.fa \
--gtf /faststorage/project/THOR/rna-seq/reference/gencode_r39_GRCh38.p13/gencode.v39.annotation.gtf \
--star_index /faststorage/project/THOR/rna-seq/reference/gencode_r39_GRCh38.p13/index/star \
--salmon_index /faststorage/project/THOR/rna-seq/reference/gencode_r39_GRCh38.p13/index/salmon \
--rseqc_modules 'bam_stat,inner_distance,infer_experiment,junction_annotation,junction_saturation,read_distribution,read_duplication' \
-profile singularity \
-resume

```

## Sequencing batch 3 (replication of several experiments with target gene knockdown)

```bash

# Check that all files exist
while IFS="," read -r col1 col2 col3 col4
do
    echo "$col1"
        if test -f "$col2"; then
            echo "EXISTS: $col2"
        fi
        if test -f "$col3"; then
            echo "EXISTS: $col3"
        fi
done < <(tail -n +2 samplesheet_kd_r3.csv)
# Run nf-core/rnaseq pipeline
nextflow run /faststorage/project/THOR/rna-seq/nf-core-rnaseq-3.5/workflow \
--input samplesheet_kd_r3.csv \
--aligner star_salmon \
--gencode \
--fasta /faststorage/project/THOR/rna-seq/reference/gencode_r39_GRCh38.p13/GRCh38.p13.genome.fa \
--gtf /faststorage/project/THOR/rna-seq/reference/gencode_r39_GRCh38.p13/gencode.v39.annotation.gtf \
--star_index /faststorage/project/THOR/rna-seq/reference/gencode_r39_GRCh38.p13/index/star \
--salmon_index /faststorage/project/THOR/rna-seq/reference/gencode_r39_GRCh38.p13/index/salmon \
--rseqc_modules 'bam_stat,inner_distance,infer_experiment,junction_annotation,junction_saturation,read_distribution,read_duplication' \
-profile singularity \
-resume

```

## Baseline experiments: effects of cholesterol loading and stretch


```bash

# Check that all files exist
while IFS="," read -r col1 col2 col3 col4
do
    echo "$col1"
        if test -f "$col2"; then
            echo "EXISTS: $col2"
        fi
        if test -f "$col3"; then
            echo "EXISTS: $col3"
        fi
done < <(tail -n +2 samplesheet.csv)
# Run nf-core/rnaseq pipeline
cd /faststorage/project/THOR/anton/rna-seq/thor_kd_zero
nextflow run /faststorage/project/THOR/rna-seq/nf-core-rnaseq-3.5/workflow \
--input samplesheet.csv \
--aligner star_salmon \
--gencode \
--fasta /faststorage/project/THOR/rna-seq/reference/gencode_r39_GRCh38.p13/GRCh38.p13.genome.fa \
--gtf /faststorage/project/THOR/rna-seq/reference/gencode_r39_GRCh38.p13/gencode.v39.annotation.gtf \
--star_index /faststorage/project/THOR/rna-seq/reference/gencode_r39_GRCh38.p13/index/star \
--salmon_index /faststorage/project/THOR/rna-seq/reference/gencode_r39_GRCh38.p13/index/salmon \
--rseqc_modules 'bam_stat,inner_distance,infer_experiment,junction_annotation,junction_saturation,read_distribution,read_duplication' \
-profile singularity

```
