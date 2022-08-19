#!/bin/bash
set -eo pipefail

# Create a Cache Dir
mkdir -p conda_cache_dir

# Get partial human ref genome
gunzip $PWD/.github/ci-data/partial_hg38_ref.fa.gz

# Run Illumina Pipeline
nextflow run ./main.nf \
    -profile conda,test \
    --cache ./conda_cache_dir \
    --directory $PWD/.github/ci-data/ \
    --human_ref $PWD/.github/ci-data/partial_hg38_ref.fa

# Check that the output is as expected
# 1. Num Reads
READS=`awk -F, '$1 == "Sample1" {print $2}' ./results/overall_sample_quality.csv`
if [[ "$READS" != "26280" ]]; then 
    echo "Incorrect output: Number of reads mapped"
    echo "  Expected: 26280, Got: $READS"
    exit 1
fi
# 2. Genome Completeness
COMPLETENESS=`awk -F, '$1 == "Sample1" {print $6}' ./results/overall_sample_quality.csv | xargs printf "%.2f"`
if [[ "$COMPLETENESS" != "0.85" ]]; then 
    echo "Incorrect output: Genome Completeness"
    echo "  Expected: 0.85, Got: $COMPLETENESS"
    exit 1
fi

# Reset and Track
mv .nextflow.log artifacts/nextflow.log
rm -rf results work/ .nextflow*

echo "Passed Tests"
