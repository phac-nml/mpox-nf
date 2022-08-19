#!/bin/bash
set -eo pipefail

# Test is run second so no need to re-create envs or unzip files

# Run Illumina Pipeline
nextflow run ./main.nf \
    -profile conda,test \
    --cache ./conda_cache_dir \
    --directory $PWD/.github/ci-data/ \
    --human_ref $PWD/.github/ci-data/partial_hg38_ref.fa
    --metadata_csv $PWD/.github/ci-data/test_metadata.csv

# Check that the output is as expected
# 1. Num Reads
if [[ `awk -F, '$1 == "Sample1" {print $2}' ./results/overall_sample_quality.csv` != 26280 ]]; then 
    echo "Incorrect output: Number of reads mapped"
    exit 1
fi
# 2. Genome Completeness
if [[ `awk -F, '$1 == "Sample1" {print $6}' ./results/overall_sample_quality.csv` != 0.845519 ]]; then 
    echo "Incorrect output: Genome Completeness"
    exit 1
fi
# 3. Check metadata put in correctly
if [[ `awk -F, '$1 == "Sample1" {print $7}' ./results/overall_sample_quality.csv` != 231 ]]; then 
    echo "Incorrect output: Metadata not parsed correctly"
    exit 1
fi

# Reset and Track
mv .nextflow.log artifacts/nextflow2.log
rm -rf results work/ .nextflow*

echo "Passed Tests"
