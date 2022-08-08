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

# Reset and Track
mv .nextflow.log artifacts/nextflow.log
ls results
rm -rf results work/ .nextflow*

echo "Done"
