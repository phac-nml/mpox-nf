#!/bin/bash
set -eo pipefail

# Test is run second so no need to re-create envs or unzip files

# Run Illumina Pipeline
nextflow run ./main.nf \
    -profile conda,test \
    --cache ./conda_cache_dir \
    --directory $PWD/.github/ci-data/ \
    --human_ref $PWD/.github/ci-data/partial_hg38_ref.fa \
    --metadata_csv $PWD/.github/ci-data/test_metadata.csv

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
# 3. Check metadata put in correctly
CELL=`awk -F, '$1 == "Sample1" {print $7}' ./results/overall_sample_quality.csv`
if [[ "$CELL" != "231" ]]; then 
    echo "Incorrect output: Metadata not parsed correctly"
    echo "  Expected: 231, Got: $CELL"
    exit 1
fi
# 4. Check that the upload creation script works
mkdir -p dehosted_fastqs
mv results/*/*_dehosted* dehosted_fastqs/
python bin/irida_upload_csv_generator.py -d dehosted_fastqs -s $PWD/.github/ci-data/test_metadata.csv --fastq
if ! diff -q dehosted_fastqs/SampleList.csv $PWD/.github/output/SampleList.csv &>/dev/null; then 
    echo "Upload SampleList Values are different"
    exit 1
fi

# Reset and Track
mv .nextflow.log artifacts/nextflow2.log
rm -rf ./results ./work/ .nextflow* ./dehosted_fastqs

echo "Passed Tests"
