## Release 0.2.0 - 2022/08/06

#### Argument Changes:
- Added argument `kraken_viral_id` with the default Kraken2 ID of `10239`

#### Output Changes:
- All samples CSV:
    - Renamed the following columns:
        - proportion_basecalled  -> genome_completeness
        - total_n_count          -> num_consensus_n
        - avg_read_depth         -> mean_sequencing_depth
    - Added `median_sequencing_depth` column
    - Added `percent_viral_reads` column for when Kraken2 is ran (not there otherwise)

#### Other Changes:
- Added changelog
- Added testing reports for host removal metrics
- Modified workflows:
    - Initial analysis workflow added (for fastqc and any other pre-checks we add)
    - The host removal workflow emits a channel called `kraken_results`
        - If kraken was run, it is populated with the sample and its kraken report
        - If no kraken was run, then it is populated with just the samples IDs
    - Added channel info to all sub-workflows
        - Shows what each input/output looks like
- Setup proper resource usage and more info on how a user can set their resources
- Added license
- Readme and help command updated

## Release 0.1.0 - Unreleased
Initial release of pipeline
