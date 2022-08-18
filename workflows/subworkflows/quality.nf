// Modules to include
include {
    assessSimpleQuality;
    runNextclade;
    concatQuality
} from '../../modules/quality_modules.nf'

// Workflow for quality metrics and other checks
workflow assess_quality {
    take:
        ch_consensus_fasta          // channel: [ val(sampleID), path(fasta) ]
        ch_filtered_bam             // channel: [ val(sampleID), path(filteredBam) ]
        ch_composite_bam            // channel: [ val(sampleID), path(compositeBam), path(compositeBamBai) ]
        ch_kraken_results           // channel: [ val(sampleID), path(krakenReport) ]
        ch_metadata                 // channel: path(metadata_csv) or [] if none

    main:
        // Join fasta with bam and get some info
        ch_qual = ch_consensus_fasta
                    .join( ch_filtered_bam, by:0 )
                    .join( ch_composite_bam, by:0 )
                    .join( ch_kraken_results, by:0 )
        assessSimpleQuality( ch_qual )

        // Nextclade with fastas
        ch_fasta_only = ch_consensus_fasta.multiMap { it ->
                            names: it[0]
                            files: it[1]
                        }
        runNextclade( ch_fasta_only.files.collect() )

        // Singular CSV output
        concatQuality(
            assessSimpleQuality.out.collect(),
            runNextclade.out,
            ch_metadata
        )

    emit:
    sequence_metrics = concatQuality.out        // channel: path(overall_sample_quality.csv)
}
