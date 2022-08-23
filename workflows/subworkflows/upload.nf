// Modules to include
include {
    uploadMetadata;
    uploadSequenceData
} from '../../modules/upload_modules.nf'

// Workflow for quality metrics and other checks
workflow upload {
    take:
        ch_consensus_fasta          // channel: [ val(sampleID), path(fasta) ]
        ch_dehosted_fastq           // channel: [ val(sampleID), path(Read1), path(Read2), val(bool) ]
        ch_sequence_metrics         // channel: path(overall_sample_quality.csv)
        ch_metadata                 // channel: path(metadata_csv) or [] if none

    main:
        // Get config setup
        ch_upload_conf = file( params.upload_config )

        // Upload
        uploadSequenceData(
            ch_consensus_fasta.collect { it[1] },
            ch_dehosted_fastq.collect { it[1] },
            ch_dehosted_fastq.collect { it[2] },
            ch_upload_conf,
            ch_metadata
        )
        uploadMetadata(
            ch_sequence_metrics,
            ch_upload_conf,
            uploadSequenceData.out // Added so that uploads are one at a time and samples don't conflict
        )

    //emit:
}
