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

    main:
        // Upload
        uploadMetadata()
        uploadSequenceData()

    //emit:
}
