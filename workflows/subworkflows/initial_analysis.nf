// Modules to include
include {
    runFastQC;
} from '../../modules/initial_analysis_modules.nf'

// Workflow for quality metrics and other checks
workflow initial_analysis {
    take:
        ch_paired_fastqs        // channel: [ val(sampleID), path(Read1), path(Read2), val(bool) ]

    main:
        // Add any other pre-analysis checks here
        runFastQC( ch_paired_fastqs )

    //emit:
}
