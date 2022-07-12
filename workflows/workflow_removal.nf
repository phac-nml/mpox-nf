// Modules to include
include {
    generateHostRemovedFastq;
    runKraken2
} from '../modules/removal_modules.nf'

// Workflow to generate host removed reads and get read quality metrics
workflow host_removal {
    // Taking in the sample, filtered bam
    take:
        ch_filtered_bam

    main:
        generateHostRemovedFastq( ch_filtered_bam )

        if ( params.kraken_db ) {
            runKraken2( generateHostRemovedFastq.out )
        }
    //emit:
}
