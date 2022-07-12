#!/usr/bin/env nextflow
/*
    Paired illumina monkeypox shotgun sequencing pipeline
*/
nextflow.enable.dsl = 2

// Validate Inputs
NextflowInitialize.validateParams(workflow, params, log)

// Include any needed workflows or modules
include { mpx_main } from './workflows/mpx_main.nf'

// Execute main workflow
workflow {
    mpx_main()
}
