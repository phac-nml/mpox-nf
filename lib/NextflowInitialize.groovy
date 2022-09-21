class NextflowInitialize {

    // Validate and initialize workflow properly
    public static void validateParams(workflow, params, log) {
        // Print help and exit
        if ( params.help ) {
            log.info """
            Usage:
                Analyze:
                nextflow run phac-nml/${workflow.manifest.name} --directory [Path/to/PairedFastqs] --human_ref [Path/to/HumanRefFasta] [workflow-options]

                Help and Exit:
                nextflow run phac-nml/${workflow.manifest.name} --help

            Profiles:
                -profile [standard,conda,nml]   Configuration profile to use. Can use both if separated by a comma
                    standard                      Standard profile, needs all tool dependencies installed on PATH and runs locally
                    conda                         Utilize conda to control tool and dependency installation (mamba for environment install)
                    nml                           NML specific profile to take advantage of NML cluster resources

            Mandatory Arguments:
                --directory [path]              Path to input directory containing paired fastq files
                --human_ref [file]              Path to human reference genome file in fasta format (ex. hg38.fa)

            Optional Unset Arguments:
                --cache [path]                      (With conda profile only) Path to cache directory containing conda environments
                --composite_bwa_index_dir [path]    Path to directory containing BWA indexed composite genome to skip indexing step
                                                      Make sure it matches the Human and MPXV reference sequences or pipeline will not work
                --metadata_csv [file]               Add metadata to output CSV. Requires column called 'sample'
                                                      If uploading to IRIDA also needs columns 'project_id' and 'sequencing_date'

            Optional Preset Arguments:
                --outdir [str]                  String name for the output results directory (Default: 'results')
                --mpx_ref [path]                Path to monkeypox reference genome (Default: 'data/NC_063383.fasta')
                --mpx_ref_id [str]              String name of reference monkeypox contig ID to keep during host removal (Default: 'NC_063383.1')
                --pcr_csv [file]                Path to CSV file containing qPCR primers for nextclade to intake (Default: 'data/nml_primers.csv')

            Optional Run Kraken2 Arguments:
                --kraken_db [path]              Path to directory containing Kraken2 database
                                                  Runs Kraken2 on the generated dehosted fastq files when used
                --kraken_viral_id [str]         String Kraken2 taxonomic identifier to use for the percent viral reads calculation. 
                                                  Default is the Integer ID for the Viruses domain - 10239

            Optional Upload to IRIDA Arguments:
                --upload_config [file]          Path to IRIDA uploader config file (example in README)
                                                  Also need to provide a metadata csv file with columns 'sample', 'project_id' and 'sequencing_date'
            """.stripIndent()
            System.exit(0)
        }

        // Check inputs
        if ( !params.directory ) {
            log.error "Please specify an input directory with --directory <Path/to/PairedFastqs>"
            System.exit(1)
        }

        if ( !params.human_ref ) {
            log.error "Please specify a human reference genome with '--human_ref <REF>'"
            System.exit(1)
        }

        if ( params.upload_config && !params.metadata_csv ) {
            log.error "IRIDA Upload requested but no metadata CSV file given. To upload, re-run with a valid metadata file adding on --metadata_csv <FILE.csv>"
            System.exit(1)
        }

        // Check CSV Metadata Headers for correct values based on if upload is wanted or not
        if ( params.metadata_csv ) {
            def firstLine = new File(params.metadata_csv).readLines().get(0)
            def entries = firstLine.split(',')
            if ( params.upload_config ) {
                if ( !entries.contains('sample') | !entries.contains('project_id') | !entries.contains('sequencing_date') ) {
                    log.error """
                    IRIDA Upload requested but input metadata CSV file does not contain the needed headers
                      Requires at minimum: ['sample', 'project_id', 'sequencing_date']
                      Got: $entries
                    """.stripIndent()
                }
            } else {
                if ( !entries.contains('sample') ) {
                    log.error "Input metadata CSV file does not contain the 'sample' header"
                }
            }
        }

        // Print Starting Niceness
        log.info """
╔═══════════════════════════╗
║   ╔════════════════════╗  ║
║   ║    MonkeyPox-NF    ║  ║
║   ╚════════════════════╝  ║
╚═══════════════════════════╝
        """.stripIndent()

        // Track Input Values
        log.info "Params:"
        params.each { key, val ->
            if ( val ) {
                log.info "  $key = $val"
            }
        }

        // Profiles
        log.info "Profile(s): $workflow.profile"
    }
}
