class NextflowInitialize {

    // Validate and initialize workflow properly
    public static void validateParams(workflow, params, log) {
        // Print help and exit
        if ( params.help ) {
            log.info """
            Usage:
                Analyze:
                nextflow run ${workflow.manifest.name} --directory [Path/to/PairedFastqs] --human_ref [Path/to/HumanRefFasta] [workflow-options]

                Help and Exit:
                nextflow run ${workflow.manifest.name} --help

            Profiles:
                -profile [conda,nml]            Configuration profile to use. Can use both if separated by a comma
                    standard                      Standard profile, needs all tool dependencies installed on PATH and runs locally
                    conda                         Utilize conda to control tool and dependency installation (uses mamba for environment install)
                    nml                           NML specific profile to take advantage of NML cluster resources

            Mandatory Arguments:
                --directory [path]              Path to input directory containing paired fastq files
                --human_ref [file]              Path to human reference genome file in fasta format (ex. hg38.fa)

            Optional Unset Arguments:
                --cache [path]                      (With conda profile only) Path to cache directory containing conda environments
                --composite_bwa_index_dir [path]    Path to directory containing BWA indexed composite genome to skip indexing step
                                                      Make sure it matches the Human and MPXV reference sequences or pipeline will not work

            Optional Preset Arguments:
                --outdir [str]                  String name for the output results directory (Default: 'results')
                --mpx_ref [path]                Path to monkeypox reference genome (Default: 'data/NC_063383.fasta')
                --mpx_ref_id [str]              String name of reference monkeypox contig ID to keep during host removal (Default: 'NC_063383.1')
                --pcr_csv [path]                Path to CSV file containing qPCR primers for nextclade to intake (Default: 'data/nml_primers.csv')

            Optional Run Kraken2 Arguments:
                --kraken_db [path]              Path to directory containing Kraken2 database
                                                  Runs Kraken2 on the generated dehosted fastq files when used
                --kraken_viral_id [str]         String Kraken2 taxonomic identifier to use for the percent viral reads calculation. 
                                                  Default is the Integer ID for the Viruses domain - 10239
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
