process assessSimpleQuality {
    publishDir "${params.outdir}/${sample}", pattern: "${sample}_quality.csv", mode: "copy"
    tag { sample }

    input:
    tuple val(sample), path(fasta), path(filteredbam), path(compositebam), path(compositebambai), path(krakenreport)

    output:
    path ("${sample}_quality.csv")

    // Script to put all wanted sample values together
    script:
    // If running kraken2, setup to calculate initial % viral
    if ( params.kraken_db ) {
        checkKraken     = true
        krakenHeader    = "initial_viral_percent,"
    } else {
        checkKraken     = false
        krakenHeader    = ""
    }
    """
    # For kraken2 if it is available/input
    if $checkKraken; then
        KRAKEN_IDS=`{ grep -P "\\s+619591\\s+" $krakenreport | cut -d\$'\t' -f 3 || :; }`
        TOTAL_READS=`samtools view -F 0x04 -c $compositebam`
        INITIAL_VIRAL_PER=`awk -v id_count="\$KRAKEN_IDS" -v reads="\$TOTAL_READS" 'BEGIN {if (id_count=="") {print 0} else {print (id_count*2)/reads * 100}}'`,
    else
        INITIAL_VIRAL_PER=""
    fi

    # Header
    printf "sample,${krakenHeader}num_reads_mapped,mean_sequencing_depth,num_consensus_n,genome_completeness\n" > ${sample}_quality.csv

    # Stats we want
    NUM_READS=`samtools view -F 0x04 -c $filteredbam`
    AVG_COV=`samtools depth -a $filteredbam | awk '{c++;s+=\$3} END {if (c=='0') {print "NA"} else {print s/c}}'`
    N_COUNT_PLUS_PROP=`seqtk comp $fasta | awk '{ OFS = "," }{ x=\$3+\$4+\$5+\$6;y=\$2; if (y=='0') {print "NA,NA"} else print y-x,1-((y-x)/y) }'`

    # Create file
    printf "${sample},\$INITIAL_VIRAL_PER\$NUM_READS,\$AVG_COV,\$N_COUNT_PLUS_PROP\n" >> ${sample}_quality.csv
    """
}
process runNextclade {
    // Running nextclade 2

    input:
    path (fastas)

    output:
    path ("nc_all_data.csv")

    script:
    """
    cat *.fa > all.fasta
    nextclade dataset get --name hMPXV --output-dir dataset
    nextclade run -D dataset/ -p ${params.pcr_csv} -t nextclade.tsv all.fasta

    cat nextclade.tsv | sed -e '1s|\\.|_|g' -e '1s|\\t|\\tnc_|g' -e 's|,|;|g' -e 's|\\t|,|g' > nc_all_data.csv
    """
}
process concatQuality {
    publishDir "${params.outdir}/", pattern: "overall_sample_quality.csv", mode: "copy"

    input:
    path(csvs)
    path(nextclade_csv)

    output:
    path ("overall_sample_quality.csv")

    script:
    """
    awk '(NR == 1) || (FNR > 1)' *quality.csv > simple_sample_quality.csv
    csvtk join -f 'sample;seqName' simple_sample_quality.csv $nextclade_csv > overall_sample_quality.csv
    """
}
