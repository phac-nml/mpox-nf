process uploadSequenceData {

    label 'upload'

    input:
    path(fastas)
    path(read1s)
    path(read2s)
    path(config)
    path(metadata)

    output:
    path("done_sequence_upload.txt")

    script:
    """
    mkdir -p irida_fastqs
    mv $read1s irida_fastqs
    mv $read2s irida_fastqs
    irida_upload_csv_generator.py --directory irida_fastqs/ --samplesheet $metadata --fastq
    irida-uploader --config $config -d irida_fastqs

    mkdir -p irida_fastas
    mv $fastas irida_fastas
    irida_upload_csv_generator.py --directory irida_fastas/ --samplesheet $metadata --fasta
    irida-uploader --config $config -d irida_fastas --upload_mode=assemblies
    touch done_sequence_upload.txt
    """
}
process uploadMetadata {
    publishDir "${params.outdir}/", pattern: "metadata_upload_status.csv", mode: "copy"

    label 'upload'

    input:
    path(metadata_csv)
    path(config)
    path(prev_upload_file)

    output:
    path("metadata_upload_status.csv")

    script:
    """
    upload.py --config $config --metadata_csv $metadata_csv
    """
}
