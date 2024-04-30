process markduplicates {

    output = params.out_dir

    tag { [run, readgroup].join(':') }
    publishDir "${output}/${readgroup}/metrics", mode:"copy"
    label "big_mem"
    degug = true

    input:
    tuple val(run), val(readgroup), path(bam_file)

    output:
    tuple val(run), val(readgroup), path("${readgroup}_dedupmk.bam*")
    

    script:
    def out_file = "${readgroup}_dedupmk.bam"
    def out_metrics = "${readgroup}.metrics.DUPLICATION"
    """
    /usr/bin/java -XX:ParallelGCThreads=8 -Xmx8g -jar /usr/local/pipeline/picard-tools/picard.jar MarkDuplicates \\
    TMP_DIR=${readgroup}_picard_tmp VALIDATION_STRINGENCY=SILENT MAX_RECORDS_IN_RAM=100000 \\
    METRICS_FILE=${out_metrics} INPUT=${bam_file[0]} OUTPUT=${out_file} && /usr/local/bin/samtools index ${out_file}
    """

}