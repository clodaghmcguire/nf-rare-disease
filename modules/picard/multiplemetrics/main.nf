process multiple_metrics {

    output = params.out_dir

    tag { [run, readgroup].join(':') }
    publishDir "${output}/${readgroup}/metrics", mode:"copy"
    label "standard_compute"
    degug = true

    input:
    tuple val(run), val(readgroup), path(clipped_bam) 
    path(ref_fasta)

    output:
    tuple val(run), val(readgroup), path("${readgroup}_*{metrics,pdf}")


    script:
    def out_file1 = "${readgroup}_multiple.metrics"
    def out_file2 = "${readgroup}_base_quality_metrics"
    def input ="${readgroup}_dedupmk.bam"

    """
    /usr/bin/java -XX:ParallelGCThreads=4 -Xmx20g -jar /usr/local/pipeline/picard-tools/picard.jar CollectMultipleMetrics TMP_DIR=${readgroup}_picard_tmp \\
        VALIDATION_STRINGENCY=SILENT MAX_RECORDS_IN_RAM=100000 INPUT=${input} OUTPUT=${out_file1} REFERENCE_SEQUENCE=${ref_fasta} \\
        PROGRAM=QualityScoreDistribution PROGRAM=MeanQualityByCycle PROGRAM=CollectInsertSizeMetrics && \\
        qxx.py ${out_file1}.quality_distribution_metrics > ${out_file2}
    """
}