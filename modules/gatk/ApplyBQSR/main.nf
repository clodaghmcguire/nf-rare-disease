process bam_recalibrate {
    output = params.out_dir

    tag { [run, readgroup].join(':') }
    publishDir "${output}/${readgroup}/align", mode:"copy"
    label "big_mem"
    degug = true

    input:
    tuple val(run), val(readgroup), path(bam_file), path(recalibration)
    path(ref_fasta)
    path(fai)
    path(ref_sequence_dict)

    output:
    tuple val(run), val(readgroup), path("${readgroup}_recalibrated.bam*")

    script:
    def input = "${readgroup}_dedupmk.bam"
    def recalibrationtable = "${readgroup}.dupemk.recalibrationtable"
    def out_file = "${readgroup}_recalibrated.bam"

    """
    /usr/bin/java -Xmx4g -Djava.io.tmpdir=${readgroup}_tmp -jar /gatk/gatk.jar ApplyBQSR \\
    -R ${ref_fasta} -bqsr ${recalibrationtable} -I ${input} -O ${out_file} && /usr/bin/samtools index ${out_file}
    """
}
