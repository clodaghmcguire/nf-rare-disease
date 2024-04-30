process recalibration_table {
    output = params.out_dir

    tag { [run, readgroup].join(':') }
    publishDir "${output}/${readgroup}/align", mode:"copy"
    label "big_mem"
    degug = true

    input:
    tuple val(run), val(readgroup), path(bam_file)
    path(ref_fasta)
    path(fai)
    path(ref_sequence_dict)
    path(ref_gatk_resources)

    output:
    tuple val(run), val(readgroup), path("${readgroup}.dupemk.recalibrationtable")

    script:
    def input = "${readgroup}_dedupmk.bam"
    def out_file = "${readgroup}.dupemk.recalibrationtable"
    def dbsnp = "dbsnp_138.b37.vcf.gz"
    def hapmap = "hapmap_3.3.b37.vcf.gz"
    def omni = "1000G_omni2.5.b37.vcf.gz"
    def high_confidence_snps = "1000G_phase1.snps.high_confidence.b37.vcf.gz"

    """
    /usr/bin/java -Xmx4g -Djava.io.tmpdir=${readgroup}_tmp -jar /gatk/gatk.jar BaseRecalibrator \\
    -R ${ref_fasta} --known-sites ${dbsnp} --known-sites ${hapmap} \\
    --known-sites ${omni} --known-sites ${high_confidence_snps} -I ${input} -O ${out_file}
    """
}