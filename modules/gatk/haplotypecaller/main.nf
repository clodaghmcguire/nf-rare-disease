process haplotype_caller {
    
    output = params.out_dir

    tag { [run, readgroup].join(':') }
    publishDir "${output}/${readgroup}/variants", mode:"copy"
    label "big_mem"
    degug = true

    input:
    tuple val(run), val(readgroup), path(recalibrated_bam)
    path(targets)
    path(ref_fasta)
    path(fai)
    path(ref_gatk_resources)
    path(ref_sequence_dict)

    output:
    tuple val(run), val(readgroup), path("${readgroup}.vcf.gz")


    script:
    def dbsnp = "dbsnp_138.b37.vcf.gz"
    def input ="${readgroup}_recalibrated.bam"
    def out_file = "${readgroup}.vcf.gz"

    """ 
    /usr/bin/java -Xmx14g -Djava.io.tmpdir=${readgroup}_tmp -jar /gatk/gatk.jar HaplotypeCaller \\
    -R ${ref_fasta} -L ${targets} --dbsnp ${dbsnp} -I ${input} -O ${out_file}   
    """
}