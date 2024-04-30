process selectVariants {

    output = params.out_dir

    tag { [run, readgroup].join(':') }
    publishDir "${output}/${readgroup}/variants", mode:"copy"
    label "big_mem"
    degug = true

    input:
    tuple val(run), val(readgroup), path(vcf)
    path(ref_fasta)
    path(fai)
    path(ref_sequence_dict)

    output:
    tuple val(run), val(readgroup), path("${readgroup}.subset*")


    script:
    def input = "${readgroup}.vcf.gz"
    def out_file1 = "${readgroup}.subset.OTHER.vcf.gz"
    def out_file2 = "${readgroup}.subset.INDEL.vcf.gz"
    def out_file3 = "${readgroup}.subset.SNP.vcf.gz"

    """ 
    /usr/bin/java -jar /gatk/gatk.jar IndexFeatureFile -F ${input}
    /usr/bin/java -Xmx4g -Djava.io.tmpdir=${readgroup}_tmp -jar /gatk/gatk.jar SelectVariants \\
    -R ${ref_fasta} -V ${input} -O ${out_file1} -xl-select-type INDEL -xl-select-type SNP && \\
    /usr/bin/java -Xmx4g -Djava.io.tmpdir=${readgroup}_tmp -jar /gatk/gatk.jar SelectVariants \\
    -R ${ref_fasta} -V ${input} -O ${out_file2} -select-type INDEL && /usr/bin/java -Xmx4g \\
    -Djava.io.tmpdir=${readgroup}_tmp -jar /gatk/gatk.jar SelectVariants -R ${ref_fasta} \\
    -V ${input} -O ${out_file3} -select-type SNP 
    """
}