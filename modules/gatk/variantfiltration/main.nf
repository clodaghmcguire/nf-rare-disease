process variantFiltration {

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
    tuple val(run), val(readgroup), path("${readgroup}.subset.*.filtered.vcf.gz")


    script:
    def input_SNP = "${readgroup}.subset.SNP.vcf.gz"
    def out_file_SNP = "${readgroup}.subset.SNP.filtered.vcf.gz"
    def input_INDEL = "${readgroup}.subset.INDEL.vcf.gz"
    def out_file_INDEL = "${readgroup}.subset.INDEL.filtered.vcf.gz"
    def input_OTHER = "${readgroup}.subset.OTHER.vcf.gz"
    def out_file_OTHER = "${readgroup}.subset.OTHER.filtered.vcf.gz"

    """  
    /usr/bin/java -Xmx4g -Djava.io.tmpdir=${readgroup}_tmp -jar /gatk/gatk.jar VariantFiltration \\
    -R ${ref_fasta} -V ${input_SNP} -O ${out_file_SNP} --filter-expression 'MQRankSum < -12.5' \\
    --filter-name MQRankSum-12.5 --filter-expression 'SOR > 3.0' --filter-name SOR3 --filter-expression 'FS > 60.0' \\
    --filter-name FS60 --filter-expression 'ReadPosRankSum < -8.0' --filter-name ReadPosRankSum-8 \\
    --filter-expression 'QUAL < 30.0' --filter-name QUAL30 --filter-expression 'QD < 2.0' --filter-name QD2 \\
    --filter-expression 'MQ < 40.0' --filter-name MQ40 && \\
    /usr/bin/java -Xmx4g -Djava.io.tmpdir=${readgroup}_tmp -jar /gatk/gatk.jar VariantFiltration \\
    -R ${ref_fasta} -V ${input_INDEL} -O ${out_file_INDEL} --filter-expression 'QD < 2.0' --filter-name QD2 \\
    --filter-expression 'QUAL < 30.0' --filter-name QUAL30 --filter-expression 'FS > 200.0' \\
    --filter-name FS200 --filter-expression 'ReadPosRankSum < -20.0' --filter-name ReadPosRankSum-20 && \\
    /usr/bin/java -Xmx4g -Djava.io.tmpdir=${readgroup}_tmp -jar /gatk/gatk.jar VariantFiltration \\
    -R ${ref_fasta} -V ${input_OTHER} -O ${out_file_OTHER} --filter-expression 'QUAL < 30.0' \\
    --filter-name QUAL30 --filter-expression 'QD < 2.0' --filter-name QD2
    """
}