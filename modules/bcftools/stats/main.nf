process variantStats {
    
    output = params.out_dir

    tag { [run, readgroup].join(':') }
    publishDir "${output}/${readgroup}/metrics", mode:"copy"
    label "big_mem"
    degug = true

    input:
    tuple val(run), val(readgroup), path(vcf)
    path(ref_fasta)
    path(ref_fai)

    output:
    tuple val(run), val(readgroup), path("${readgroup}.metrics.variation*")

    script:
    def input ="${readgroup}.vcf.gz"
    def out_file = "${readgroup}.metrics.variation"
    """
    /usr/local/bin/bcftools stats -F ${ref_fasta} -s - ${input} > ${out_file}
    """
}