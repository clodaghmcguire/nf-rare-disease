process verifyBamID {
    
    output = params.out_dir

    tag { [run, readgroup].join(':') }
    publishDir "${output}/${readgroup}/metrics", mode:"copy"
    label "big_mem"
    degug = true

    input:
    tuple val(run), val(readgroup), path(dedupmk_bam)
    path(variation)

    output:
    tuple val(run), val(readgroup), path("${readgroup}.metrics.variation*")

    script:
    def input ="${readgroup}_dedupmk.bam"
    def out_file = "${readgroup}.metrics.variation"
    """
    /usr/local/bin/verifyBamID --noPhoneHome --vcf ${variation} --bam ${input} --out ${out_file} --self --maxDepth 1000 --precise --ignoreRG
    """
}