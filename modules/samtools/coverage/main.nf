process samtools_coverage {

    output = params.out_dir

    tag { [run, readgroup].join(':') }
    publishDir "${output}/${readgroup}/align", mode:"copy"
    label "standard_compute"
    degug = true

    input:
    tuple val(run), val(readgroup), path(dedupmk_bam)
    path(ref_chromsizes)
    path(targets)

    output:
    tuple val(run), val(readgroup), path("${readgroup}_coverage.bedgraph")

    script:
    def out_file = "${readgroup}_coverage.bedgraph"
    def input ="${readgroup}_dedupmk.bam"

    """
    /usr/local/bin/samtools view -b -q 10 ${input} | /usr/local/bin/bedtools bamtobed -i - | awk '(\$3>\$2) { print \$0 }' | \\
    /usr/local/bin/bedtools genomecov -bg -i - -g ${ref_chromsizes} | /usr/local/bin/bedtools intersect -u -a - -b ${targets} > ${out_file}

    """
}