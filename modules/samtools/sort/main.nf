process bamSort {
     
    output = params.out_dir
    
    tag { [run, readgroup].join(':') }
    publishDir "${output}/${readgroup}/align", mode:"copy"
    label "big_mem"
    degug = true

    input:
    tuple val(run), val(readgroup), path(bam_file)

    output:
    tuple val(run), val(readgroup), path("${readgroup}_sorted.bam*")

    script:
    def out_file = "${readgroup}_sorted.bam"

    """
    /usr/local/bin/samtools sort -@ 2 -O bam -T tmpfile_${readgroup} -m 2G "${bam_file}" > ${out_file} \\
        && /usr/local/bin/samtools index ${out_file}
    """
}