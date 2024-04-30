process exomedepth {
    output = params.out_dir

    tag { [run, readgroup].join(':') }
    publishDir "${output}/${readgroup}/sv", mode:"copy"
    label "big_mem"
    degug = true
    
    input:
    tuple val(run), val(readgroup), path(bam)
    path(readCount)
    path(exomedepth_targets)
    path(ref_fasta)
    path(ref_fai)

    output:
    tuple val(run), val(readgroup), path("${readgroup}.exomedepth*")

    script:
    def sample_name = "${readgroup}"
    def input = "${readgroup}_dedupmk.bam".replaceAll(/-/, ".")
    def threshold = "${params.exomedepth_threshold}"
    def out_file = "${readgroup}.exomedepth.pdf"
    """
    exomeDepth.R v1 ${out_file} ${exomedepth_targets}:llgp2 \\
    ${readCount} ${input}:${sample_name}:${threshold} ${params.exomedepth_plot_pon} ${params.exomedepth_plot_batch}
    ls
    """
}