process readCount {
    output = params.out_dir

    publishDir "${output}", mode:"copy"
    label "big_mem"
    degug = true

    input:
    path(bam_files)
    path(exomedepth_targets)
    path(ref_fasta)
    path(ref_fai)
    path(panel_normals)

    output:
    path("readCount.RData")

    script:
    def out_file = "readCount.RData"
    """
    echo ${bam_files}
    readCount.R ${out_file} ${ref_fasta} ${exomedepth_targets} ${bam_files} ${panel_normals}
    """
}