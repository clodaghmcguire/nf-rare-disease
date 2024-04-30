process collateMetrics {
    output = params.out_dir

    tag { [run, readgroup].join(':') }
    publishDir "${output}/${readgroup}/metrics", mode:"copy"
    label "standard_compute"
    degug = true

    input:
    tuple val(run), val(readgroup), path(fastqc), path(multiplemetrics), path(coverage), path(verifyBamID) 

    output:
    tuple val(run), val(readgroup), path("${readgroup}.json")

    script:
    def outfile = "${readgroup}.json"
    """
    collateMetrics.py ${outfile}
    """
}