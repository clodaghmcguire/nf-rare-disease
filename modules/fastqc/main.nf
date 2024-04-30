process fastqc {

    output = params.out_dir

    tag { [run, readgroup].join(':') }
    publishDir "${output}/${readgroup}/metrics", mode:"copy"
    label "standard_compute"

    input:
    tuple val(run), val(readgroup), path(fastq_files)
    val(previous_std_out)

    output:
    tuple val(run), val(readgroup), path("${readgroup}*")

    script:

    """

    /usr/local/pipeline/FastQC/fastqc -o . "${fastq_files[0]}" "${fastq_files[2]}"

    """
}