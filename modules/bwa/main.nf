process bwa_mem {

    output = params.out_dir
    
    tag { [run, readgroup].join(':') }
    publishDir "${output}/${readgroup}/align", mode:"copy"
    label "bwa_mem"
    degug = true

    input:
    tuple val(run), val(readgroup), path(fastq_files)
    val(previous_std_out)
    path bwa_indexes

    output:
    tuple val(run), val(readgroup), path("${readgroup}.bam")

    script:
    def idxbase = bwa_indexes[0].baseName
    def out_file = "${readgroup}.bam"
    def config_settings ="${params.bwa_mem_config}"
    def cpu ="${params.resources.bwa_mem.cpu}"

    """
    /usr/local/pipeline/bwa-0.7.16a/bwa mem \\
        -R '@RG\\tID:${readgroup}\\tLB:TAS\\tPL:illumina\\tPU:runname\\tSM:${run}' -t ${cpu} \\
        ${config_settings} \\
        "${idxbase}"\\
        "${fastq_files[0]}" \\
        "${fastq_files[2]}" | \\
        /usr/local/bin/samtools view -b -o ${out_file} -
    """
}
