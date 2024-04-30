process trimmomatic{

    output = params.out_dir

    tag { [run, readgroup].join(':') }
    publishDir "${output}/${readgroup}/fastq", mode:"copy"
    label "standard_compute"

    input:
    tuple val(run), val(readgroup), path(reads)
    path adapters

    output:
    tuple val(run), val(readgroup), path("${readgroup}_*{filtered,metrics}.{fq.gz,READS}")
    stdout emit: trimlogging

    script:
    def out_fastq_files = "${readgroup}_1_filtered.fq.gz ${readgroup}_1_rubbish_filtered.fq.gz ${readgroup}_2_filtered.fq.gz ${readgroup}_2_rubbish_filtered.fq.gz"
    def out_metrics_file = "${readgroup}_metrics.READS"
    def (r1, r2) = reads

    """
    java -Xmx4G -jar /usr/local/pipeline/Trimmomatic-0.32/trimmomatic-0.32.jar \\
        PE -trimlog /dev/stdout  \\
        -threads 4 \\
        ${r1} \\
        ${r2} \\
        ${out_fastq_files} \\
        ILLUMINACLIP:${adapters}:2:30:10:1:true SLIDINGWINDOW:4:15 LEADING:3 TRAILING:3 MINLEN:36 \\
        |parseTrimlog.py > ${out_metrics_file}
    """
}