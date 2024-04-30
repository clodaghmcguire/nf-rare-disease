process sortVCF {

    output = params.out_dir
    
    tag { [run, readgroup].join(':') }
    publishDir "${output}/${readgroup}/variants", mode:"copy"
    label "big_mem"
    degug = true

    input:
    tuple val(run), val(readgroup), path(vcf)
    path(ref_sequence_dict)

    output:
    tuple val(run), val(readgroup), path("${readgroup}.filtered.vcf.gz")


    script:
    def input_SNP = "${readgroup}.subset.SNP.filtered.vcf.gz"
    def input_INDEL = "${readgroup}.subset.INDEL.filtered.vcf.gz"
    def input_OTHER = "${readgroup}.subset.OTHER.filtered.vcf.gz"
    def out_file = "${readgroup}.filtered.vcf.gz"

    """  
    /usr/bin/java -XX:ParallelGCThreads=8 -Xmx8g -jar /usr/local/pipeline/picard-tools/picard.jar SortVcf \\
    TMP_DIR=${readgroup}_tmp VALIDATION_STRINGENCY=SILENT MAX_RECORDS_IN_RAM=100000 \\
    SEQUENCE_DICTIONARY=${ref_sequence_dict} O=${out_file} I=${input_INDEL} I=${input_OTHER} I=${input_SNP}
    """
}