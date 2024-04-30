process roi_coverage {

    output = params.out_dir

    tag { [run, readgroup].join(':') }
    publishDir "${output}/${readgroup}/metrics", mode:"copy"
    label "standard_compute"
    degug = true

    input:
    tuple val(run), val(readgroup), path(bedgraph) 
    path(ref_fasta)
    path(ref_gene)
    path(roi_bed)

    output:
    tuple val(run), val(readgroup), path("${readgroup}_*{metrics,tab}")

    script:
    def inputfile = "${readgroup}_coverage.bedgraph"
    def out_file1 = "${readgroup}_exoncoverage_metrics"
    def out_file2 = "${readgroup}_lowcoverage_metrics"
    def out_file3 = "${readgroup}_hgvslowcoverage_metrics.tab"
    def out_file4 = "${readgroup}_genecoverage_metrics.tab"

    """
    bed4.awk ${roi_bed} | /usr/local/bin/bedtools intersect -wao -split -a - -b ${inputfile} | \\
        bedCoverage.py -s 1 5 10 15 20 25 30 40 50 60 100 200 300 400 500 600 -t 20 | \\
        tee ${out_file1} | \\
        Metrics.py -d -b 2 -n /dev/stdin | \\
        sort -k1,1 -k2n,3n | tee ${out_file2} | \\
        bed2hgvs.py -g ${ref_fasta} -a ${ref_gene} -c -e n /dev/stdin | \\
        sort -k6,7 -k1,1 -k2n,3n | cut -f1-3,5-9 | \\
        sed -e "1ichrom\\tchromStart\\tchromEnd\\tcoverage\\tgene\ttranscript\\thgvsStart\\thgvsEnd" > ${out_file3} && \\
        Metrics.py  -b 3 -r ${out_file1}  > ${out_file4}
    """
}