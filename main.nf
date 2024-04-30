/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include {trimmomatic} from './modules/trimmomatic/main'
include {fastqc} from './modules/fastqc/main'
include {bwa_mem} from './modules/bwa/main'
include {bamSort} from './modules/samtools/sort/main'
include {markduplicates} from './modules/picard/markduplicates/main'
include {recalibration_table} from './modules/gatk/baserecalibrator/main'
include {bam_recalibrate} from './modules/gatk/ApplyBQSR/main'
include {multiple_metrics} from './modules/picard/multiplemetrics/main'
include {collateMetrics} from './modules/collatemetrics/main'
include {verifyBamID} from './modules/verifyBam/main'
include {samtools_coverage} from './modules/samtools/coverage/main'
include {haplotype_caller} from './modules/gatk/haplotypecaller/main'
include {roi_coverage} from './modules/bedtools/coverage/main'
include {selectVariants} from './modules/gatk/selectVariants/main'
include {variantFiltration} from './modules/gatk/variantfiltration/main'
include {sortVCF} from './modules/picard/sort/main'
include {variantStats} from './modules/bcftools/stats/main'
include {readCount} from './modules/readCount/main'
include {exomedepth} from './modules/exomedepth/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    output = params.out_dir    

    Channel.fromFilePairs( params.fastq_files ) \
        | map { readgroup, reads ->
            def (run_name) = reads*.parent.baseName as Set

            tuple( run_name, readgroup, reads )
        } \
        | groupTuple() \
        | map { run, readgroups, reads ->
            tuple( groupKey(run, readgroups.size()), readgroups, reads )
        } \
        | transpose() \
        | set { sample_readgroups }

    gatk_resources = file("${params.ref_gatk_resources}")

    trim1_ch = trimmomatic(sample_readgroups, params.sureselect_adapters)
    trimmomatic.out.trimlogging.collectFile(name: "butt.txt", newLine: true)
    fastqc_ch = fastqc(trim1_ch) 
    bwa_index = file( "${params.ref_bwa_index}" )
    align_ch = bwa_mem(trim1_ch, bwa_index)
    bamSort_ch = bamSort(align_ch)
    dedup_ch = markduplicates(bamSort_ch)
    multiple_metrics_ch = multiple_metrics(dedup_ch, params.ref_fasta)
    samtools_coverage_ch = samtools_coverage(dedup_ch, params.ref_chromsizes, params.target_file)
    coverage_ch = roi_coverage(samtools_coverage_ch, params.ref_fasta, params.ref_gene, params.roi_file)
    bamID_ch = verifyBamID(dedup_ch, params.variation)
    all_metrics_ch = fastqc_ch.join(multiple_metrics_ch, by: 1).map { readgroup, run1, fastqc, run2, multiplemetrics -> tuple(run1, readgroup, fastqc, multiplemetrics) }
    .join(coverage_ch, by: 1).map {readgroup, run1, fastqc, multiplemetrics, run2, coverage -> tuple(run1, readgroup, fastqc, multiplemetrics, coverage) }
    .join(bamID_ch, by: 1).map {readgroup, run1, fastqc, multiplemetrics, coverage, run2, verifyBam -> tuple(run1, readgroup, fastqc, multiplemetrics, coverage, verifyBam) }
    collate_metrics_ch = collateMetrics(all_metrics_ch)
    recalibrationtable_ch = recalibration_table(dedup_ch, params.ref_fasta, params.ref_fasta_fai, params.ref_sequence_dict, gatk_resources)
    join_ch = dedup_ch.join(recalibrationtable_ch, by: 1).map { readgroup, run1, bam, run2, table -> tuple(run1, readgroup, bam, table) }
    recalibrateBam_ch = bam_recalibrate(join_ch, params.ref_fasta, params.ref_fasta_fai, params.ref_sequence_dict)
    haplotype_ch = haplotype_caller(recalibrateBam_ch, params.target_file, params.ref_fasta, params.ref_fasta_fai, gatk_resources, params.ref_sequence_dict)
    variantStats_ch = variantStats(haplotype_ch, params.ref_fasta, params.ref_fasta_fai)
    selectVariants_ch = selectVariants(haplotype_ch, params.ref_fasta, params.ref_fasta_fai, params.ref_sequence_dict)
    variantFiltration_ch = variantFiltration(selectVariants_ch, params.ref_fasta, params.ref_fasta_fai, params.ref_sequence_dict)
    sortVCF_ch = sortVCF(variantFiltration_ch, params.ref_sequence_dict)
    all_bams_ch = dedup_ch.map { run, sample, bam_file ->

        tuple( bam_file )
    }.collect()
    readcount_ch = readCount(all_bams_ch, params.exomedepth, params.ref_fasta, params.ref_fasta_fai, params.normals)
    exomedepth_ch = exomedepth(dedup_ch, readcount_ch, params.exomedepth, params.ref_fasta, params.ref_fasta_fai)

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
