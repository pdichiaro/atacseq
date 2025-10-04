include { SAMTOOLS_SORT            } from '../../modules/nf-core/modules/samtools/sort/main'
include { SAMTOOLS_INDEX           } from '../../modules/nf-core/modules/samtools/index/main'
include { SAMTOOLS_FLAGSTAT        } from '../../modules/nf-core/modules/samtools/flagstat/main'
include { DEEPTOOLS_ALIGNMENTSIEVE } from '../../modules/nf-core/deeptools/alignmentsieve/main'

workflow BAM_SHIFT_READS {
    take:
    ch_bam_bai                   // channel: [ val(meta), [ bam ], [bai] ]
    minFragmentLength            // val: minimum fragment length
    maxFragmentLength            // val: maximum fragment length

    main:
    def ch_versions = Channel.empty()

    //
    // Shift reads using deepTools alignmentSieve with --ATACshift
    // This applies the standard +4/-5 bp offset for Tn5 binding sites
    // Fragment length filtering is applied via minFragmentLength and maxFragmentLength
    //
    DEEPTOOLS_ALIGNMENTSIEVE (
        ch_bam_bai
    )
    ch_versions = ch_versions.mix(DEEPTOOLS_ALIGNMENTSIEVE.out.versions)

    //
    // Sort shifted reads
    //
    SAMTOOLS_SORT (
        DEEPTOOLS_ALIGNMENTSIEVE.out.bam,
        [[],[]]  // No reference fasta needed for coordinate sorting
    )
    ch_versions = ch_versions.mix(SAMTOOLS_SORT.out.versions)

    //
    // Index shifted and sorted BAM files
    //
    SAMTOOLS_INDEX (
        SAMTOOLS_SORT.out.bam
    )
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions)

    //
    // Generate flagstat for QC metrics
    //
    SAMTOOLS_FLAGSTAT (
        SAMTOOLS_SORT.out.bam.join(SAMTOOLS_INDEX.out.bai, by: [0])
    )
    ch_versions = ch_versions.mix(SAMTOOLS_FLAGSTAT.out.versions)

    emit:
    bam      = SAMTOOLS_SORT.out.bam                // channel: [ val(meta), [ bam ] ]
    bai      = SAMTOOLS_INDEX.out.bai               // channel: [ val(meta), [ bai ] ]
    csi      = SAMTOOLS_INDEX.out.csi               // channel: [ val(meta), [ csi ] ]
    flagstat = SAMTOOLS_FLAGSTAT.out.flagstat       // channel: [ val(meta), [ flagstat ] ]
    versions = ch_versions                          // channel: [ versions.yml ]
}
