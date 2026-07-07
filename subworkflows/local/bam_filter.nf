//
// Filter BAM file using standard ATAC-seq approach (removes multimappers)
//

include { BAM_FILTER as BAM_FILTER_PROCESS } from '../../modules/local/bam_filter'
include { BAM_SORT_SAMTOOLS  } from '../nf-core/bam_sort_samtools'
 
workflow BAM_FILTER {
    take:
    ch_bam_bai   // channel: [ val(meta), [ bam ] ]
    ch_bed       // channel: [ bed ]

    main:

    ch_versions = Channel.empty()

    //
    // STANDARD ATAC-seq FILTERING
    // Removes multimappers (NH:i:1 filter), duplicates, blacklist regions,
    // unmapped reads, and filters by fragment size
    //
    
    BAM_FILTER_PROCESS( ch_bam_bai, ch_bed )
    ch_versions = ch_versions.mix(BAM_FILTER_PROCESS.out.versions.first())

    // Use filtered BAM for downstream processing
    BAM_SORT_SAMTOOLS(BAM_FILTER_PROCESS.out.bam)
    ch_versions = ch_versions.mix(BAM_SORT_SAMTOOLS.out.versions.first())

    emit:

    bam      = BAM_SORT_SAMTOOLS.out.bam      // channel: [ val(meta), [ bam ] ]
    bai      = BAM_SORT_SAMTOOLS.out.bai      // channel: [ val(meta), [ bai ] ]
    stats    = BAM_SORT_SAMTOOLS.out.stats    // channel: [ val(meta), [ stats ] ]
    flagstat = BAM_SORT_SAMTOOLS.out.flagstat // channel: [ val(meta), [ flagstat ] ]
    idxstats = BAM_SORT_SAMTOOLS.out.idxstats // channel: [ val(meta), [ idxstats ] ]

    versions = ch_versions                    // channel: [ versions.yml ]
}

