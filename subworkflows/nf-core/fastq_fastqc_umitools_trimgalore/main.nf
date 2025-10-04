//
// Read QC, UMI extraction and trimming
//

include { FASTQC           } from '../../../modules/nf-core/modules/fastqc/main'
include { UMITOOLS_EXTRACT } from '../../../modules/nf-core/modules/umitools/extract/main'
include { TRIMGALORE       } from '../../../modules/nf-core/modules/trimgalore/main'

//
// Function that parses TrimGalore log output file to get total number of reads after trimming
//
def getTrimGaloreReadsAfterFiltering(log_file) {
    def total_reads = 0
    def filtered_reads = 0
    log_file.eachLine { line ->
        def total_reads_matcher = line =~ /([\d\.]+)\ssequences processed in total/
        def filtered_reads_matcher = line =~ /shorter than the length cutoff[^:]+:\s([\d\.]+)/
        if (total_reads_matcher) total_reads = total_reads_matcher[0][1].toFloat()
        if (filtered_reads_matcher) filtered_reads = filtered_reads_matcher[0][1].toFloat()
    }
    return total_reads - filtered_reads
}

workflow FASTQ_FASTQC_UMITOOLS_TRIMGALORE {
    take:
    reads             // channel: [ val(meta), [ reads ] ]
    skip_fastqc       // boolean: true/false
    with_umi          // boolean: true/false
    skip_umi_extract  // boolean: true/false
    skip_trimming     // boolean: true/false
    umi_discard_read  // integer: 0, 1 or 2
    min_trimmed_reads // integer: > 0

    main:
    ch_versions = channel.empty()

    fastqc_html = channel.empty()
    fastqc_zip  = channel.empty()
    if (!skip_fastqc) {
        FASTQC (reads)
        fastqc_html = FASTQC.out.html
        fastqc_zip  = FASTQC.out.zip
        ch_versions = ch_versions.mix(FASTQC.out.versions.first())
    }

    umi_reads = reads
    umi_log   = channel.empty()
    if (with_umi && !skip_umi_extract) {
            UMITOOLS_EXTRACT (reads)
            umi_reads   = UMITOOLS_EXTRACT.out.reads
            umi_log     = UMITOOLS_EXTRACT.out.log
            ch_versions = ch_versions.mix(UMITOOLS_EXTRACT.out.versions.first())

            // Discard R1 / R2 if required
            if (umi_discard_read in [1,2]) {
                UMITOOLS_EXTRACT
                    .out
                    .reads
                    .map { meta_in, reads_in ->
                        meta_in.single_end ? [ meta_in, reads_in ] : [ meta_in + ['single_end': true], reads_in[umi_discard_read % 2] ]
                    }
                    .set { umi_reads }
            }
    }

    trim_reads      = umi_reads
    trim_unpaired   = channel.empty()
    trim_html       = channel.empty()
    trim_zip        = channel.empty()
    trim_log        = channel.empty()
    trim_read_count = channel.empty()
    if (!skip_trimming) {
        TRIMGALORE (umi_reads)
        trim_unpaired = TRIMGALORE.out.unpaired
        trim_html     = TRIMGALORE.out.html
        trim_zip      = TRIMGALORE.out.zip
        trim_log      = TRIMGALORE.out.log
        ch_versions   = ch_versions.mix(TRIMGALORE.out.versions.first())

        //
        // Filter FastQ files based on minimum trimmed read count after adapter trimming
        //
        TRIMGALORE
            .out
            .reads
            .join(trim_log, remainder: true)
            .map { meta2, reads2, trim_log2 ->
                if (trim_log2) {
                    def num_reads = getTrimGaloreReadsAfterFiltering(meta2.single_end ? trim_log2 : trim_log2[-1])
                    [ meta2, reads2, num_reads ]
                } else {
                    [ meta2, reads2, min_trimmed_reads.toFloat() + 1 ]
                }
            }
            .set { ch_num_trimmed_reads }

        ch_num_trimmed_reads
            .filter { _meta3, _reads3, num_reads3 -> num_reads3 >= min_trimmed_reads.toFloat() }
            .map { meta3, reads3, _num_reads3 -> [ meta3, reads3 ] }
            .set { trim_reads }

        ch_num_trimmed_reads
            .map { meta, _reads, num_reads -> [ meta, num_reads ] }
            .set { trim_read_count }
    }

    emit:
    reads = trim_reads // channel: [ val(meta), [ reads ] ]

    fastqc_html        // channel: [ val(meta), [ html ] ]
    fastqc_zip         // channel: [ val(meta), [ zip ] ]

    umi_log            // channel: [ val(meta), [ log ] ]

    trim_unpaired      // channel: [ val(meta), [ reads ] ]
    trim_html          // channel: [ val(meta), [ html ] ]
    trim_zip           // channel: [ val(meta), [ zip ] ]
    trim_log           // channel: [ val(meta), [ txt ] ]
    trim_read_count    // channel: [ val(meta), val(count) ]

    versions = ch_versions.ifEmpty(null) // channel: [ versions.yml ]
}
