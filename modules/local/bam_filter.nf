/*
 * Filter BAM file using standard ATAC-seq filtering approach (nf-core compatible)
 * 
 * Primary/Secondary alignment handling:
 * - ALWAYS removes secondary (0x100) and supplementary (0x800) alignments
 * - Only primary alignments are processed downstream
 * 
 * Multi-mapper handling (aligned with nf-core/atacseq):
 * - Bowtie2 with NO -k flag: Reports single best primary alignment (MAPQ=0 if ambiguous)
 * - Bowtie2 with -k 100: Reports primary (MAPQ=0) + 99 secondary alignments (flag 0x100)
 * - Secondary alignments are ALWAYS filtered out (via -F 0x100)
 * - keep_multi_map = false (default): MAPQ >= 1 filter removes primary with MAPQ=0
 * - keep_multi_map = true: Keeps primary alignment even with MAPQ=0 (ambiguous mapping)
 * 
 * Fragment size filtering (PE reads only):
 * - Bowtie2 uses -X 1000 (permissive search limit, ensures no valid pairs are missed)
 * - BAM_FILTER applies params.insert_size (default: 500bp) for biological quality filtering
 * - Removes fragments >500bp to filter chimeras, artifacts, and unusually long fragments
 * - Recommended values: 400-500bp (narrow peaks), 500-600bp (broad marks)
 * 
 * Other filters:
 * - Removes duplicates (unless keep_dups = true)
 * - Removes blacklisted regions
 * - Ensures proper paired-end reads (when applicable)
 */
process BAM_FILTER {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::bedtools=2.30.0 bioconda::samtools=1.15.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-8186960447c5cb2faa697666dc1e6d919ad23f3e:3127fcae6b6bdaf8181e21a26ae61231030a9fcb-0':
        'quay.io/biocontainers/mulled-v2-8186960447c5cb2faa697666dc1e6d919ad23f3e:3127fcae6b6bdaf8181e21a26ae61231030a9fcb-0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    path bed

    output:
    tuple val(meta), path("*.filter2.bam"), emit: bam
    path "versions.yml"           , emit: versions

    script:

    def prefix           = task.ext.prefix ?: "${meta.id}"
    // ALWAYS exclude secondary (0x100) and supplementary (0x800) alignments
    // This ensures only primary alignments are processed, even when keep_multi_map=true
    def base_filter      = '-F 0x0100 -F 0x0800'  // Exclude secondary + supplementary
    def filter_params    = meta.single_end ? 
        "${base_filter} -F 0x004" : 
        "${base_filter} -F 0x004 -F 0x0008 -f 0x001 -f 0x002"  // proper pair selection
    def dup_params       = params.keep_dups ? '' : '-F 0x0400'
    def blacklist_params = params.blacklist ? "-L $bed" : ''
    def max_frag = params.insert_size ? params.insert_size.toInteger() : 500

    if(params.keep_multi_map == false ) {

        """
        # General filtering of the bam
        samtools view \\
            $filter_params \\
            $dup_params \\
            $blacklist_params \\
            -b $bam > ${prefix}.filter1.bam
        
        # Remove multi-mappers (MAPQ < 1) and filter by fragment size
        # -q 1: Keep only reads with MAPQ >= 1 (removes primary alignments with MAPQ=0)
        # Note: Secondary alignments already removed by -F 0x100 in filter_params
        # awk: Filter pairs with insert size <= max_frag (default: params.insert_size = 500bp)

        samtools view -q 1 -h ${prefix}.filter1.bam | \\
            awk -v var="$max_frag" '{if(substr(\$0,1,1)=="@" || ((\$9>=0?\$9:-\$9)<=var)) print \$0}' | \\
            samtools view -b > ${prefix}.filter2.bam

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            samtools: \$(samtools --version 2>&1 | head -n1 | sed 's/^samtools //')
        END_VERSIONS
        """
        
    } else {

        """
        # General filtering of the bam
        samtools view \\
            $filter_params \\
            $dup_params \\
            $blacklist_params \\
            -b $bam > ${prefix}.filter1.bam

        # Filter pairs by insert size (keep_multi_map=true mode)
        # Keep primary alignments with MAPQ=0 (multi-mappers)
        # Note: Secondary alignments already removed by -F 0x100 in filter_params
        # awk: Filter pairs with insert size <= max_frag (default: params.insert_size = 500bp)
        
        samtools view -h ${prefix}.filter1.bam | \\
            awk -v var="$max_frag" '{if(substr(\$0,1,1)=="@" || ((\$9>=0?\$9:-\$9)<=var)) print \$0}' | \\
            samtools view -b > ${prefix}.filter2.bam

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            samtools: \$(samtools --version 2>&1 | head -n1 | sed 's/^samtools //')
        END_VERSIONS
        """

    }


}
