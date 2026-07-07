/*
 * Generate a detailed log showing reads removed by all filtering steps
 * Specifically calculates reads removed by blacklist filtering separately from other filters
 */
process BLACKLIST_LOG {
    tag "$meta.id"
    label 'process_low'
    publishDir path: { "${params.outdir}/${params.aligner}/mergedLibrary/filtering_metrics" }, mode: params.publish_dir_mode

    conda (params.enable_conda ? "bioconda::bedtools=2.30.0 bioconda::samtools=1.15.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-8186960447c5cb2faa697666dc1e6d919ad23f3e:3127fcae6b6bdaf8181e21a26ae61231030a9fcb-0':
        'quay.io/biocontainers/mulled-v2-8186960447c5cb2faa697666dc1e6d919ad23f3e:3127fcae6b6bdaf8181e21a26ae61231030a9fcb-0' }"

    input:
    tuple val(meta), path(bam_before), path(bai_before), path(bam_after), path(bai_after)
    path filtered_bed
    path blacklist_bed

    output:
    path "*.filtering.log", emit: log
    path "versions.yml"   , emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Count reads BEFORE any filtering (from MARK_DUPLICATES output)
    READS_BEFORE=\$(samtools view -c ${bam_before})
    
    # Count reads AFTER all filtering (duplicates + blacklist + MAPQ + fragment size)
    READS_AFTER=\$(samtools view -c ${bam_after})
    
    # Count reads that overlap blacklist regions in the BEFORE BAM
    READS_IN_BLACKLIST=\$(samtools view -c -L ${blacklist_bed} ${bam_before})
    
    # Count duplicates marked in BEFORE BAM
    DUPLICATES_MARKED=\$(samtools view -c -f 0x0400 ${bam_before})
    
    # Calculate total reads removed by ALL filters
    TOTAL_REMOVED=\$((READS_BEFORE - READS_AFTER))
    
    # Calculate percentages
    PERCENT_BLACKLIST=\$(awk "BEGIN {printf \\"%.2f\\", (\$READS_IN_BLACKLIST / \$READS_BEFORE) * 100}")
    PERCENT_DUPLICATES=\$(awk "BEGIN {printf \\"%.2f\\", (\$DUPLICATES_MARKED / \$READS_BEFORE) * 100}")
    PERCENT_TOTAL_REMOVED=\$(awk "BEGIN {printf \\"%.2f\\", (\$TOTAL_REMOVED / \$READS_BEFORE) * 100}")
    PERCENT_RETAINED=\$(awk "BEGIN {printf \\"%.2f\\", (\$READS_AFTER / \$READS_BEFORE) * 100}")
    
    # Calculate OTHER filters (MAPQ, fragment size, etc.)
    OTHER_FILTERS=\$((TOTAL_REMOVED - DUPLICATES_MARKED - READS_IN_BLACKLIST))
    PERCENT_OTHER=\$(awk "BEGIN {printf \\"%.2f\\", (\$OTHER_FILTERS / \$READS_BEFORE) * 100}")
    
    # Number of blacklist regions
    NUM_BL_REGIONS=\$(wc -l < ${blacklist_bed})
    
    # Generate log
    cat > ${prefix}.filtering.log <<EOF
========================================================================
BAM FILTERING LOG - Sample: ${meta.id}
========================================================================

Date: \$(date '+%Y-%m-%d %H:%M:%S')
Input BAM (MARK_DUPLICATES):   ${bam_before}
Output BAM (after filtering):  ${bam_after}
Blacklist file:                ${blacklist_bed}

------------------------------------------------------------------------
FILTERING STATISTICS
------------------------------------------------------------------------

Total reads (input):                      \$(printf "%15s" "\$(printf "%'d" \$READS_BEFORE)")

Reads overlapping blacklist regions:      \$(printf "%15s" "\$(printf "%'d" \$READS_IN_BLACKLIST)")  (\${PERCENT_BLACKLIST}%)
Duplicate reads (marked by Picard):       \$(printf "%15s" "\$(printf "%'d" \$DUPLICATES_MARKED)")  (\${PERCENT_DUPLICATES}%)
Reads removed by other filters*:          \$(printf "%15s" "\$(printf "%'d" \$OTHER_FILTERS)")  (\${PERCENT_OTHER}%)
  (*MAPQ < 1, fragment size > 500bp, secondary/supplementary alignments)

------------------------------------------------------------------------
TOTAL FILTERING IMPACT
------------------------------------------------------------------------

Total reads REMOVED (all filters):        \$(printf "%15s" "\$(printf "%'d" \$TOTAL_REMOVED)")  (\${PERCENT_TOTAL_REMOVED}%)
Total reads RETAINED:                     \$(printf "%15s" "\$(printf "%'d" \$READS_AFTER)")  (\${PERCENT_RETAINED}%)

Number of blacklist regions:              \$(printf "%15s" "\$(printf "%'d" \$NUM_BL_REGIONS)")

------------------------------------------------------------------------
NOTE
------------------------------------------------------------------------
- Blacklist count shows reads overlapping blacklist regions
- Duplicate count shows reads marked by Picard MarkDuplicates
- Other filters include: multi-mappers (MAPQ<1), large fragments (>500bp),
  secondary/supplementary alignments, unmapped reads
- Some reads may be counted in multiple categories (e.g., a duplicate
  read in a blacklist region contributes to both counts)

========================================================================
EOF

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version 2>&1 | head -n1 | sed 's/^samtools //')
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
END_VERSIONS
    """
}
