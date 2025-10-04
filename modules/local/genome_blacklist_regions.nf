/*
 * Prepare genome intervals for filtering by removing blacklist regions and optionally chrM
 * 
 * Filters applied:
 * 1. Blacklist regions (if keep_blacklist = false): Removes ENCODE blacklist regions
 * 2. Mitochondrial chromosome (if keep_mito = false): Removes chrM/MT/Mt based on mito_name
 * 
 * Output: BED file with genomic regions to INCLUDE in analysis
 */
process GENOME_BLACKLIST_REGIONS {
    tag "$sizes"

    conda (params.enable_conda ? "bioconda::bedtools=2.30.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedtools:2.30.0--hc088bd4_0':
        'quay.io/biocontainers/bedtools:2.30.0--hc088bd4_0' }"

    input:
    path sizes
    path blacklist
    val mito_name
    val keep_mito
    val keep_blacklist

    output:
    path '*.bed'       , emit: bed
    path "versions.yml", emit: versions

    script:
    def file_out = "${sizes.simpleName}.include_regions.bed"
    
    // Prepare mito filter for awk
    def mito_filter = (mito_name && !keep_mito) ? 
        "| awk '\$1 !~ /${mito_name}/' " : 
        ''
    
    if (blacklist && !keep_blacklist) {
        """
        # Remove blacklist regions AND optionally chrM
        sortBed -i $blacklist -g $sizes | \\
            complementBed -i stdin -g $sizes $mito_filter > $file_out

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
        END_VERSIONS
        """
    } else {
        """
        # Create full genome BED, optionally excluding chrM
        awk '{print \$1, "0" , \$2}' OFS='\\t' $sizes $mito_filter > $file_out

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
        END_VERSIONS
        """
    }
}
