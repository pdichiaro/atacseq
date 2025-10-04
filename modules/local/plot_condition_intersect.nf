/*
 * Plot peak intersections across CONDITIONS (not individual replicates)
 * This creates an UpSet plot showing overlaps between condition-level consensus peaks
 * Only runs if there are >= 2 conditions for an antibody
 */
process PLOT_CONDITION_INTERSECT {
    tag "$antibody"
    label 'process_medium'

    conda (params.enable_conda ? "conda-forge::biopython conda-forge::r-optparse=1.7.1 conda-forge::r-upsetr=1.4.0 bioconda::bedtools=2.30.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-2f48cc59b03027e31ead6d383fe1b8057785dd24:5d182f583f4696f4c4d9f3be93052811b383341f-0':
        'quay.io/biocontainers/mulled-v2-2f48cc59b03027e31ead6d383fe1b8057785dd24:5d182f583f4696f4c4d9f3be93052811b383341f-0' }"

    input:
    tuple val(antibody), val(condition_ids), path(peaks)

    output:
    tuple val(antibody), path("*.conditions.intersect.txt"), emit: intersect_txt
    tuple val(antibody), path("*.conditions.intersect.plot.pdf"), emit: pdf
    path "versions.yml", emit: versions

    when:
    (task.ext.when == null || task.ext.when) && condition_ids.size() >= 2

    script:
    def prefix = task.ext.prefix ?: "${antibody}"
    def peak_type = params.narrow_peak ? 'narrowPeak' : 'broadPeak'
    def mergecols = params.narrow_peak ? (2..10).join(',') : (2..9).join(',')
    def collapsecols = params.narrow_peak ? (['collapse']*9).join(',') : (['collapse']*8).join(',')
    def expandparam = params.narrow_peak ? '--is_narrow_peak' : ''
    // AWK print statement depends on peak type (narrowPeak=10 cols, broadPeak=9 cols)
    def awk_print = params.narrow_peak ? 
        '{print $1,$2,$3,cond"_"$4,$5,$6,$7,$8,$9,$10}' : 
        '{print $1,$2,$3,cond"_"$4,$5,$6,$7,$8,$9}'
    // Create pairs of (condition, peak_file) and sort by condition
    def condition_peak_pairs = [condition_ids, peaks].transpose().sort { a, b -> a[0] <=> b[0] }
    def sorted_conditions = condition_peak_pairs.collect { it[0] }
    def sorted_peaks = condition_peak_pairs.collect { it[1].toString() }
    
    """
    # Debug: Print what we received
    echo "DEBUG: Antibody = ${antibody}"
    echo "DEBUG: Number of conditions = ${condition_ids.size()}"
    echo "DEBUG: Condition IDs (sorted) = ${sorted_conditions.join(', ')}"
    echo "DEBUG: Peak files (sorted):"
    ls -lh ${sorted_peaks.join(' ')}
    
    # Use the condition IDs directly - they are already sorted to match peak files
    PEAK_FILES=(${sorted_peaks.join(' ')})
    CONDITION_NAMES=(${sorted_conditions.collect{"'$it'"}.join(' ')})
    
    echo "DEBUG: PEAK_FILES array: \${PEAK_FILES[@]}"
    echo "DEBUG: CONDITION_NAMES array: \${CONDITION_NAMES[@]}"
    
    # Tag each peak file with its corresponding condition name
    for i in "\${!PEAK_FILES[@]}"; do
        PEAK_FILE="\${PEAK_FILES[\$i]}"
        COND_NAME="\${CONDITION_NAMES[\$i]}"
        # Add condition name as prefix to peak name (column 4)
        # Number of columns depends on peak type (narrowPeak=10, broadPeak=9)
        awk -v cond="\$COND_NAME" 'BEGIN{OFS="\\t"} $awk_print' "\$PEAK_FILE" > "\${COND_NAME}.tagged.bed"
    done
    
    # Sort and merge all tagged peaks
    cat *.tagged.bed | sort -k1,1 -k2,2n | \\
        mergeBed -c $mergecols -o $collapsecols > ${prefix}.merged.txt
    
    # Use custom script to generate intersection data directly
    # Create comma-separated list from bash array
    CONDITION_NAMES_CSV=\$(IFS=, ; echo "\${CONDITION_NAMES[*]}")
    
    plot_condition_intersect_custom.py \\
        ${prefix}.merged.txt \\
        "\$CONDITION_NAMES_CSV" \\
        ${prefix}.conditions.intersect.txt

    # Generate UpSet plot showing peak overlaps between CONDITIONS
    if [ -s ${prefix}.conditions.intersect.txt ]; then
        plot_peak_intersect.r \\
            -i ${prefix}.conditions.intersect.txt \\
            -o ${prefix}.conditions.intersect.plot.pdf
    else
        echo "No intersect data available for conditions" > ${prefix}.conditions.intersect.plot.pdf
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        bedtools: \$(bedtools --version | sed 's/bedtools v//')
    END_VERSIONS
    """
}
