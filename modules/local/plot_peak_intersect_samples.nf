/*
 * Plot peak intersections for individual samples
 * This creates an UpSet plot showing overlaps between all individual sample peaks
 */
process PLOT_PEAK_INTERSECT_SAMPLES {
    tag "$antibody"
    label 'process_medium'

    conda (params.enable_conda ? "conda-forge::biopython conda-forge::r-optparse=1.7.1 conda-forge::r-upsetr=1.4.0 bioconda::bedtools=2.30.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-2f48cc59b03027e31ead6d383fe1b8057785dd24:5d182f583f4696f4c4d9f3be93052811b383341f-0':
        'quay.io/biocontainers/mulled-v2-2f48cc59b03027e31ead6d383fe1b8057785dd24:5d182f583f4696f4c4d9f3be93052811b383341f-0' }"

    input:
    tuple val(antibody), path(peaks)

    output:
    tuple val(antibody), path("*.intersect.txt"), emit: intersect_txt
    tuple val(antibody), path("*.intersect.plot.pdf"), emit: pdf
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${antibody}"
    def peak_type = params.narrow_peak ? 'narrowPeak' : 'broadPeak'
    def mergecols = params.narrow_peak ? (2..10).join(',') : (2..9).join(',')
    def collapsecols = params.narrow_peak ? (['collapse']*9).join(',') : (['collapse']*8).join(',')
    def expandparam = params.narrow_peak ? '--is_narrow_peak' : ''
    """
    # Sort and merge all individual sample peaks
    sort -T '.' -k1,1 -k2,2n ${peaks.collect{it.toString()}.sort().join(' ')} \\
        | mergeBed -c $mergecols -o $collapsecols > ${prefix}.merged.txt
  
    # Create boolean matrix to show which samples have peaks in each region
    macs2_merged_expand.py \\
        ${prefix}.merged.txt \\
        ${peaks.collect{it.toString()}.sort().join(',').replaceAll("_peaks.${peak_type}","").replaceAll("_chr[^,]*","")} \\
        ${prefix}.boolean.txt \\
        --min_replicates 1 \\
        $expandparam

    # Generate UpSet plot showing peak overlaps between all samples
    if [ -s ${prefix}.boolean.intersect.txt ]; then
        plot_peak_intersect.r \\
            -i ${prefix}.boolean.intersect.txt \\
            -o ${prefix}.samples.intersect.plot.pdf
    else
        echo "No intersect data available for individual samples" > ${prefix}.samples.intersect.plot.pdf
    fi

    # Keep the intersect file for output
    cp ${prefix}.boolean.intersect.txt ${prefix}.samples.intersect.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        bedtools: \$(bedtools --version | sed 's/bedtools v//')
    END_VERSIONS
    """
}
