process DEEPTOOLS_BIGWIG {
    tag "$meta.id"
    label 'process_high'

    conda (params.enable_conda ? 'bioconda::deeptools=3.5.1' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/deeptools:3.5.1--py_0' :
        'quay.io/biocontainers/deeptools:3.5.1--py_0' }"

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("*.extend.bw"), emit: bigwig
    tuple val(meta), path("*.extend.center.bw"), emit: center_bigwig
    path "versions.yml"                , emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def extend = (meta.single_end && params.fragment_size > 0) ? "--extendReads ${params.fragment_size}" : '--extendReads'
    """

    bamCoverage \\
        --bam $bam \\
        --binSize 1 \\
        --numberOfProcessors $task.cpus \\
        --normalizeUsing CPM \\
        $extend \\
        --maxFragmentLength 10000 \\
        -o ${prefix}.extend.bw

    bamCoverage \\
        --bam $bam \\
        --binSize 1 \\
        --numberOfProcessors $task.cpus \\
        --normalizeUsing CPM \\
        $extend \\
        --maxFragmentLength 10000 \\
        --centerReads \\
        -o ${prefix}.extend.center.bw


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deeptools: \$(computeMatrix --version | sed -e "s/computeMatrix //g")
    END_VERSIONS
    """
}
