process PREPARE_HOMER_ANNOTATION {
    tag "homer_annotation"
    label 'process_low'

    conda (params.enable_conda ? 'bioconda::homer=4.11' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/homer:4.11--pl526hc9558a2_3' :
        'quay.io/biocontainers/homer:4.11--pl526hc9558a2_3' }"

    input:
    path gtf
    path fasta

    output:
    path "gene_annotation.txt", emit: annotation
    path "versions.yml"       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def VERSION = '4.11'
    """
    # Create a BED file from GTF (gene level)
    awk 'BEGIN{OFS="\\t"} \$3=="gene" {
        split(\$9, a, ";");
        gene_id = "";
        gene_name = "";
        gene_type = "";
        for (i in a) {
            if (a[i] ~ /gene_id/) {
                gsub(/.*gene_id "/, "", a[i]);
                gsub(/".*/, "", a[i]);
                gene_id = a[i];
            }
            if (a[i] ~ /gene_name/) {
                gsub(/.*gene_name "/, "", a[i]);
                gsub(/".*/, "", a[i]);
                gene_name = a[i];
            }
            if (a[i] ~ /gene_type/ || a[i] ~ /gene_biotype/) {
                if (a[i] ~ /gene_type/) {
                    gsub(/.*gene_type "/, "", a[i]);
                } else {
                    gsub(/.*gene_biotype "/, "", a[i]);
                }
                gsub(/".*/, "", a[i]);
                gene_type = a[i];
            }
        }
        if (gene_name == "") gene_name = gene_id;
        if (gene_type == "") gene_type = "unknown";
        print \$1, \$4, \$5, gene_id"|"gene_name"|"gene_type, ".", \$7;
    }' ${gtf} > genes.bed

    # Use HOMER annotatePeaks.pl to annotate genes
    annotatePeaks.pl \\
        genes.bed \\
        ${fasta} \\
        -gtf ${gtf} \\
        -cpu ${task.cpus} \\
        > gene_annotation.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        homer: $VERSION
    END_VERSIONS
    """

    stub:
    def VERSION = '4.11'
    """
    touch gene_annotation.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        homer: $VERSION
    END_VERSIONS
    """
}
