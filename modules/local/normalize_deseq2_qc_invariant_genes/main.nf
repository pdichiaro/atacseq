process NORMALIZE_DESEQ2_QC_INVARIANT_GENES {
    tag "deseq2_qc_invariant_genes"
    label "process_medium"

    container 'docker://pdichiaro/env_r_ngs'

    input:
    path counts
    val quantifier
    path homer_annotation, stageAs: 'homer_annotation.txt'

    output:
    path "normalization"        , emit: normalization_results
    path "scaling_dat.txt"      , emit: scaling_factors
    path "scaling_factors/*_scaling_factor.txt", emit: scaling_factors_individual
    path "*_normalized_counts.txt", emit: normalized_counts
    path "*_rlog_counts.txt"    , optional:true, emit: rlog_counts
    path "*.RData"              , optional:true, emit: rdata
    path "*.pca.vals.txt"       , optional:true, emit: pca_all_genes_txt
    path "*.pca.top*.vals.txt"  , optional:true, emit: pca_top_genes_txt
    path "*.sample.dists.txt"   , optional:true, emit: sample_distances_txt
    path "*.read.distribution.normalized.txt", optional:true, emit: read_dist_norm_txt
    path "Quality_Control"      , optional:true, emit: quality_control_plots
    path "Read_Distribution"    , optional:true, emit: read_dist_plots
    path "*.log"                , optional:true, emit: log
    path "versions.yml"         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args  = task.ext.args  ?: ''
    def sigma_times = params.sigma_times ?: '1'
    def n_pop = params.n_pop ?: '1'
    prefix = task.ext.prefix ?: "deseq2_invariant"
    """
    normalize_deseq2_qc_invariant_genes.r \\
        --count_file $counts \\
        --count_col 7 \\
        --outdir ./ \\
        --outprefix "$prefix" \\
        --quantifier "$quantifier" \\
        --sigma_times $sigma_times \\
        --n_pop $n_pop \\
        --cores $task.cpus \\
        --homer_annotation homer_annotation.txt \\
        $args

    echo "=== FILES GENERATED FOR MULTIQC ==="
    ls -la *.txt 2>/dev/null || echo "No .txt files found"
    echo "===================================="

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        bioconductor-deseq2: \$(Rscript -e "library(DESeq2); cat(as.character(packageVersion('DESeq2')))")
        r-generalnormalizer: \$(Rscript -e "library(GeneralNormalizer); cat(as.character(packageVersion('GeneralNormalizer')))")
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "deseq2_invariant"
    """
    touch ${prefix}.dds.RData
    touch ${prefix}_rlog_counts.txt
    touch ${prefix}.normalized_filt.txt
    touch ${prefix}.pca.vals.txt
    touch ${prefix}.pca.top500.vals.txt
    touch ${prefix}.sample.dists.txt
    touch ${prefix}.sample.correlations.txt
    touch ${prefix}.read.distribution.normalized.txt
    touch ${prefix}.read.distribution.histogram.normalized.txt
    touch ${prefix}.read.distribution.density.normalized.txt
    touch ${prefix}_normalized_counts.txt
    touch scaling_dat.txt
    touch R_sessionInfo.log
    
    mkdir -p Quality_Control
    touch Quality_Control/Sample_Distance_Heatmap.pdf
    touch Quality_Control/Sample_Distance_Heatmap.png
    touch Quality_Control/Sample_Correlation_Heatmap.pdf
    touch Quality_Control/Sample_Correlation_Heatmap.png
    touch Quality_Control/PCA_Plot_All_Genes.pdf
    touch Quality_Control/PCA_Plot_Top500_Genes.pdf
    
    mkdir -p Read_Distribution
    touch Read_Distribution/Read_Distribution_Invariant_Norm_boxplot.pdf
    touch Read_Distribution/Read_Distribution_Invariant_Norm_histogram.pdf
    touch Read_Distribution/Read_Distribution_Invariant_Norm_density.pdf

    mkdir -p normalization
    touch normalization/invariant_genes.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        bioconductor-deseq2: \$(Rscript -e "library(DESeq2); cat(as.character(packageVersion('DESeq2')))")
        r-generalnormalizer: \$(Rscript -e "library(GeneralNormalizer); cat(as.character(packageVersion('GeneralNormalizer')))")
    END_VERSIONS
    """
}