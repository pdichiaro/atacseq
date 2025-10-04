process DESEQ2_TRANSFORM {
    label 'process_single'
    tag "$deseq2_file"

    conda "conda-forge::sed=4.7"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'nf-core/ubuntu:20.04' }"

    input:
    path deseq2_file
    path pca_header
    path clustering_header
    path read_dist_header

    output:
    path "*_mqc.txt", emit: multiqc_files, optional: true
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def file_name = deseq2_file.getName()
    // Create output name with _mqc.txt suffix (replace .txt with _mqc.txt)
    def output_name = file_name.replaceAll(/\.txt$/, '_mqc.txt')
    """
    echo "=== DESEQ2_TRANSFORM DEBUG ==="
    echo "Input file: ${file_name}"
    echo "Output name: ${output_name}"
    echo "File exists: \$(test -f ${deseq2_file} && echo YES || echo NO)"
    echo "File size: \$(wc -l ${deseq2_file} 2>/dev/null || echo 'Cannot read')"
    echo "=============================="
    
    # Initialize SKIP_FILE flag (used for Read Distribution plots)
    SKIP_FILE=false
    
    # Detect quantifier and level from filename for unique IDs and section anchors
    # Filename patterns: featureCounts.deseq2.all_genes.*, featureCounts.deseq2.invariant_genes.*, etc.
    QUANTIFIER=""
    QUANTIFIER_SHORT=""
    LEVEL=""
    SECTION_NAME=""
    
    # Extract quantifier - ATAC-seq specific (featureCounts is primary quantifier)
    # Use case-insensitive matching for featurecounts/featureCounts
    shopt -s nocasematch
    if [[ "${file_name}" == featurecounts.* || "${file_name}" == featureCounts.* || "${file_name}" == subread.* ]]; then
        QUANTIFIER="deseq2-featurecounts-qc"
        QUANTIFIER_SHORT="featurecounts"
        PARENT_NAME="DESeq2 FeatureCounts QC"
    else
        echo "Warning: Could not determine quantifier from filename: ${file_name}"
        QUANTIFIER="deseq2-qc"
        QUANTIFIER_SHORT="unknown"
        PARENT_NAME="DESeq2 QC"
    fi
    shopt -u nocasematch
    
    # Extract level (all_genes or invariant_genes)
    if [[ "${file_name}" == *.all_genes.* ]]; then
        LEVEL="all_genes"
        SECTION_NAME="All Genes"
    elif [[ "${file_name}" == *.invariant_genes.* ]]; then
        LEVEL="invariant_genes"
        SECTION_NAME="Invariant Genes"
    else
        LEVEL="unknown"
        SECTION_NAME="Unknown Level"
    fi
    
    echo "Detected quantifier: \${QUANTIFIER_SHORT}, level: \${LEVEL}, section: \${QUANTIFIER}"
    echo "Output file will be: ${output_name}"

    # Determine number prefix based on gene set (01-03 for all_genes, 04-06 for invariant_genes)
    # Read Distribution plots are EXCLUDED from MultiQC
    OFFSET=0
    if [[ "\${LEVEL}" == "invariant_genes" ]]; then
        OFFSET=3
    fi
    
    # Add appropriate header to each file type for MultiQC custom content module
    # Each plot gets nested under parent section with unique ID
    # Plot ordering controlled by section_order field in metadata
    # Number ranges: 01-03 for All Genes, 04-06 for Invariant Genes
    # NOTE: Read Distribution plots are skipped (not included in MultiQC)
    # IMPORTANT: Check .pca.top*.vals.txt BEFORE .pca.vals.txt to avoid false matches
    if [[ "${file_name}" == *".pca.top"*".vals.txt" ]]; then
        # PCA top variable genes (pattern: *.pca.top500.vals.txt) - ORDER: 03 or 06
        PLOT_NUM=\$((3 + OFFSET))
        PLOT_ID="\$(printf '%02d' \$PLOT_NUM)_deseq2_pca_top500_\${QUANTIFIER_SHORT}_\${LEVEL}"
        SECTION_TITLE="PCA Top 500 (\${SECTION_NAME})"
        PLOT_TITLE="PCA Top 500 (\${SECTION_NAME})"
        output_file="\${PLOT_ID}_mqc.txt"
        sed "s|#section_anchor:.*|#parent_id: '\${QUANTIFIER}'\\n#parent_name: '\${PARENT_NAME}'\\n#section_order: \${PLOT_NUM}|; s|#section_name:.*|#section_name: '\${SECTION_TITLE}'|; s|#id:.*|#id: '\${PLOT_ID}'|; s|title:.*|title: '\${PLOT_TITLE}'|" ${pca_header} > temp_header.txt
        cat temp_header.txt ${deseq2_file} > temp_output.txt
        mv temp_output.txt "\${output_file}"
        echo "Created \${output_file} with PCA-500 header (ID: \${PLOT_ID}, section: \${SECTION_TITLE}, parent: \${QUANTIFIER})"
    elif [[ "${file_name}" == *".pca.vals.txt" ]]; then
        # PCA all genes (pattern: *.pca.vals.txt) - ORDER: 02 or 05
        PLOT_NUM=\$((2 + OFFSET))
        PLOT_ID="\$(printf '%02d' \$PLOT_NUM)_deseq2_pca_\${QUANTIFIER_SHORT}_\${LEVEL}"
        SECTION_TITLE="PCA (\${SECTION_NAME})"
        PLOT_TITLE="PCA (\${SECTION_NAME})"
        output_file="\${PLOT_ID}_mqc.txt"
        sed "s|#section_anchor:.*|#parent_id: '\${QUANTIFIER}'\\n#parent_name: '\${PARENT_NAME}'\\n#section_order: \${PLOT_NUM}|; s|#section_name:.*|#section_name: '\${SECTION_TITLE}'|; s|#id:.*|#id: '\${PLOT_ID}'|; s|title:.*|title: '\${PLOT_TITLE}'|" ${pca_header} > temp_header.txt
        cat temp_header.txt ${deseq2_file} > temp_output.txt
        mv temp_output.txt "\${output_file}"
        echo "Created \${output_file} with PCA header (ID: \${PLOT_ID}, section: \${SECTION_TITLE}, parent: \${QUANTIFIER})"
    elif [[ "${file_name}" == *".sample.dists.txt" ]]; then
        # Sample distance - ORDER: 01 or 04
        PLOT_NUM=\$((1 + OFFSET))
        PLOT_ID="\$(printf '%02d' \$PLOT_NUM)_deseq2_sample_distance_\${QUANTIFIER_SHORT}_\${LEVEL}"
        SECTION_TITLE="Sample Distances (\${SECTION_NAME})"
        PLOT_TITLE="Sample Distances (\${SECTION_NAME})"
        output_file="\${PLOT_ID}_mqc.txt"
        sed "s|#section_anchor:.*|#parent_id: '\${QUANTIFIER}'\\n#parent_name: '\${PARENT_NAME}'\\n#section_order: \${PLOT_NUM}|; s|#section_name:.*|#section_name: '\${SECTION_TITLE}'|; s|#id:.*|#id: '\${PLOT_ID}'|; s|title:.*|title: '\${PLOT_TITLE}'|" ${clustering_header} > temp_header.txt
        cat temp_header.txt ${deseq2_file} > temp_output.txt
        mv temp_output.txt "\${output_file}"
        echo "Created \${output_file} with sample distance header (ID: \${PLOT_ID}, section: \${SECTION_TITLE}, parent: \${QUANTIFIER})"
    elif [[ "${file_name}" == *".read.distribution.normalized.txt" ]]; then
        # Read distribution - SKIPPED (not included in MultiQC report)
        echo "⏭️  Skipping Read Distribution plot for MultiQC (file: ${file_name})"
        # Set a flag to skip the file existence check
        output_file=""
        SKIP_FILE=true
    else
        # Unknown file type - ERROR: All file types should be explicitly handled
        echo "❌ ERROR: Unknown file type for: ${file_name}"
        echo "   This file does not match any expected pattern:"
        echo "   - *.pca.top*.vals.txt"
        echo "   - *.pca.vals.txt"
        echo "   - *.sample.dists.txt"
        echo "   - *.read.distribution.normalized.txt"
        echo ""
        echo "   Please add explicit handling for this file type in deseq2_transform/main.nf"
        exit 1
    fi
    
    # Verify output file was created (unless explicitly skipped)
    if [[ "\${SKIP_FILE}" != "true" ]]; then
        if [[ ! -f "\${output_file}" ]]; then
            echo "❌ ERROR: Failed to create output file: \${output_file}"
            echo "All files in directory:"
            ls -lah
            exit 1
        fi
        
        echo "=== OUTPUT FILE CREATED ==="
        ls -lh "\${output_file}"
        echo "✅ Successfully created MultiQC file"
        echo "============================"
    else
        echo "✅ File skipped as requested (Read Distribution excluded from MultiQC)"
    fi

    cat <<END_VERSIONS > versions.yml
"${task.process}":
    bash: \$(bash --version | head -n1 | awk '{print \$4}')
END_VERSIONS
    """

    stub:
    """
    touch stub_file.txt

    cat <<END_VERSIONS > versions.yml
"${task.process}":
    bash: \$(bash --version | head -n1 | awk '{print \$4}')
END_VERSIONS
    """
}
