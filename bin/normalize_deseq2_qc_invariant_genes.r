#!/usr/bin/env Rscript

################################################
################################################
## INVARIANT GENES NORMALIZATION + DESEQ2 QC  ##
################################################
################################################

## This script performs advanced normalization using GeneralNormalizer
## with invariant gene identification, followed by comprehensive DESeq2 QC analysis
## Adapted to work with count matrices from Kallisto, RSEM, and Salmon

################################################
################################################
## LOAD LIBRARIES                             ##
################################################
################################################

library(optparse)
library(DESeq2)
library(dplyr)
library(tidyr)
library(tibble)
library(purrr)
library(GeneralNormalizer)
library(ComplexHeatmap)
library(ggplot2)
library(RColorBrewer)
library(pheatmap)


suppressPackageStartupMessages(library(ComplexHeatmap))

################################################
################################################
## PARSE COMMAND-LINE PARAMETERS              ##
################################################
################################################

option_list <- list(
    make_option(c("-i", "--count_file"    ), type="character", default=NULL    , metavar="path"   , help="Count file matrix where rows are genes and columns are samples."                        ),
    make_option(c("-f", "--count_col"     ), type="integer"  , default=7       , metavar="integer", help="First column containing sample count data."                                             ),
    make_option(c("-d", "--id_col"        ), type="integer"  , default=1       , metavar="integer", help="Column containing identifiers to be used."                                              ),
    make_option(c("-r", "--sample_suffix" ), type="character", default=''      , metavar="string" , help="Suffix to remove after sample name in columns e.g. '.rmDup.bam' if 'DRUG_R1.rmDup.bam'."),
    make_option(c("-o", "--outdir"        ), type="character", default='./'    , metavar="path"   , help="Output directory."                                                                      ),
    make_option(c("-p", "--outprefix"     ), type="character", default='invariant_deseq2_qc', metavar="string" , help="Output prefix."                                                                         ),
    make_option(c("-v", "--vst"           ), type="logical"  , default=FALSE   , metavar="boolean", help="Run vst transform instead of rlog."                                                     ),
    make_option(c("-s", "--sigma_times"   ), type="integer"  , default=1       , metavar="integer", help="Number of sigma times for invariant gene selection."                                   ),
    make_option(c("-z", "--n_pop"         ), type="integer"  , default=1       , metavar="integer", help="Number of populations for normalization."                                               ),
    make_option(c("-c", "--cores"         ), type="integer"  , default=1       , metavar="integer", help="Number of cores."                                                                       ),
    make_option(c("-q", "--quantifier"    ), type="character", default=''      , metavar="string" , help="Quantification method for output directory organization."                              ),
    make_option(c("-a", "--homer_annotation"), type="character", default=NULL  , metavar="path"   , help="HOMER annotatepeaks.txt file for peak annotations."                                     )
)

opt_parser <- OptionParser(option_list=option_list)
opt        <- parse_args(opt_parser)

if (is.null(opt$count_file)){
    print_help(opt_parser)
    stop("Please provide a counts file.", call.=FALSE)
}

# Create file prefix based on quantifier for MultiQC compatibility
# MultiQC config expects files like: kallisto.deseq2.invariant_genes.*.txt
if (opt$quantifier != '' && opt$quantifier != 'null') {
    quantifier_lower <- tolower(opt$quantifier)
    if (quantifier_lower == "star_rsem") {
        file_prefix <- "star.rsem.deseq2.invariant_genes"
    } else if (quantifier_lower == "star_salmon") {
        file_prefix <- "star.salmon.deseq2.invariant_genes"
    } else if (quantifier_lower == "star_genome") {
        file_prefix <- "star.genome.deseq2.invariant_genes"
    } else {
        file_prefix <- paste0(quantifier_lower, ".deseq2.invariant_genes")
    }
} else {
    file_prefix <- opt$outprefix
}

cat("=== INVARIANT GENES NORMALIZATION + DESEQ2 QC ===\n")
cat("Input file:", opt$count_file, "\n")
cat("File prefix for outputs:", file_prefix, "\n")

################################################
################################################
## READ IN COUNTS FILE                        ##
################################################
################################################

# featureCounts files start with comment lines (starting with #)
# We need to skip these lines before reading the data
count.table           <- read.delim(file=opt$count_file, header=TRUE, row.names=NULL, check.names=FALSE, comment.char="#")
rownames(count.table) <- count.table[,opt$id_col]

# Identify annotation columns vs sample count columns before subsetting
full_count_table <- count.table  # Keep the full table with annotations
# Ensure count_col is valid for this dataset
if (opt$count_col > ncol(count.table)) {
    cat("Warning: count_col (", opt$count_col, ") exceeds number of columns (", ncol(count.table), "). Using all columns as samples.\n")
    opt$count_col <- 2  # Default to column 2 as first sample column
}

sample_cols <- colnames(count.table)[opt$count_col:ncol(count.table)]
annotation_cols <- colnames(count.table)[1:(opt$count_col-1)]

cat("Total columns:", ncol(count.table), "\n")
cat("Count column start:", opt$count_col, "\n") 
cat("Detected annotation columns:", paste(annotation_cols, collapse = ", "), "\n")
cat("Detected sample columns:", paste(head(sample_cols, 5), collapse = ", "), "...\n")

# Extract annotation data from the input file (safely)
annotation_from_input <- NULL
if (length(annotation_cols) > 0 && opt$count_col > 1) {
    annotation_from_input <- count.table[, annotation_cols, drop = FALSE]
    cat("Extracted", ncol(annotation_from_input), "annotation columns from input file\n")
} else {
    cat("No annotation columns in input file\n")
    annotation_from_input <- data.frame(gene_id = rownames(count.table))
}

# Now subset count.table to only sample columns for DESeq2 processing
count.table           <- count.table[,opt$count_col:ncol(count.table),drop=FALSE]
colnames(count.table) <- gsub(opt$sample_suffix,"",colnames(count.table))
colnames(count.table) <- gsub(pattern='\\.$', replacement='', colnames(count.table))

cat("Count matrix dimensions:", nrow(count.table), "genes x", ncol(count.table), "samples\n")



################################################
################################################
## RUN INVARIANT GENES NORMALIZATION          ##
################################################
################################################

cpus=opt$cores
if(!is.null(cpus)){
    param <- BiocParallel::SnowParam(workers=cpus,tasks=0,stop.on.error=TRUE,progressbar = TRUE,type="SOCK")
}else{
    param <- NULL
}

# Create temporary normalization folder for GeneralNormalizer outputs
# The organized directory structure will be created later
temp_norm_dir <- if (opt$quantifier != '') {
    paste0(opt$quantifier, "_normalization_temp/")
} else {
    "normalization_temp/"
}

save_folder_norm <- temp_norm_dir
if (!dir.exists(save_folder_norm)) {
    dir.create(save_folder_norm, recursive=TRUE)
}

cat("\n=== STARTING INVARIANT GENES NORMALIZATION ===\n")
cat("Temporary normalization folder:", save_folder_norm, "\n")

# Assemble the experimental design using sample name components
cn =  colnames(count.table)
name_components <- strsplit(cn, "_")
n_components    <- length(unlist(name_components))
decompose       <- n_components!=1 && all(sapply(name_components, length)==n_components)
coldata         <- data.frame(cn, Sample_ID=cn, row.names=1)
if (decompose) {
    groupings        <- as.data.frame(lapply(1:n_components, function(i) sapply(name_components, "[[", i)))
    names(groupings) <- paste0("Group", 1:n_components)
    n_distinct       <- sapply(groupings, function(grp) length(unique(grp)))
    groupings        <- groupings[n_distinct!=1 & n_distinct!=length(cn)]
    if (ncol(groupings)!=0) {
        coldata <- cbind(coldata, groupings)
    } else {
        decompose <- FALSE
    }
}else{
    coldata$Group1 = gsub("_R[0-9].*","",cn)
}

# Take only the R[0-9] part so remove text before any _R[0-9]
coldata$Sample_Replicate = paste0("R",gsub(".*_R","",colnames(count.table)))

design = coldata

x = count.table
x_depth_p = x/colSums(x)*1000000

cat("Original gene count:", nrow(x), "\n")

# Remove outlier genes (top 0.1% and bottom 0.1% by mean expression)
# This enriches for stably expressed genes
x = x[order(rowMeans(x_depth_p), decreasing = T),]
x = x[round(nrow(x)*0.001):round(nrow(x)*0.999),]

cat("Genes after outlier removal:", nrow(x), "\n")
cat("Removed", nrow(count.table) - nrow(x), "outlier genes (top/bottom 0.1%)\n")

# Compute the normalization using GeneralNormalizer
sigma_times = opt$sigma_times
n_pop = opt$n_pop

cat("Running GeneralNormalizer with parameters:\n")
cat("  - sigma_times:", sigma_times, "\n")
cat("  - n_pop:", n_pop, "\n")
cat("  - cores:", cpus, "\n")

Result = RunNorm(   x,
                    design,
                    sigma_times = sigma_times,
                    fix_reference="random",
                    row_name_index=rownames(x),
                    saving_path=save_folder_norm,
                    n_pop=n_pop,
                    n_pop_reference=n_pop,
                    BiocParam=param
                    )

colData_norm = Result$scaling_factors
mat = Result$norm_mat

cat("Invariant genes normalization completed!\n")
cat("Scaling factors calculated for", nrow(colData_norm), "samples\n")
                               
################################################
################################################
## SAVE SCALING FACTORS                       ##
################################################
################################################

# Define base_dir for output directory structure
if (opt$quantifier != '') {
    base_dir <- paste0(opt$quantifier, "/deseq2/invariant_genes/")
} else {
    base_dir <- "deseq2/invariant_genes/"
}

cat("Quantifier parameter:", opt$quantifier, "\n")
cat("Creating output directory structure:", base_dir, "\n")

# Create early output directory for scaling factors
early_save_folder = paste0(base_dir, "scaling_factors/")
if (file.exists(early_save_folder) == FALSE) {
    dir.create(early_save_folder, recursive=TRUE)
}

# Save GeneralNormalizer scaling factors in text format (no RData needed)
cat("Saving scaling factors...\n")

# Save individual scaling factor files for each sample
for (sample_name in rownames(colData_norm)) {
    scalingFactorFile <- paste0(early_save_folder, sample_name, ".txt")
    write(as.numeric(colData_norm[sample_name, "scaling"]), file=scalingFactorFile)
}

# Save scaling factors summary table (more useful than RData)
scaling_factors_summary <- data.frame(
    sample = rownames(colData_norm),
    scaling_factor = colData_norm$scaling,
    stringsAsFactors = FALSE
)
write.table(scaling_factors_summary, 
           file = paste0(early_save_folder, opt$outprefix, ".scaling_factors_summary.txt"),
           sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

# Create scaling_dat.txt (scaling factors in simple format)
scaling_dat_early <- data.frame(
    sample = rownames(colData_norm),
    scaling_factor = colData_norm$scaling,
    stringsAsFactors = FALSE
)
write.table(scaling_dat_early, file = paste0(early_save_folder, "scaling_dat.txt"), sep = "\t", 
           row.names = FALSE, col.names = TRUE, quote = FALSE)
# Also write to root for Nextflow
write.table(scaling_dat_early, file = "scaling_dat.txt", sep = "\t", 
           row.names = FALSE, col.names = TRUE, quote = FALSE)

# Create per-sample scaling factor files for safe DeepTools normalization
# This avoids grep parsing issues and is more robust
cat("Creating per-sample scaling factor files...\n")
scaling_factors_dir <- "scaling_factors"
dir.create(scaling_factors_dir, showWarnings = FALSE, recursive = TRUE)

for (i in 1:nrow(scaling_dat_early)) {
    sample_name <- scaling_dat_early$sample[i]
    scaling_factor <- scaling_dat_early$scaling_factor[i]
    
    # Write individual scaling factor file (just the numeric value)
    sample_file <- file.path(scaling_factors_dir, paste0(sample_name, "_scaling_factor.txt"))
    writeLines(as.character(scaling_factor), sample_file)
    
    cat("  - Created:", sample_file, "with scaling factor:", scaling_factor, "\n")
}



cat("Scaling factors saved to:", early_save_folder, "\n")

################################################
################################################
## RUN DESEQ2 QC WITH CUSTOM SCALING FACTORS  ##
################################################
################################################

if (file.exists(opt$outdir) == FALSE) {
    dir.create(opt$outdir, recursive=TRUE)
}
setwd(opt$outdir)

cat("\n=== STARTING DESEQ2 QC ANALYSIS ===\n")

samples.vec     <- colnames(count.table)
name_components <- strsplit(samples.vec, "_")
n_components    <- length(name_components[[1]])
decompose       <- n_components!=1 && all(sapply(name_components, length)==n_components)
coldata         <- data.frame(samples.vec, sample=samples.vec, row.names=1)
if (decompose) {
    groupings        <- as.data.frame(lapply(1:n_components, function(i) sapply(name_components, "[[", i)))
    n_distinct       <- sapply(groupings, function(grp) length(unique(grp)))
    groupings        <- groupings[n_distinct!=1 & n_distinct!=length(samples.vec)]
    if (ncol(groupings)!=0) {
        names(groupings) <- paste0("Group", 1:ncol(groupings))
        coldata <- cbind(coldata, groupings)
    } else {
        decompose <- FALSE
    }
}

DDSFile <- paste(opt$outprefix,".dds.RData",sep="")

################################################
## ZERO COUNT REMOVAL FROM NORMALIZED MATRIX ##
################################################

# Apply zero count filtering to the GeneralNormalizer normalized matrix
cat("Applying zero count filtering to GeneralNormalizer normalized matrix...\n")
normalized_matrix <- Result$norm_mat
cat("Starting with", nrow(normalized_matrix), "genes from GeneralNormalizer\n")

# Remove genes with zero counts across all samples from normalized data
keep_genes <- rowSums(normalized_matrix) > 0
normalized_matrix_filtered <- normalized_matrix[keep_genes, ]

cat("Genes after zero count removal from normalized matrix:", nrow(normalized_matrix_filtered), "\n")
cat("Filtered out", sum(!keep_genes), "genes with zero counts from normalized data\n")

# Ensure proper column and row names for DESeq2 
# Keep the original sample names from the normalized matrix
original_sample_names <- colnames(normalized_matrix)
colnames(normalized_matrix_filtered) <- original_sample_names
rownames(normalized_matrix_filtered) <- rownames(normalized_matrix)[keep_genes]

cat("Normalized matrix sample names:", paste(head(colnames(normalized_matrix_filtered), 3), collapse = ", "), "...\n")


# Create base directory if it doesn't exist
if (file.exists(base_dir) == FALSE) {
    dir.create(base_dir, recursive=TRUE)
}

# Create final normalization folder and move temp files if needed
final_norm_dir = paste0(base_dir, "normalization/")
if (!dir.exists(final_norm_dir)) {
    dir.create(final_norm_dir, recursive=TRUE)
}

# Check for normalization files and move them to organized location
if (dir.exists(save_folder_norm)) {
    temp_files <- list.files(save_folder_norm, full.names = TRUE)
    cat("Found", length(temp_files), "files in temp normalization folder:", save_folder_norm, "\n")
    if (length(temp_files) > 0) {
        # List the files for debugging
        cat("Files found:", paste(basename(temp_files), collapse = ", "), "\n")
        file.copy(temp_files, final_norm_dir, recursive = TRUE)
        cat("Moved", length(temp_files), "normalization files to organized location\n")
        # Clean up temp directory only if different from final
        if (save_folder_norm != final_norm_dir) {
            unlink(save_folder_norm, recursive = TRUE)
        }
    }
} else {
    cat("Warning: temp normalization folder", save_folder_norm, "does not exist\n")
}

# Create normalization directory at root level for Nextflow module output
if (!dir.exists("normalization")) {
    dir.create("normalization", recursive = TRUE)
    cat("Created root normalization directory for Nextflow\n")
}

# Copy contents from final_norm_dir to root level normalization for Nextflow
if (dir.exists(final_norm_dir)) {
    norm_files <- list.files(final_norm_dir, full.names = TRUE, recursive = TRUE)
    cat("Found", length(norm_files), "files in final_norm_dir:", final_norm_dir, "\n")
    if (length(norm_files) > 0) {
        cat("Copying files:", paste(basename(norm_files), collapse = ", "), "\n")
        file.copy(norm_files, "normalization/", recursive = TRUE)
        cat("Copied", length(norm_files), "files to root normalization directory for Nextflow\n")
    } else {
        cat("Warning: No files found in final_norm_dir to copy to root\n")
    }
} else {
    cat("Warning: final_norm_dir", final_norm_dir, "does not exist\n")
    # Create empty normalization directory with a placeholder file
    if (!dir.exists("normalization")) {
        dir.create("normalization", recursive = TRUE)
    }
    write("Normalization files not available", file = "normalization/placeholder.txt")
    cat("Created placeholder in normalization directory\n")
}

                               
################################################
## SAVE NORMALIZED COUNT MATRIX               ##
################################################

# Use the filtered GeneralNormalizer normalized matrix (zero-count filtered matrix_normalized.txt)
deseq2_normalized_counts <- normalized_matrix_filtered
deseq2_normalized_df <- data.frame(
    gene_id = rownames(deseq2_normalized_counts),
    deseq2_normalized_counts,
    check.names = FALSE,
    stringsAsFactors = FALSE
)

# Add annotation merge by gene_id for normalized_counts.txt (same as rlog_counts.txt approach)
cat("Adding annotation to normalized counts using merge by gene_id...\n")

# Get the original count table with all columns (annotation + counts) for annotation merge
original_table <- read.delim(file=opt$count_file, header=TRUE, row.names=NULL, comment.char="#")
original_annotation <- original_table[, 1:(opt$count_col-1), drop = FALSE]

# Filter annotation to match normalized counts
genes_to_keep <- deseq2_normalized_df$gene_id
original_annotation_filtered <- original_annotation %>% filter(Geneid %in% genes_to_keep)                          
                               
cat("Original annotation rows:", nrow(original_annotation), "\n")
cat("Original annotation filtered rows:", nrow(original_annotation_filtered), "\n")
cat("Normalized counts rows:", nrow(deseq2_normalized_df), "\n")

# Load HOMER annotations if provided
homer_annotation <- NULL
if (!is.null(opt$homer_annotation) && file.exists(opt$homer_annotation)) {
    cat("Loading HOMER annotations from:", opt$homer_annotation, "\n")
    homer_annotation <- read.delim(opt$homer_annotation, header=TRUE, stringsAsFactors=FALSE, comment.char="#")
    # Rename the first column to "Gwneid" (HOMER often has malformed column name)
    colnames(homer_annotation)[1] <- "Geneid"
    cat("HOMER annotation column names:", paste(colnames(homer_annotation), collapse=", "), "\n")
    cat("HOMER annotation dimensions:", nrow(homer_annotation), "x", ncol(homer_annotation), "\n")
}

# Merge: first with original annotation, then with HOMER if available
if (!is.null(homer_annotation)) {
    # First left_join original annotation with HOMER (keeps original_annotation columns, adds HOMER columns)
    cat("Merging original annotation with HOMER annotations...\n")
    combined_annotation <- left_join(original_annotation_filtered, homer_annotation, by = c("Geneid","Chr","Start","End"))
    cat("Combined annotation dimensions:", nrow(combined_annotation), "x", ncol(combined_annotation), "\n")
    
    # Then left_join with normalized counts
    final_normalized_counts <- left_join(combined_annotation, deseq2_normalized_df, by = c("Geneid" = "gene_id"))
    
} else {
    # Simple left_join by gene_id
    final_normalized_counts <- left_join(original_annotation_filtered, deseq2_normalized_df,  by = c("Geneid" = "gene_id"))
}

cat("Final normalized counts dimensions:", nrow(final_normalized_counts), "x", ncol(final_normalized_counts), "\n")

write.table(final_normalized_counts, file = paste0(base_dir, "normalized_counts.txt"), sep = "\t",
           row.names = FALSE, col.names = TRUE, quote = FALSE)
cat("Normalized count matrix (with annotation) saved as:", paste0(base_dir, "normalized_counts.txt"), "\n")

# Also create normalized_counts.txt at root level for Nextflow module output (with prefix)
normalized_counts_root <- paste0(opt$outprefix, "_normalized_counts.txt")
write.table(final_normalized_counts, file = normalized_counts_root, sep = "\t",
           row.names = FALSE, col.names = TRUE, quote = FALSE)
cat("Normalized count matrix (with annotation) also saved to root as:", normalized_counts_root, "\n")

################################################
## READ DISTRIBUTION QC PLOTS                ##
################################################

# Create directories for read distribution plots
ReadDistDir <- paste0(base_dir, "Read_Distribution/")
if (file.exists(ReadDistDir) == FALSE) {
    dir.create(ReadDistDir, recursive=TRUE)
}

# Also create Read_Distribution at root level for Nextflow module output
ReadDistDirRoot <- "Read_Distribution/"
if (file.exists(ReadDistDirRoot) == FALSE) {
    dir.create(ReadDistDirRoot, recursive=TRUE)
}

# Plot read distribution using the normalized count data
cat("Creating read distribution plots for invariant gene normalized data...\n")

# Use the final_normalized_counts data for plotting (without annotation columns)
# Extract only the count columns (skip annotation columns)
# Use sample names directly from the normalized matrix to avoid annotation column issues
sample_columns <- colnames(deseq2_normalized_counts)
cat("Sample columns for plotting:", paste(sample_columns, collapse=", "), "\n")

# Verify all sample columns exist in final_normalized_counts
missing_cols <- sample_columns[!sample_columns %in% colnames(final_normalized_counts)]
if (length(missing_cols) > 0) {
    cat("WARNING: Missing sample columns:", paste(missing_cols, collapse=", "), "\n")
    sample_columns <- sample_columns[sample_columns %in% colnames(final_normalized_counts)]
}

plot_data <- final_normalized_counts[, c("Geneid", sample_columns), drop = FALSE]
cat("Plot data dimensions:", nrow(plot_data), "x", ncol(plot_data), "\n")

# Check if we have grouping information for coldata setup
samples.vec <- sample_columns  # Use the actual sample column names
name_components <- strsplit(samples.vec, "_")
n_components <- length(name_components[[1]])
decompose_plots <- n_components != 1 && all(sapply(name_components, length) == n_components)
coldata_plots <- data.frame(samples.vec, sample = samples.vec, row.names = 1)

if (decompose_plots) {
    groupings <- as.data.frame(lapply(1:n_components, function(i) sapply(name_components, "[[", i)))
    n_distinct <- sapply(groupings, function(grp) length(unique(grp)))
    groupings <- groupings[n_distinct != 1 & n_distinct != length(samples.vec)]
    if (ncol(groupings) != 0) {
        names(groupings) <- paste0("Group", 1:ncol(groupings))
        coldata_plots <- cbind(coldata_plots, groupings)
    } else {
        decompose_plots <- FALSE
    }
}

# Set up Bio_replicates for plotting
if (decompose_plots && ncol(coldata_plots) > 1) {
    group_col <- names(coldata_plots)[2]  # First column after 'sample'
    coldata_plots$Bio_replicates <- coldata_plots[[group_col]]
} else {
    coldata_plots$Bio_replicates <- "Sample"
}

# Prepare data for plotting - use log-transformed normalized counts
plot_matrix <- as.matrix(plot_data[, sample_columns])
rownames(plot_matrix) <- plot_data$Geneid

# Verify matrix is numeric
if (!is.numeric(plot_matrix)) {
    cat("ERROR: plot_matrix contains non-numeric values!\n")
    cat("Class of plot_matrix:", class(plot_matrix), "\n")
    cat("Sample of plot_matrix (first 5 rows, first 3 cols):\n")
    print(plot_matrix[1:min(5, nrow(plot_matrix)), 1:min(3, ncol(plot_matrix))])
    stop("Cannot proceed with non-numeric matrix")
}

all.reads_z <- as.data.frame(log(plot_matrix + 1))
all.reads_zz <- all.reads_z[, rownames(coldata_plots)] %>% rownames_to_column(var="gene") %>% as_tibble()            
gathered_all.reads <- gather(all.reads_zz, key = "samplename", value = "normalized_counts", -gene)
gathered_all.reads <- gathered_all.reads %>% 
  left_join(coldata_plots %>% rownames_to_column("samplename") %>% dplyr::select(samplename, Bio_replicates), by = c("samplename")) %>%
  arrange(Bio_replicates) %>% 
  mutate(samplename = factor(samplename, levels = unique(samplename)))

# Create PNG plots for read distribution (invariant gene normalized)
p.1 <- ggplot(gathered_all.reads, aes(x = samplename, y = normalized_counts, fill = Bio_replicates)) + 
  geom_boxplot(show.legend = FALSE) + xlab("") +
  ylab(expression(ln(count + 1))) + theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p.2 <- ggplot(gathered_all.reads, aes(x = normalized_counts, colour = Bio_replicates, fill = Bio_replicates)) +
  geom_histogram(binwidth = 1) + xlab(expression(ln(count + 1))) + ylab("frequency") + 
  ylim(c(0,200000)) + theme(legend.position = "top") + theme_classic()

p.3 <- ggplot(gathered_all.reads, aes(x = normalized_counts, colour = Bio_replicates, fill = Bio_replicates)) + 
  geom_density(alpha = 0.2, linewidth = 1.25) + xlab(expression(ln(count))) + ylim(c(0, 0.5)) +
  theme(legend.position = "top") + theme_classic()

# Save PDF plots to Read_Distribution folder (root level for easy access)
nm_pdf1 <- file.path(ReadDistDirRoot, "Read_Distribution_Invariant_Norm_boxplot.pdf")
nm_pdf2 <- file.path(ReadDistDirRoot, "Read_Distribution_Invariant_Norm_histogram.pdf")
nm_pdf3 <- file.path(ReadDistDirRoot, "Read_Distribution_Invariant_Norm_density.pdf")

ggsave(nm_pdf1, p.1, width=8, height=6, dpi=300, units="in", device="pdf")
ggsave(nm_pdf2, p.2, width=8, height=6, dpi=300, units="in", device="pdf")
ggsave(nm_pdf3, p.3, width=8, height=6, dpi=300, units="in", device="pdf")
cat("Invariant gene normalized read distribution PDF plots saved:", nm_pdf1, nm_pdf2, nm_pdf3, "\n")

# Generate boxplot data for MultiQC interactive plot (invariant genes normalized data)
# MultiQC box/violin plots expect wide format: one column per sample
# Each column contains all gene expression values for that sample
sample_cols <- rownames(coldata_plots)
boxplot_data_norm_wide <- all.reads_z[, sample_cols]
colnames(boxplot_data_norm_wide) <- sample_cols

# Write boxplot data for MultiQC (wide format)
read_dist_norm_file <- paste0(file_prefix, ".read.distribution.normalized.txt")
write.table(boxplot_data_norm_wide, file = read_dist_norm_file, 
           sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)
cat("Read distribution boxplot data (normalized) for interactive MultiQC saved:", read_dist_norm_file, "\n")

cat("Invariant gene normalized read distribution plots saved in:", ReadDistDir, "\n")

# Create DESeq2 object using the filtered normalized matrix
cat("Creating DESeq2 object using filtered normalized matrix from GeneralNormalizer...\n")

# Ensure coldata matches the filtered normalized matrix samples
cat("Normalized matrix samples:", paste(head(colnames(normalized_matrix_filtered), 3), collapse = ", "), "...\n")
cat("Coldata samples:", paste(head(rownames(coldata), 3), collapse = ", "), "...\n")

# Verify that coldata and normalized matrix have matching samples
if (nrow(coldata) != ncol(normalized_matrix_filtered)) {
    stop("Mismatch between coldata rows (", nrow(coldata), ") and normalized matrix columns (", ncol(normalized_matrix_filtered), ")")
}

# Ensure sample names match between coldata and normalized matrix
if (!all(rownames(coldata) %in% colnames(normalized_matrix_filtered))) {
    cat("Warning: Some coldata samples not found in normalized matrix\n")
    cat("Missing samples:", paste(setdiff(rownames(coldata), colnames(normalized_matrix_filtered)), collapse = ", "), "\n")
}

# Try to add grouping information if sample names follow standard format
name_components_filtered <- strsplit(rownames(coldata), "_")
if (length(name_components_filtered) > 0 && all(sapply(name_components_filtered, length) > 1)) {
    # Add Group1 column based on sample name prefix (remove _R[0-9] part)
    coldata$Group1 <- gsub("_R[0-9].*", "", rownames(coldata))
    cat("Added Group1 column to coldata\n")
}

cat("Final coldata dimensions:", nrow(coldata), "x", ncol(coldata), "\n")

# Create DESeq2 object with the filtered normalized data
# Since the data is already normalized, we set the size factors to 1 (no additional normalization needed)
# Use the original coldata that matches the sample structure (like all_genes module does)
dds <- DESeqDataSetFromMatrix(countData=round(normalized_matrix_filtered), colData=coldata, design=~1)

# Set size factors to 1 since data is already normalized by GeneralNormalizer
sizeFactors(dds) <- rep(1, ncol(normalized_matrix_filtered))
cat("DESeq2 object created using filtered, pre-normalized data from GeneralNormalizer\n")

if (min(dim(count.table))<=1)  { # No point if only one sample, or one gene
    save(dds,file=DDSFile)
    saveRDS(dds, file=sub("\\.dds\\.RData$", ".rds", DDSFile))
    warning("Not enough samples or genes in counts file for PCA.", call.=FALSE)
    quit(save = "no", status = 0, runLast = FALSE)
}

# Apply variance stabilization
if (!opt$vst) {
    vst_name <- "rlog"
    cat("Applying rlog transformation...\n")
    rld      <- rlog(dds)
} else {
    vst_name <- "vst"
    cat("Applying VST transformation...\n")
    rld      <- varianceStabilizingTransformation(dds)
}

assay(dds, vst_name) <- assay(rld)
                               
################################################
################################################
## SAVE TRANSFORMED COUNTS                    ##
################################################
################################################

# Create VST/rlog transformed rlog_counts.txt (required by module)
# This uses the transformed data
transformed_counts <- assay(dds, vst_name)

cat("Creating VST/rlog transformed counts...\n")
cat("Original annotation rows:", nrow(original_annotation), "\n")
cat("Original annotation filtered rows:", nrow(original_annotation_filtered), "\n")
cat("VST/rlog counts dimensions:", nrow(transformed_counts), "x", ncol(transformed_counts), "\n")

# Create dataframe with gene_id from rownames
transformed_counts_df <- data.frame(
    gene_id = rownames(transformed_counts),
    transformed_counts,
    stringsAsFactors = FALSE,
    check.names = FALSE
)
                         
# Simple merge by gene_id for rlog counts too
cat("Adding annotation to VST/rlog transformed counts using merge by gene_id...\n")

# Load HOMER annotations if provided (reuse from previous section)
homer_annotation_rlog <- NULL
if (!is.null(opt$homer_annotation) && file.exists(opt$homer_annotation)) {
    cat("Loading HOMER annotations for rlog counts from:", opt$homer_annotation, "\n")
    homer_annotation_rlog <- read.delim(opt$homer_annotation, header=TRUE, stringsAsFactors=FALSE, comment.char="#")
    # Rename the first column to "Geneid" (HOMER often has malformed column name)
    colnames(homer_annotation_rlog)[1] <- "Geneid"
    cat("HOMER annotation column names:", paste(colnames(homer_annotation_rlog), collapse=", "), "\n")
    cat("HOMER annotation dimensions:", nrow(homer_annotation_rlog), "x", ncol(homer_annotation_rlog), "\n")
}

# Merge: first with original annotation, then with HOMER if available
if (!is.null(homer_annotation_rlog)) {
    # First left_join original annotation with HOMER
    cat("Merging original annotation with HOMER for rlog counts...\n")
    combined_annotation_rlog <- left_join(original_annotation_filtered, homer_annotation_rlog, by = c("Geneid","Chr","Start","End"))
    cat("Combined annotation for rlog dimensions:", nrow(combined_annotation_rlog), "x", ncol(combined_annotation_rlog), "\n")
    
    # Then left_join with rlog counts
    rlog_counts_output <- left_join(combined_annotation_rlog, transformed_counts_df, by = c("Geneid" = "gene_id"))

} else {
    # Simple left_join by gene_id
    rlog_counts_output <- left_join(original_annotation_filtered, transformed_counts_df, by = c("Geneid" = "gene_id"))
}

cat("Final rlog counts dimensions:", nrow(rlog_counts_output), "x", ncol(rlog_counts_output), "\n")
                 
# Write to organized location (base_dir)
rlog_file <- paste0(base_dir, "rlog_counts.txt")
write.table(rlog_counts_output, 
           file = rlog_file, 
           sep = "\t", 
           quote = FALSE, 
           row.names = FALSE, 
           col.names = TRUE)

# Also write to root level for Nextflow module output (with prefix)
rlog_file_root <- paste0(opt$outprefix, "_rlog_counts.txt")
write.table(rlog_counts_output, 
           file = rlog_file_root, 
           sep = "\t", 
           quote = FALSE, 
           row.names = FALSE, 
           col.names = TRUE)

# Verify files were created
if (file.exists(rlog_file)) {
    cat("VST/rlog transformed counts written to:", rlog_file, "\n")
    cat("File size:", file.size(rlog_file), "bytes\n")
} else {
    cat("ERROR: rlog_counts.txt file was not created in base_dir\n")
}

if (file.exists(rlog_file_root)) {
    cat("VST/rlog transformed counts written to root:", rlog_file_root, "\n")
    cat("File size:", file.size(rlog_file_root), "bytes\n")
    cat("Dimensions:", nrow(transformed_counts), "x", ncol(transformed_counts), "\n")
} else {
    cat("ERROR: root level rlog_counts.txt file was not created\n")
}
save(dds,file=DDSFile)
saveRDS(dds, file=sub("\\.dds\\.RData$", ".rds", DDSFile))

################################################
################################################
## PLOT QC                                    ##
################################################
################################################

##' PCA pre-processor
##'
##' Generate all the necessary information to plot PCA from a DESeq2 object
##' in which an assay containing a variance-stabilised matrix of counts is
##' stored. Copied from DESeq2::plotPCA, but with additional ability to
##' say which assay to run the PCA on.
##'
##' @param object The DESeq2DataSet object.
##' @param ntop number of top genes to use for principal components, selected by highest row variance.
##' @param assay the name or index of the assay that stores the variance-stabilised data.
##' @return A data.frame containing the projected data alongside the grouping columns.
##' A 'percentVar' attribute is set which includes the percentage of variation each PC explains,
##' and additionally how much the variation within that PC is explained by the grouping variable.
##' @author Gavin Kelly
plotPCA_vst <- function (object,  ntop = 500, assay=length(assays(object))) {
    rv         <- rowVars(assay(object, assay))
    select     <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
    pca        <- prcomp(t(assay(object, assay)[select, ]), center=TRUE, scale=FALSE)
    percentVar <- pca$sdev^2/sum(pca$sdev^2)
    df         <- cbind( as.data.frame(colData(object)), pca$x)
    

    
    #Order points so extreme samples are more likely to get label
    ord        <- order(abs(rank(df$PC1)-median(df$PC1)), abs(rank(df$PC2)-median(df$PC2)))
    df         <- df[ord,]
    attr(df, "percentVar") <- data.frame(PC=seq(along=percentVar), percentVar=100*percentVar)
    return(df)
}

cat("\n=== GENERATING QC PLOTS ===\n")

cat("Creating PNG plots only (no PDF generation) for invariant genes normalization\n")

# Suppress automatic Rplots.pdf generation
options(device = function() png(tempfile(), width=1, height=1))

# Ensure base directory exists
if (!dir.exists(base_dir)) {
    dir.create(base_dir, recursive = TRUE)
    cat("Created directory:", base_dir, "\n")
}

## SAMPLE CORRELATION HEATMAP
sampleDists      <- dist(t(assay(dds, vst_name)))
sampleDistMatrix <- as.matrix(sampleDists)
# Blues scale: dark blue (low distance/similar) to light blue (high distance/different)
colors           <- colorRampPalette( brewer.pal(9, "Blues") )(255)

# Create directory for quality control plots
QualityControlDir <- "Quality_Control/"
if (file.exists(QualityControlDir) == FALSE) {
    dir.create(QualityControlDir, recursive=TRUE)
}

# Generate PDF for sample distance heatmap
distance_pdf <- file.path(QualityControlDir, "Sample_Distance_Heatmap.pdf")
pdf(distance_pdf, width = 10, height = 8)
pheatmap(
    sampleDistMatrix,
    clustering_distance_rows=sampleDists,
    clustering_distance_cols=sampleDists,
    col=colors,
    main=paste("Euclidean distance between", vst_name, "of samples")
)
dev.off()
cat("Sample distance PDF saved:", distance_pdf, "\n")

## WRITE SAMPLE DISTANCES TO FILE (for MultiQC)
write.table(cbind(sample = rownames(sampleDistMatrix), sampleDistMatrix),file=paste0(file_prefix, ".sample.dists.txt"),
            row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)

## PEARSON CORRELATION HEATMAP
cat("Generating Pearson correlation heatmap...\n")
Pearson.corr <- cor(assay(dds, vst_name), method="pearson")  

# Create annotation for samples (use existing coldata groupings)
# Create annotation for samples (use coldata to match DESeq2 object)
annot_c <- coldata

# Debug information
cat("coldata dimensions:", nrow(coldata), "x", ncol(coldata), "\n")
cat("coldata row names:", paste(head(rownames(coldata), 3), collapse = ", "), "...\n")
cat("Available columns in coldata:", paste(colnames(coldata), collapse = ", "), "\n")

if ("Bio_replicates" %in% colnames(annot_c)) {
    annot_c <- annot_c[, c("Bio_replicates"), drop=FALSE]
    cat("Using Bio_replicates column for annotation\n")
} else {
    # Use first available grouping column or create simple annotation
    if (ncol(coldata) > 1) {
        annot_c <- coldata[, 2, drop=FALSE]  # Use first group column
        cat("Using first group column:", colnames(coldata)[2], "\n")
    } else {
        # Fix: Ensure the row names match the coldata structure
        # Create annotation that matches the actual samples in the DESeq2 object
        cat("Creating simple Sample_Type annotation\n")
        cat("Number of samples for annotation:", nrow(coldata), "\n")
        
        # Verify row names exist and are valid
        if (is.null(rownames(coldata)) || length(rownames(coldata)) == 0) {
            stop("coldata has no row names - this indicates a data structure issue")
        }
        
        annot_c <- data.frame(
            Sample_Type = rep("Sample", nrow(coldata)), 
            row.names = rownames(coldata),
            stringsAsFactors = FALSE
        )
        cat("Created annotation with dimensions:", nrow(annot_c), "x", ncol(annot_c), "\n")
    }
}

# Use consolidated Quality_Control directory (already created above)
# Generate PDF for storage (TSV not needed for MultiQC)
correlation_pdf <- file.path(QualityControlDir, "Sample_Correlation_Heatmap.pdf")

# Correlation data not saved to TSV (not needed for MultiQC)

# Generate PDF for organized storage
pdf(correlation_pdf, width = 10, height = 8)
pheatmap(Pearson.corr,
         annotation_col = annot_c,
         fontsize = 8,
         fontsize_row = 8,
         fontsize_col = 8,
         border_color = FALSE,
         col = rev(brewer.pal(8, "RdBu")),
         main = "Sample Correlation Heatmap")
dev.off()
cat("Sample correlation PDF saved:", correlation_pdf, "\n")

# Correlation data available in PNG plot and Pearson.corr matrix - no TSV needed

## PCA

# Create directory for PCA plots
# Use consolidated Quality_Control directory (already created above)

ntop <- c(500, Inf)
for (n_top_var in ntop) {
    pca.data      <- plotPCA_vst(dds, assay=vst_name, ntop=n_top_var)
    percentVar    <- round(attr(pca.data, "percentVar")$percentVar)
    plot_subtitle <- ifelse(n_top_var==Inf, "All genes", paste("Top", n_top_var, "genes"))
    pl <- ggplot(pca.data, aes(PC1, PC2, label=paste0(" ", sample, " "))) +
        geom_point() +
        geom_text(check_overlap=TRUE, vjust=0.5, hjust="inward") +
        xlab(paste0("PC1: ",percentVar[1],"% variance")) +
        ylab(paste0("PC2: ",percentVar[2],"% variance")) +
        labs(title = paste0("PCA Plot - ", plot_subtitle)) +
        theme(legend.position="top",
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            panel.border = element_rect(colour = "black", fill=NA, linewidth=1))
    
    # Save PCA data for MultiQC interactive plots (nf-core pattern) + PDF for storage
    if (n_top_var == Inf) {
        # TSV data for MultiQC interactive plot
        pca_data_file <- paste0(file_prefix, ".pca.vals.txt")
        pca_vals_mqc <- pca.data[,c("PC1","PC2")]
        colnames(pca_vals_mqc) <- paste0(colnames(pca_vals_mqc), ": ", percentVar[1:2], '% variance')
        pca_vals_mqc <- cbind(sample = rownames(pca_vals_mqc), pca_vals_mqc)
        write.table(pca_vals_mqc, file = pca_data_file, row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
        cat("PCA data for interactive MultiQC plot saved:", pca_data_file, "\n")
        
        # PDF for organized storage (Quality_Control folder)
        pca_pdf <- file.path(QualityControlDir, "PCA_Plot_All_Genes.pdf")
        ggsave(pca_pdf, pl, width = 10, height = 8, dpi = 300, units = "in", device = "pdf")
        cat("PCA plot PDF saved:", pca_pdf, "\n")
    } else {
        # TSV data for MultiQC interactive plot (500 genes)
        pca_data_file <- paste0(file_prefix, ".pca.top", n_top_var, ".vals.txt")
        pca_vals_mqc <- pca.data[,c("PC1","PC2")]
        colnames(pca_vals_mqc) <- paste0(colnames(pca_vals_mqc), ": ", percentVar[1:2], '% variance')
        pca_vals_mqc <- cbind(sample = rownames(pca_vals_mqc), pca_vals_mqc)
        write.table(pca_vals_mqc, file = pca_data_file, row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
        cat("PCA data for interactive MultiQC plot saved:", pca_data_file, "\n")
        
        # PDF for organized storage (Quality_Control folder)
        pca_pdf <- file.path(QualityControlDir, paste0("PCA_Plot_Top", n_top_var, "_Genes.pdf"))
        ggsave(pca_pdf, pl, width = 10, height = 8, dpi = 300, units = "in", device = "pdf")
        cat("PCA plot PDF saved:", pca_pdf, "\n")
    }

    if (decompose) {
        pc_names <- paste0("PC", attr(pca.data, "percentVar")$PC)
        long_pc <- reshape(pca.data, varying=pc_names, direction="long", sep="", timevar="component", idvar="pcrow")
        long_pc <- subset(long_pc, component<=5)
        long_pc_grp <- reshape(long_pc, varying=names(groupings), direction="long", sep="", timevar="grouper")
        long_pc_grp <- subset(long_pc_grp, grouper<=5)
        long_pc_grp$component <- paste("PC", long_pc_grp$component)
        long_pc_grp$grouper <- paste0(long_pc_grp$grouper, c("st","nd","rd","th","th")[long_pc_grp$grouper], " prefix")
        pl_decomp <- ggplot(long_pc_grp, aes(x=Group, y=PC)) +
            geom_point() +
            stat_summary(fun=mean, geom="line", aes(group = 1)) +
            labs(x=NULL, y=NULL, title=paste0("PCA Decomposed - ", plot_subtitle)) +
            facet_grid(component~grouper, scales="free_x") +
            scale_x_discrete(guide = guide_axis(n.dodge = 3))
        
        # Save decomposed PCA plots: PNG for MultiQC (root level), PDF for storage (organized)
        # PNG for MultiQC (root level only)
        decomp_png_name <- if(n_top_var==Inf) "pca_decomposed_all_genes_mqc.png" else paste0("pca_decomposed_top", n_top_var, "_genes_mqc.png")
        decomp_png <- paste0(opt$outprefix, "_", decomp_png_name)
        ggsave(decomp_png, pl_decomp, width = 12, height = 10, dpi = 300, units = "in")
        cat("Decomposed PCA plot PNG saved for MultiQC:", decomp_png, "\n")
        
        # PDF for organized storage (Quality_Control folder)
        decomp_pdf_name <- if(n_top_var==Inf) "PCA_Decomposed_All_Genes.pdf" else paste0("PCA_Decomposed_Top", n_top_var, "_Genes.pdf")
        decomp_pdf <- file.path(QualityControlDir, decomp_pdf_name)
        ggsave(decomp_pdf, pl_decomp, width = 12, height = 10, dpi = 300, units = "in", device = "pdf")
        cat("Decomposed PCA plot PDF saved:", decomp_pdf, "\n")
    }
} # at end of loop, we'll be using the user-defined ntop if any, else all genes

# Log information about decomposed plots
if (!decompose) {
    cat("Note: Decomposed PCA plots not generated (samples lack sufficient grouping structure)\n")
}

# PCA data already written above in the loop for MultiQC interactive plots
# PNG plots have been saved individually - no PDF device to close
cat("All individual PNG plots have been created for invariant genes normalization\n")






cat("\n=== INVARIANT GENES NORMALIZATION + DESEQ2 QC COMPLETE ===\n")
cat("Generated PNG files:\n")
cat("  - Sample distance heatmaps, correlation heatmaps, and PCA plots (individual PNG files)\n")
cat("  - Read distribution plots as PNG files\n")
cat("  - ", paste0(base_dir, file_prefix, ".pca.vals.txt"), " (PCA values)\n")
cat("  - ", paste0(base_dir, file_prefix, ".sample.dists.txt"), " (Sample distances)\n")

cat("  - ", final_norm_dir, " folder (GeneralNormalizer invariant gene analysis)\n")

cat("  - ", paste0(base_dir, "scaling_dat.txt"), " (Size factors table)\n")
cat("  - ", paste0(base_dir, "normalized_counts.txt"), " (Zero-count filtered GeneralNormalizer normalized matrix with annotation)\n")
cat("  - ", normalized_counts_root, " (Zero-count filtered GeneralNormalizer normalized matrix with annotation - root level)\n")
cat("  - ", paste0(base_dir, "rlog_counts.txt"), " (VST/rlog transformed matrix with annotation)\n")
cat("  - ", paste0(base_dir, "Read_Distribution/"), " folder (Read distribution QC plots)\n")

################################################
################################################
## R SESSION INFO                             ##
################################################
################################################

RLogFile <- "R_sessionInfo.log"

sink(RLogFile)
a <- sessionInfo()
print(a)
sink()

################################################
################################################
## TSV DATA FILES FOR INTERACTIVE MULTIQC    ##
################################################
################################################

cat("\n=== TSV DATA FILES GENERATED FOR INTERACTIVE MULTIQC ===\n")

# TSV data files are generated for MultiQC to create interactive plots
# PDF files provide high-quality static versions for publications and analysis

# Look for TSV files with the standard file prefix pattern
tsv_pattern <- paste0("^", gsub("\\.", "\\\\.", file_prefix), "\\.(pca|sample|read\\.distribution)\\.(vals|dists|normalized)\\.txt$")
tsv_files <- list.files(".", pattern = tsv_pattern, full.names = TRUE)

if (length(tsv_files) > 0) {
    cat("TSV data files generated for interactive MultiQC plots:\n")
    for (tsv_file in tsv_files) {
        if (file.exists(tsv_file)) {
            cat("  Generated:", tsv_file, "\n")
        }
    }
    cat("\nThese TSV files will be processed by MultiQC to create interactive plots.\n")
    cat("PDF files are organized in folders for high-quality analysis and publications.\n")
} else {
    cat("Warning: No TSV data files found - check script execution for errors.\n")
}

################################################
################################################
################################################
################################################
