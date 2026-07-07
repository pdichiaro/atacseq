#!/usr/bin/env Rscript

################################################
################################################
## ALL GENES NORMALIZATION + DESEQ2 QC        ##
################################################
################################################

## This script performs standard DESeq2 normalization using median-of-ratios method
## followed by comprehensive DESeq2 QC analysis
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
    make_option(c("-p", "--outprefix"     ), type="character", default='all_genes_deseq2_qc', metavar="string" , help="Output prefix."                                                                         ),
    make_option(c("-v", "--vst"           ), type="logical"  , default=FALSE   , metavar="boolean", help="Run vst transform instead of rlog."                                                     ),
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
# MultiQC config expects files like: kallisto.deseq2.all_genes.*.txt
if (opt$quantifier != '' && opt$quantifier != 'null') {
    # Normalize quantifier name to lowercase for consistency with MultiQC config
    quantifier_lower <- tolower(opt$quantifier)
    # Special handling for STAR-based quantifiers
    if (quantifier_lower == "star_rsem") {
        file_prefix <- "star.rsem.deseq2.all_genes"
    } else if (quantifier_lower == "star_salmon") {
        file_prefix <- "star.salmon.deseq2.all_genes"
    } else if (quantifier_lower == "star_genome") {
        file_prefix <- "star.genome.deseq2.all_genes"
    } else {
        file_prefix <- paste0(quantifier_lower, ".deseq2.all_genes")
    }
} else {
    # Fallback to outprefix if no quantifier specified
    file_prefix <- opt$outprefix
}

cat("=== ALL GENES NORMALIZATION + DESEQ2 QC ===\n")
cat("File prefix for outputs:", file_prefix, "\n")
cat("Input file:", opt$count_file, "\n")
cat("Working directory:", getwd(), "\n")
cat("File exists:", file.exists(opt$count_file), "\n")
if (file.exists(opt$count_file)) {
    cat("File size:", file.size(opt$count_file), "bytes\n")
} else {
    cat("Available files in current directory:\n")
    files <- list.files(".", full.names = FALSE)
    for (f in head(files, 10)) {
        cat("  ", f, "\n")
    }
}

################################################
################################################
## READ IN COUNTS FILE                        ##
################################################
################################################

# Read featureCounts output with skip=1 to skip the command line header
# Also use comment.char to skip any lines starting with #
count.table           <- read.delim(file=opt$count_file, header=TRUE, row.names=NULL, skip=1, check.names=FALSE, comment.char="#")
rownames(count.table) <- count.table[,opt$id_col]
# Extract only count columns (skip annotation columns from featureCounts)
count.table           <- count.table[,opt$count_col:ncol(count.table),drop=FALSE]
colnames(count.table) <- gsub(opt$sample_suffix,"",colnames(count.table))
colnames(count.table) <- gsub(pattern='\\.$', replacement='', colnames(count.table))

cat("Count matrix dimensions:", nrow(count.table), "genes x", ncol(count.table), "samples\n")

################################################
################################################
## RUN DESEQ2 STANDARD NORMALIZATION + QC     ##
################################################
################################################

if (file.exists(opt$outdir) == FALSE) {
    dir.create(opt$outdir, recursive=TRUE)
}
setwd(opt$outdir)

cat("\n=== STARTING ALL GENES NORMALIZATION + DESEQ2 QC ===\n")

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

counts  <- count.table[,samples.vec,drop=FALSE]

################################################
## ZERO COUNT REMOVAL                        ##
################################################

# Remove genes with zero counts across all samples
cat("Removing genes with zero counts across all samples...\n")
cat("Starting with", nrow(counts), "genes\n")

# Simple zero count filter
keep_genes <- rowSums(counts) > 0
counts <- counts[keep_genes, ]

cat("Genes after zero count removal:", nrow(counts), "\n")
cat("Filtered out", sum(!keep_genes), "genes with zero counts\n")

# `design=~1` creates intercept-only model, equivalent to setting `blind=TRUE` for transformation.
dds     <- DESeqDataSetFromMatrix(countData=round(counts), colData=coldata, design=~1)

# Use DESeq2's standard median-of-ratios normalization
cat("Applying DESeq2 standard normalization (median-of-ratios method)...\n")
dds     <- estimateSizeFactors(dds)

cat("DESeq2 size factors calculated for", length(sizeFactors(dds)), "samples\n")
cat("Size factor range:", round(min(sizeFactors(dds)), 3), "to", round(max(sizeFactors(dds)), 3), "\n")

################################################
## SAVE SIZE FACTORS                          ##
################################################

# Create hierarchical output structure: quantifier/deseq2/normalization_method/
if (opt$quantifier != '') {
    base_dir <- paste0(opt$quantifier, "/deseq2/all_genes/")
} else {
    base_dir <- "deseq2/all_genes/"
}

cat("Quantifier parameter:", opt$quantifier, "\n")
cat("Creating output directory structure:", base_dir, "\n")
if (file.exists(base_dir) == FALSE) {
    dir.create(base_dir, recursive=TRUE)
}

# Save size factors summary table (text format for scaling_dat.txt)
size_factors_summary <- data.frame(
    sample = names(sizeFactors(dds)),
    size_factor = as.numeric(sizeFactors(dds)),
    stringsAsFactors = FALSE
)
SizeFactorsDir <- paste0(base_dir, "size_factors/")
if (file.exists(SizeFactorsDir) == FALSE) {
    dir.create(SizeFactorsDir, recursive=TRUE)
}
write.table(size_factors_summary, 
           file = paste0(SizeFactorsDir, opt$outprefix, ".size_factors_summary.txt"),
           sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

# Create required output files for Nextflow module
cat("Creating additional output files for module compatibility...\n")

# Create scaling_dat.txt (size factors in simple format)
scaling_dat <- data.frame(
    sample = names(sizeFactors(dds)),
    size_factor = as.numeric(sizeFactors(dds)),
    stringsAsFactors = FALSE
)
write.table(scaling_dat, file = paste0(base_dir, "scaling_dat.txt"), sep = "\t", 
           row.names = FALSE, col.names = TRUE, quote = FALSE)
# Also write to root for Nextflow
write.table(scaling_dat, file = "scaling_dat.txt", sep = "\t", 
           row.names = FALSE, col.names = TRUE, quote = FALSE)

# Create per-sample scaling factor files for safe DeepTools normalization
# This avoids grep parsing issues and is more robust
cat("Creating per-sample scaling factor files...\n")
scaling_factors_dir <- "scaling_factors"
dir.create(scaling_factors_dir, showWarnings = FALSE, recursive = TRUE)

for (i in 1:nrow(scaling_dat)) {
    sample_name <- scaling_dat$sample[i]
    scaling_factor <- scaling_dat$size_factor[i]
    
    # Write individual scaling factor file (just the numeric value)
    sample_file <- file.path(scaling_factors_dir, paste0(sample_name, "_scaling_factor.txt"))
    writeLines(as.character(scaling_factor), sample_file)
    
    cat("  - Created:", sample_file, "with scaling factor:", scaling_factor, "\n")
}

################################################
## SAVE NORMALIZED COUNT MATRIX               ##
################################################
                               
# Create normalized_counts.txt (normalized count matrix WITH annotation)
cat("Creating normalized_counts.txt...\n")
normalized_counts <- counts(dds, normalized = TRUE)

# Create dataframe with gene_id from rownames
deseq2_normalized_df <- data.frame(
    gene_id = rownames(normalized_counts),
    normalized_counts,
    stringsAsFactors = FALSE,
    check.names = FALSE
)
                               
# Read original featureCounts file to get annotation columns
cat("Reading annotation from featureCounts file...\n")
original_table <- read.delim(file=opt$count_file, header=TRUE, row.names=NULL, comment.char="#")
original_annotation <- original_table[, 1:(opt$count_col-1), drop = FALSE]

# Filter annotation to match normalized counts
genes_to_keep <- deseq2_normalized_df$gene_id
original_annotation_filtered <- original_annotation %>% filter(Geneid %in% genes_to_keep)
                               
cat("Original annotation rows:", nrow(original_annotation), "\n")
cat("Original annotation filtered rows:", nrow(original_annotation_filtered), "\n")
cat("Normalized counts dimensions:", nrow(normalized_counts), "x", ncol(normalized_counts), "\n")
                               
# Load HOMER annotations (*consensus_peaks.annotatePeaks.txt) if provided
homer_annotation <- NULL
if (!is.null(opt$homer_annotation) && file.exists(opt$homer_annotation)) {
    cat("Loading HOMER annotations from:", opt$homer_annotation, "\n")
    homer_annotation <- read.delim(opt$homer_annotation, header=TRUE, stringsAsFactors=FALSE, comment.char="#")
    # Rename the first column to "GeneID" (HOMER often has malformed column name)
    colnames(homer_annotation)[1] <- "Geneid"
    cat("HOMER annotation column names:", paste(colnames(homer_annotation), collapse=", "), "\n")
    cat("HOMER annotation dimensions:", nrow(homer_annotation), "x", ncol(homer_annotation), "\n")
}

# Merge: first with original annotation filtered, then with HOMER if available
if (!is.null(homer_annotation)) {
    # First merge original annotation with HOMER
    cat("Merging original annotation with HOMER annotations...\n")
    combined_annotation <- left_join(original_annotation_filtered, homer_annotation, by = c("Geneid","Chr","Start","End"))
    cat("Combined annotation dimensions:", nrow(combined_annotation), "x", ncol(combined_annotation), "\n")

    # Then merge with normalized counts
    normalized_counts_output <- left_join(combined_annotation, deseq2_normalized_df, by = c("Geneid" = "gene_id"))
    
} else {
    # Merge annotation with normalized counts
    normalized_counts_output <- left_join(original_annotation_filtered, deseq2_normalized_df,  by = c("Geneid" = "gene_id"))
}

cat("Merged normalized counts dimensions:", nrow(normalized_counts_output), "x", ncol(normalized_counts_output), "\n")

write.table(normalized_counts_output, file = paste0(base_dir, "normalized_counts.txt"), sep = "\t",
           row.names = FALSE, col.names = TRUE, quote = FALSE)
# Also write to root for Nextflow - with prefix
normalized_counts_file_root <- paste0(opt$outprefix, "_normalized_counts.txt")
write.table(normalized_counts_output, file = normalized_counts_file_root, sep = "\t",
           row.names = FALSE, col.names = TRUE, quote = FALSE)

cat("Size factors and normalized counts written immediately after estimateSizeFactors\n")

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

# Plot read distribution - Raw data
cat("Creating read distribution plots...\n")

# Check if we have grouping information, otherwise create simple sample groups
if (decompose && ncol(coldata) > 1) {
    # Use the first grouping column if available
    group_col <- names(coldata)[2]  # First column after 'sample'
    coldata$Bio_replicates <- coldata[[group_col]]
} else {
    # Create simple grouping based on sample names or just use "Sample" for all
    coldata$Bio_replicates <- "Sample"
}

# Raw count data (log transformed)
all.reads_t <- as.data.frame(log(count.table[,samples.vec]+1))
all.reads_tt <- all.reads_t[,rownames(coldata)] %>% rownames_to_column(var="gene") %>% as_tibble()            
gathered_all.reads <- gather(all.reads_tt, key = "samplename", value = "normalized_counts", -gene)
gathered_all.reads <- gathered_all.reads %>% 
  left_join(coldata %>% rownames_to_column("samplename") %>% dplyr::select(samplename, Bio_replicates), by = c("samplename")) %>%
  arrange(Bio_replicates) %>% 
  mutate(samplename = factor(samplename, levels = unique(samplename)))

# Create PNG plots for read distribution (raw data)
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
nm_pdf1 <- file.path(ReadDistDirRoot, "Read_Distribution_Raw_boxplot.pdf")
nm_pdf2 <- file.path(ReadDistDirRoot, "Read_Distribution_Raw_histogram.pdf") 
nm_pdf3 <- file.path(ReadDistDirRoot, "Read_Distribution_Raw_density.pdf")

ggsave(nm_pdf1, p.1, width=8, height=6, dpi=300, units="in", device="pdf")
ggsave(nm_pdf2, p.2, width=8, height=6, dpi=300, units="in", device="pdf")
ggsave(nm_pdf3, p.3, width=8, height=6, dpi=300, units="in", device="pdf")
cat("Read distribution raw PDF plots saved:", nm_pdf1, nm_pdf2, nm_pdf3, "\n")

# Generate boxplot data for MultiQC interactive plot (raw data)
# MultiQC box/violin plots expect wide format: one column per sample
# Each column contains all gene expression values for that sample
boxplot_data_raw_wide <- all.reads_t[, samples.vec]
colnames(boxplot_data_raw_wide) <- samples.vec

# Raw read distribution data not saved to file (not needed for MultiQC)

# Plot read distribution - Normalized data  
all.reads_d <- counts(dds, normalized = TRUE)  # Get DESeq2 normalized counts
all.reads_z <- as.data.frame(log(all.reads_d[,samples.vec]+1))
all.reads_zz <- all.reads_z[,rownames(coldata)] %>% rownames_to_column(var="gene") %>% as_tibble()            
gathered_all.reads <- gather(all.reads_zz, key = "samplename", value = "normalized_counts", -gene)
gathered_all.reads <- gathered_all.reads %>% 
  left_join(coldata %>% rownames_to_column("samplename") %>% dplyr::select(samplename, Bio_replicates), by = c("samplename")) %>%
  arrange(Bio_replicates) %>% 
  mutate(samplename = factor(samplename, levels = unique(samplename)))

# Create PNG plots for normalized read distribution
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
nm_pdf1 <- file.path(ReadDistDirRoot, "Read_Distribution_Norm_Filt_boxplot.pdf")
nm_pdf2 <- file.path(ReadDistDirRoot, "Read_Distribution_Norm_Filt_histogram.pdf")
nm_pdf3 <- file.path(ReadDistDirRoot, "Read_Distribution_Norm_Filt_density.pdf")

ggsave(nm_pdf1, p.1, width=8, height=6, dpi=300, units="in", device="pdf")
ggsave(nm_pdf2, p.2, width=8, height=6, dpi=300, units="in", device="pdf")
ggsave(nm_pdf3, p.3, width=8, height=6, dpi=300, units="in", device="pdf")
cat("Normalized read distribution PDF plots saved:", nm_pdf1, nm_pdf2, nm_pdf3, "\n")

# Generate boxplot data for MultiQC interactive plot (normalized data)
# MultiQC box/violin plots expect wide format: one column per sample
# Each column contains all gene expression values for that sample
boxplot_data_norm_wide <- all.reads_z[, samples.vec]
colnames(boxplot_data_norm_wide) <- samples.vec

# Write boxplot data for MultiQC (wide format)
read_dist_norm_file <- paste0(file_prefix, ".read.distribution.normalized.txt")
write.table(boxplot_data_norm_wide, file = read_dist_norm_file, 
           sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)
cat("Read distribution boxplot data (normalized) for interactive MultiQC saved:", read_dist_norm_file, "\n")

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
# This uses the transformed data WITH annotation
transformed_counts <- assay(dds, vst_name)

cat("Creating VST/rlog transformed counts...\n")
cat("Original annotation filtered rows:", nrow(original_annotation_filtered), "\n")
cat("VST/rlog counts dimensions:", nrow(transformed_counts), "x", ncol(transformed_counts), "\n")

# Create dataframe with gene_id from rownames
transformed_counts_df <- data.frame(
    gene_id = rownames(transformed_counts),
    transformed_counts,
    stringsAsFactors = FALSE,
    check.names = FALSE
)

# Load HOMER annotations for rlog if provided (reuse from earlier)
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
    # First merge original annotation with HOMER
    cat("Merging original annotation with HOMER for rlog counts...\n")
    combined_annotation_rlog <- left_join(original_annotation_filtered, homer_annotation_rlog, by = c("Geneid","Chr","Start","End"))
    cat("Combined annotation for rlog dimensions:", nrow(combined_annotation_rlog), "x", ncol(combined_annotation_rlog), "\n")
    
    # Then merge with rlog counts
    rlog_counts_output <- left_join(combined_annotation_rlog, transformed_counts_df, by = c("Geneid" = "gene_id"))

} else {
    # Merge annotation with transformed counts (reuse original_annotation filtered from earlier)
    rlog_counts_output <- left_join(original_annotation_filtered, transformed_counts_df, by = c("Geneid" = "gene_id"))
}

cat("Merged rlog counts dimensions:", nrow(rlog_counts_output), "x", ncol(rlog_counts_output), "\n")

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

cat("Creating PNG plots only (no PDF generation)\n")

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
write.table(cbind(sample = rownames(sampleDistMatrix), sampleDistMatrix),
            file=paste0(file_prefix, ".sample.dists.txt"),
            row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)
cat("Sample distance data (text format for MultiQC) saved:", paste0(file_prefix, ".sample.dists.txt"), "\n")

## PEARSON CORRELATION HEATMAP
cat("Generating Pearson correlation heatmap...\n")
Pearson.corr <- cor(assay(dds, vst_name), method="pearson")  

# Create annotation for samples (use existing coldata groupings)
annot_c <- coldata
if ("Bio_replicates" %in% colnames(annot_c)) {
    annot_c <- annot_c[, c("Bio_replicates"), drop=FALSE]
} else {
    # Use first available grouping column or create simple annotation
    if (decompose && ncol(coldata) > 1) {
        annot_c <- coldata[, 2, drop=FALSE]  # Use first group column
    } else {
        annot_c <- data.frame(Sample_Type = "Sample", row.names = rownames(coldata))
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

# PCA data already written above in the loop for MultiQC interactive plots





cat("\n=== ALL GENES NORMALIZATION + DESEQ2 QC COMPLETE ===\n")
cat("Generated files:\n")
cat("Interactive MultiQC data files (root level):\n")
cat("  - TSV data files for interactive plots (.pca.vals.txt, .sample.dists.txt)\n")
cat("PDF files for publications (organized folders):\n")
cat("  - Read_Distribution/ folder (read distribution QC plots)\n")
cat("  - Quality_Control/ folder (sample distance, correlation, and PCA plots)\n")
cat("Data files:\n")
cat("  - ", paste0(base_dir, "scaling_dat.txt"), " (Size factors table)\n")
cat("  - ", paste0(base_dir, "normalized_counts.txt"), " (Normalized count matrix with annotation)\n")
cat("  - ", normalized_counts_file_root, " (Normalized count matrix with annotation - root level)\n")
cat("  - ", paste0(base_dir, "rlog_counts.txt"), " (VST/rlog transformed matrix with annotation)\n")

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
