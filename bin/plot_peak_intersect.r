#!/usr/bin/env Rscript

################################################
################################################
## LOAD LIBRARIES                             ##
################################################
################################################

library(optparse)
library(UpSetR)

################################################
################################################
## PARSE COMMAND-LINE PARAMETERS              ##
################################################
################################################

option_list <- list(make_option(c("-i", "--input_file"), type="character", default=NULL, help="Path to tab-delimited file containing two columns i.e sample1&sample2&sample3 indicating intersect between samples <TAB> set size.", metavar="path"),
                    make_option(c("-o", "--output_file"), type="character", default=NULL, help="Path to output file with '.pdf' extension.", metavar="path"))

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$input_file)){
    print_help(opt_parser)
    stop("Input file must be supplied.", call.=FALSE)
}
if (is.null(opt$output_file)){
    print_help(opt_parser)
    stop("Output pdf file must be supplied.", call.=FALSE)
}

OutDir <- dirname(opt$output_file)
if (file.exists(OutDir) == FALSE) {
    dir.create(OutDir,recursive=TRUE)
}

################################################
################################################
## PLOT DATA                                  ##
################################################
################################################

comb.dat <- read.table(opt$input_file,sep="\t",header=FALSE)
comb.vec <- comb.dat[,2]
comb.vec <- setNames(comb.vec,comb.dat[,1])
sets <- sort(unique(unlist(strsplit(names(comb.vec),split='&'))), decreasing = TRUE)

# Check if we have enough data for UpSet plot
# Need at least 2 sets AND at least 2 different intersections
num_sets <- length(sets)
num_intersects <- length(comb.vec)

# Debug output
cat("DEBUG: Input file:", opt$input_file, "\n")
cat("DEBUG: Number of sets:", num_sets, "\n")
cat("DEBUG: Sets:", paste(sets, collapse=", "), "\n")
cat("DEBUG: Number of intersections:", num_intersects, "\n")

if (num_sets < 2 || num_intersects < 2) {
    warning("Not enough data for UpSet plot (need at least 2 sets and 2 intersections, found ", num_sets, " sets and ", num_intersects, " intersections). Creating placeholder PDF.")
    pdf(opt$output_file,onefile=F,height=10,width=20)
    plot.new()
    text(0.5, 0.5, paste0("Not enough data for UpSet plot\n", 
                          "Found ", num_sets, " sets and ", num_intersects, " intersection(s)\n",
                          "Need at least 2 sets and 2 different intersections"), 
         cex=1.5)
    dev.off()
    quit(save="no", status=0)
}

nintersects = length(names(comb.vec))
if (nintersects > 70) {
    nintersects <- 70
    comb.vec <- sort(comb.vec, decreasing = TRUE)[1:70]
    sets <- sort(unique(unlist(strsplit(names(comb.vec),split='&'))), decreasing = TRUE)
}

pdf(opt$output_file,onefile=F,height=10,width=20)

upset(
    fromExpression(comb.vec),
    nsets = length(sets),
    nintersects = nintersects,
    sets = sets,
    keep.order = TRUE,
    sets.bar.color = "#56B4E9",
    point.size = 3,
    line.size = 1,
    mb.ratio = c(0.55, 0.45),
    order.by = "freq",
    number.angles = 30,
    text.scale = c(1.5, 1.5, 1.5, 1.5, 1.5, 1.2)
)

dev.off()

################################################
################################################
################################################
################################################
