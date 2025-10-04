#!/bin/sh
CONDA_BASE_PATH=$(conda info --base)
source $CONDA_BASE_PATH/etc/profile.d/conda.sh

#####################################
# The main script will run on the frontend (it is not consuming many resources).
# All processes will run in parellel.
####################################


# paths
NF_FOLDER=/mnt/ngs_ricerca/NEXTFLOW/
UsefulData=/mnt/ngs_ricerca/Software/
work_dir=/mnt/ngs_ricerca/NEXTFLOW/nextflow_temp/
conf_file=/mnt/ngs_ricerca/NEXTFLOW/local.config

######## PARAMS #########
sample_file=/mnt/ngs_ricerca/processed_data/Chipseq_test/scripts/samplesheet.csv
outdir=/mnt/ngs_ricerca/processed_data/Chipseq_test/output/
bowtie2_index=$UsefulData/reference_genome/hg38_UCSC/Bowtie2Index/
GTF=$UsefulData/reference_genome/gencode.v47.annotation.gtf  
bed_file=$UsefulData/reference_genome/gencode.v47.bed
genome_fasta=$UsefulData/reference_genome/GRCh38.primary_assembly.genome.fa
transcriptome_fasta=$UsefulData/reference_transcriptome/Gencode_v47_GRCh38_p14/gencode.v47.transcripts.fa
blacklist=$UsefulData/blacklist_regions/hg38-blacklist.v2.bed
FDR=0.05   #1e-05  #1e-15  #0.01  
########################

# Where will be stored work files, it can be deleted once the script is finished
mkdir -p $work_dir

project_name=Chipseq_test
work_dir_project=$work_dir/$project_name/

mkdir -p $work_dir_project
chmod -R 777 $work_dir_project
cd $work_dir_project


###--- Run the script ---###
conda activate nextflow

export NXF_ASSETS=$NF_FOLDER/Nextflow_pipeline

NXF_VER=25.04.7 nextflow run pdichiaro/chipseq -r main \
            --input $sample_file \
            --outdir $outdir \
            --fasta $genome_fasta \
            --gene_bed $bed_file \
            --gtf $GTF \
            --bowtie2_index $bowtie2_index \
            --blacklist $blacklist \
            --trimmer trimgalore \
            --extra_trimgalore_args '--quality 20 --stringency 3 --length 20' \
            --min_trimmed_reads 100000 \
            --read_length 50 \
            --aligner bowtie2 \
            --keep_dups False \
            --keep_multi_map False \
            --keep_blacklist False \
            --with_inputs False \
            --narrow_peak False \
            --macs_model False \
            --fragment_size 300 \
            --broad_cutoff 0.1 \
            --macs_fdr $FDR \
            --min_reps_consensus 2 \
            --normalization_method 'all_genes,invariant_genes' \
            --deseq2_vst True \
            --skip_fastqc False \
            --skip_trimming False \
            --skip_picard_metrics False \
            --skip_peak_qc False \
            --skip_peak_annotation False \
            --skip_consensus_peaks False \
            --skip_plot_profile False \
            --skip_deseq2_qc False \
            --skip_deeptools_norm False \
            --save_merged_fastq False \
            --save_macs_pileup False \
            --save_reference False \
            --save_align_intermeds False \
            --multiqc_title $project_name \
            -profile singularity \
            -c $conf_file \
            -process.echo \
            -resume
            
if [ $? -eq 0 ]; then
    echo "Nextflow finished successfully"
    scp -r {./report*,./timeline*,./trace*,./flowchart*} $outdir
    chmod -R 777 $work_dir_project 
    #rm -rf $work_dir_project 
else
    echo "Nextflow encountered an error"
fi

conda deactivate
