/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def valid_params = [
    aligners       : [ 'bowtie2' ]
]

// Validate input parameters
// WorkflowAtacseq.initialise(params, log, valid_params)  // Disabled: validation now handled by PIPELINE_INITIALISATION

// Check input path parameters to see if they exist
def checkPathParamList = [
    params.input, params.multiqc_config,
    params.fasta,
    params.gtf, params.gff, params.gene_bed,
    params.bowtie2_index,
    params.blacklist
]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters - moved inside workflow to avoid early exit on --help
// if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }
// if (params.rerpmsk) { ch_rerpmsk = file(params.rerpmsk) } else { exit 1, 'rerpmsk must be provided!' }

// Save AWS IGenomes file containing annotation version
def anno_readme = params.genomes[ params.genome ]?.readme
if (anno_readme && file(anno_readme).exists()) {
    file("${params.outdir}/genome/").mkdirs()
    file(anno_readme).copyTo("${params.outdir}/genome/")
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()

// Header files for MultiQC
ch_spp_nsc_header           = file("$projectDir/assets/multiqc/spp_nsc_header.txt", checkIfExists: true)
ch_spp_rsc_header           = file("$projectDir/assets/multiqc/spp_rsc_header.txt", checkIfExists: true)
ch_spp_correlation_header   = file("$projectDir/assets/multiqc/spp_correlation_header.txt", checkIfExists: true)
ch_peak_count_header        = file("$projectDir/assets/multiqc/peak_count_header.txt", checkIfExists: true)
ch_frip_score_header        = file("$projectDir/assets/multiqc/frip_score_header.txt", checkIfExists: true)
ch_peak_annotation_header   = file("$projectDir/assets/multiqc/peak_annotation_header.txt", checkIfExists: true)
ch_deseq2_pca_header        = file("$projectDir/assets/multiqc/deseq2_pca_header.txt", checkIfExists: true)
ch_deseq2_clustering_header = file("$projectDir/assets/multiqc/deseq2_clustering_header.txt", checkIfExists: true)
ch_deseq2_read_dist_header  = file("$projectDir/assets/multiqc/read_distribution_normalized_header.txt", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


include { FRIP_SCORE                                               } from '../modules/local/frip_score'
include { PLOT_MACS2_QC                                            } from '../modules/local/plot_macs2_qc'
include { PLOT_MACS2_QC_CONSENSUS                                  } from '../modules/local/plot_macs2_qc_consensus'
include { PLOT_HOMER_ANNOTATEPEAKS                                 } from '../modules/local/plot_homer_annotatepeaks'
include { PLOT_HOMER_ANNOTATEPEAKS_CONSENSUS                       } from '../modules/local/plot_homer_annotatepeaks_consensus'
include { PLOT_PEAK_INTERSECT_SAMPLES          } from '../modules/local/plot_peak_intersect_samples'
include { PLOT_CONDITION_INTERSECT             } from '../modules/local/plot_condition_intersect'
include { MACS2_CONSENSUS                     } from '../modules/local/macs2_consensus'
include { MACS2_CONSENSUS_BY_CONDITION        } from '../modules/local/macs2_consensus_by_condition'
include { ANNOTATE_BOOLEAN_PEAKS              } from '../modules/local/annotate_boolean_peaks'
// include { COUNT_NORM                          } from '../modules/local/count_normalization'  // Module not found
include { NORMALIZE_DESEQ2_QC_INVARIANT_GENES } from '../modules/local/normalize_deseq2_qc_invariant_genes'
include { NORMALIZE_DESEQ2_QC_ALL_GENES       } from '../modules/local/normalize_deseq2_qc_all_genes'
include { DESEQ2_TRANSFORM                    } from '../modules/local/deseq2_transform'
include { MULTIQC                             } from '../modules/local/multiqc'
include { MULTIQC_CUSTOM_PHANTOMPEAKQUALTOOLS } from '../modules/local/multiqc_custom_phantompeakqualtools'
include { MULTIQC_CUSTOM_PEAKS                } from '../modules/local/multiqc_custom_peaks'
include { BLACKLIST_LOG                       } from '../modules/local/blacklist_log'

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK         } from '../subworkflows/local/input_check'
include { PREPARE_GENOME      } from '../subworkflows/local/prepare_genome'
include { BAM_FILTER as BAM_FILTER_SUBWF } from '../subworkflows/local/bam_filter'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//

include { PICARD_MERGESAMFILES          } from '../modules/nf-core/modules/picard/mergesamfiles/main'
include { PICARD_COLLECTMULTIPLEMETRICS } from '../modules/nf-core/modules/picard/collectmultiplemetrics/main'
include { PHANTOMPEAKQUALTOOLS          } from '../modules/nf-core/modules/phantompeakqualtools/main'
include { DEEPTOOLS_BIGWIG              } from '../modules/local/deeptools_bw'
include { DEEPTOOLS_BIGWIG_NORM         } from '../modules/local/deeptools_bw_norm'
include { DEEPTOOLS_BIGWIG_NORM as DEEPTOOLS_BIGWIG_NORM_INVARIANT } from '../modules/local/deeptools_bw_norm'
include { DEEPTOOLS_BIGWIG_NORM as DEEPTOOLS_BIGWIG_NORM_ALL_GENES } from '../modules/local/deeptools_bw_norm'
include { DEEPTOOLS_COMPUTEMATRIX       } from '../modules/nf-core/modules/deeptools/computematrix/main'
include { DEEPTOOLS_PLOTPROFILE         } from '../modules/nf-core/modules/deeptools/plotprofile/main'
include { DEEPTOOLS_PLOTHEATMAP         } from '../modules/nf-core/modules/deeptools/plotheatmap/main'
include { DEEPTOOLS_PLOTFINGERPRINT     } from '../modules/nf-core/modules/deeptools/plotfingerprint/main'
include { KHMER_UNIQUEKMERS             } from '../modules/nf-core/modules/khmer/uniquekmers/main'
include { MACS2_CALLPEAK as MACS2_CALLPEAK_SINGLE          } from '../modules/nf-core/modules/macs2/callpeak/main'
include { MACS2_CALLPEAK as MACS2_CALLPEAK_MERGED          } from '../modules/nf-core/modules/macs2/callpeak/main'
include { SUBREAD_FEATURECOUNTS         } from '../modules/nf-core/modules/subread/featurecounts/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS   } from '../modules/nf-core/modules/custom/dumpsoftwareversions/main'

include { HOMER_ANNOTATEPEAKS as HOMER_ANNOTATEPEAKS_MACS2              } from '../modules/nf-core/modules/homer/annotatepeaks/main'
include { HOMER_ANNOTATEPEAKS as HOMER_ANNOTATEPEAKS_CONSENSUS          } from '../modules/nf-core/modules/homer/annotatepeaks/main'
include { HOMER_ANNOTATEPEAKS as HOMER_ANNOTATEPEAKS_CONSENSUS_CONDITION} from '../modules/nf-core/modules/homer/annotatepeaks/main'

//
// SUBWORKFLOW: Consisting entirely of nf-core/modules
//

include { FASTQ_FASTQC_UMITOOLS_TRIMGALORE } from '../subworkflows/nf-core/fastq_fastqc_umitools_trimgalore/main'

include { FASTQC_TRIMGALORE      } from '../subworkflows/nf-core/fastqc_trimgalore'
include { FASTQ_ALIGN_BOWTIE2    } from '../subworkflows/nf-core/fastq_align_bowtie2'
include { MARK_DUPLICATES_PICARD } from '../subworkflows/nf-core/mark_duplicates_picard'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow ATACSEQ {

    ch_versions = Channel.empty()
    
    // Check mandatory input samplesheet parameter
    if (!params.input) {
        error('Input samplesheet not specified! Please provide --input <samplesheet.csv>')
    }
    ch_input = Channel.fromPath(params.input, checkIfExists: true)
    
    // Create summary_params for MultiQC
    def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

    //
    // Get mitochondrial chromosome name from genome config
    //
    def mito_name = params.genomes && params.genome && params.genomes[params.genome] ? 
        params.genomes[params.genome].mito_name : null
    
    //
    // SUBWORKFLOW: Uncompress and prepare reference genome files
    //
    PREPARE_GENOME (
        params.aligner,
        mito_name
    )
    ch_versions = ch_versions.mix(PREPARE_GENOME.out.versions)
    
    //
    // Validate bowtie2 index is available (using view to avoid consuming the channel)
    //
    PREPARE_GENOME.out.bowtie2_index
        .ifEmpty {
            error "ERROR: Bowtie2 index not generated or provided. " +
                  "Please provide --bowtie2_index or ensure --fasta is provided for index generation."
        }
        .view { index -> "✓ Bowtie2 index available: ${index}" }


    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        ch_input
    )
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

    // Set the reads channel from INPUT_CHECK output
    INPUT_CHECK.out.reads
        .set { ch_reads }

    //
    // SUBWORKFLOW: Read QC, extract UMI and trim adapters with TrimGalore!
    //
    FASTQ_FASTQC_UMITOOLS_TRIMGALORE (
        ch_reads,
        params.skip_fastqc || params.skip_qc,
        false,
        false,
        params.skip_trimming,
        0,
        params.min_trimmed_reads
    )
    ch_filtered_reads_raw  = FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.reads
    ch_fastqc_raw_multiqc  = FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.fastqc_zip
    ch_fastqc_trim_multiqc = FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.trim_zip
    ch_trim_log_multiqc    = FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.trim_log
    ch_trim_read_count     = FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.trim_read_count
    ch_versions = ch_versions.mix(FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.versions)

    //
    // Validate and debug filtered reads channel
    //
    ch_filtered_reads_raw
        .map { meta, reads ->
            // Validate meta object has required fields
            if (!meta) {
                error "ERROR: meta object is null after trimming for reads: ${reads}"
            }
            if (!meta.containsKey('single_end')) {
                error "ERROR: meta object missing 'single_end' field for sample ${meta.id}"
            }
            if (!meta.containsKey('id')) {
                error "ERROR: meta object missing 'id' field"
            }
            
            log.info "✓ Sample ${meta.id} passed filtering: single_end=${meta.single_end}, reads=${reads.size()} file(s)"
            
            return [meta, reads]
        }
        .set { ch_filtered_reads }

    //
    // SUBWORKFLOW: Alignment with Bowtie2 & BAM QC
    //
    FASTQ_ALIGN_BOWTIE2 (
        ch_filtered_reads,
        PREPARE_GENOME.out.bowtie2_index,  // Now emits path directly, not [meta, path]
        false,  // save_unaligned
        false,  // sort_bam
        PREPARE_GENOME.out.fasta  // Pass real fasta reference
    )
    ch_genome_bam        = FASTQ_ALIGN_BOWTIE2.out.bam
    ch_genome_bam_index  = FASTQ_ALIGN_BOWTIE2.out.csi
    
    ch_samtools_stats    = FASTQ_ALIGN_BOWTIE2.out.stats
    ch_samtools_flagstat = FASTQ_ALIGN_BOWTIE2.out.flagstat
    ch_samtools_idxstats = FASTQ_ALIGN_BOWTIE2.out.idxstats
    ch_bowtie2_multiqc   = FASTQ_ALIGN_BOWTIE2.out.log_out  // Capture Bowtie2 logs for MultiQC

    ch_versions = ch_versions.mix(FASTQ_ALIGN_BOWTIE2.out.versions)

    //
    // MODULE: Merge resequenced BAM files - 
    // Sligtly off.. would be better to identify before hand which ones are to be merged.. this can be done from the sample_sheet
    // This is a point of collection to is stops before proceeding:
    // It removes the "T" and merges all of them later will remove the "R[0-9]" bit
    //
    // ch_genome_bam
    //     .map {
    //         meta, bam ->
    //             new_id = meta.id - ~/_T\d+/
    //             [  meta + [id: new_id], bam ] 
    //     }
    
    ch_genome_bam
    .map { meta, bam ->
        // Use regex to find the last underscore and remove any text from that point onwards
        def new_id = meta.id.replaceAll(/_[^_]+$/, "")
        [meta + [id: new_id], bam]
    }
        .groupTuple(by: [0])
        .map { 
            it ->
                [ it[0], it[1].flatten() ] 
        }
        .set { ch_sort_bam }
    
    ch_sort_bam.view()

    PICARD_MERGESAMFILES (
        ch_sort_bam
    )
    ch_versions = ch_versions.mix(PICARD_MERGESAMFILES.out.versions.first().ifEmpty(null))

    //
    // SUBWORKFLOW: Mark duplicates & filter BAM files after merging
    //
    MARK_DUPLICATES_PICARD (
        PICARD_MERGESAMFILES.out.bam
    )
    ch_versions = ch_versions.mix(MARK_DUPLICATES_PICARD.out.versions)

    //
    // SUBWORKFLOW: Filter BAM file with BamTools 
    //
    
    BAM_FILTER_SUBWF (
        MARK_DUPLICATES_PICARD.out.bam.join(MARK_DUPLICATES_PICARD.out.bai, by: [0]),
        PREPARE_GENOME.out.filtered_bed.first()
    )
    ch_versions = ch_versions.mix(BAM_FILTER_SUBWF.out.versions.first().ifEmpty(null))

    //
    // MODULE: Generate filtering metrics log (blacklist + duplicates + other filters)
    //
    if (params.blacklist) {
        // Join all channels by meta to ensure correct pairing
        ch_blacklist_input = MARK_DUPLICATES_PICARD.out.bam
            .join(MARK_DUPLICATES_PICARD.out.bai, by: [0])
            .join(BAM_FILTER_SUBWF.out.bam, by: [0])
            .join(BAM_FILTER_SUBWF.out.bai, by: [0])
            .map { meta, bam_before, bai_before, bam_after, bai_after ->
                tuple(meta, bam_before, bai_before, bam_after, bai_after)
            }
        
        BLACKLIST_LOG (
            ch_blacklist_input,
            PREPARE_GENOME.out.filtered_bed.first(),
            PREPARE_GENOME.out.blacklist.first()
        )
        ch_versions = ch_versions.mix(BLACKLIST_LOG.out.versions.first().ifEmpty(null))
    }

    //
    // MODULE: Picard post alignment QC
    //
    ch_picardcollectmultiplemetrics_multiqc = Channel.empty()
    if (!params.skip_picard_metrics) {
        PICARD_COLLECTMULTIPLEMETRICS (
            BAM_FILTER_SUBWF.out.bam,
            PREPARE_GENOME.out.fasta,
            []
        )
        ch_picardcollectmultiplemetrics_multiqc = PICARD_COLLECTMULTIPLEMETRICS.out.metrics
        ch_versions = ch_versions.mix(PICARD_COLLECTMULTIPLEMETRICS.out.versions.first())
    }

    //
    // MODULE: Phantompeaktools strand cross-correlation and QC metrics
    //
    PHANTOMPEAKQUALTOOLS (
        BAM_FILTER_SUBWF.out.bam
    )
    ch_versions = ch_versions.mix(PHANTOMPEAKQUALTOOLS.out.versions.first())

    //
    // MODULE: MultiQC custom content for Phantompeaktools
    //
    MULTIQC_CUSTOM_PHANTOMPEAKQUALTOOLS (
        PHANTOMPEAKQUALTOOLS.out.spp.join(PHANTOMPEAKQUALTOOLS.out.rdata, by: [0]),
        ch_spp_nsc_header,
        ch_spp_rsc_header,
        ch_spp_correlation_header
    )


    //
    // Create channels: [ meta, [ ip_bam, control_bam ] [ ip_bai, control_bai ] ]
    // Differently from standard nf-core atacseq we can evaluate the possibility to run the chip-seq w/o inputs
    // This needs to be assessed on the fly i.e. check if there ar
    
    BAM_FILTER_SUBWF
        .out
        .bam
        .join(BAM_FILTER_SUBWF.out.bai, by: [0])
        .set { ch_genome_bam_bai }

    //
    // MODULE: deepTools plotFingerprint thi will assess in sample only
    //
    ch_deeptoolsplotfingerprint_multiqc = Channel.empty()
    if (!params.skip_plot_fingerprint ) {
        DEEPTOOLS_PLOTFINGERPRINT (
            ch_genome_bam_bai
        ) 
        ch_deeptoolsplotfingerprint_multiqc = DEEPTOOLS_PLOTFINGERPRINT.out.matrix
        ch_versions = ch_versions.mix(DEEPTOOLS_PLOTFINGERPRINT.out.versions.first())
    }
    
    // Create channels: [ meta, ip_bam, ([] for control_bam) ]
    // Always pass empty array as control (no input controls support)
    ch_genome_bam_bai
        .map {
            meta, bam, bai -> 
                [ meta , bam, [] ]
        }
        .set { ch_ip_control_bam }
    
    // Merge all bams by antibody: from meta,bam,control 
    ch_ip_control_bam
        .map {
            meta, bam1, bam2 ->
            def new_meta = meta.clone()
            new_meta.id =  meta.antibody
            [new_meta, bam1, bam2]
        }
        .groupTuple(by: 0)
        .map {
            meta, bam1, bam2 ->
                [ meta , bam1, [] ]
        }
        .set { ch_antibody_bam }
    
    // 
    // MODULE: Calculute genome size with khmer
    //
    ch_macs_gsize                     = Channel.empty()
    ch_custompeaks_frip_multiqc       = Channel.empty()
    ch_custompeaks_count_multiqc      = Channel.empty()
    ch_plothomerannotatepeaks_multiqc = Channel.empty()
    ch_subreadfeaturecounts_multiqc   = Channel.empty()
    ch_macs_gsize = params.macs_gsize
     
    if (!params.macs_gsize) {
        KHMER_UNIQUEKMERS (
            PREPARE_GENOME.out.fasta,
            params.read_length
        )
        ch_macs_gsize = KHMER_UNIQUEKMERS.out.kmers.map { it.text.trim() }
    }

    //
    // MODULE: Call peaks with MACS2
    // MACS2_CALLPEAK_SINGLE runs on individual samples (ch_ip_control_bam)
    //
    MACS2_CALLPEAK_SINGLE (
         ch_ip_control_bam,
         ch_macs_gsize
    )
    ch_versions = ch_versions.mix(MACS2_CALLPEAK_SINGLE.out.versions.first())

    //
    // Filter out samples with 0 MACS2 peaks called with warning system
    //
    MACS2_CALLPEAK_SINGLE
        .out
        .peak
        .branch { meta, peaks ->
            passed: peaks.size() > 0
                return [meta, peaks]
            failed: true
                return [meta, peaks]
        }
        .set { ch_macs2_branched }

    // Emit warning for each sample with zero peaks
    ch_macs2_branched
        .failed
        .subscribe { meta, peaks ->
            log.warn """
            ╔════════════════════════════════════════════════════════════════════════════════╗
            ║                          ⚠️  MACS2 ZERO PEAKS WARNING                          ║
            ╚════════════════════════════════════════════════════════════════════════════════╝
            
            Sample '${meta.id}' produced 0 peaks from MACS2 peak calling.
            This sample will be excluded from downstream analysis.
            
            Possible causes and solutions:
            
            1. Poor ChIP enrichment
               → Check ATAC-seq quality metrics (FRiP, NSC, RSC scores)
               → Verify antibody quality and ChIP protocol
            
            2. Insufficient sequencing depth
               → Current depth may be too low for peak detection
               → Consider deeper sequencing (≥20M reads for TF, ≥40M for histone marks)
            
            3. Overly stringent MACS2 parameters
               → Try adjusting q-value threshold (default: 0.05)
               → Use --broad flag for broad histone marks
               → Adjust --mfold parameter
            
            4. Poor quality control/input sample
               → Check if control sample has issues
               → Verify control matches treatment conditions
            
            5. Wrong genome size parameter
               → Verify --macs_gsize matches your genome
            
            6. Biological factors
               → Low/absent protein binding in experimental conditions
               → Check positive controls
            
            ════════════════════════════════════════════════════════════════════════════════
            """.stripIndent()
        }

    // Check if ALL samples failed - this is critical error
    ch_macs2_branched
        .passed
        .count()
        .subscribe { count ->
            if (count == 0) {
                log.error """
                ╔════════════════════════════════════════════════════════════════════════════════╗
                ║                      🔴 CRITICAL: ALL SAMPLES FAILED                           ║
                ╚════════════════════════════════════════════════════════════════════════════════╝
                
                ALL samples produced 0 peaks from MACS2 peak calling!
                The pipeline cannot continue with downstream analysis.
                
                IMMEDIATE ACTIONS REQUIRED:
                
                1. Review MACS2 parameters:
                   → Check --macs_gsize parameter
                   → Review q-value/p-value thresholds
                   → Consider --broad flag for histone marks
                
                2. Verify input data quality:
                   → Check sequencing depth (FastQC reports)
                   → Review alignment rates
                   → Inspect control samples
                
                3. Check ATAC-seq quality:
                   → Review PhantomPeakQualTools metrics
                   → Check cross-correlation plots
                   → Verify FRiP scores (if available)
                
                4. Review experimental design:
                   → Verify antibody specificity
                   → Check ChIP protocol efficiency
                   → Ensure proper controls
                
                Please address these issues before re-running the pipeline.
                ════════════════════════════════════════════════════════════════════════════════
                """.stripIndent()
            } else {
                log.info "✅ MACS2 peak calling successful for ${count} sample(s)"
            }
        }

    // Use only samples that passed
    ch_macs2_branched
        .passed
        .set { ch_macs2_peaks }

    // If is narrow we call high conf summits by merging all BAMS:

    MACS2_CALLPEAK_MERGED(
        ch_antibody_bam,
        ch_macs_gsize
    )

    // Create channels: [ meta, ip_bam, peaks ]
    // Join individual samples with their peaks
    ch_ip_control_bam
        .join(ch_macs2_peaks, by: [0])
        .map { 
            it -> 
                [ it[0], it[1], it[3] ] 
        }
        .set { ch_ip_bam_peaks }


    //
    // MODULE: Calculate FRiP score
    //
    FRIP_SCORE (
        ch_ip_bam_peaks
    )
    ch_versions = ch_versions.mix(FRIP_SCORE.out.versions.first())

    // Create channels: [ meta, peaks, frip ]
    ch_ip_bam_peaks
        .join(FRIP_SCORE.out.txt, by: [0])
        .map { 
            it -> 
                [ it[0], it[2], it[3] ] 
        }
        .set { ch_ip_peaks_frip }

    //
    // MODULE: FRiP score custom content for MultiQC
    //
    MULTIQC_CUSTOM_PEAKS (
        ch_ip_peaks_frip,
        ch_peak_count_header,
        ch_frip_score_header
    )
    ch_custompeaks_frip_multiqc  = MULTIQC_CUSTOM_PEAKS.out.frip
    ch_custompeaks_count_multiqc = MULTIQC_CUSTOM_PEAKS.out.count

    if (!params.skip_peak_annotation) {
        //
        // MODULE: Annotate peaks with MACS2
        //
        HOMER_ANNOTATEPEAKS_MACS2 (
            ch_macs2_peaks,
            PREPARE_GENOME.out.fasta,
            PREPARE_GENOME.out.gtf
        )
        ch_versions = ch_versions.mix(HOMER_ANNOTATEPEAKS_MACS2.out.versions.first())

        if (!params.skip_peak_qc) {
            //
            // MODULE: MACS2 QC plots with R
            //
            PLOT_MACS2_QC (
                ch_macs2_peaks
                    .map { meta, peaks -> peaks }
                    .collect()
            )
            ch_versions = ch_versions.mix(PLOT_MACS2_QC.out.versions)

            //
            // MODULE: Peak annotation QC plots with R
            //
            PLOT_HOMER_ANNOTATEPEAKS (
                HOMER_ANNOTATEPEAKS_MACS2.out.txt
                    .map { meta, txt -> txt }
                    .collect(),
                ch_peak_annotation_header,
                "_peaks.annotatePeaks.txt"
            )
            ch_plothomerannotatepeaks_multiqc = PLOT_HOMER_ANNOTATEPEAKS.out.tsv
            ch_versions = ch_versions.mix(PLOT_HOMER_ANNOTATEPEAKS.out.versions)
        }
    }

    //
    // MODULE: Plot peak intersections for all individual samples (before consensus)
    // This creates an UpSet plot showing overlaps between all individual sample peaks
    // grouped by antibody
    //
    if (!params.skip_peak_qc) {
        ch_macs2_peaks
            .map { 
                meta, peak ->
                    [ meta.antibody, peak ]
            }
            .groupTuple()
            .set { ch_antibody_individual_peaks }

        PLOT_PEAK_INTERSECT_SAMPLES (
            ch_antibody_individual_peaks
        )
        ch_versions = ch_versions.mix(PLOT_PEAK_INTERSECT_SAMPLES.out.versions.first())
    }

    //
    //  Consensus peaks analysis
    //  Here the aim is to generate a global Consensus and a "By_Condition" consensus
    //  Consider selecting by IDR score as best ENCODE practice:
    //  

    ch_macs2_consensus_bed_lib   = Channel.empty()
    ch_macs2_consensus_txt_lib   = Channel.empty()
    ch_deseq2_all_multiqc        = Channel.empty()

    // It makes by default a consensus - this is used to quantify and compute scaling FACTORS:

    // STEP 1: Create consensus peaks BY CONDITION (intermediate files)
    // Extract condition from sample ID and group peaks by condition+antibody
    // Example: TLBR2_shMCM5_CT_REP1 -> group_id = TLBR2_shMCM5_CT
    // Example: TLBR2_shMCM5_CT_RNAseA_REP1 -> group_id = TLBR2_shMCM5_CT_RNAseA
    // Example: TLBR2_shMCM5_DOX_RNAseA_REP2 -> group_id = TLBR2_shMCM5_DOX_RNAseA
    
    ch_macs2_peaks
        .map { 
            meta, peak ->
                // Extract group identifier by removing ONLY the replicate suffix (_REP{N})
                // This preserves all condition information including treatment suffixes like _RNAseA
                // Handles both patterns: _REP\d+_T\d+ (old format) and _REP\d+ (current format)
                def group_id = meta.id.replaceAll(/_REP\d+(_T\d+)?$/, '')
                
                // Debug: verify grouping is correct
                println "════════════════════════════════════════════════════════════════"
                println "🔍 CONDITION GROUPING DEBUG:"
                println "   Original sample ID: ${meta.id}"
                println "   → Condition group:  ${group_id}"
                println "   → Antibody:         ${meta.antibody}"
                println "   → Peak file:        ${peak.getName()}"
                println "════════════════════════════════════════════════════════════════"
                
                [ group_id, meta.antibody, peak ] 
        }
        .groupTuple(by: 0)
        .map {
            group_id, antibodies, peaks ->
                def meta_new = [:]
                meta_new.id = group_id
                meta_new.antibody = antibodies[0]  // All should have same antibody
                
                // Debug: verify final grouping with all replicate names
                println "════════════════════════════════════════════════════════════════"
                println "✅ FINAL CONDITION GROUP:"
                println "   Condition:     ${group_id}"
                println "   Antibody:      ${antibodies[0]}"
                println "   Num replicates: ${peaks.size()}"
                println "   Peak files:"
                peaks.eachWithIndex { peak, idx ->
                    println "      ${idx+1}. ${peak.getName()}"
                }
                println "════════════════════════════════════════════════════════════════"
                
                [ meta_new, peaks ] 
        }
        .set { ch_condition_peaks }
    
    //
    // MODULE: Generate consensus peaks BY CONDITION (e.g., WT_BCATENIN, NAIVE_BCATENIN)
    // These are intermediate files published to consensus_peaks/{antibody}/by_condition/
    //
    MACS2_CONSENSUS_BY_CONDITION (
        ch_condition_peaks
    )
    ch_versions = ch_versions.mix(MACS2_CONSENSUS_BY_CONDITION.out.versions)

    //
    // MODULE: Annotate and plot QC for consensus peaks BY CONDITION
    //
    if (!params.skip_peak_annotation) {
        HOMER_ANNOTATEPEAKS_CONSENSUS_CONDITION (
            MACS2_CONSENSUS_BY_CONDITION.out.bed,
            PREPARE_GENOME.out.fasta,
            PREPARE_GENOME.out.gtf
        )
        ch_versions = ch_versions.mix(HOMER_ANNOTATEPEAKS_CONSENSUS_CONDITION.out.versions)

        if (!params.skip_peak_qc) {
            //
            // MODULE: MACS2 QC plots for consensus peaks by condition
            //
            PLOT_MACS2_QC_CONSENSUS (
                MACS2_CONSENSUS_BY_CONDITION.out.peaks
                    .map { meta, peaks -> [ meta.antibody, meta, peaks ] }
                    .groupTuple(by: 0)
                    .map { antibody, metas, peaks -> 
                        def meta_new = metas[0].clone()
                        meta_new.id = antibody
                        [ meta_new, peaks ]
                    }
            )
            ch_versions = ch_versions.mix(PLOT_MACS2_QC_CONSENSUS.out.versions)

            //
            // MODULE: Peak annotation QC plots for consensus peaks by condition
            //
            PLOT_HOMER_ANNOTATEPEAKS_CONSENSUS (
                HOMER_ANNOTATEPEAKS_CONSENSUS_CONDITION.out.txt
                    .map { meta, txt -> [ meta.antibody, meta, txt ] }
                    .groupTuple(by: 0)
                    .map { antibody, metas, txts -> 
                        def meta_new = metas[0].clone()
                        meta_new.id = antibody
                        [ meta_new, txts ]
                    },
                ch_peak_annotation_header,
                ".condition_consensus_peaks.annotatePeaks.txt"
            )
            ch_versions = ch_versions.mix(PLOT_HOMER_ANNOTATEPEAKS_CONSENSUS.out.versions)

            //
            // MODULE: Plot condition intersections (UpSet plot across conditions, not replicates)
            // Only generates plot if antibody has >= 2 conditions
            //
            MACS2_CONSENSUS_BY_CONDITION
                .out
                .peaks
                .map { meta, peak -> [ meta.antibody, meta.id, peak ] }
                .groupTuple(by: 0)
                .set { ch_condition_peaks_for_plot }
            
            PLOT_CONDITION_INTERSECT (
                ch_condition_peaks_for_plot
            )
            ch_versions = ch_versions.mix(PLOT_CONDITION_INTERSECT.out.versions.first())
        }
    }

    // STEP 2: Merge condition consensus peaks by ANTIBODY for final analysis
    // Group the condition-level peaks by antibody for the final merge
    // Example: WT_BCATENIN + NAIVE_BCATENIN -> BCATENIN
    
    MACS2_CONSENSUS_BY_CONDITION
        .out
        .peaks  // Use .peaks (narrowPeak/broadPeak format) instead of .bed for proper format
        .map { 
            meta, peak_file ->
                // Use antibody from metadata instead of extracting from ID
                [ meta.antibody, meta.id, peak_file ] 
        }
        .groupTuple(by: 0)
        .map {
            antibody, group_ids, peak_files ->
                def meta_new = [:]
                meta_new.id = antibody
                meta_new.multiple_groups = group_ids.size() > 1
                meta_new.replicates_exist = true  // Conditions already have consensus from replicates
                [ meta_new, peak_files ] 
        }
        .set { ch_antibody_peaks }
    
    ch_antibody_peaks.view()
    //
    //  MODULE: Generate final consensus peaks by ANTIBODY across all conditions
    //  This merges condition-level consensus (e.g., WT_BCATENIN + NAIVE_BCATENIN -> BCATENIN)
    //  NOTE: We set min_reps_consensus = 1 for this step because each input is already
    //        a consensus from biological replicates (filtered at condition level)
    //

    MACS2_CONSENSUS ( 
        ch_antibody_peaks
    )
    ch_macs2_consensus_bed_lib = MACS2_CONSENSUS.out.bed
    ch_macs2_consensus_txt_lib = MACS2_CONSENSUS.out.txt
    ch_versions = ch_versions.mix(MACS2_CONSENSUS.out.versions)

    if (!params.skip_peak_annotation) {
        //
        // MODULE: Annotate consensus peaks
        //
        HOMER_ANNOTATEPEAKS_CONSENSUS (
            MACS2_CONSENSUS.out.bed,
            PREPARE_GENOME.out.fasta,
            PREPARE_GENOME.out.gtf
        )
        ch_versions = ch_versions.mix(HOMER_ANNOTATEPEAKS_CONSENSUS.out.versions)
        //
        // MODULE: Add boolean fields to annotated consensus peaks to aid filtering
        //
        ANNOTATE_BOOLEAN_PEAKS (
            MACS2_CONSENSUS.out.boolean_txt.join(HOMER_ANNOTATEPEAKS_CONSENSUS.out.txt, by: [0]),
        )
        ch_versions = ch_versions.mix(ANNOTATE_BOOLEAN_PEAKS.out.versions)
    }

    // Create channels: [ antibody, [ ip_bams ] ]
    ch_ip_control_bam
        .map { 
            meta, ip_bam, control_bam ->
                [ meta.antibody, ip_bam ]
        }
        .groupTuple()
        .set { ch_antibody_bams }
    

    // Create channels: [ meta, [ ip_bams ], saf ]
    MACS2_CONSENSUS
        .out
        .saf
        .map { 
            meta, saf -> 
                [ meta.id, meta, saf ] 
        }
        .join(ch_antibody_bams)
        .map {
            antibody, meta, saf, bams ->
                [ meta, bams.flatten().sort(), saf ]
        }
        .set { ch_saf_bams }

    //
    // MODULE: Quantify peaks across samples with featureCounts
    //
    SUBREAD_FEATURECOUNTS (
        ch_saf_bams
    )
    ch_subreadfeaturecounts_multiqc = SUBREAD_FEATURECOUNTS.out.summary
    ch_versions = ch_versions.mix(SUBREAD_FEATURECOUNTS.out.versions.first())

    //
    // Normalize samples with Gualdrini et al. 2016 Method compute scaling factor and generate a channel - [meta , bam, scaling]
    // Here we compute the scaling factor used to normalize the counts and the bigWigs
    // Generate counts as SUBREAD so norm counts rounded for DESEQ2 to work
    // out has: counts, scalings, 
    //
    // MODULE: Compute normalization and quality plots with conditional methods:
    //
    def normalization_methods = params.normalization_method instanceof List ? 
        params.normalization_method : params.normalization_method.split(',').collect{it.trim()}
    
    ch_deseq2_all_multiqc        = Channel.empty()
    ch_size_factors              = Channel.empty()
    ch_scaling_factors_all       = Channel.empty()
    ch_deseq2_raw_files          = Channel.empty()
    ch_normalization_versions    = Channel.empty()
    ch_normalization_scaling_factors = Channel.empty()
    
    //
    // Use consensus peaks annotation instead of gene annotation for normalization
    // Extract the first consensus peaks annotation file (e.g., pRPA.consensus_peaks.annotatePeaks.txt)
    // Only available if peak annotation was not skipped
    //
    def ch_consensus_annotation = !params.skip_peak_annotation 
        ? HOMER_ANNOTATEPEAKS_CONSENSUS.out.txt
            .map { meta, txt -> txt }
            .first()
        : Channel.empty()
    
    //
    // MODULE: Invariant genes normalization (stable genes only)
    //
    if (normalization_methods.contains('invariant_genes')) {
        NORMALIZE_DESEQ2_QC_INVARIANT_GENES (
            SUBREAD_FEATURECOUNTS.out.counts.map { meta, counts -> counts },
            "featureCounts",
            ch_consensus_annotation
        )
        
        // Collect raw files for DESEQ2_TRANSFORM
        ch_deseq2_raw_files = ch_deseq2_raw_files.mix(NORMALIZE_DESEQ2_QC_INVARIANT_GENES.out.read_dist_norm_txt)
        ch_deseq2_raw_files = ch_deseq2_raw_files.mix(NORMALIZE_DESEQ2_QC_INVARIANT_GENES.out.sample_distances_txt)
        ch_deseq2_raw_files = ch_deseq2_raw_files.mix(NORMALIZE_DESEQ2_QC_INVARIANT_GENES.out.pca_all_genes_txt)
        ch_deseq2_raw_files = ch_deseq2_raw_files.mix(NORMALIZE_DESEQ2_QC_INVARIANT_GENES.out.pca_top_genes_txt)
        
        ch_normalization_versions = ch_normalization_versions.mix(NORMALIZE_DESEQ2_QC_INVARIANT_GENES.out.versions)
        ch_normalization_scaling_factors = ch_normalization_scaling_factors.mix(NORMALIZE_DESEQ2_QC_INVARIANT_GENES.out.scaling_factors)
        ch_scaling_factors_all = ch_scaling_factors_all.mix(NORMALIZE_DESEQ2_QC_INVARIANT_GENES.out.scaling_factors)
        
        NORMALIZE_DESEQ2_QC_INVARIANT_GENES
            .out
            .scaling_factors
            .splitCsv ( header:true, sep:'\t' )
            .map { row -> 
                // Support both 'sample' and 'Sample_ID' column names
                def id = row.Sample_ID ?: row.sample
                // Remove .bam suffix and processing suffixes to match meta.id
                def clean_id = id.replaceAll(/\.mLb\.clN\.sorted\.bam$/, '').replaceAll(/\.bam$/, '')
                // Support multiple column name variations: size_factor, scaling_factor, scaling
                def value = row.size_factor ?: (row.scaling_factor ?: row.scaling)
                log.info "🔍 SCALING PARSED: sample='${clean_id}', value='${value}'"
                [ clean_id, value, 'invariant_genes' ]
            }
            .set { ch_size_factors_invariant }
        
        ch_size_factors = ch_size_factors.mix(ch_size_factors_invariant)
    }
    
    //
    // MODULE: All genes normalization (default DESeq2 method)
    //
    if (normalization_methods.contains('all_genes')) {
        NORMALIZE_DESEQ2_QC_ALL_GENES (
            SUBREAD_FEATURECOUNTS.out.counts.map { meta, counts -> counts },
            "featureCounts",
            ch_consensus_annotation
        )
        
        // Collect raw files for DESEQ2_TRANSFORM
        ch_deseq2_raw_files = ch_deseq2_raw_files.mix(NORMALIZE_DESEQ2_QC_ALL_GENES.out.read_dist_norm_txt)
        ch_deseq2_raw_files = ch_deseq2_raw_files.mix(NORMALIZE_DESEQ2_QC_ALL_GENES.out.sample_distances_txt)
        ch_deseq2_raw_files = ch_deseq2_raw_files.mix(NORMALIZE_DESEQ2_QC_ALL_GENES.out.pca_all_genes_txt)
        ch_deseq2_raw_files = ch_deseq2_raw_files.mix(NORMALIZE_DESEQ2_QC_ALL_GENES.out.pca_top_genes_txt)
        
        ch_normalization_versions = ch_normalization_versions.mix(NORMALIZE_DESEQ2_QC_ALL_GENES.out.versions)
        ch_normalization_scaling_factors = ch_normalization_scaling_factors.mix(NORMALIZE_DESEQ2_QC_ALL_GENES.out.scaling_factors)
        ch_scaling_factors_all = ch_scaling_factors_all.mix(NORMALIZE_DESEQ2_QC_ALL_GENES.out.scaling_factors)
        
        NORMALIZE_DESEQ2_QC_ALL_GENES
            .out
            .scaling_factors
            .splitCsv ( header:true, sep:'\t' )
            .map { row -> 
                // Support both 'sample' and 'Sample_ID' column names
                def id = row.Sample_ID ?: row.sample
                // Remove .bam suffix and processing suffixes to match meta.id
                def clean_id = id.replaceAll(/\.mLb\.clN\.sorted\.bam$/, '').replaceAll(/\.bam$/, '')
                // Support multiple column name variations: size_factor, scaling_factor, scaling
                def value = row.size_factor ?: (row.scaling_factor ?: row.scaling)
                log.info "🔍 SCALING PARSED (all_genes): sample='${clean_id}', value='${value}'"
                [ clean_id, value, 'all_genes' ]
            }
            .set { ch_size_factors_all_genes }
        
        ch_size_factors = ch_size_factors.mix(ch_size_factors_all_genes)
    }
    
    //
    // MODULE: Transform DESeq2 files for MultiQC with proper headers
    //
    DESEQ2_TRANSFORM (
        ch_deseq2_raw_files.flatten(),
        ch_deseq2_pca_header,
        ch_deseq2_clustering_header,
        ch_deseq2_read_dist_header
    )
    ch_versions = ch_versions.mix(DESEQ2_TRANSFORM.out.versions.first())
    
    // Populate MultiQC channels with transformed DESeq2 files
    // Separate files into 3 channels, then combine and sort at the end
    // Pass all DESeq2 files to MultiQC - they are already numbered (01_, 02_, etc.) by DESEQ2_TRANSFORM
    // MultiQC will sort them alphabetically by filename, preserving the numeric order
    ch_deseq2_all_multiqc = DESEQ2_TRANSFORM.out.multiqc_files.flatten()

    ch_versions = ch_versions.mix(ch_normalization_versions)

    // Assemble the channel 
    // Given a tab separated matrix with the first column : Sample_id, Scaling_factor convert the matrix to a channel with [Sample_id, Scaling_factor] pairs
    // Consider that the first line is the header - in principle Sample_id must match the meta.id from BAM_FILTER_SUBWF.out.bam    

    // Split size factors by normalization method
    ch_size_factors_invariant = ch_size_factors
        .filter { id, scaling, method -> method == 'invariant_genes' }
        .map { id, scaling, method -> 
            log.info "📊 SCALING FACTOR (invariant): id='${id}', scaling=${scaling}"
            [id, scaling] 
        }
    
    ch_size_factors_all_genes_only = ch_size_factors
        .filter { id, scaling, method -> method == 'all_genes' }
        .map { id, scaling, method -> 
            log.info "📊 SCALING FACTOR (all_genes): id='${id}', scaling=${scaling}"
            [id, scaling] 
        }
    
    // Create combined channel for invariant genes
    // Prepare BAM channel for deeptools
    ch_bam_for_deeptools = ch_genome_bam_bai

    // CHANNEL OPERATION: Combine BAM files with scaling factors for invariant genes
    // Use .combine() to create cartesian product, then filter by matching sample IDs
    // Strategy from rnaseq: use .map with null + .filter for cleaner separation
    ch_bam_bai_scale_invariant = ch_bam_for_deeptools
        .combine(ch_size_factors_invariant)
        .map { meta, bam, bai, sample_id, scaling -> 
            if (meta.id == sample_id) {
                def new_meta = meta.clone()
                new_meta.id = "${meta.id}.invariant_genes"
                new_meta.norm_method = 'invariant_genes'
                log.info "✅ MATCHED sample for invariant normalization: ${meta.id} (scaling=${scaling})"
                [new_meta, bam, bai, scaling]
            } else {
                null
            }
        }
        .filter { it != null }
    
    // CHANNEL OPERATION: Combine BAM files with scaling factors for all genes
    // Use same strategy as invariant genes for consistency
    ch_bam_bai_scale_all_genes = ch_bam_for_deeptools
        .combine(ch_size_factors_all_genes_only)
        .map { meta, bam, bai, sample_id, scaling -> 
            if (meta.id == sample_id) {
                def new_meta = meta.clone()
                new_meta.id = "${meta.id}.all_genes"
                new_meta.norm_method = 'all_genes'
                log.info "✅ MATCHED sample for all_genes normalization: ${meta.id} (scaling=${scaling})"
                [new_meta, bam, bai, scaling]
            } else {
                null
            }
        }
        .filter { it != null }

    ch_deeptoolsplotprofile_multiqc = Channel.empty()
    //
    // Scale to depth of sequencing using Deeptools:
    // 
    DEEPTOOLS_BIGWIG (
        ch_genome_bam_bai
    )
    ch_versions = ch_versions.mix(DEEPTOOLS_BIGWIG.out.versions.first())
    ch_big_wig = DEEPTOOLS_BIGWIG.out.bigwig

    if ( !params.skip_deeptools_norm ) {
        //
        // MODULE: DESeq2 normalized BigWig coverage tracks - invariant genes
        //
        if (normalization_methods.contains('invariant_genes')) {
            DEEPTOOLS_BIGWIG_NORM_INVARIANT (
                ch_bam_bai_scale_invariant
            )
            ch_versions = ch_versions.mix(DEEPTOOLS_BIGWIG_NORM_INVARIANT.out.versions.first())
            ch_big_wig = ch_big_wig.mix(DEEPTOOLS_BIGWIG_NORM_INVARIANT.out.bigwig)
        }
        
        //
        // MODULE: DESeq2 normalized BigWig coverage tracks - all genes
        //
        if (normalization_methods.contains('all_genes')) {
            DEEPTOOLS_BIGWIG_NORM_ALL_GENES (
                ch_bam_bai_scale_all_genes
            )
            ch_versions = ch_versions.mix(DEEPTOOLS_BIGWIG_NORM_ALL_GENES.out.versions.first())
            ch_big_wig = ch_big_wig.mix(DEEPTOOLS_BIGWIG_NORM_ALL_GENES.out.bigwig)
        }
    } 
    
    if (!params.skip_plot_profile ) {

        // Add an if so that if DEEPTOOLS_BIGWIG_NORM.out.bigwig is empty it will use DEEPTOOLS_BIGWIG.out.bigwig
        // MODULE: deepTools matrix generation for plotting
        //
            
        DEEPTOOLS_COMPUTEMATRIX (
            ch_big_wig,
            PREPARE_GENOME.out.gene_bed
        )

        ch_versions = ch_versions.mix(DEEPTOOLS_COMPUTEMATRIX.out.versions.first())

        //
        // MODULE: deepTools profile plots
        //
        DEEPTOOLS_PLOTPROFILE (
            DEEPTOOLS_COMPUTEMATRIX.out.matrix
        )
        ch_deeptoolsplotprofile_multiqc = DEEPTOOLS_PLOTPROFILE.out.table
        ch_versions = ch_versions.mix(DEEPTOOLS_PLOTPROFILE.out.versions.first())

        //
        // MODULE: deepTools heatmaps
        //
        DEEPTOOLS_PLOTHEATMAP (
            DEEPTOOLS_COMPUTEMATRIX.out.matrix
        )
        ch_versions = ch_versions.mix(DEEPTOOLS_PLOTHEATMAP.out.versions.first())
    } 

    //
    // MODULE: Pipeline reporting
    //
    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    //
    // MODULE: MultiQC
    //
    if (!params.skip_multiqc) {
        workflow_summary    = WorkflowAtacseq.paramsSummaryMultiqc(workflow, summary_params)
        ch_workflow_summary = Channel.value(workflow_summary)

        MULTIQC (
            ch_multiqc_config,
            ch_multiqc_custom_config.collect().ifEmpty([]),
            CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect(),
            ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'),

            ch_fastqc_raw_multiqc.collect{it[1]}.ifEmpty([]),
            ch_fastqc_trim_multiqc.collect{it[1]}.ifEmpty([]),
            ch_trim_log_multiqc.collect{it[1]}.ifEmpty([]),

            ch_bowtie2_multiqc.collect{it[1]}.ifEmpty([]),

            ch_samtools_stats.collect{it[1]}.ifEmpty([]),
            ch_samtools_flagstat.collect{it[1]}.ifEmpty([]),
            ch_samtools_idxstats.collect{it[1]}.ifEmpty([]),

            MARK_DUPLICATES_PICARD.out.stats.collect{it[1]}.ifEmpty([]),
            MARK_DUPLICATES_PICARD.out.flagstat.collect{it[1]}.ifEmpty([]),
            MARK_DUPLICATES_PICARD.out.idxstats.collect{it[1]}.ifEmpty([]),
            MARK_DUPLICATES_PICARD.out.metrics.collect{it[1]}.ifEmpty([]),

            BAM_FILTER_SUBWF.out.stats.collect{it[1]}.ifEmpty([]),
            BAM_FILTER_SUBWF.out.flagstat.collect{it[1]}.ifEmpty([]),
            BAM_FILTER_SUBWF.out.idxstats.collect{it[1]}.ifEmpty([]),
            ch_picardcollectmultiplemetrics_multiqc.collect{it[1]}.ifEmpty([]),
    
            ch_deeptoolsplotprofile_multiqc.collect{it[1]}.ifEmpty([]),
            ch_deeptoolsplotfingerprint_multiqc.collect{it[1]}.ifEmpty([]),
    
            PHANTOMPEAKQUALTOOLS.out.spp.collect{it[1]}.ifEmpty([]),
            MULTIQC_CUSTOM_PHANTOMPEAKQUALTOOLS.out.nsc.collect{it[1]}.ifEmpty([]),
            MULTIQC_CUSTOM_PHANTOMPEAKQUALTOOLS.out.rsc.collect{it[1]}.ifEmpty([]),
            MULTIQC_CUSTOM_PHANTOMPEAKQUALTOOLS.out.correlation.collect{it[1]}.ifEmpty([]),

            ch_custompeaks_frip_multiqc.collect{it[1]}.ifEmpty([]),
            ch_custompeaks_count_multiqc.collect{it[1]}.ifEmpty([]),
            ch_plothomerannotatepeaks_multiqc.collect().ifEmpty([]),
            ch_subreadfeaturecounts_multiqc.collect{it[1]}.ifEmpty([]),

            ch_deseq2_all_multiqc.collect().ifEmpty([])
        )
        multiqc_report = MULTIQC.out.report.toList()
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// NOTE: Completion handling is now managed by the PIPELINE_COMPLETION subworkflow
// in subworkflows/local/utils_nfcore_atacseq_pipeline/main.nf
// This includes email notifications, completion summary, and webhook notifications (imNotification)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
