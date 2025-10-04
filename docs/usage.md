# pdichiaro/atacseq: Usage

## Pipeline Overview

This ATAC-seq pipeline uses **Bowtie2** for alignment and **MACS2** for peak calling:

- **Bowtie2 Alignment**: 
  - Local alignment mode (`--local --very-sensitive-local`) for adapter soft-clipping
  - Two-stage fragment filtering for paired-end data (see [BOWTIE2_AND_BAM_FILTERING.md](BOWTIE2_AND_BAM_FILTERING.md))
  - Configurable fragment size filtering with `--insert_size` (default: 500bp)
  
- **MACS2 Peak Calling**:
  - Default: `--nomodel` mode (no model building, fixed extension)
  - Optional: Model building enabled with `--macs_gsize` parameter
  - For PE data: Uses actual fragment sizes from BAM files

## Quick Start

### Standard ATAC-seq Analysis with Input Controls
```bash
nextflow run pdichiaro/atacseq \
    --input samplesheet.csv \
    --outdir results \
    --genome GRCh38 \
    --with_inputs true \
    -profile singularity
```

### ATAC-seq Analysis without Input Controls
```bash
nextflow run pdichiaro/atacseq \
    --input samplesheet.csv \
    --outdir results \
    --genome GRCh38 \
    --with_inputs false \
    -profile singularity
```

| Parameter | Status | Description | Default |
|-----------|--------|-------------|---------|
| `--input` | **MANDATORY** | Path to samplesheet CSV file | `null` |
| `--outdir` | **MANDATORY** | Output directory path | `null` |
| `-profile` | **MANDATORY** | Execution environment (docker/singularity/conda) | None |


## Samplesheet input

You will need to create a samplesheet with information about the samples you would like to analyse before running the pipeline. Use this parameter to specify its location. It has to be a comma-separated file with 7 columns, and a header row as shown in the examples below.

```bash
--input '[path to samplesheet file]'
```

### Required Columns

| Column | Required | Description |
|--------|----------|-------------|
| `sample` | Yes | Sample name (alphanumeric, underscores, dots, dashes only) |
| `fastq_1` | Yes | Path to Read 1 FASTQ file (`.fastq.gz` or `.fq.gz`) |
| `fastq_2` | For PE | Path to Read 2 FASTQ file (leave empty for single-end) |
| `replicate` | Yes | Replicate number (integer: 1, 2, 3, ...) |
| `antibody` | For ChIP | Antibody/target name (empty for input samples) |
| `control` | For ChIP | Name of the control/input sample to use |
| `control_replicate` | For ChIP | Replicate number of the control sample to use |

### With Input Controls (Paired-End)

This is the standard format for ATAC-seq experiments with input controls:

```csv title="samplesheet.csv"
sample,fastq_1,fastq_2,replicate,antibody,control,control_replicate
WT_BCATENIN_IP,IP_rep1_R1.fastq.gz,IP_rep1_R2.fastq.gz,1,BCATENIN,WT_INPUT,1
WT_BCATENIN_IP,IP_rep2_R1.fastq.gz,IP_rep2_R2.fastq.gz,2,BCATENIN,WT_INPUT,2
WT_INPUT,input_rep1_R1.fastq.gz,input_rep1_R2.fastq.gz,1,,,
WT_INPUT,input_rep2_R1.fastq.gz,input_rep2_R2.fastq.gz,2,,,
```

### With Input Controls (Single-End)

For single-end sequencing, leave the `fastq_2` column empty:

```csv title="samplesheet.csv"
sample,fastq_1,fastq_2,replicate,antibody,control,control_replicate
WT_BCATENIN_IP,IP_rep1_R1.fastq.gz,,1,BCATENIN,WT_INPUT,1
WT_BCATENIN_IP,IP_rep2_R1.fastq.gz,,2,BCATENIN,WT_INPUT,2
WT_INPUT,input_rep1_R1.fastq.gz,,1,,,
WT_INPUT,input_rep2_R1.fastq.gz,,2,,,
```

### Without Input Controls

If you don't have input controls, leave `antibody`, `control`, and `control_replicate` empty:

```csv title="samplesheet.csv"
sample,fastq_1,fastq_2,replicate,antibody,control,control_replicate
SAMPLE1_H3K27ac,sample1_rep1_R1.fastq.gz,sample1_rep1_R2.fastq.gz,1,,,
SAMPLE1_H3K27ac,sample1_rep2_R1.fastq.gz,sample1_rep2_R2.fastq.gz,2,,,
SAMPLE2_H3K4me3,sample2_rep1_R1.fastq.gz,sample2_rep1_R2.fastq.gz,1,,,
```

> **Note:** Set `--with_inputs false` when running without input controls.

### Multiple Technical Replicates

The pipeline automatically merges multiple technical replicates (same sample, different sequencing runs). Technical replicates are identified by the `_T` suffix added automatically:

```csv title="samplesheet.csv"
sample,fastq_1,fastq_2,replicate,antibody,control,control_replicate
SAMPLE1_BCATENIN,run1_R1.fastq.gz,run1_R2.fastq.gz,1,BCATENIN,INPUT,1
SAMPLE1_BCATENIN,run2_R1.fastq.gz,run2_R2.fastq.gz,1,BCATENIN,INPUT,1
SAMPLE1_BCATENIN,run3_R1.fastq.gz,run3_R2.fastq.gz,2,BCATENIN,INPUT,2
INPUT,input_rep1_R1.fastq.gz,input_rep1_R2.fastq.gz,1,,,
INPUT,input_rep2_R1.fastq.gz,input_rep2_R2.fastq.gz,2,,,
```

In this example:
- Replicate 1 has 2 technical replicates (run1 and run2) → merged as `SAMPLE1_BCATENIN_REP1`
- Replicate 2 has 1 technical replicate (run3) → becomes `SAMPLE1_BCATENIN_REP2`

### Important Notes

- **Replicate IDs** must be consecutive integers starting from 1 (1, 2, 3, ...)
- **Control matching** is validated: the `control` + `control_replicate` combination must exist in the samplesheet
- **Sample names** cannot contain spaces (will be replaced with underscores)
- **File paths** can be absolute or relative to the launch directory
- **Mixed data types** are not allowed: all replicates of a sample must be either paired-end or single-end


## Complete Parameter Reference

### Mandatory Parameters

| Parameter | Type | Description | Example |
|-----------|------|-------------|---------|
| `--input` | path | Path to samplesheet CSV | `samplesheet.csv` |
| `--outdir` | path | Output directory | `results/` |

### Conditionally Mandatory Parameters

| Parameter | Condition | Type | Description |
|-----------|-----------|------|-------------|
| `--genome` | If no custom refs | string | iGenomes reference ID |
| `--fasta` | If no --genome | path | Genome FASTA file |
| `--gtf` | If no --genome | path | Gene annotation GTF |

### Key Optional Parameters

| Category | Parameter | Default | Description |
|----------|-----------|---------|-------------|
| **ATAC-seq** | `--with_inputs` | `true` | Use input control samples |
| | `--aligner` | `star` | Alignment method (star only) |
| | `--read_length` | `50` | Read length for MACS2 gsize |
| | `--fragment_size` | `200` | Estimated fragment size (SE) |
| **Peak Calling** | `--macs_gsize` | `null` | MACS2 genome size (auto-calculated) |
| | `--blacklist` | `null` | Regions to exclude from analysis |
| **Normalization** | `--skip_deeptools_norm` | `false` | Skip DESeq2 normalization |
| | `--normalization_method` | `invariant_genes` | Normalization method |
| **Quality** | `--skip_trimming` | `false` | Skip read trimming |
| | `--skip_fastqc` | `false` | Skip FastQC reports |
| | `--skip_qc` | `false` | Skip all QC steps |

### DESeq2 Normalization

By default, the pipeline generates **two types of BigWig coverage tracks**:

1. **Standard CPM normalization** (`DEEPTOOLS_BIGWIG`)
   - Always generated for all samples
   - Uses `--normalizeUsing CPM` (Counts Per Million)
   - Output: `*.extend.bw` and `*.extend.center.bw`

2. **DESeq2 size factor normalization** (`DEEPTOOLS_BIGWIG_NORM`) 
   - Generated by default (`--skip_deeptools_norm false`)
   - Uses DESeq2-calculated scaling factors
   - Better for differential binding analysis
   - To skip: `--skip_deeptools_norm true`

The `--normalization_method` parameter controls DESeq2 normalization:
- `invariant_genes` - Normalization using stable genes (default)
- `all_genes` - Standard DESeq2 normalization
- `all_genes,invariant_genes` - Run both methods

### Skip Options (All default to false)

- `--skip_fastqc` - Skip FastQC reports
- `--skip_trimming` - Skip read trimming with TrimGalore
- `--skip_picard_metrics` - Skip Picard QC metrics
- `--skip_plot_fingerprint` - Skip deepTools fingerprint plot
- `--skip_plot_profile` - Skip deepTools profile plots
- `--skip_peak_annotation` - Skip HOMER peak annotation
- `--skip_peak_qc` - Skip peak QC plots
- `--skip_multiqc` - Skip MultiQC report

### Reference Genome Options

| Parameter | Type | Description |
|-----------|------|-------------|
| `--fasta` | path | Reference genome FASTA |
| `--gtf` | path | Gene annotation GTF file |
| `--gff` | path | Gene annotation GFF file (alternative to GTF) |
| `--gene_bed` | path | Gene BED file (auto-generated if not provided) |
| `--star_index` | path | Pre-built STAR index |
| `--blacklist` | path | Blacklist regions BED file |
| `--save_reference` | boolean | Save generated indices (default: false) |

### Advanced Options

#### Trimming Options

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--clip_r1` | `null` | Remove bp from 5' end of R1 |
| `--clip_r2` | `null` | Remove bp from 5' end of R2 |
| `--three_prime_clip_r1` | `null` | Remove bp from 3' end of R1 |
| `--three_prime_clip_r2` | `null` | Remove bp from 3' end of R2 |
| `--trim_nextseq` | `null` | NextSeq/NovaSeq poly-G trimming |
| `--save_trimmed` | `false` | Save trimmed FastQ files |
| `--min_trimmed_reads` | `10000` | Min reads after trimming |
| `--extra_trimgalore_args` | `null` | Additional TrimGalore arguments |

#### Alignment Options (Bowtie2)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--aligner` | `bowtie2` | Alignment tool (fixed to bowtie2) |
| `--bowtie2_mode` | `--local` | Bowtie2 alignment mode (local vs end-to-end) |
| `--bowtie2_sensitivity` | `--very-sensitive-local` | Bowtie2 sensitivity preset |
| `--insert_size` | `500` | **Max fragment size for BAM filtering (PE only)** |
| `--keep_dups` | `false` | Keep duplicate reads |
| `--keep_multi_map` | `false` | Keep multimapping reads |
| `--keep_blacklist` | `false` | Keep blacklist regions (false = remove) |
| `--save_align_intermeds` | `false` | Save intermediate BAM files |
| `--save_unaligned` | `false` | Save unaligned reads |

**Note:** For Paired-End data, Bowtie2 uses a fixed `-X 1000` during alignment to search for fragments up to 1000bp. Post-alignment filtering with `--insert_size` (default 500bp) removes artifacts while keeping biologically relevant fragments. See [BOWTIE2_AND_BAM_FILTERING.md](BOWTIE2_AND_BAM_FILTERING.md) for details.

#### Peak Calling Options (MACS2)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--narrow_peak` | `false` | Call narrow peaks (vs broad) |
| `--broad_cutoff` | `0.1` | Broad peak FDR cutoff |
| `--macs_fdr` | `null` | MACS2 FDR threshold (q-value) |
| `--macs_pvalue` | `null` | MACS2 p-value threshold |
| `--macs_gsize` | `null` | **Effective genome size for MACS2 model building** |
| `--min_reps_consensus` | `1` | Min replicates for consensus peaks |
| `--save_macs_pileup` | `false` | Save MACS2 pileup tracks |

**MACS2 Fragment Length Strategy:**
- **Default behavior** (`--macs_gsize` not set):
  - Uses `--nomodel` (no model building)
  - Fixed extension size (200bp default)
  - Faster and more consistent
  - Recommended for most experiments
  
- **With `--macs_gsize` provided** (e.g., `--macs_gsize 2.7e9` for human):
  - MACS2 builds predictive model from data
  - Estimates fragment length automatically
  - Generates model plots (`*_model.r`)
  - Recommended for narrow peaks (transcription factors)
  
- **For Paired-End data**: Actual fragment sizes are used from BAM (TLEN field)

#### DESeq2 Normalization Options

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--sigma_times` | `1` | Sigma multiplier for invariant genes |
| `--n_pop` | `1` | Min samples for DESeq2 analysis |
| `--deseq2_vst` | `true` | Use VST transformation |
| `--skip_deseq2_qc` | `false` | Skip DESeq2 QC plots |

#### Quality Control Options

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--fingerprint_bins` | `500000` | Bins for deepTools fingerprint |
| `--skip_preseq` | `false` | Skip library complexity analysis |
| `--skip_spp` | `false` | Skip strand cross-correlation QC |

#### Other Advanced Options

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--seq_center` | `null` | Sequencing center for BAM RG |
| `--multiqc_title` | `null` | Custom MultiQC report title |
| `--email` | `null` | Email for completion summary |
| `--email_on_fail` | `null` | Email for failure notification |
| `--bamtools_filter_pe_config` | `assets/bamtools_filter_pe.json` | Paired-end BAM filtering config |
| `--bamtools_filter_se_config` | `assets/bamtools_filter_se.json` | Single-end BAM filtering config |

