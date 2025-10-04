# BAM Shift Reads Integration Summary

## Overview
Successfully integrated the BAM_SHIFT_READS subworkflow into the pdichiaro/atacseq pipeline to apply Tn5 transposase binding site offset corrections (+4 bp forward strand, -5 bp reverse strand) using deepTools alignmentSieve.

## Files Modified

### 1. Configuration Files

#### nextflow.config (Lines 50-53)
Added new parameters after `keep_mito`:
```groovy
shift_reads             = false  // Apply Tn5 transposase binding site offset (+4/-5 bp)
minFragmentLength       = 0      // Min fragment length for shifted reads (0 = no minimum)
maxFragmentLength       = 120    // Max fragment length for shifted reads (120 bp default)
```

#### nextflow_schema.json (Lines 234-254)
Added parameter definitions in the "alignment_options" section:
- `shift_reads`: Boolean flag to enable/disable read shifting (default: false)
- `minFragmentLength`: Integer for minimum fragment length filtering (default: 0)
- `maxFragmentLength`: Integer for maximum fragment length filtering (default: 120)

### 2. Subworkflow Created

#### subworkflows/local/bam_shift_reads.nf
New subworkflow containing:
- **DEEPTOOLS_ALIGNMENTSIEVE**: Applies standard ATAC-seq Tn5 offset (+4/-5 bp)
- **SAMTOOLS_SORT**: Sorts shifted reads by coordinate
- **SAMTOOLS_INDEX**: Creates BAM index (BAI/CSI)
- **SAMTOOLS_FLAGSTAT**: Generates alignment statistics

Outputs:
- `bam`: Shifted and sorted BAM files
- `bai`: BAM index files
- `csi`: CSI index files (if applicable)
- `flagstat`: Alignment statistics
- `versions`: Software versions

### 3. Main Workflow Integration

#### workflows/atacseq.nf

**Line 87**: Added include statement
```groovy
include { BAM_SHIFT_READS     } from '../subworkflows/local/bam_shift_reads'
```

**Lines 303-322**: Added conditional logic after BAM_FILTER_SUBWF
```groovy
if (params.shift_reads) {
    BAM_SHIFT_READS (
        BAM_FILTER_SUBWF.out.bam.join(BAM_FILTER_SUBWF.out.bai, by: [0]),
        params.minFragmentLength,
        params.maxFragmentLength
    )
    ch_final_bam = BAM_SHIFT_READS.out.bam
    ch_final_bai = BAM_SHIFT_READS.out.bai
    ch_final_csi = BAM_SHIFT_READS.out.csi
    ch_final_flagstat = BAM_SHIFT_READS.out.flagstat
    ch_versions = ch_versions.mix(BAM_SHIFT_READS.out.versions.first().ifEmpty(null))
} else {
    ch_final_bam = BAM_FILTER_SUBWF.out.bam
    ch_final_bai = BAM_FILTER_SUBWF.out.bai
    ch_final_csi = BAM_FILTER_SUBWF.out.csi
    ch_final_flagstat = BAM_FILTER_SUBWF.out.flagstat
}
```

**Downstream Channel Updates**: Replaced all `BAM_FILTER_SUBWF.out.bam` and `BAM_FILTER_SUBWF.out.bai` references with `ch_final_bam` and `ch_final_bai` in:
- BLACKLIST_LOG input (line 328-331)
- PICARD_COLLECTMULTIPLEMETRICS input (line 350)
- PHANTOMPEAKQUALTOOLS input (line 362)
- ch_genome_bam_bai channel creation (line 378-380)

## Usage

### Enable Read Shifting
```bash
nextflow run pdichiaro/atacseq \
    --input samplesheet.csv \
    --genome GRCh38 \
    --shift_reads \
    --minFragmentLength 0 \
    --maxFragmentLength 120 \
    --outdir results
```

### Disable Read Shifting (Default)
```bash
nextflow run pdichiaro/atacseq \
    --input samplesheet.csv \
    --genome GRCh38 \
    --outdir results
```

## Technical Details

### Tn5 Offset Correction
The Tn5 transposase creates a 9 bp duplication at the insertion site:
- **Forward strand reads**: Shift +4 bp (to the actual binding site)
- **Reverse strand reads**: Shift -5 bp (to the actual binding site)

This correction is standard practice in ATAC-seq analysis to accurately represent the true Tn5 binding locations.

### Fragment Length Filtering
The `minFragmentLength` and `maxFragmentLength` parameters control which reads are retained:
- **Default**: 0-120 bp (nucleosome-free fragments)
- **Nucleosome-free**: 0-100 bp
- **Mononucleosome**: 180-247 bp
- **All fragments**: Set min=0, max=1000 (or higher)

### deepTools alignmentSieve
Uses the following command structure:
```bash
alignmentSieve \
    --bam input.bam \
    --outFile shifted.bam \
    --ATACshift \
    --minFragmentLength ${min} \
    --maxFragmentLength ${max}
```

## Validation Steps

1. ✅ Configuration parameters added to nextflow.config
2. ✅ Schema parameters added to nextflow_schema.json
3. ✅ Subworkflow created with proper module dependencies
4. ✅ Include statement added to main workflow
5. ✅ Conditional logic implemented (if params.shift_reads)
6. ✅ Channel routing updated for all downstream processes
7. ⚠️  Nextflow lint skipped (SSL certificate issue in sandbox)
8. ⏳ Integration testing pending (requires sample data)

## Next Steps

### Testing
1. Run pipeline with `--shift_reads false` (default behavior)
2. Run pipeline with `--shift_reads true` on test dataset
3. Compare BAM file statistics between shifted and non-shifted reads
4. Verify downstream analyses (peak calling, QC metrics)

### Documentation
1. Update pipeline README with new parameters
2. Add usage examples to documentation
3. Include technical notes on Tn5 offset correction

## Dependencies
- deepTools ≥3.5.0 (for alignmentSieve)
- samtools ≥1.10
- Existing nf-core modules (samtools/sort, samtools/index, samtools/flagstat)

## Author
Seqera AI - Implementation completed on 2025
