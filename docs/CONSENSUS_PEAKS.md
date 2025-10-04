# Consensus Peaks Generation

## Overview

The ATAC-seq pipeline implements a **two-level process** for consensus peak generation, intelligently combining peaks identified across biological replicates and multiple experimental conditions.

## Workflow

```
Individual peaks (MACS2)
         ↓
    [LEVEL 1]
         ↓
Consensus per CONDITION (e.g., WT_BCATENIN, NAIVE_BCATENIN)
         ↓
    [LEVEL 2]
         ↓
Final consensus per ANTIBODY (e.g., BCATENIN)
```

---

## LEVEL 1: Consensus per Condition

### Process: `MACS2_CONSENSUS_BY_CONDITION`

**Objective**: Create consensus peaks for each condition + antibody combination (e.g., `WT_BCATENIN`, `NAIVE_BCATENIN`).

### Input
- **Individual peaks** from biological replicates for each condition
- Format: narrowPeak or broadPeak (files `*_peaks.narrowPeak` or `*_peaks.broadPeak`)

### Algorithm

#### Step 1: Peak Merging
```bash
# Sort and merge overlapping genomic regions
sort -T '.' -k1,1 -k2,2n [peak_files] | mergeBed -c 2,3,4,5,6,7,8,9,10 -o collapse,collapse,...
```

- **Sorts** peaks by chromosome and position
- **mergeBed** merges overlapping peaks while preserving all original columns as collapsed data

#### Step 2: Expansion and Filtering
```bash
macs2_merged_expand.py \
    merged.txt \
    sample_names \
    output.boolean.txt \
    --min_replicates N
```

The `macs2_merged_expand.py` script:

1. **Analyzes merged regions** identifying which samples contribute to each region
2. **Groups by biological sample** (removing `_R1`, `_R2`, etc. suffixes)
3. **Applies replicate threshold**: retains only regions present in ≥ N replicates
   - `--min_replicates` defined by `params.min_reps_consensus` (default: 1)
4. **Generates boolean matrix** indicating which samples contribute to each peak

#### Step 3: Output Generation

**Files generated for each condition:**

| File | Description |
|------|-------------|
| `*.bed` | BED6 format (chr, start, end, name, score, strand) |
| `*_peaks.narrowPeak` / `*_peaks.broadPeak` | Complete MACS2 format for downstream analysis |
| `*.saf` | SAF format for featureCounts (quantification) |
| `*.condition.txt` | File → path mapping |
| `*.boolean.txt` | Sample presence/absence matrix + statistics |
| `*.intersect.txt` | Intersections between samples (for UpSet plots) |

**Publishing directory:**
```
{outdir}/{aligner}/mergedLibrary/macs2/{narrowPeak|broadPeak}/consensus/{antibody}/by_condition/
```

### Practical Example

**Input**: 3 replicates for condition `WT_BCATENIN`
```
WT_BCATENIN_REP1_peaks.narrowPeak
WT_BCATENIN_REP2_peaks.narrowPeak
WT_BCATENIN_REP3_peaks.narrowPeak
```

**Output**: Consensus peaks `WT_BCATENIN`
- Only peaks present in ≥ `min_reps_consensus` replicates (e.g., 2/3)
- Unified coordinates for overlapping regions
- Aggregated metadata (fold-change, p-value, q-value for each replicate)

---

## LEVEL 2: Final Consensus per Antibody

### Process: `MACS2_CONSENSUS`

**Objective**: Merge consensus peaks from all conditions for the same antibody into a final set.

### Input
- **Consensus peaks per condition** (output from Level 1)
- E.g., `WT_BCATENIN_peaks.narrowPeak` + `NAIVE_BCATENIN_peaks.narrowPeak`

### Key Differences from Level 1

#### Parameter `min_replicates`
```groovy
def min_reps = meta.replicates_exist ? 1 : params.min_reps_consensus
```

- **If `replicates_exist = true`** (Level 2): `min_reps = 1`
  - Each input is already a filtered consensus → accepts all conditions
- **If `replicates_exist = false`** (Level 1): uses `params.min_reps_consensus`

This prevents **double filtering**: since peaks have already been filtered by replicates at Level 1, Level 2 accepts any contributing condition.

### Algorithm

Same as Level 1, but with adapted filtering logic:

1. **Merge** peaks from all conditions
2. **Expansion** with `macs2_merged_expand.py` (min_replicates=1)
3. **Generate** final output for the antibody

### Output

**Files generated for each antibody:**

| File | Description |
|------|-------------|
| `{antibody}.bed` | Consensus peaks in BED6 format |
| `{antibody}.saf` | For quantification with featureCounts |
| `{antibody}.antibody.txt` | File mapping |
| `{antibody}.boolean.txt` | Presence/absence matrix for conditions |
| `{antibody}.intersect.txt` | Intersections between conditions |

**Publishing directory:**
```
{outdir}/{aligner}/mergedLibrary/macs2/{narrowPeak|broadPeak}/consensus/{antibody}/
```

### Practical Example

**Input**: Consensus from 2 conditions for antibody `BCATENIN`
```
WT_BCATENIN_peaks.narrowPeak      (from 3 replicates)
NAIVE_BCATENIN_peaks.narrowPeak   (from 3 replicates)
```

**Output**: Final consensus `BCATENIN`
- Merges peaks present in any condition (≥1)
- Unified coordinates for overlapping regions between conditions
- Metadata for each condition

---

## Python Script: `macs2_merged_expand.py`

### Main Features

#### 1. Merged Peaks Parsing
```python
# Input: file from mergeBed with collapsed columns
chromID, mstart, mend, starts, ends, names, fcs, pvals, qvals, summits
```

#### 2. Sample Grouping
```python
# Removes replicate suffixes (_R1, _R2) to identify biological samples
sample = "_".join(names[idx].split("_")[:-2])
```

#### 3. Replicate Threshold
```python
if len(sample_replicates) >= minReplicates:
    passRepThreshList.append(sample)
```

#### 4. Boolean Matrix
For each genomic region:
```python
boolList = ["TRUE" if sample in region else "FALSE" for sample in all_samples]
```

#### 5. Intersection File
```python
# Output for UpSet plot
"Sample1&Sample2&Sample3\t42"  # 42 common peaks
"Sample1&Sample2\t150"          # 150 peaks in these 2
```

---

## Output File Formats

### 1. `.boolean.txt` File

Table with columns:

| Column | Description |
|---------|-------------|
| `chr`, `start`, `end` | Genomic coordinates |
| `interval_id` | Unique ID (Interval_1, Interval_2, ...) |
| `num_peaks` | Total number of original peaks in the region |
| `num_samples` | Number of samples passing threshold |
| `{sample}.bool` | TRUE/FALSE for presence in each sample |
| `{sample}.fc` | Fold-change (multiple values separated by `;`) |
| `{sample}.qval` | Q-value (multiple values separated by `;`) |
| `{sample}.pval` | P-value (multiple values separated by `;`) |
| `{sample}.start` | Original start coordinates (multiple separated by `;`) |
| `{sample}.end` | Original end coordinates (multiple separated by `;`) |
| `{sample}.summit` | (narrowPeak only) Summit offset (multiple separated by `;`) |

### 2. `.saf` File (Simplified Annotation Format)

Format for featureCounts:
```
GeneID    Chr    Start    End    Strand
Interval_1    chr1    1000    1500    +
Interval_2    chr1    2000    2300    +
```

### 3. `.intersect.txt` File

For UpSet plots:
```
WT_BCATENIN&NAIVE_BCATENIN    523
WT_BCATENIN    127
NAIVE_BCATENIN    89
```

---

## Configurable Parameters

### `params.min_reps_consensus`

**Default**: `1` (but recommended: `2`)

**Description**: Minimum number of biological replicates required to consider a peak in the consensus.

**Example**:
```groovy
// nextflow.config
params.min_reps_consensus = 2  // Requires ≥2 replicates
```

**When applied**:
- ✅ **Level 1** (consensus per condition): applies threshold to replicates
- ❌ **Level 2** (consensus per antibody): ignores (uses min_reps=1)

### `params.narrow_peak`

**Default**: `true`

**Description**: Determines peak format and columns to process.

| Value | Format | Columns |
|--------|---------|---------|
| `true` | narrowPeak | 10 columns (includes summit) |
| `false` | broadPeak | 9 columns (no summit) |

---

## Quantification with featureCounts

The generated `.saf` files are used to quantify reads in each consensus peak:

```groovy
SUBREAD_FEATURECOUNTS (
    [ meta, bams, saf ]
)
```

**Output**: Read count matrix per peak × sample, used for:
- Normalization (DESeq2)
- Differential analysis
- QC (PCA, clustering)

---

