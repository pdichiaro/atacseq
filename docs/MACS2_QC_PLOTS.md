# QC Plots Timeline - When MACS2 and Consensus Plots Are Generated

## 📊 QC Plot Execution Order

### 🔵 LEVEL 1: Individual Peak QC (MACS2_CALLPEAK_SINGLE)

**When they start:**
1. ✅ **MACS2_CALLPEAK_SINGLE** completes peak calling on each individual sample
2. ✅ Samples with 0 peaks are filtered out and a warning is emitted
3. ✅ **FRIP_SCORE** calculates FRiP (Fraction of Reads in Peaks) for each sample
4. ✅ **MULTIQC_CUSTOM_PEAKS** prepares data for MultiQC

**Then, IF `!params.skip_peak_annotation`:**
5. ✅ **HOMER_ANNOTATEPEAKS_MACS2** annotates individual peaks

**Then, IF `!params.skip_peak_qc`:**

#### 🎨 PLOT_MACS2_QC (QC plots for individual peaks)
```groovy
PLOT_MACS2_QC (
    ch_macs2_peaks.collect{it[1]}  // Collects ALL individual peak files
)
```
**Input**: All `.narrowPeak` or `.broadPeak` files from individual samples

**Output** (directory: `results/macs2/qc/`):
- `macs2_peak.counts_plot.pdf` - Bar plot of peak counts per sample
- `macs2_peak.widths_plot.pdf` - Distribution of peak widths
- Other statistical plots for individual peaks

#### 🎨 PLOT_HOMER_ANNOTATEPEAKS (Individual annotation plots)
```groovy
PLOT_HOMER_ANNOTATEPEAKS (
    HOMER_ANNOTATEPEAKS_MACS2.out.txt.collect{it[1]},
    ch_peak_annotation_header,
    "_peaks.annotatePeaks.txt"
)
```
**Input**: All Homer annotation files from individual samples

**Output** (directory: `results/macs2/qc/`):
- `peaks.annotatePeaks.summary.pdf` - Genomic distribution of peaks (promoter, intron, exon, etc.)
- `peaks.annotatePeaks.summary.tsv` - Data for MultiQC

---

### 🟢 LEVEL 2: Consensus Peaks QC BY CONDITION

**When they start:**
1. ✅ **MACS2_CONSENSUS_BY_CONDITION** completes consensus calling for each condition
   - E.g., `WT_BCATENIN` (from WT_BCATENIN_REP1 + WT_BCATENIN_REP2)
   - E.g., `NAIVE_BCATENIN` (from NAIVE_BCATENIN_REP1 + NAIVE_BCATENIN_REP2)

**Then, IF `!params.skip_peak_annotation`:**
2. ✅ **HOMER_ANNOTATEPEAKS_CONSENSUS_CONDITION** annotates consensus peaks by condition

**Then, IF `!params.skip_peak_qc`:**

#### 🎨 PLOT_MACS2_QC_CONSENSUS_CONDITION (QC plots for consensus by condition)
```groovy
PLOT_MACS2_QC_CONSENSUS_CONDITION (
    MACS2_CONSENSUS_BY_CONDITION.out.peaks.collect{it[1]}
)
```
**Input**: All consensus `.narrowPeak` or `.broadPeak` files for each condition

**Output** (directory: `results/macs2/consensus_peaks/{antibody}/by_condition/qc/`):
- `macs2_peak.condition.counts_plot.pdf` - Peak counts per condition
- `macs2_peak.condition.widths_plot.pdf` - Width distribution per condition
- Statistics of filtered consensus peaks (min_reps_consensus = 2)

#### 🎨 PLOT_HOMER_ANNOTATEPEAKS_CONSENSUS_CONDITION (Consensus annotation plots by condition)
```groovy
PLOT_HOMER_ANNOTATEPEAKS_CONSENSUS_CONDITION (
    HOMER_ANNOTATEPEAKS_CONSENSUS_CONDITION.out.txt.collect{it[1]},
    ch_peak_annotation_header,
    "_peaks.condition.annotatePeaks.txt"
)
```
**Input**: Homer annotations for each condition

**Output** (directory: `results/macs2/consensus_peaks/{antibody}/by_condition/qc/`):
- `peaks.condition.annotatePeaks.summary.pdf` - Genomic distribution per condition
- `peaks.condition.annotatePeaks.summary.tsv` - Data for MultiQC

---

### 🟣 LEVEL 3: Final Consensus BY ANTIBODY (MACS2_CONSENSUS)

**When it starts:**
1. ✅ **MACS2_CONSENSUS** completes final merge by antibody
   - E.g., `BCATENIN` (from WT_BCATENIN + NAIVE_BCATENIN)

**Then, IF `!params.skip_peak_annotation`:**
2. ✅ **HOMER_ANNOTATEPEAKS_CONSENSUS** annotates final consensus
3. ✅ **ANNOTATE_BOOLEAN_PEAKS** adds boolean columns

**⚠️ IMPORTANT NOTE**: 
There are no specific QC plots for the final by-antibody level in the current workflow.
Plots are only generated for:
- Individual peaks (level 1)
- Consensus by condition (level 2)

---

## 🕐 Complete Timeline

```
Start
  │
  ├─► MACS2_CALLPEAK_SINGLE (on each sample)
  │     ↓
  │   FRIP_SCORE
  │     ↓
  │   HOMER_ANNOTATEPEAKS_MACS2
  │     ↓
  │   🎨 PLOT_MACS2_QC (LEVEL 1)
  │     ↓
  │   🎨 PLOT_HOMER_ANNOTATEPEAKS (LEVEL 1)
  │
  ├─► MACS2_CONSENSUS_BY_CONDITION (per condition)
  │     ↓
  │   HOMER_ANNOTATEPEAKS_CONSENSUS_CONDITION
  │     ↓
  │   🎨 PLOT_MACS2_QC_CONSENSUS_CONDITION (LEVEL 2)
  │     ↓
  │   🎨 PLOT_HOMER_ANNOTATEPEAKS_CONSENSUS_CONDITION (LEVEL 2)
  │
  └─► MACS2_CONSENSUS (per antibody - final)
        ↓
      HOMER_ANNOTATEPEAKS_CONSENSUS
        ↓
      ANNOTATE_BOOLEAN_PEAKS
        ↓
      (NO QC PLOTS - final output ready)
```

---

## ⚙️ Parameters to Control Plots

### Disable all QC plots:
```bash
--skip_peak_qc
```
This skips **ALL** QC plots (level 1 and 2)

### Disable annotations (and therefore plots):
```bash
--skip_peak_annotation
```
This skips Homer annotations and related plots

---

## 📂 Output Directories

```
results/
└── macs2/
    ├── qc/                                    # QC LEVEL 1 (individual peaks)
    │   ├── macs2_peak.counts_plot.pdf
    │   ├── macs2_peak.widths_plot.pdf
    │   └── peaks.annotatePeaks.summary.pdf
    │
    └── consensus_peaks/
        └── {ANTIBODY}/                        # E.g., BCATENIN/
            ├── by_condition/                  # QC LEVEL 2 (consensus by condition)
            │   ├── qc/
            │   │   ├── macs2_peak.condition.counts_plot.pdf
            │   │   ├── macs2_peak.condition.widths_plot.pdf
            │   │   └── peaks.condition.annotatePeaks.summary.pdf
            │   ├── WT_BCATENIN_peaks.narrowPeak
            │   └── NAIVE_BCATENIN_peaks.narrowPeak
            │
            └── {ANTIBODY}.bed                 # LEVEL 3 (final - no QC plots)
                {ANTIBODY}.boolean.txt
                {ANTIBODY}.saf
```

---

## 🎯 What to Expect During Execution

1. **First plots to start**: PLOT_MACS2_QC and PLOT_HOMER_ANNOTATEPEAKS
   - Start immediately after MACS2_CALLPEAK_SINGLE
   - Use individual peak files

2. **Intermediate plots**: PLOT_MACS2_QC_CONSENSUS_CONDITION and PLOT_HOMER_ANNOTATEPEAKS_CONSENSUS_CONDITION
   - Start after consensus by condition
   - Show statistics on filtered peaks (min_reps = 2)

3. **Final output**: No automatic QC plots
   - The `{ANTIBODY}.boolean.txt` file contains all info for custom analysis
   - You can create custom plots using R/Python on final files

---

## 💡 Suggestion for Custom Level 3 Plots

If you want to generate plots for the final consensus by antibody, you can use the boolean.txt file:

```R
# R example for custom plots
library(ggplot2)

# Load boolean file
peaks <- read.table("BCATENIN.boolean.txt", header=TRUE, sep="\t")

# Plot: number of peaks per condition
conditions <- colnames(peaks)[5:ncol(peaks)]
peak_counts <- colSums(peaks[,5:ncol(peaks)])

ggplot(data.frame(Condition=conditions, Peaks=peak_counts), 
       aes(x=Condition, y=Peaks)) +
  geom_bar(stat="identity") +
  theme_minimal() +
  labs(title="BCATENIN Peaks by Condition")
```

---

## 📋 Summary Table

| Level | Stage | QC Plots Generated | Output Directory |
|-------|-------|-------------------|------------------|
| **1** | Individual samples | ✅ PLOT_MACS2_QC<br>✅ PLOT_HOMER_ANNOTATEPEAKS | `results/macs2/qc/` |
| **2** | Consensus by condition | ✅ PLOT_MACS2_QC_CONSENSUS_CONDITION<br>✅ PLOT_HOMER_ANNOTATEPEAKS_CONSENSUS_CONDITION | `results/macs2/consensus_peaks/{antibody}/by_condition/qc/` |
| **3** | Final consensus by antibody | ❌ No automatic plots | `results/macs2/consensus_peaks/{antibody}/` |

---

## 🔍 Key Differences Between Levels

### Level 1 (Individual Peaks)
- **Input**: Raw peak calls from each replicate
- **Purpose**: QC individual samples, identify outliers
- **Min peaks filter**: Samples with 0 peaks are excluded

### Level 2 (Consensus by Condition)
- **Input**: Merged peaks across replicates of same condition
- **Purpose**: Assess consistency between replicates
- **Filter**: `min_reps_consensus` parameter (default: 2)
- **Example**: WT_BCATENIN consensus requires peak present in ≥2 WT_BCATENIN replicates

### Level 3 (Final Consensus by Antibody)
- **Input**: Merged peaks across all conditions for same antibody
- **Purpose**: Create final peak set for downstream analysis
- **Output**: Boolean matrix showing which conditions support each peak
- **No automatic plots**: Use boolean.txt for custom visualization

---

## 🚀 Workflow Optimization Tips

1. **Quick QC check**: Run with `--skip_peak_qc` to skip plot generation during initial testing
2. **Full QC report**: Include all plots for final analysis and publication
3. **Custom analysis**: Use boolean.txt files from Level 3 for condition-specific comparisons
4. **MultiQC integration**: All QC plots contribute data to the final MultiQC report
