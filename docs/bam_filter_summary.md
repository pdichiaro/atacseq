# BAM Filtering Pipeline Summary

## Complete Filtering Workflow (Bowtie2 → Analysis-Ready BAM)

This document provides a synthetic overview of the complete BAM filtering pipeline used in the ATAC-seq workflow, from raw alignment to the final analysis-ready BAM file.

---

## 📊 Complete Pipeline Flowchart (Paired-End)

```
┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
┃         COMPLETE FILTERING PIPELINE (PAIRED-END)                ┃
┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

📁 INPUT: sample_R1.fq.gz + sample_R2.fq.gz
         │
         ▼
┌───────────────────────────────────────────────────────────────────┐
│  STEP 1: BOWTIE2 ALIGNMENT                                        │
│  ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━     │
│                                                                   │
│  bowtie2 --very-sensitive -X 1000 -x genome                       │
│                                                                   │
│  Output: All concordant pairs (0-1000bp search space)             │
└───────────────────────────────────────────────────────────────────┘
         │
         ▼
    raw.sam (all alignments, 0-1000bp)
         │
         ▼
┌───────────────────────────────────────────────────────────────────┐
│  STEP 2: PICARD MARK DUPLICATES                                   │
│  ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━     │
│                                                                   │
│  picard MarkDuplicates \                                          │
│    REMOVE_DUPLICATES=false \  ← Only MARKS, doesn't remove        │
│    ASSUME_SORT_ORDER=coordinate                                   │
│                                                                   │
│  Sets FLAG 0x0400 for duplicate reads                             │
│                                                                   │
│  Output: sample.mkD.bam (duplicates MARKED)                       │
└───────────────────────────────────────────────────────────────────┘
         │
         ▼
    sample.mkD.bam
    ├─ All aligned reads
    ├─ FLAG 0x0400 on duplicates
    └─ Metrics: sample.mkD.MarkDuplicates.metrics.txt
         │
         ▼
┌───────────────────────────────────────────────────────────────────┐
│  STEP 3: BAM_FILTER (2-PASS FILTERING)                            │
│  ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━     │
│                                                                   │
│  ┌─────────────────────────────────────────────────────────────┐  │
│  │ PASS 1: Apply Most Filters                                  │  │
│  │ ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━    │  │
│  │                                                              │ │
│  │ samtools view -b -h \                                        │ │
│  │   -F 0x0100          ← Remove secondary alignments           │ │
│  │   -F 0x0800          ← Remove supplementary alignments       │ │
│  │   -F 0x0004          ← Remove unmapped reads                 │ │
│  │   -F 0x0008          ← Remove reads with unmapped mate       │ │
│  │   -f 0x0001          ← Keep only paired reads                │ │
│  │   -f 0x0002          ← Keep only proper pairs                │ │
│  │   -F 0x0400          ← Remove duplicates (if !keep_dups)     │ │
│  │   -L include_regions.bed  ← Keep only NON-blacklist reads    │ │
│  │   sample.mkD.bam > sample.filter1.bam                        │ │
│  │                                                              │ │
│  │ Where include_regions.bed is created by:                     │ │
│  │   bedtools complement -i blacklist.bed -g genome.sizes       │ │
│  │                                                              │ │
│  │ Output: sample.filter1.bam                                   │ │
│  └─────────────────────────────────────────────────────────────┘  │
│                                                                   │
│         │                                                         │
│         ▼                                                         │
│    sample.filter1.bam                                             │
│    ├─ Proper pairs only                                           │
│    ├─ No secondary/supplementary alignments                       │
│    ├─ No unmapped reads/mates                                     │
│    ├─ No blacklist regions                                        │
│    └─ No duplicates (if keep_dups=false)                          │
│         │                                                         │
│         ▼                                                         │
│  ┌─────────────────────────────────────────────────────────────┐  │
│  │ PASS 2: MAPQ + Fragment Size Filtering                      │  │
│  │ ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━     │ │
│  │                                                              │ │
│  │ samtools view -q 1 -h sample.filter1.bam | \                 │ │
│  │   awk -v max="500" \                                         │ │
│  │     '{if($0~/^@/ || (($9>=0?$9:-$9)<=max)) print}' | \       │ │
│  │   samtools view -b > sample.filter2.bam                      │ │
│  │                                                              │ │
│  │ Filters applied:                                             │ │
│  │   • -q 1: Remove MAPQ < 1 (multi-mappers)                    │ │
│  │   • awk: Keep fragments ≤ 500bp (params.insert_size)         │ │
│  │                                                              │ │
│  │ Output: sample.filter2.bam (FINAL)                           │ │
│  └─────────────────────────────────────────────────────────────┘  │
└───────────────────────────────────────────────────────────────────┘
         │
         ▼
    sample.filter2.bam ✅ FINAL OUTPUT
    ├─ Proper pairs only
    ├─ Fragments 0-500bp
    ├─ MAPQ ≥ 1
    ├─ No secondary/supplementary alignments
    ├─ No duplicates (if params.keep_dups=false)
    └─ No blacklist regions
         │
         ▼
┌───────────────────────────────────────────────────────────────────┐
│  STEP 4: BAM_SORT_SAMTOOLS                                        │
│  ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━     │
│                                                                   │
│  samtools sort sample.filter2.bam > sample.filter2.sorted.bam     │
│  samtools index sample.filter2.sorted.bam                         │
│                                                                   │       
└───────────────────────────────────────────────────────────────────┘
```

---

## 🔑 Key SAMtools Commands

### BAM_FILTER Implementation: 2-Pass Filtering

The BAM_FILTER module uses a **2-pass approach** to efficiently filter BAM files:

```bash
# ────────────────────────────────────────────────────────────────────
# PASS 1: Apply most filters (creates filter1.bam)
# ────────────────────────────────────────────────────────────────────

# First, create include regions (inverse of blacklist)
bedtools complement \
    -i blacklist.bed \
    -g genome.sizes > include_regions.bed

# Then apply all flags-based and blacklist filters
samtools view -b -h \
    -F 0x0100          # Remove secondary alignments (flag 256)
    -F 0x0800          # Remove supplementary alignments (flag 2048)
    -F 0x0004          # Remove unmapped reads (flag 4)
    -F 0x0008          # Remove reads with unmapped mate (flag 8)
    -f 0x0001          # Keep only paired reads (flag 1)
    -f 0x0002          # Keep only proper pairs (flag 2)
    -F 0x0400          # Remove duplicates (flag 1024) [if keep_dups=false]
    -L include_regions.bed  # Keep only reads NOT in blacklist regions
    sample.mkD.bam > sample.filter1.bam

# ────────────────────────────────────────────────────────────────────
# PASS 2: Apply MAPQ + fragment size filters (creates filter2.bam)
# ────────────────────────────────────────────────────────────────────

samtools view -q 1 -h sample.filter1.bam | \
    awk -v max="500" '{if($0~/^@/ || (($9>=0?$9:-$9)<=max)) print}' | \
    samtools view -b > sample.filter2.bam

# Where:
# -q 1: Keep only reads with MAPQ >= 1 (removes multi-mappers)
# awk: Keep only fragments with |TLEN| <= 500bp (params.insert_size)
```

**Why 2 passes instead of 1?**
- **Pass 1** applies all flag-based filters + blacklist filtering in one `samtools view` command
- **Pass 2** adds MAPQ filtering (`-q 1`) + fragment size filtering (via `awk`)
- This separation allows for efficient processing while maintaining clarity

**Why use include regions instead of excluding blacklist?**
- SAMtools `-L` (include) is faster and more reliable than `-U` (exclude)
- `bedtools complement` inverts the blacklist → creates "allowed" regions
- Result: Only reads in allowed regions are kept

---

## 📋 SAM Flags Reference

| Flag | Hex | Meaning | Action in Pipeline |
|------|-----|---------|-------------------|
| **KEEP flags (-f)** |
| 0x0001 | 1 | Read paired | ✅ KEEP paired reads |
| 0x0002 | 2 | Proper pair | ✅ KEEP proper pairs |
| **REMOVE flags (-F)** |
| 0x0004 | 4 | Unmapped | ❌ REMOVE unmapped reads |
| 0x0008 | 8 | Mate unmapped | ❌ REMOVE reads with unmapped mate |
| 0x0100 | 256 | Secondary alignment | ❌ REMOVE secondary alignments |
| 0x0400 | 1024 | PCR/optical duplicate | ❌ REMOVE if keep_dups=false |
| 0x0800 | 2048 | Supplementary alignment | ❌ REMOVE supplementary alignments |

---


## 🔍 Single-End vs Paired-End Differences

| Aspect | Paired-End | Single-End |
|--------|-----------|-----------|
| **Fragment size filter** | ✅ Yes (TLEN ≤ 500bp) | ❌ No (no TLEN field) |
| **Mate filters** | ✅ Yes (-F 0x0008, -f 0x0002) | ❌ No mate information |
| **Blacklist filter** | ✅ Yes | ✅ Yes |
| **Duplicate removal** | ✅ Yes | ✅ Yes |
| **MAPQ filter** | ✅ Yes | ✅ Yes |

**Single-End filtering (simplified):**

```bash
# SE: No fragment size or mate filtering
samtools view -b -h \
    -F 0x0004              # Remove unmapped
    -F 0x0100              # Remove secondary
    -F 0x0800              # Remove supplementary
    -F 0x0400              # Remove duplicates (if keep_dups=false)
    -q 1                   # Remove MAPQ < 1
    -L include_regions.bed # Keep only NON-blacklist
    sample.mkD.bam > sample.filter2.bam
```

---

## 📈 Expected Filtering Rates

For typical high-quality ATAC-seq data:

| Filter Category | Expected % | Concern if > |
|----------------|-----------|-------------|
| **Blacklist** | 5-15% | 20% |
| **Duplicates** | 20-40% | 60% |
| **Other filters** | 5-15% | 25% |
| **Total removed** | 30-60% | 75% |
| **Retained** | 40-70% | < 25% |

**Interpreting high removal rates:**

| High Category | Likely Cause | Recommended Action |
|---------------|-------------|-------------------|
| **Blacklist > 20%** | Non-specific antibody, wrong blacklist file | Check antibody specificity; verify genome |
| **Duplicates > 60%** | Over-amplification, low input material | Reduce PCR cycles, increase starting material |
| **Other > 25%** | Quality issues, multi-mappers | Check FastQC reports, verify alignment settings |
| **Retained < 25%** | Multiple issues | Review entire library prep protocol |

---

