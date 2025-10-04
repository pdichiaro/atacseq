# Bowtie2 Alignment and BAM Filtering Strategy

## 📘 Overview

This document explains the complete alignment and filtering strategy used in this ATAC-seq pipeline for both **Single-End (SE)** and **Paired-End (PE)** sequencing data.

---

## 🗺️ Pipeline Workflow Overview

```
┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
┃                     PAIRED-END WORKFLOW                           ┃
┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

📁 INPUT: sample_R1.fq.gz + sample_R2.fq.gz
         │
         ▼
┌────────────────────────────────────────────────────────────────┐
│  STEP 1: BOWTIE2 ALIGNMENT                                     │
│  ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━   │
│                                                                │
│  bowtie2 --very-sensitive \                                    │
│    -X 1000           ← Max fragment search (FIXED)             │
│    -x genome -1 R1.fq.gz -2 R2.fq.gz                           │
│                                                                │
│  Output: Unfiltered SAM                                        │
│  ├─ Concordant pairs (0-1000bp)                                │
│  ├─ Both reads mapped                                          │
│  └─ Correct orientation                                        │
└────────────────────────────────────────────────────────────────┘
         │
         ▼
    raw.sam (all fragments 0-1000bp)
         │
         ▼
┌────────────────────────────────────────────────────────────────┐
│  STEP 2: CONVERT TO BAM + INITIAL FILTER                       │
│  ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━  │
│                                                                 │
│  samtools view -h -b \                                         │
│    -F 0x004   ← Remove unmapped reads                          │
│    -F 0x0008  ← Remove reads with unmapped mate                │
│    -f 0x001   ← Keep only paired reads                         │
│    -q 1       ← Remove MAPQ < 1                                │
│                                                                 │
└────────────────────────────────────────────────────────────────┘
         │
         ▼
    temp1.bam
         │
         ▼
┌────────────────────────────────────────────────────────────────┐
│  STEP 3: SORT BY NAME                                          │
│  ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━  │
│                                                                 │
│  samtools sort -n    ← Group read pairs together               │
│                                                                 │
└────────────────────────────────────────────────────────────────┘
         │
         ▼
    temp2.bam (name-sorted)
         │
         ▼
┌────────────────────────────────────────────────────────────────┐
│  STEP 4: FIX MATE INFORMATION                                  │
│  ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━  │
│                                                                 │
│  samtools fixmate -r  ← Fix mate info, remove secondary        │
│                                                                 │
└────────────────────────────────────────────────────────────────┘
         │
         ▼
    temp3.bam
         │
         ▼
┌────────────────────────────────────────────────────────────────┐
│  STEP 5: FRAGMENT SIZE FILTER (AWK)                            │
│  ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━  │
│                                                                 │
│  awk -v var=500 '                                              │
│    if (header line)           → KEEP                           │
│    if (|TLEN| <= 500)         → KEEP                          │
│    else                       → DISCARD                        │
│  '                                                             │
│                                                                 │
│  Actual AWK command:                                           │
│  awk -v var="$max_frag" '{                                     │
│    if(substr($0,1,1)=="@" || (($9>=0?$9:-$9)<=var))           │
│      print $0                                                  │
│  }'                                                            │
│                                                                 │
│  ⚙️  params.insert_size = 500 (DEFAULT, USER-CONFIGURABLE)     │
│                                                                 │
│  Removes:                                                      │
│  ❌ Fragments > 500bp (chimeras, artifacts)                    │
│  ❌ Unusually long inserts                                     │
│                                                                 │
└────────────────────────────────────────────────────────────────┘
         │
         ▼
    temp4.sam (fragments 0-500bp only)
         │
         ▼
┌────────────────────────────────────────────────────────────────┐
│  STEP 6: CONVERT BACK TO BAM                                   │
│  ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━  │
│                                                                 │
│  samtools view -h -b                                           │
│                                                                 │
└────────────────────────────────────────────────────────────────┘
         │
         ▼
    temp5.bam
         │
         ▼
┌────────────────────────────────────────────────────────────────┐
│  STEP 7: SORT BY COORDINATE                                    │
│  ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━  │
│                                                                 │
│  samtools sort  ← Position-sorted for downstream tools         │
│                                                                 │
└────────────────────────────────────────────────────────────────┘
         │
         ▼
    sample.filtered.bam
         │
         ▼
┌────────────────────────────────────────────────────────────────┐
│  STEP 8: INDEX BAM                                             │
│  ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━  │
│                                                                 │
│  samtools index                                                │
│                                                                 │
└────────────────────────────────────────────────────────────────┘
         │
         ▼
📁 OUTPUT: sample.filtered.bam + sample.filtered.bam.bai
         │
         ▼
   [Next: Picard MarkDuplicates → MACS2 Peak Calling]


┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
┃                      SINGLE-END WORKFLOW                          ┃
┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

📁 INPUT: sample.fq.gz
         │
         ▼
┌────────────────────────────────────────────────────────────────┐
│  STEP 1: BOWTIE2 ALIGNMENT                                     │
│  ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━  │
│                                                                 │
│  bowtie2 --very-sensitive \                                    │
│    -U sample.fq.gz -x genome                                   │
│                                                                 │
└────────────────────────────────────────────────────────────────┘
         │
         ▼
    raw.sam
         │
         ▼
┌────────────────────────────────────────────────────────────────┐
│  STEP 2: CONVERT TO BAM + FILTER                               │
│  ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━  │
│                                                                 │
│  samtools view -h -b \                                         │
│    -F 0x004   ← Remove unmapped reads                          │
│    -q 1       ← Remove MAPQ < 1                                │
│                                                                 │
│  ⚠️  NO FRAGMENT SIZE FILTER (no TLEN in SE)                   │
│                                                                 │
└────────────────────────────────────────────────────────────────┘
         │
         ▼
    temp.bam
         │
         ▼
┌────────────────────────────────────────────────────────────────┐
│  STEP 3: SORT BY COORDINATE                                    │
│  ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━  │
│                                                                 │
│  samtools sort                                                 │
│                                                                 │
└────────────────────────────────────────────────────────────────┘
         │
         ▼
    sample.filtered.bam
         │
         ▼
┌────────────────────────────────────────────────────────────────┐
│  STEP 4: INDEX BAM                                             │
│  ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━  │
│                                                                 │
│  samtools index                                                │
│                                                                 │
└────────────────────────────────────────────────────────────────┘
         │
         ▼
📁 OUTPUT: sample.filtered.bam + sample.filtered.bam.bai
         │
         ▼
   [Next: Picard MarkDuplicates → MACS2 Peak Calling]
```

---

## 🎯 Key Differences: PE vs SE

| Aspect | Paired-End (PE) | Single-End (SE) |
|--------|----------------|----------------|
| **Bowtie2 -X** | Fixed at 1000bp | N/A (ignored) |
| **Fragment info** | Yes (TLEN field) | No (TLEN = 0) |
| **Fragment filter** | ✅ Yes (AWK by TLEN) | ❌ No (not applicable) |
| **Mate filtering** | ✅ Yes (fixmate step) | ❌ No mate info |
| **Complexity** | 8 steps | 4 steps |
| **Filter control** | `--insert_size` param | None (only MAPQ) |

---

## 🔧 Two-Stage Filtering (PE Only)

```
┌─────────────────────────────────────────────────────────────────┐
│                   STAGE 1: BOWTIE2 ALIGNMENT                    │
│                   ═══════════════════════════                   │
│                                                                 │
│  Parameter: -X 1000 (FIXED - never changes)                    │
│  Purpose:   Permissive search for concordant pairs             │
│  Result:    All valid fragments 0-1000bp aligned               │
│                                                                 │
│  ✅ Ensures no valid pairs are missed during alignment          │
│  ✅ Covers biological range + potential artifacts               │
└─────────────────────────────────────────────────────────────────┘
                              │
                              ▼
                    Unfiltered BAM (0-1000bp)
                              │
                              ▼
┌─────────────────────────────────────────────────────────────────┐
│                    STAGE 2: BAM_FILTER (AWK)                    │
│                    ════════════════════════                     │
│                                                                 │
│  Parameter: params.insert_size = 500 (DEFAULT, configurable)   │
│  Purpose:   Biological quality control filtering               │
│  Result:    Only high-quality fragments 0-500bp retained       │
│                                                                 │
│  ✅ Removes chimeras and sequencing artifacts (500-1000bp)      │
│  ✅ User can adjust based on experiment type                    │
│  ✅ No re-alignment needed to change filtering                  │
└─────────────────────────────────────────────────────────────────┘
                              │
                              ▼
                    Filtered BAM (0-500bp)
                              │
                              ▼
                     Ready for peak calling
```

---

## 🧬 Part 1: Bowtie2 Alignment

### What is Bowtie2?

**Bowtie2** is a fast and memory-efficient alignment tool that maps sequencing reads to a reference genome. It uses the **Burrows-Wheeler Transform (BWT)** algorithm for efficient alignment.

### Key Bowtie2 Parameters Used in This Pipeline

#### Common Parameters (SE and PE)

```bash
--local                    # Local alignment mode (soft-clips ends)
--very-sensitive-local     # Preset for high sensitivity
--no-mixed                 # Suppress unpaired alignments for paired reads (PE only)
--no-discordant           # Suppress discordant alignments (PE only)
--phred33                  # Input quality scores are Phred+33
--minins 0                 # Minimum fragment length (PE only)
-I 0                       # Minimum fragment length (PE only)
-X 1000                    # Maximum fragment length (PE only) - FIXED VALUE
--no-unal                  # Suppress unaligned reads in SAM output
```

#### Alignment Mode: --local vs --end-to-end

| Mode | Description | When to Use |
|------|-------------|-------------|
| **--local** | Soft-clips read ends that don't match well | ATAC-seq (used in this pipeline) |
| **--end-to-end** | Requires entire read to align | RNA-seq, strict matching scenarios |

**Why --local for ATAC-seq?**
- ATAC-seq reads may contain adapter sequences or low-quality ends
- Soft-clipping improves alignment rate without sacrificing accuracy
- Focuses on the high-quality core of each read

#### Sensitivity Preset: --very-sensitive-local

This preset configures multiple alignment parameters for high sensitivity:

```bash
# Equivalent to:
-D 20      # Number of consecutive seed extension attempts
-R 3       # Number of re-seeding rounds
-N 0       # Number of mismatches allowed in seed (0 = exact match)
-L 20      # Seed length
-i S,1,0.50 # Seed interval function
```

**Trade-off:**
- ✅ Higher sensitivity → detects more alignments
- ⚠️ Slower runtime, more memory usage

---

## 🔍 Part 2: Paired-End (PE) Specific Behavior

### Fragment Length Filtering (-X parameter)

```bash
-X 1000  # Maximum insert size (FIXED in this pipeline)
```

**What does -X control?**
- Maximum distance between read pairs to be considered concordant
- Only affects **alignment search** during Bowtie2 execution
- Does NOT filter the final BAM file

### Understanding Insert Size vs Fragment Length

```
5' ═══R1═══>           <═══R2═══ 3'
   |←─────────────────────────→|
        Fragment Length
   
   |←────────────→|
    Insert Size (sequenced region)
```

**In practice:**
- **Fragment length** = Full DNA fragment from library prep
- **Insert size** = Distance between read pairs (TLEN in BAM)
- For 150bp PE reads with 500bp fragments: insert ≈ 200bp

### Concordant vs Discordant Alignments

#### Concordant Pairs (KEPT)
```
─────────→ R1         R2 ←─────────
    |←───── ≤ 1000bp ─────→|
```
- Correctly oriented (R1 forward, R2 reverse)
- Distance ≤ 1000bp

#### Discordant Pairs (REMOVED with --no-discordant)
```
─────────→ R1   R2 ─────────→  (same orientation)
                                
R1 ←─────────   ─────────→ R2  (reversed)

R1 ───→                   ←─── R2  (too far apart: > 1000bp)
```

### Mixed Alignments (REMOVED with --no-mixed)

**Mixed alignment** = One read aligns, its mate doesn't

```
R1 ───────→ (aligned to chr1)
R2          (unaligned)
```

**Why suppress mixed alignments?**
- In ATAC-seq PE, we want both reads to support the fragment
- Single-read alignment is less reliable for fragment inference

---

## 🔬 Part 3: Single-End (SE) Specific Behavior

### No Fragment Inference

For SE data:
- Only one read per fragment
- No insert size information
- No concordance concept
- **-X, --no-mixed, --no-discordant are IGNORED**

### SE Alignment Strategy

```bash
bowtie2 \
    --local \
    --very-sensitive-local \
    --phred33 \
    --no-unal \
    -U reads.fastq.gz \
    -x genome_index \
    -S output.sam
```

**Key differences from PE:**
- `-U` (unpaired) instead of `-1/-2` (paired)
- No fragment length constraints
- Alignment quality depends on single read quality

---

## 🧹 Part 4: Post-Alignment BAM Filtering

After Bowtie2 alignment, the pipeline applies additional filtering in the **BAM_FILTER** module.

### Paired-End (PE) Filtering

#### Stage 1: SAMtools Initial Filter
```bash
samtools view -h -b \
    -F 0x004 \     # Remove unmapped reads
    -F 0x0008 \    # Remove reads with unmapped mate
    -f 0x001 \     # Require read paired
    -q 1 \         # Remove MAPQ < 1
```

**SAM Flags Explained:**

| Flag | Hex | Meaning |
|------|-----|---------|
| 0x001 | 1 | Read paired in sequencing |
| 0x002 | 2 | Read mapped in proper pair |
| 0x004 | 4 | Read unmapped |
| 0x008 | 8 | Mate unmapped |
| 0x010 | 16 | Read reverse strand |
| 0x020 | 32 | Mate reverse strand |
| 0x040 | 64 | First in pair |
| 0x080 | 128 | Second in pair |
| 0x100 | 256 | Secondary alignment |
| 0x200 | 512 | Failed quality check |
| 0x400 | 1024 | PCR or optical duplicate |
| 0x800 | 2048 | Supplementary alignment |

#### Stage 2: Coordinate Sorting
```bash
samtools sort -n  # Sort by read name (pairs together)
```

#### Stage 3: Fix Mate Information
```bash
samtools fixmate -r  # -r removes secondary/unmapped reads
```

#### Stage 4: Fragment Size Filtering (AWK)

This is where **params.insert_size** is applied:

```bash
# Filter fragments by insert size
# Actual command used in the pipeline:
awk -v var="$max_frag" '{
    if(substr($0,1,1)=="@" || (($9>=0?$9:-$9)<=var)) 
        print $0
}'
```

**How this AWK filter works:**

1. **`-v var="$max_frag"`**: Sets AWK variable `var` to the maximum fragment size (default: 500bp from `params.insert_size`)

2. **`substr($0,1,1)=="@"`**: Checks if the line is a SAM header (starts with `@`)
   - If true: **KEEP** the header line (required for valid SAM/BAM files)

3. **`(($9>=0?$9:-$9)<=var)`**: Checks the absolute value of TLEN (field 9)
   - `$9>=0?$9:-$9`: Computes |TLEN| (absolute value)
   - `<=var`: Checks if |TLEN| ≤ max_frag
   - If true: **KEEP** the read pair

4. **Otherwise**: **DISCARD** the read pair (fragment too long)

**TLEN (Template Length) Field:**
- Column 9 in SAM format
- Signed integer: + for forward read, - for reverse read
- Represents insert size (distance between read pairs)

**Example:**
```
Read1  99   chr1  1000  60  100M  =  1200  300  ...  (TLEN = +300)
Read2  147  chr1  1200  60  100M  =  1000  -300 ...  (TLEN = -300)
```

For `params.insert_size = 500`:
- ✅ Keeps fragments with |TLEN| ≤ 500bp (both reads of the pair)
- ❌ Removes fragments with |TLEN| > 500bp (potential chimeras, artifacts)

**Why use absolute value?**
- TLEN is positive for the forward read (+300) and negative for the reverse read (-300)
- We want to filter based on fragment **size**, not direction
- `($9>=0?$9:-$9)` computes the absolute value in AWK (equivalent to `abs($9)`)

#### Stage 5: Position Sorting and Indexing
```bash
samtools sort        # Sort by genomic coordinate
samtools index       # Create BAM index
```

### Single-End (SE) Filtering

SE filtering is simpler since there's no fragment information:

#### Stage 1: SAMtools Filter
```bash
samtools view -h -b \
    -F 0x004 \     # Remove unmapped reads
    -q 1           # Remove MAPQ < 1
```

**No fragment filtering:**
- No TLEN field (always 0)
- No mate information
- Only read quality matters

#### Stage 2: Sorting and Indexing
```bash
samtools sort        # Sort by genomic coordinate
samtools index       # Create BAM index
```

---

## 📊 Part 5: Two-Stage Filtering Strategy (PE Only)

### Why Two Stages?

```
┌─────────────────────────────────────────────────┐
│ STAGE 1: Bowtie2 Alignment                     │
│ Parameter: -X 1000 (FIXED)                     │
│                                                 │
│ Purpose: Permissive alignment search           │
│ Result: All valid pairs up to 1000bp aligned   │
└──────────────────┬──────────────────────────────┘
                   │
                   ▼
            Unfiltered BAM
         (fragments 0-1000bp)
                   │
                   ▼
┌─────────────────────────────────────────────────┐
│ STAGE 2: BAM_FILTER                            │
│ Parameter: params.insert_size = 500 (default)  │
│                                                 │
│ Purpose: Biological quality control            │
│ Result: Only high-quality fragments kept       │
└──────────────────┬──────────────────────────────┘
                   │
                   ▼
            Filtered BAM
         (fragments 0-500bp)
                   │
                   ▼
          Peak calling (MACS2)
```

### Benefits of Decoupling

| Aspect | Bowtie2 -X | BAM_FILTER insert_size |
|--------|-----------|------------------------|
| **Purpose** | Find valid alignments | Remove biological artifacts |
| **Fixed vs Variable** | Fixed at 1000 | User-configurable |
| **When Applied** | During alignment | Post-alignment |
| **Can Change** | No (requires re-alignment) | Yes (just re-filter BAM) |

### Recommended Values by Experiment Type

```bash
# Narrow peaks (transcription factors)
nextflow run pdichiaro/atacseq --insert_size 400

# Standard ATAC-seq (H3K4me3, H3K27ac)
nextflow run pdichiaro/atacseq --insert_size 500  # DEFAULT

# Broad marks (H3K27me3, H3K36me3)
nextflow run pdichiaro/atacseq --insert_size 600

# Very permissive (keep longer fragments)
nextflow run pdichiaro/atacseq --insert_size 800

# No post-alignment filtering (use only Bowtie2)
nextflow run pdichiaro/atacseq --insert_size 1000
```

---

## 🧪 Part 6: Quality Metrics Impact

### Fragment Size Distribution Analysis

The pipeline generates fragment size distribution plots showing:
- Distribution before BAM_FILTER (0-1000bp)
- Distribution after BAM_FILTER (0-500bp)

**Expected patterns:**

#### Good ATAC-seq Library
```
  Count
    ▲
    │     ╱╲
    │    ╱  ╲
    │   ╱    ╲___
    │  ╱         ╲___
    │ ╱              ╲___
    └─────────────────────→ Fragment size (bp)
      0   200  500  800 1000
      
Peak around 200-300bp (nucleosome-sized)
```

#### Library with Artifacts
```
  Count
    ▲
    │     ╱╲
    │    ╱  ╲  ← Nucleosomal
    │   ╱    ╲
    │  ╱      ╲╱╲ ← Chimeras/artifacts
    │ ╱           ╲╱╲
    └─────────────────────→ Fragment size (bp)
      0   200  500  800 1000
      
Secondary peaks > 500bp indicate problems
```

### Impact on Alignment Statistics

| Metric | Before BAM_FILTER | After BAM_FILTER |
|--------|------------------|------------------|
| **Aligned pairs** | 100% | 95-98% |
| **Fragments retained** | 0-1000bp | 0-500bp |
| **Avg fragment size** | 350bp | 280bp |
| **Potential chimeras** | Included | Removed |

---

## 🔍 Part 7: Troubleshooting

### Issue 1: Low Alignment Rate (< 70%)

**Possible causes:**
- Wrong reference genome
- Poor read quality
- High adapter contamination
- Wrong species

**Solutions:**
```bash
# Check FastQC reports
# Run adapter trimming:
nextflow run pdichiaro/atacseq --trim_nextseq 20
```

### Issue 2: High Discordant Rate (> 10%)

**Possible causes:**
- Library prep issues
- Contamination
- Incorrect -X value (too restrictive)

**Our solution:**
- Fixed -X 1000 → prevents false discordant calls

### Issue 3: Unusual Fragment Size Distribution

**Pattern: Bimodal distribution with peak > 500bp**

```bash
# Investigate with:
samtools view filtered.bam | awk '{print sqrt($9*$9)}' | \
  sort | uniq -c | sort -rn | head -50
```

**Action:**
- If real biology → increase `--insert_size`
- If artifacts → keep default (500bp)

### Issue 4: Too Few Fragments After Filtering

**Symptom:** BAM_FILTER removes > 20% of fragments

**Diagnosis:**
```bash
# Compare unfiltered vs filtered
samtools flagstat unfiltered.bam
samtools flagstat filtered.bam
```

**Solutions:**
- Increase `--insert_size` to 600 or 700
- Check if library has unusually long fragments
- Verify library prep protocol

---

## 📚 Part 8: Technical References

### SAM/BAM Format Fields

```
QNAME  FLAG  RNAME  POS  MAPQ  CIGAR  RNEXT  PNEXT  TLEN  SEQ  QUAL
  1      2     3     4     5      6      7      8     9    10   11
```

**Key fields:**
- **TLEN (9)**: Template length (fragment size)
- **FLAG (2)**: Bit flag encoding read properties
- **MAPQ (5)**: Mapping quality (Phred-scaled)

### Bowtie2 Exit Codes and Messages

```
# Normal completion
100 reads; of these:
  100 (100.00%) were paired; of these:
    10 (10.00%) aligned concordantly 0 times
    80 (80.00%) aligned concordantly exactly 1 time
    10 (10.00%) aligned concordantly >1 times
90.00% overall alignment rate
```

**Target metrics:**
- Concordant alignment rate > 70%
- Multiple alignments < 20%

### Pipeline File Flow

```
input_fastq.gz
    ↓ [FASTQC]
quality_reports/
    ↓ [BOWTIE2]
unfiltered.bam (0-1000bp fragments)
    ↓ [BAM_FILTER]
filtered.bam (0-500bp fragments)
    ↓ [PICARD]
deduped.bam
    ↓ [MACS2]
peaks.narrowPeak
```

---

## ✅ Part 9: Best Practices Summary

### For Paired-End ATAC-seq

1. ✅ **Use --local alignment mode** (handles adapter contamination)
2. ✅ **Keep Bowtie2 -X at 1000** (permissive search)
3. ✅ **Filter BAM with insert_size = 500** (remove artifacts)
4. ✅ **Remove discordant and mixed alignments** (clean pairs only)
5. ✅ **Check fragment size distribution** (should peak ~250bp)

### For Single-End ATAC-seq

1. ✅ **Use --local alignment mode**
2. ✅ **Apply MAPQ filtering** (q ≥ 1)
3. ✅ **No fragment filtering needed** (no TLEN available)
4. ✅ **Extend reads in MACS2** (--extsize based on expected fragment size)

### General

1. ✅ **Always run FastQC first** (detect quality issues early)
2. ✅ **Use deduplication** (remove PCR duplicates)
3. ✅ **Generate alignment statistics** (monitor pipeline health)
4. ✅ **Compare biological replicates** (ensure consistency)

---

## 🎯 Conclusion

This pipeline implements a **two-stage filtering strategy**:

1. **Bowtie2 (-X 1000)**: Permissive alignment to avoid missing valid fragments
2. **BAM_FILTER (insert_size = 500)**: Stringent biological QC to remove artifacts

This approach balances:
- ✅ **Sensitivity** (don't miss real alignments)
- ✅ **Specificity** (remove chimeras and artifacts)
- ✅ **Flexibility** (users can tune filtering without re-alignment)

For most ATAC-seq experiments, the default configuration provides optimal results. Adjust `--insert_size` based on your specific biological system and experimental design.

---

**Document Version:** 1.0  
**Last Updated:** 2026-03-11  
**Pipeline:** pdichiaro/atacseq  
**Compatible with:** Nextflow 25.04+, Bowtie2 2.3+, SAMtools 1.9+
