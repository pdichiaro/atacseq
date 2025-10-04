# 📊 Esempio pratico: rimozione reads mitocondriali

## Scenario: Campione ATAC-seq umano (GRCh38)

### Input iniziale
```
Sample: ATAC_heart_rep1
Total raw reads: 50,000,000 paired-end reads
Genome: GRCh38 (chrM present)
```

---

## 📍 Tracciamento reads attraverso il pipeline

### FASE 1: Preparazione Genoma

**Processo: GENOME_BLACKLIST_REGIONS**

Input:
```bash
# genome.sizes
chr1    248956422
chr2    242193529
...
chrX    156040895
chrY    57227415
chrM    16569       ← Cromosoma mitocondriale presente
```

Con `params.keep_mito = false`:
```bash
awk '$1 !~ /chrM/' genome.sizes > genome.include_regions.bed
```

Output: `genome.include_regions.bed`
```bash
chr1    0    248956422
chr2    0    242193529
...
chrX    0    156040895
chrY    0    57227415
# chrM NON è presente → verrà usato per filtrare le reads
```

---

### FASE 2: Allineamento

**Processo: BOWTIE2_ALIGN**

```bash
bowtie2 -x genome -1 R1.fq.gz -2 R2.fq.gz | samtools view -bS - > aligned.bam
samtools index aligned.bam
```

Statistiche `aligned.bam`:
```bash
samtools idxstats aligned.bam

chr1    248956422    8500000    50000
chr2    242193529    7200000    45000
...
chrX    156040895    2100000    15000
chrY    57227415     500000     3000
chrM    16569        12000000   80000  ← 12 milioni di reads su chrM!
*       0            1200000    0
# Total: 50,000,000 reads
```

**Breakdown:**
- Autosomal chromosomes: 35,500,000 reads (71%)
- Sex chromosomes: 2,600,000 reads (5.2%)
- **Mitochondrial (chrM): 12,000,000 reads (24%)**  ← Molto comune in ATAC-seq!
- Unmapped: 1,200,000 reads (2.4%)

---

**Processo: MARK_DUPLICATES_PICARD**

```bash
picard MarkDuplicates I=aligned.bam O=aligned.marked.bam M=metrics.txt
```

Statistiche `aligned.marked.bam`:
```bash
samtools idxstats aligned.marked.bam

chr1    248956422    8500000    50000    (800k duplicates marked)
chr2    242193529    7200000    45000    (650k duplicates marked)
...
chrM    16569        12000000   80000    (5M duplicates marked!)
*       0            1200000    0
# Total ancora: 50,000,000 reads (duplicati solo marcati, non rimossi)
```

**Note:** 
- chrM ha MOLTISSIMI duplicati (5M su 12M = 42%)
- Normale per mitocondrio (alta copia, piccolo genoma)
- Duplicati ancora presenti ma marcati con flag 0x400

---

### FASE 3: Filtraggio BAM ⚠️ RIMOZIONE chrM

**Processo: BAM_FILTER**

```bash
# Step 1: General filtering
samtools view \
    -F 0x0100 \      # Rimuove secondary alignments
    -F 0x0800 \      # Rimuove supplementary alignments  
    -F 0x0400 \      # Rimuove duplicates (keep_dups=false)
    -L genome.include_regions.bed \  # 🎯 APPLICA FILTRO REGIONI (senza chrM!)
    -b aligned.marked.bam \
    > filtered.step1.bam

# Step 2: MAPQ and fragment size filtering
samtools view -q 1 -h filtered.step1.bam | \
    awk -v var="500" '{if(substr($0,1,1)=="@" || (($9>=0?$9:-$9)<=var)) print $0}' | \
    samtools view -b > filtered.bam
```

**Reads rimosse in ogni step:**

| Step | Reads Removed | Reason |
|------|---------------|--------|
| `-L genome.include_regions.bed` | 12,000,000 | **Mitochondrial reads (chrM non nel BED)** |
| `-F 0x0100` | 500,000 | Secondary alignments |
| `-F 0x0800` | 200,000 | Supplementary alignments |
| `-F 0x0400` | 8,500,000 | PCR duplicates (genomic + mitochondrial) |
| `-q 1` (MAPQ) | 3,200,000 | Multi-mappers |
| Fragment size >500bp | 1,100,000 | Long fragments/artifacts |
| **TOTAL REMOVED** | **25,500,000** | **51% of original reads** |

Statistiche `filtered.bam`:
```bash
samtools idxstats filtered.bam

chr1    248956422    4200000    0
chr2    242193529    3600000    0
...
chrX    156040895    1050000    0
chrY    57227415     250000     0
chrM    16569        0          0    ← ZERO reads su chrM!
*       0            0          0
# Total: 24,500,000 reads (49% retained)
```

**Breakdown finale:**
- Autosomal chromosomes: 22,800,000 reads (93.1%)
- Sex chromosomes: 1,300,000 reads (5.3%)
- **Mitochondrial (chrM): 0 reads (0%)**  ✅ Completamente rimosso!
- Unmapped: 0 reads (0%)

---

### FASE 4: Analisi Downstream

Tutti i moduli successivi usano `filtered.bam` (SENZA chrM):

**Peak Calling (MACS2):**
```bash
macs2 callpeak -t filtered.bam -n sample -f BAM --shift -75 --extsize 150
```

Output: `sample_peaks.narrowPeak`
```bash
chr1    1245000    1245500    peak_1    100
chr2    5678000    5678300    peak_2    95
...
# NO peaks su chrM (non ci sono reads!)
```

**BigWig generation:**
```bash
bamCoverage -b filtered.bam -o sample.bigWig --normalizeUsing CPM
```

Coverage su chrM:
```bash
# Visualizzando in IGV su chrM → coverage = 0
```

**MultiQC metriche:**
```
Sample: ATAC_heart_rep1
Total reads: 50,000,000
Aligned reads: 48,800,000 (97.6%)
Mitochondrial reads: 12,000,000 (24.0%)  ← Reportato per QC
Duplicates: 8,500,000 (17.4%)
Filtered reads (final): 24,500,000 (49.0%)
```

---

## 🎚️ Confronto: keep_mito = true vs false

### Con `--keep_mito false` (default):

```
Reads in filtered.bam: 24,500,000
  ├─ Nuclear genome: 24,500,000 (100%)
  └─ Mitochondrial: 0 (0%)

Peaks called: 45,230
  ├─ Nuclear: 45,230 (100%)
  └─ Mitochondrial: 0 (0%)
```

### Con `--keep_mito true`:

```
Reads in filtered.bam: 31,500,000
  ├─ Nuclear genome: 24,500,000 (77.8%)
  └─ Mitochondrial: 7,000,000 (22.2%)  ← chrM reads retained!

Peaks called: 45,891
  ├─ Nuclear: 45,230 (98.6%)
  └─ Mitochondrial: 661 (1.4%)  ← Peaks called on chrM!

⚠️ Issues:
- Normalization skewed by high chrM coverage
- Nuclear peaks might be under-called
- TSS enrichment calculation affected
```

---

## 📈 Impatto sulle metriche QC

### FRiP (Fraction of Reads in Peaks):

**With keep_mito = false:**
```
FRiP = (reads in peaks) / (total filtered reads)
     = 18,375,000 / 24,500,000
     = 0.75 (75%)  ← Excellent ATAC-seq quality!
```

**With keep_mito = true:**
```
FRiP = 18,375,000 / 31,500,000
     = 0.58 (58%)  ← Artificially lowered by chrM reads
```

### TSS Enrichment:

**With keep_mito = false:**
```
TSS enrichment = 15.2  ← Good signal-to-noise
```

**With keep_mito = true:**
```
TSS enrichment = 11.8  ← Diluted by chrM background
```

---

## ✅ Conclusione

Il cromosoma mitocondriale viene rimosso:

1. **Dove**: Nel processo `BAM_FILTER`
2. **Quando**: Dopo allineamento e rimozione duplicati, prima del peak calling
3. **Come**: Via `samtools view -L genome.include_regions.bed`
4. **Perché**: Migliora qualità analisi cromatina nucleare

**File con chrM:**
- ✅ `aligned.bam` (dopo BOWTIE2)
- ✅ `aligned.marked.bam` (dopo PICARD)

**File senza chrM:**
- ❌ `filtered.bam` (dopo BAM_FILTER)
- ❌ Tutti gli output downstream

**Comando per verificare:**
```bash
# Check reads on chrM before filtering
samtools view -c aligned.marked.bam chrM
# Output: 12000000

# Check reads on chrM after filtering
samtools view -c filtered.bam chrM
# Output: 0
```
