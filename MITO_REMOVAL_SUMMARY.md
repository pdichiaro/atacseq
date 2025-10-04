# 🧬 Rimozione del Genoma Mitocondriale - Sommario Rapido

## 🎯 Risposta breve

Il cromosoma mitocondriale viene rimosso **durante il filtraggio dei file BAM**, nella fase 3 del pipeline, DOPO l'allineamento e la rimozione dei duplicati.

---

## 📋 Flusso semplificato

```
FASE 1: Preparazione Genoma
   └─ GENOME_BLACKLIST_REGIONS
      └─ Crea: genome.include_regions.bed
         └─ Contiene: tutti i cromosomi TRANNE chrM
         
         
FASE 2: Allineamento
   ├─ BOWTIE2_ALIGN → aligned.bam (con chrM ✓)
   └─ MARK_DUPLICATES → aligned.marked.bam (con chrM ✓)
   

FASE 3: Filtraggio ⚠️ QUI VIENE RIMOSSO chrM
   └─ BAM_FILTER
      ├─ samtools view -L genome.include_regions.bed
      └─ Output: filtered.bam (SENZA chrM ✗)


FASE 4: Analisi downstream
   └─ Peak calling, BigWig, QC, ecc.
      └─ Usano tutti filtered.bam (SENZA chrM)
```

---

## 🔧 Meccanismo tecnico

1. **Creazione del file BED di filtraggio** (`GENOME_BLACKLIST_REGIONS`):
   ```bash
   # Se keep_mito = false
   awk '$1 !~ /chrM/' genome.sizes > genome.include_regions.bed
   ```
   
   Risultato:
   ```
   chr1    0    248956422
   chr2    0    242193529
   ...
   chrX    0    156040895
   chrY    0    57227415
   # chrM NON è presente!
   ```

2. **Applicazione del filtro** (`BAM_FILTER`):
   ```bash
   samtools view -L genome.include_regions.bed aligned.marked.bam
   ```
   
   Il flag `-L` mantiene SOLO reads nelle regioni del BED.
   
   chrM non è nel BED → reads mitocondriali vengono scartate.

---

## 📊 Confronto file BAM

| File BAM | Fase | Contiene chrM? |
|----------|------|----------------|
| `aligned.bam` | Dopo BOWTIE2 | ✅ SÌ |
| `aligned.marked.bam` | Dopo PICARD | ✅ SÌ |
| `filtered.bam` | Dopo BAM_FILTER | ❌ NO (se keep_mito=false) |
| `*.bigWig` | Visualizzazione | ❌ NO |
| `*.peaks.bed` | Peak calling | ❌ NO |

---

## 🎚️ Controllo utente

```bash
# Default: rimuove chrM
nextflow run pdichiaro/atacseq --input samples.csv --genome GRCh38

# Mantieni chrM
nextflow run pdichiaro/atacseq --input samples.csv --genome GRCh38 --keep_mito
```

---

## 🧪 Come verificare

```bash
# Prima del filtraggio (file con chrM)
samtools idxstats aligned.marked.bam | grep chrM
# Output: chrM    16569    123456    0

# Dopo il filtraggio (file senza chrM)
samtools idxstats filtered.bam | grep chrM
# Output: chrM    16569    0    0    ← ZERO reads su chrM!
```

---

## ❓ Domande frequenti

### Perché rimuovere il mitocondrio?

- Reads mitocondriali rappresentano 5-50% del totale in ATAC-seq
- Distorcono normalizzazione e peak calling
- Standard practice per analisi cromatina nucleare

### Perché rimuoverlo DOPO l'allineamento?

- Permette di quantificare contaminazione mitocondriale (metrica QC)
- Rimuovere duplicati prima ottimizza performance
- Mantiene flessibilità per analisi alternative

### Cosa succede se --keep_mito è attivo?

Il file `genome.include_regions.bed` INCLUDE chrM:
```
chr1    0    248956422
...
chrM    0    16569      ← chrM presente!
```

Reads mitocondriali vengono mantenute nell'analisi.

---

## 📚 File da consultare

- **Implementazione**: `modules/local/bam_filter.nf`
- **Creazione BED**: `modules/local/genome_blacklist_regions.nf`
- **Dettagli completi**: `MITO_REMOVAL_TIMELINE.md`
- **Configurazione genomi**: `conf/igenomes.config`
