# Timeline della rimozione del genoma mitocondriale

## 📍 Quando viene rimosso il cromosoma mitocondriale?

Il cromosoma mitocondriale viene rimosso **durante il filtraggio dei file BAM**, dopo l'allineamento e la rimozione dei duplicati.

---

## 🔄 Flusso completo del pipeline

```
┌─────────────────────────────────────────────────────────────────┐
│ FASE 1: PREPARAZIONE DEL GENOMA                                 │
│ (Subworkflow PREPARE_GENOME)                                    │
└─────────────────────────────────────────────────────────────────┘

1. CUSTOM_GETCHROMSIZES
   ├─ Input: genome.fasta
   └─ Output: genome.sizes (lista cromosomi con dimensioni)

2. GENOME_BLACKLIST_REGIONS ⚠️ [QUI VIENE CREATO IL FILE DI FILTRAGGIO]
   ├─ Input: 
   │  ├─ genome.sizes
   │  ├─ blacklist.bed (opzionale)
   │  ├─ mito_name (es. "chrM", "MT", "Mt")
   │  ├─ keep_mito (true/false)
   │  └─ keep_blacklist (true/false)
   │
   ├─ Logica:
   │  │
   │  ├─ Se keep_mito = false E mito_name è definito:
   │  │  └─ Rimuove cromosoma mitocondriale con:
   │  │     awk '$1 !~ /chrM/' (o MT, Mt, ecc.)
   │  │
   │  └─ Se keep_blacklist = false:
   │     └─ Rimuove regioni blacklist ENCODE
   │
   └─ Output: genome.include_regions.bed
      └─ Contiene SOLO le regioni da INCLUDERE nell'analisi
         (tutto il genoma ESCLUSO chrM e blacklist)


┌─────────────────────────────────────────────────────────────────┐
│ FASE 2: ALLINEAMENTO                                            │
└─────────────────────────────────────────────────────────────────┘

3. BOWTIE2_ALIGN
   ├─ Input: reads FASTQ
   └─ Output: aligned.bam (contiene TUTTE le reads allineate)
      └─ Include reads allineate a chrM ✓

4. MARK_DUPLICATES_PICARD
   ├─ Input: aligned.bam
   └─ Output: aligned.marked.bam (con duplicati marcati)
      └─ Include ancora reads allineate a chrM ✓


┌─────────────────────────────────────────────────────────────────┐
│ FASE 3: FILTRAGGIO BAM (QUI VIENE RIMOSSO chrM) 🎯              │
│ (Subworkflow BAM_FILTER_SUBWF)                                  │
└─────────────────────────────────────────────────────────────────┘

5. BAM_FILTER (Process)
   ├─ Input: 
   │  ├─ aligned.marked.bam (con tutte le reads)
   │  └─ genome.include_regions.bed (SENZA chrM)
   │
   ├─ Operazioni:
   │  │
   │  ├─ samtools view -L genome.include_regions.bed
   │  │  └─ Mantiene SOLO reads nelle regioni del BED
   │  │     ⚠️ chrM NON è nel BED → reads mitocondriali RIMOSSE
   │  │
   │  ├─ Altri filtri applicati:
   │  │  ├─ Rimuove secondary alignments (-F 0x100)
   │  │  ├─ Rimuove supplementary alignments (-F 0x800)
   │  │  ├─ Rimuove duplicati se keep_dups=false (-F 0x0400)
   │  │  ├─ Rimuove multi-mappers se keep_multi_map=false (-q 1)
   │  │  └─ Filtra fragment size (default ≤500bp)
   │  │
   │  └─ Output: filtered.bam
   │     └─ SENZA reads mitocondriali ✗
   │
   └─ Output: sample.filter2.bam
      └─ BAM pulito, pronto per analisi downstream


┌─────────────────────────────────────────────────────────────────┐
│ FASE 4: ANALISI DOWNSTREAM                                      │
└─────────────────────────────────────────────────────────────────┘

6. Tutti i moduli successivi usano filtered.bam:
   ├─ Peak calling (MACS2)
   ├─ BigWig generation
   ├─ Fragment size distribution
   ├─ TSS enrichment
   ├─ FRiP score calculation
   └─ MultiQC report
   
   → Nessuno di questi vede reads mitocondriali
```

---

## 🔍 Dettaglio tecnico del filtraggio

### Contenuto di `genome.include_regions.bed` (esempio GRCh38):

**Con keep_mito = false (default):**
```
chr1    0    248956422
chr2    0    242193529
chr3    0    198295559
...
chrX    0    156040895
chrY    0    57227415
# chrM NON è presente → reads su chrM vengono rimosse
```

**Con keep_mito = true:**
```
chr1    0    248956422
chr2    0    242193529
chr3    0    198295559
...
chrX    0    156040895
chrY    0    57227415
chrM    0    16569      ← chrM è incluso
```

### Comando samtools che applica il filtro:

```bash
samtools view \
    -F 0x0100 \           # Rimuove secondary alignments
    -F 0x0800 \           # Rimuove supplementary alignments
    -F 0x0400 \           # Rimuove duplicati (se keep_dups=false)
    -L genome.include_regions.bed \  # 🎯 FILTRA PER REGIONI NEL BED
    -b aligned.marked.bam \
    > filtered.step1.bam

# Poi applica altri filtri (MAPQ, fragment size)
```

Il flag `-L` di samtools dice: **"Mantieni SOLO le reads che intersecano le regioni nel file BED"**

Se chrM non è nel BED → reads mitocondriali vengono scartate.

---

## 📊 Statistiche di filtraggio

Il modulo `BLACKLIST_LOG` genera statistiche che mostrano:

```
Reads before filtering: 10,000,000
Reads after filtering:   8,500,000
Reads removed:          1,500,000

Breakdown:
- Mitochondrial reads:    800,000  (53%)  ← se keep_mito=false
- Blacklist regions:      300,000  (20%)
- Duplicates:             200,000  (13%)
- Multi-mappers:          150,000  (10%)
- Fragment size filter:    50,000  (4%)
```

---

## ⚙️ Parametri che controllano il comportamento

| Parametro | Default | Effetto |
|-----------|---------|---------|
| `--keep_mito` | `false` | Se `false`, rimuove reads mitocondriali |
| `--keep_blacklist` | `false` | Se `false`, rimuove regioni blacklist |
| `--keep_dups` | `false` | Se `false`, rimuove duplicati PCR |
| `--keep_multi_map` | `false` | Se `false`, rimuove multi-mappers (MAPQ=0) |
| `--insert_size` | `500` | Massima dimensione frammento da mantenere (bp) |

---

## 🎯 Risposta diretta alla domanda

**Il genoma mitocondriale viene rimosso:**

1. **Fase**: Durante il filtraggio dei file BAM (dopo allineamento e rimozione duplicati)
2. **Modulo**: `BAM_FILTER` (process in `modules/local/bam_filter.nf`)
3. **Meccanismo**: `samtools view -L genome.include_regions.bed`
4. **File chiave**: `genome.include_regions.bed` (creato da `GENOME_BLACKLIST_REGIONS`)
5. **Timing**: Dopo `MARK_DUPLICATES_PICARD`, prima di peak calling
6. **Controllo utente**: Parametro `--keep_mito` (default: false)

**Le reads mitocondriali sono presenti in:**
- ✓ File BAM iniziali dopo allineamento (`BOWTIE2_ALIGN.out.bam`)
- ✓ File BAM dopo rimozione duplicati (`MARK_DUPLICATES_PICARD.out.bam`)

**Le reads mitocondriali sono ASSENTI in:**
- ✗ File BAM filtrati (`BAM_FILTER.out.bam`)
- ✗ Tutti i file downstream (peaks, BigWig, metriche QC)

---

## 💡 Perché questo approccio?

1. **Standard ATAC-seq**: Reads mitocondriali (5-50% del totale) distorcono normalizzazione
2. **Ottimizzazione**: Meglio rimuovere dopo duplicati (evita duplicati su chrM)
3. **Flessibilità**: Utente può ispezionare BAM pre-filtro se necessario
4. **Performance**: Un solo passaggio di filtraggio per tutti i criteri

---

## 🧪 Come verificare

```bash
# Conta reads mitocondriali PRIMA del filtraggio
samtools view -c aligned.marked.bam chrM

# Conta reads mitocondriali DOPO il filtraggio
samtools view -c filtered.bam chrM
# Output: 0 (se keep_mito=false)

# Verifica contenuto del BED di filtraggio
grep "chrM\|MT" genome.include_regions.bed
# Nessun output se keep_mito=false
```
