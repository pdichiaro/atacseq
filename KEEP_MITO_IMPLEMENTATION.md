# Implementazione del parametro keep_mito

## Sommario
Il parametro `keep_mito` è stato implementato nel pipeline **pdichiaro/atacseq** per permettere agli utenti di controllare se mantenere o rimuovere le reads mitocondriali dall'analisi ATAC-seq.

## Modifiche implementate

### 1. Configurazione (`nextflow.config`)
- ✅ Parametro `keep_mito = false` già definito (default: rimuove reads mitocondriali)

### 2. Schema JSON (`nextflow_schema.json`)
- ✅ Documentazione completa del parametro con descrizione e icona

### 3. Configurazione genomi (`conf/igenomes.config`)
- ✅ Aggiunto `mito_name` ai genomi mancanti:
  - **Gm01** (Glycine max): `mito_name = "Mt"`
  - **Sbi1** (Sorghum bicolor): `mito_name = "Mt"`
- ✅ Tutti i genomi eucarioti ora hanno `mito_name` configurato
- ℹ️ Solo batteri (EB1, EB2) non hanno `mito_name` (corretto, non hanno mitocondri)

### 4. Modulo di filtraggio (`modules/local/genome_blacklist_regions.nf`)
- ✅ Implementata logica per filtrare cromosoma mitocondriale:
  - Accetta parametri `mito_name` e `keep_mito`
  - Se `keep_mito = false` E `mito_name` è definito → rimuove cromosoma mitocondriale
  - Utilizza filtro AWK: `awk '$1 !~ /MT/'` (o chrM, Mt, M, MtDNA, ecc.)
- ✅ Funziona sia con blacklist che senza blacklist

### 5. Subworkflow preparazione genoma (`subworkflows/local/prepare_genome.nf`)
- ✅ Aggiunto parametro `mito_name` come input del workflow
- ✅ Passa `mito_name` e `params.keep_mito` al modulo `GENOME_BLACKLIST_REGIONS`

### 6. Workflow principale (`workflows/atacseq.nf`)
- ✅ Estrae `mito_name` dalla configurazione del genoma selezionato
- ✅ Passa `mito_name` al subworkflow `PREPARE_GENOME`

## Comportamento

### Default (keep_mito = false)
```bash
nextflow run pdichiaro/atacseq --input samplesheet.csv --genome GRCh38
```
- Rimuove automaticamente reads da cromosoma MT (in GRCh38)
- Ottimale per analisi ATAC-seq standard

### Mantenere reads mitocondriali (keep_mito = true)
```bash
nextflow run pdichiaro/atacseq --input samplesheet.csv --genome GRCh38 --keep_mito
```
- Mantiene tutte le reads, incluse quelle mitocondriali
- Utile per analisi specifiche o troubleshooting

## Nomi cromosomi mitocondriali per genoma

| Organismo | Genoma | mito_name |
|-----------|--------|-----------|
| Homo sapiens | GRCh37 | MT |
| Homo sapiens | GRCh38 | chrM |
| Mus musculus | GRCm38/GRCm39 | MT |
| Arabidopsis thaliana | TAIR10 | Mt |
| Bos taurus | UMD3.1 | MT |
| Caenorhabditis elegans | WBcel235 | MtDNA |
| Canis familiaris | CanFam3.1 | MT |
| Danio rerio | GRCz10/GRCz11 | MT |
| Drosophila melanogaster | BDGP6 | M |
| Equus caballus | EquCab2 | MT |
| Gallus gallus | Galgal4 | MT |
| Glycine max | Gm01 | Mt |
| Macaca mulatta | Mmul_1 | MT |
| Oryza sativa | IRGSP-1.0 | Mt |
| Pan troglodytes | CHIMP2.1.4 | MT |
| Rattus norvegicus | Rnor_5.0/6.0 | MT |
| Saccharomyces cerevisiae | R64-1-1 | MT |
| Schizosaccharomyces pombe | EF2 | MT |
| Sorghum bicolor | Sbi1 | Mt |
| Sus scrofa | Sscrofa10.2/11.1 | MT |
| Zea mays | AGPv3/AGPv4 | Mt |

## Vantaggi dell'implementazione

1. **Compatibilità nf-core**: Approccio allineato con nf-core/atacseq
2. **Flessibilità**: Utente può scegliere basandosi sulle esigenze dell'analisi
3. **Automatico**: Rimozione gestita automaticamente dal pipeline
4. **Genoma-specifico**: Usa il nome corretto del cromosoma per ogni organismo
5. **Documentato**: Parametro documentato nello schema JSON

## Test consigliati

```bash
# Test 1: Default (rimozione mitocondriale)
nextflow run pdichiaro/atacseq --input test.csv --genome GRCh38

# Test 2: Mantenimento mitocondriale
nextflow run pdichiaro/atacseq --input test.csv --genome GRCh38 --keep_mito

# Verifica risultati
samtools idxstats sample.bam | grep -E "(chrM|MT|Mt)"
```

## File modificati

1. `conf/igenomes.config` - Aggiunto mito_name a Gm01 e Sbi1
2. `modules/local/genome_blacklist_regions.nf` - Implementata logica filtraggio
3. `subworkflows/local/prepare_genome.nf` - Aggiunto parametro mito_name
4. `workflows/atacseq.nf` - Estrazione mito_name da config genoma

## Note

- **Batteri (EB1, EB2)**: Non hanno `mito_name` perché non hanno mitocondri
- **Performance**: La rimozione avviene a livello di regioni genomiche, molto efficiente
- **Compatibilità**: Funziona con tutti i genomi configurati in igenomes.config
