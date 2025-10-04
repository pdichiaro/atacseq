#!/bin/bash
# Verifica l'efficacia della blacklist filtering

SAMPLE=$1
OUTDIR=${2:-results}
BLACKLIST=${3:-/path/to/blacklist.bed}

echo "=== VERIFICA BLACKLIST FILTERING ==="
echo "Sample: $SAMPLE"
echo "Blacklist: $BLACKLIST"
echo ""

INPUT_BAM="${OUTDIR}/bowtie2/mergedLibrary/mdup/${SAMPLE}.mLb.mkD.sorted.bam"
FILTERED_BAM="${OUTDIR}/bowtie2/mergedLibrary/Final_BAM/${SAMPLE}.mLb.clN.sorted.bam"

# 1. CONTA READS IN BLACKLIST REGIONS
echo "### 1. READS IN BLACKLIST REGIONS ###"
echo ""

READS_IN_BL_INPUT=$(bedtools intersect -a "$INPUT_BAM" -b "$BLACKLIST" -u | samtools view -c -)
READS_IN_BL_FILTERED=$(bedtools intersect -a "$FILTERED_BAM" -b "$BLACKLIST" -u | samtools view -c -)

TOTAL_INPUT=$(samtools view -c "$INPUT_BAM")
TOTAL_FILTERED=$(samtools view -c "$FILTERED_BAM")

PCT_BL_INPUT=$(awk "BEGIN {printf \"%.2f\", ($READS_IN_BL_INPUT / $TOTAL_INPUT) * 100}")
PCT_BL_FILTERED=$(awk "BEGIN {printf \"%.2f\", ($READS_IN_BL_FILTERED / $TOTAL_FILTERED) * 100}")

echo "INPUT BAM:"
echo "  Total reads:              $TOTAL_INPUT"
echo "  Reads in blacklist:       $READS_IN_BL_INPUT"
echo "  Percentage in blacklist:  ${PCT_BL_INPUT}%"
echo ""

echo "FILTERED BAM:"
echo "  Total reads:              $TOTAL_FILTERED"
echo "  Reads in blacklist:       $READS_IN_BL_FILTERED"
echo "  Percentage in blacklist:  ${PCT_BL_FILTERED}%"
echo ""

REMOVED=$(($READS_IN_BL_INPUT - $READS_IN_BL_FILTERED))
echo "Reads removed from blacklist: $REMOVED"
echo ""

# 2. TOP BLACKLIST REGIONS CON PIÙ READS
echo "### 2. TOP 10 BLACKLIST REGIONS (INPUT BAM) ###"
echo ""

bedtools coverage -a "$BLACKLIST" -b "$INPUT_BAM" | \
    sort -k4 -rn | \
    head -10 | \
    awk '{printf "%-10s %10d-%10d | %8d reads | %6.2f%% coverage\n", 
          $1, $2, $3, $4, $7*100}'
echo ""

# 3. DISTRIBUZIONE PER CROMOSOMA
echo "### 3. BLACKLIST READS PER CROMOSOMA ###"
echo ""

echo "--- INPUT BAM ---"
bedtools intersect -a "$INPUT_BAM" -b "$BLACKLIST" -u | \
    samtools view - | \
    awk '{print $3}' | \
    sort | \
    uniq -c | \
    sort -rn | \
    head -10 | \
    awk '{printf "%-15s %10d reads\n", $2, $1}'
echo ""

# 4. COVERAGE PLOT DATA (per visualizzazione)
echo "### 4. GENERAZIONE DATI PER PLOT ###"
echo ""

OUTPUT_COV="${OUTDIR}/qc/blacklist_coverage.txt"
mkdir -p "$(dirname $OUTPUT_COV)"

echo -e "chrom\tstart\tend\tcoverage_input\tcoverage_filtered" > "$OUTPUT_COV"
bedtools coverage -a "$BLACKLIST" -b "$INPUT_BAM" | \
    awk '{print $1"\t"$2"\t"$3"\t"$4}' | \
    paste - <(bedtools coverage -a "$BLACKLIST" -b "$FILTERED_BAM" | awk '{print $4}') \
    >> "$OUTPUT_COV"

echo "Coverage data saved to: $OUTPUT_COV"
echo ""

# 5. STATISTICHE RIASSUNTIVE
echo "### 5. RIASSUNTO ###"
echo ""

NUM_BL_REGIONS=$(wc -l < "$BLACKLIST")
TOTAL_BL_SIZE=$(awk '{sum+=$3-$2} END {print sum}' "$BLACKLIST")
AVG_BL_SIZE=$(awk '{sum+=$3-$2; count++} END {printf "%.0f", sum/count}' "$BLACKLIST")

echo "Blacklist statistics:"
echo "  Total regions:       $NUM_BL_REGIONS"
echo "  Total size:          $TOTAL_BL_SIZE bp"
echo "  Average region size: $AVG_BL_SIZE bp"
echo ""

echo "Filtering effectiveness:"
if [ "$READS_IN_BL_FILTERED" -eq 0 ]; then
    echo "  ✅ PERFETTO: 100% dei reads in blacklist rimossi"
elif [ "$READS_IN_BL_FILTERED" -lt 100 ]; then
    echo "  ✅ OTTIMO: Solo $READS_IN_BL_FILTERED reads residui in blacklist"
else
    REMOVAL_PCT=$(awk "BEGIN {printf \"%.2f\", (($READS_IN_BL_INPUT - $READS_IN_BL_FILTERED) / $READS_IN_BL_INPUT) * 100}")
    echo "  ⚠️  ATTENZIONE: ${REMOVAL_PCT}% reads blacklist rimossi"
    echo "     (potrebbero esserci reads overlapping parziale)"
fi

echo ""
echo "=== FINE VERIFICA BLACKLIST ==="
