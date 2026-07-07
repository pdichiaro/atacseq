#!/bin/bash
# Script semplice per generare un log della rimozione delle reads in blacklist regions
# Usage: ./blacklist_removal_log.sh INPUT_BAM BLACKLIST OUTPUT_LOG

set -euo pipefail

INPUT_BAM=$1
BLACKLIST=$2
OUTPUT_LOG=${3:-"blacklist_removal.log"}

# Verifica che i file esistano
if [ ! -f "$INPUT_BAM" ]; then
    echo "ERROR: Input BAM not found: $INPUT_BAM"
    exit 1
fi

if [ ! -f "$BLACKLIST" ]; then
    echo "ERROR: Blacklist file not found: $BLACKLIST"
    exit 1
fi

# Estrai nome sample dal path del BAM
SAMPLE=$(basename "$INPUT_BAM" | sed 's/\.bam$//')

# Crea directory output se non esiste
mkdir -p "$(dirname "$OUTPUT_LOG")"

# Calcola statistiche
echo "Calculating blacklist filtering statistics for: $SAMPLE"
echo "Input BAM: $INPUT_BAM"
echo "Blacklist: $BLACKLIST"
echo ""

TOTAL_READS=$(samtools view -c "$INPUT_BAM")
echo "Total reads in input BAM: $(printf "%'d" $TOTAL_READS)"

# Reads in blacklist regions
READS_IN_BL=$(bedtools intersect -a "$INPUT_BAM" -b "$BLACKLIST" -u | samtools view -c -)
echo "Reads in blacklist regions: $(printf "%'d" $READS_IN_BL)"

# Reads after blacklist removal
READS_AFTER_BL=$((TOTAL_READS - READS_IN_BL))
echo "Reads after blacklist removal: $(printf "%'d" $READS_AFTER_BL)"

# Calcola percentuali
PERCENT_REMOVED=$(awk "BEGIN {printf \"%.2f\", ($READS_IN_BL / $TOTAL_READS) * 100}")
PERCENT_RETAINED=$(awk "BEGIN {printf \"%.2f\", ($READS_AFTER_BL / $TOTAL_READS) * 100}")

# Numero regioni blacklist
NUM_BL_REGIONS=$(wc -l < "$BLACKLIST")

# Genera log
cat > "$OUTPUT_LOG" << EOF
========================================================================
BLACKLIST FILTERING LOG - Sample: $SAMPLE
========================================================================

Date: $(date '+%Y-%m-%d %H:%M:%S')
Input BAM: $INPUT_BAM
Blacklist: $BLACKLIST

------------------------------------------------------------------------
STATISTICS
------------------------------------------------------------------------

Total reads in input BAM:           $(printf "%15s" "$(printf "%'d" $TOTAL_READS)")

Reads IN blacklist regions:         $(printf "%15s" "$(printf "%'d" $READS_IN_BL)")  (${PERCENT_REMOVED}%)
Reads OUTSIDE blacklist (retained): $(printf "%15s" "$(printf "%'d" $READS_AFTER_BL)")  (${PERCENT_RETAINED}%)

Number of blacklist regions:        $(printf "%15s" "$(printf "%'d" $NUM_BL_REGIONS)")

------------------------------------------------------------------------
SUMMARY
------------------------------------------------------------------------

Reads removed:     $(printf "%'d" $READS_IN_BL) (${PERCENT_REMOVED}%)
Reads retained:    $(printf "%'d" $READS_AFTER_BL) (${PERCENT_RETAINED}%)

========================================================================
EOF

# Output su stdout e salva su file
cat "$OUTPUT_LOG"

echo ""
echo "✅ Log saved to: $OUTPUT_LOG"
