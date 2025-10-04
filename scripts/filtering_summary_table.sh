#!/bin/bash
# Script per generare una tabella riassuntiva dei filtri post-alignment
# Usage: ./filtering_summary_table.sh SAMPLE OUTDIR [BLACKLIST]

set -euo pipefail

SAMPLE=$1
OUTDIR=${2:-results}
BLACKLIST=${3:-""}

# Percorsi dei file BAM
INPUT_BAM="${OUTDIR}/bowtie2/mergedLibrary/mdup/${SAMPLE}.mLb.mkD.sorted.bam"
FILTERED_BAM="${OUTDIR}/bowtie2/mergedLibrary/Final_BAM/${SAMPLE}.mLb.clN.sorted.bam"

# Verifica che i file esistano
if [ ! -f "$INPUT_BAM" ]; then
    echo "ERROR: Input BAM not found: $INPUT_BAM"
    exit 1
fi

if [ ! -f "$FILTERED_BAM" ]; then
    echo "ERROR: Filtered BAM not found: $FILTERED_BAM"
    exit 1
fi

# Output file
OUTPUT_TABLE="${OUTDIR}/qc/${SAMPLE}_filtering_summary.txt"
mkdir -p "$(dirname $OUTPUT_TABLE)"

# Funzione per formattare numeri con separatori di migliaia
format_number() {
    printf "%'d" "$1" 2>/dev/null || echo "$1"
}

# Funzione per calcolare percentuale
calc_pct() {
    local numerator=$1
    local denominator=${2:-1}
    if [ "$denominator" -eq 0 ]; then
        echo "0.00"
    else
        awk "BEGIN {printf \"%.2f\", ($numerator / $denominator) * 100}"
    fi
}

# Header
echo "======================================================================" > "$OUTPUT_TABLE"
echo "              FILTERING SUMMARY - SAMPLE: $SAMPLE" >> "$OUTPUT_TABLE"
echo "======================================================================" >> "$OUTPUT_TABLE"
echo "" >> "$OUTPUT_TABLE"
echo "Generated: $(date '+%Y-%m-%d %H:%M:%S')" >> "$OUTPUT_TABLE"
echo "" >> "$OUTPUT_TABLE"

# ============================================================================
# 1. METRICHE GENERALI
# ============================================================================
echo "┌────────────────────────────────────────────────────────────────────┐" >> "$OUTPUT_TABLE"
echo "│  1. GENERAL STATISTICS                                             │" >> "$OUTPUT_TABLE"
echo "└────────────────────────────────────────────────────────────────────┘" >> "$OUTPUT_TABLE"
echo "" >> "$OUTPUT_TABLE"

TOTAL_INPUT=$(samtools view -c "$INPUT_BAM")
TOTAL_FILTERED=$(samtools view -c "$FILTERED_BAM")
RETENTION_PCT=$(calc_pct $TOTAL_FILTERED $TOTAL_INPUT)

printf "%-40s %15s %15s %10s\n" "Metric" "Input BAM" "Filtered BAM" "Change" >> "$OUTPUT_TABLE"
printf "%-40s %15s %15s %10s\n" "------" "---------" "------------" "------" >> "$OUTPUT_TABLE"
printf "%-40s %15s %15s %9s%%\n" "Total Reads" "$(format_number $TOTAL_INPUT)" "$(format_number $TOTAL_FILTERED)" "$RETENTION_PCT" >> "$OUTPUT_TABLE"
echo "" >> "$OUTPUT_TABLE"

# Valutazione retention rate
if (( $(echo "$RETENTION_PCT >= 60" | bc -l) )) && (( $(echo "$RETENTION_PCT <= 80" | bc -l) )); then
    VERDICT="✅ OPTIMAL"
elif (( $(echo "$RETENTION_PCT >= 40" | bc -l) )); then
    VERDICT="⚠️  BORDERLINE"
else
    VERDICT="❌ LOW"
fi
printf "%-40s %s\n" "Retention Rate Assessment:" "$VERDICT" >> "$OUTPUT_TABLE"
echo "" >> "$OUTPUT_TABLE"

# ============================================================================
# 2. FILTRI SAM FLAGS
# ============================================================================
echo "┌────────────────────────────────────────────────────────────────────┐" >> "$OUTPUT_TABLE"
echo "│  2. SAM FLAGS FILTERING                                            │" >> "$OUTPUT_TABLE"
echo "└────────────────────────────────────────────────────────────────────┘" >> "$OUTPUT_TABLE"
echo "" >> "$OUTPUT_TABLE"

printf "%-40s %15s %15s %12s\n" "Flag Type" "Input BAM" "Filtered BAM" "Status" >> "$OUTPUT_TABLE"
printf "%-40s %15s %15s %12s\n" "---------" "---------" "------------" "------" >> "$OUTPUT_TABLE"

# Primary alignments
PRIMARY_INPUT=$(samtools view -c -F 0x100 -F 0x800 "$INPUT_BAM")
PRIMARY_FILTERED=$(samtools view -c -F 0x100 -F 0x800 "$FILTERED_BAM")
PRIMARY_PCT_INPUT=$(calc_pct $PRIMARY_INPUT $TOTAL_INPUT)
PRIMARY_PCT_FILTERED=$(calc_pct $PRIMARY_FILTERED $TOTAL_FILTERED)
printf "%-40s %15s %15s %12s\n" "Primary alignments" "$(format_number $PRIMARY_INPUT) (${PRIMARY_PCT_INPUT}%)" "$(format_number $PRIMARY_FILTERED) (${PRIMARY_PCT_FILTERED}%)" "✅" >> "$OUTPUT_TABLE"

# Secondary alignments
SECONDARY_INPUT=$(samtools view -c -f 0x100 "$INPUT_BAM")
SECONDARY_FILTERED=$(samtools view -c -f 0x100 "$FILTERED_BAM")
if [ "$SECONDARY_FILTERED" -eq 0 ]; then
    SECONDARY_STATUS="✅ PASS"
else
    SECONDARY_STATUS="❌ FAIL"
fi
printf "%-40s %15s %15s %12s\n" "Secondary alignments (0x100)" "$(format_number $SECONDARY_INPUT)" "$(format_number $SECONDARY_FILTERED)" "$SECONDARY_STATUS" >> "$OUTPUT_TABLE"

# Supplementary alignments
SUPP_INPUT=$(samtools view -c -f 0x800 "$INPUT_BAM")
SUPP_FILTERED=$(samtools view -c -f 0x800 "$FILTERED_BAM")
if [ "$SUPP_FILTERED" -eq 0 ]; then
    SUPP_STATUS="✅ PASS"
else
    SUPP_STATUS="❌ FAIL"
fi
printf "%-40s %15s %15s %12s\n" "Supplementary alignments (0x800)" "$(format_number $SUPP_INPUT)" "$(format_number $SUPP_FILTERED)" "$SUPP_STATUS" >> "$OUTPUT_TABLE"

# Unmapped reads
UNMAPPED_INPUT=$(samtools view -c -f 0x004 "$INPUT_BAM")
UNMAPPED_FILTERED=$(samtools view -c -f 0x004 "$FILTERED_BAM")
if [ "$UNMAPPED_FILTERED" -eq 0 ]; then
    UNMAPPED_STATUS="✅ PASS"
else
    UNMAPPED_STATUS="❌ FAIL"
fi
printf "%-40s %15s %15s %12s\n" "Unmapped reads (0x004)" "$(format_number $UNMAPPED_INPUT)" "$(format_number $UNMAPPED_FILTERED)" "$UNMAPPED_STATUS" >> "$OUTPUT_TABLE"

# Duplicates
DUPS_INPUT=$(samtools view -c -f 0x400 "$INPUT_BAM")
DUPS_FILTERED=$(samtools view -c -f 0x400 "$FILTERED_BAM")
DUPS_PCT_INPUT=$(calc_pct $DUPS_INPUT $TOTAL_INPUT)
if [ "$DUPS_FILTERED" -eq 0 ]; then
    DUPS_STATUS="✅ REMOVED"
else
    DUPS_STATUS="⚠️  KEPT"
fi
printf "%-40s %15s %15s %12s\n" "Duplicates (0x400)" "$(format_number $DUPS_INPUT) (${DUPS_PCT_INPUT}%)" "$(format_number $DUPS_FILTERED)" "$DUPS_STATUS" >> "$OUTPUT_TABLE"

echo "" >> "$OUTPUT_TABLE"

# ============================================================================
# 3. PAIRED-END SPECIFIC (se applicabile)
# ============================================================================
PAIRED_INPUT=$(samtools view -c -f 0x001 "$INPUT_BAM")
if [ "$PAIRED_INPUT" -gt 0 ]; then
    echo "┌────────────────────────────────────────────────────────────────────┐" >> "$OUTPUT_TABLE"
    echo "│  3. PAIRED-END METRICS                                             │" >> "$OUTPUT_TABLE"
    echo "└────────────────────────────────────────────────────────────────────┘" >> "$OUTPUT_TABLE"
    echo "" >> "$OUTPUT_TABLE"
    
    printf "%-40s %15s %15s %12s\n" "Metric" "Input BAM" "Filtered BAM" "Status" >> "$OUTPUT_TABLE"
    printf "%-40s %15s %15s %12s\n" "------" "---------" "------------" "------" >> "$OUTPUT_TABLE"
    
    # Properly paired
    PROPER_INPUT=$(samtools view -c -f 0x002 "$INPUT_BAM")
    PROPER_FILTERED=$(samtools view -c -f 0x002 "$FILTERED_BAM")
    PROPER_PCT_INPUT=$(calc_pct $PROPER_INPUT $PAIRED_INPUT)
    PROPER_PCT_FILTERED=$(calc_pct $PROPER_FILTERED $TOTAL_FILTERED)
    
    if [ "$PROPER_PCT_FILTERED" == "100.00" ] || (( $(echo "$PROPER_PCT_FILTERED >= 99.5" | bc -l) )); then
        PROPER_STATUS="✅ PASS"
    else
        PROPER_STATUS="⚠️  CHECK"
    fi
    printf "%-40s %15s %15s %12s\n" "Properly paired (0x002)" "$(format_number $PROPER_INPUT) (${PROPER_PCT_INPUT}%)" "$(format_number $PROPER_FILTERED) (${PROPER_PCT_FILTERED}%)" "$PROPER_STATUS" >> "$OUTPUT_TABLE"
    
    # Singletons
    SINGLETON_INPUT=$(samtools view -c -f 0x008 "$INPUT_BAM")
    SINGLETON_FILTERED=$(samtools view -c -f 0x008 "$FILTERED_BAM")
    if [ "$SINGLETON_FILTERED" -eq 0 ]; then
        SINGLETON_STATUS="✅ REMOVED"
    else
        SINGLETON_STATUS="⚠️  PRESENT"
    fi
    printf "%-40s %15s %15s %12s\n" "Singletons (mate unmapped)" "$(format_number $SINGLETON_INPUT)" "$(format_number $SINGLETON_FILTERED)" "$SINGLETON_STATUS" >> "$OUTPUT_TABLE"
    
    echo "" >> "$OUTPUT_TABLE"
fi

# ============================================================================
# 4. MAPQ FILTERING
# ============================================================================
echo "┌────────────────────────────────────────────────────────────────────┐" >> "$OUTPUT_TABLE"
echo "│  4. MAPPING QUALITY (MAPQ)                                         │" >> "$OUTPUT_TABLE"
echo "└────────────────────────────────────────────────────────────────────┘" >> "$OUTPUT_TABLE"
echo "" >> "$OUTPUT_TABLE"

printf "%-40s %15s %15s %12s\n" "MAPQ Category" "Input BAM" "Filtered BAM" "Status" >> "$OUTPUT_TABLE"
printf "%-40s %15s %15s %12s\n" "-------------" "---------" "------------" "------" >> "$OUTPUT_TABLE"

# Multi-mappers (MAPQ=0)
MAPQ0_INPUT=$(samtools view "$INPUT_BAM" | awk '$5==0' | wc -l)
MAPQ0_FILTERED=$(samtools view "$FILTERED_BAM" | awk '$5==0' | wc -l)
MAPQ0_PCT_INPUT=$(calc_pct $MAPQ0_INPUT $TOTAL_INPUT)

if [ "$MAPQ0_FILTERED" -eq 0 ]; then
    MAPQ0_STATUS="✅ REMOVED"
else
    MAPQ0_PCT_FILTERED=$(calc_pct $MAPQ0_FILTERED $TOTAL_FILTERED)
    MAPQ0_STATUS="⚠️  KEPT (${MAPQ0_PCT_FILTERED}%)"
fi
printf "%-40s %15s %15s %12s\n" "Multi-mappers (MAPQ=0)" "$(format_number $MAPQ0_INPUT) (${MAPQ0_PCT_INPUT}%)" "$(format_number $MAPQ0_FILTERED)" "$MAPQ0_STATUS" >> "$OUTPUT_TABLE"

# High quality (MAPQ>=20)
MAPQ20_INPUT=$(samtools view "$INPUT_BAM" | awk '$5>=20' | wc -l)
MAPQ20_FILTERED=$(samtools view "$FILTERED_BAM" | awk '$5>=20' | wc -l)
MAPQ20_PCT_INPUT=$(calc_pct $MAPQ20_INPUT $TOTAL_INPUT)
MAPQ20_PCT_FILTERED=$(calc_pct $MAPQ20_FILTERED $TOTAL_FILTERED)
printf "%-40s %15s %15s %12s\n" "High quality (MAPQ≥20)" "$(format_number $MAPQ20_INPUT) (${MAPQ20_PCT_INPUT}%)" "$(format_number $MAPQ20_FILTERED) (${MAPQ20_PCT_FILTERED}%)" "✅" >> "$OUTPUT_TABLE"

echo "" >> "$OUTPUT_TABLE"

# ============================================================================
# 5. INSERT SIZE FILTERING (PE only)
# ============================================================================
if [ "$PAIRED_INPUT" -gt 0 ]; then
    echo "┌────────────────────────────────────────────────────────────────────┐" >> "$OUTPUT_TABLE"
    echo "│  5. INSERT SIZE DISTRIBUTION (PE)                                 │" >> "$OUTPUT_TABLE"
    echo "└────────────────────────────────────────────────────────────────────┘" >> "$OUTPUT_TABLE"
    echo "" >> "$OUTPUT_TABLE"
    
    printf "%-40s %15s %15s\n" "Insert Size Range" "Input BAM" "Filtered BAM" >> "$OUTPUT_TABLE"
    printf "%-40s %15s %15s\n" "-----------------" "---------" "------------" >> "$OUTPUT_TABLE"
    
    # Calcola distribuzione insert size
    INSERT_INPUT_0_100=$(samtools view -f 0x002 "$INPUT_BAM" | awk '{if($9>0 && $9<=100) count++} END {print count+0}')
    INSERT_INPUT_100_200=$(samtools view -f 0x002 "$INPUT_BAM" | awk '{if($9>100 && $9<=200) count++} END {print count+0}')
    INSERT_INPUT_200_300=$(samtools view -f 0x002 "$INPUT_BAM" | awk '{if($9>200 && $9<=300) count++} END {print count+0}')
    INSERT_INPUT_300_500=$(samtools view -f 0x002 "$INPUT_BAM" | awk '{if($9>300 && $9<=500) count++} END {print count+0}')
    INSERT_INPUT_500_PLUS=$(samtools view -f 0x002 "$INPUT_BAM" | awk '{if($9>500) count++} END {print count+0}')
    
    INSERT_FILT_0_100=$(samtools view -f 0x002 "$FILTERED_BAM" | awk '{if($9>0 && $9<=100) count++} END {print count+0}')
    INSERT_FILT_100_200=$(samtools view -f 0x002 "$FILTERED_BAM" | awk '{if($9>100 && $9<=200) count++} END {print count+0}')
    INSERT_FILT_200_300=$(samtools view -f 0x002 "$FILTERED_BAM" | awk '{if($9>200 && $9<=300) count++} END {print count+0}')
    INSERT_FILT_300_500=$(samtools view -f 0x002 "$FILTERED_BAM" | awk '{if($9>300 && $9<=500) count++} END {print count+0}')
    INSERT_FILT_500_PLUS=$(samtools view -f 0x002 "$FILTERED_BAM" | awk '{if($9>500) count++} END {print count+0}')
    
    printf "%-40s %15s %15s\n" "0-100 bp" "$(format_number $INSERT_INPUT_0_100)" "$(format_number $INSERT_FILT_0_100)" >> "$OUTPUT_TABLE"
    printf "%-40s %15s %15s\n" "100-200 bp" "$(format_number $INSERT_INPUT_100_200)" "$(format_number $INSERT_FILT_100_200)" >> "$OUTPUT_TABLE"
    printf "%-40s %15s %15s\n" "200-300 bp" "$(format_number $INSERT_INPUT_200_300)" "$(format_number $INSERT_FILT_200_300)" >> "$OUTPUT_TABLE"
    printf "%-40s %15s %15s\n" "300-500 bp" "$(format_number $INSERT_INPUT_300_500)" "$(format_number $INSERT_FILT_300_500)" >> "$OUTPUT_TABLE"
    
    if [ "$INSERT_FILT_500_PLUS" -eq 0 ]; then
        INSERTSIZE_STATUS="✅ REMOVED"
    else
        INSERTSIZE_STATUS="❌ PRESENT"
    fi
    printf "%-40s %15s %15s %12s\n" ">500 bp (should be removed)" "$(format_number $INSERT_INPUT_500_PLUS)" "$(format_number $INSERT_FILT_500_PLUS)" "$INSERTSIZE_STATUS" >> "$OUTPUT_TABLE"
    
    # Media e mediana insert size
    MEDIAN_INPUT=$(samtools view -f 0x002 "$INPUT_BAM" | awk '{if($9>0) print $9}' | sort -n | awk '{a[NR]=$1} END {print a[int(NR/2)]}')
    MEDIAN_FILTERED=$(samtools view -f 0x002 "$FILTERED_BAM" | awk '{if($9>0) print $9}' | sort -n | awk '{a[NR]=$1} END {print a[int(NR/2)]}')
    
    echo "" >> "$OUTPUT_TABLE"
    printf "%-40s %15s %15s\n" "Median insert size:" "${MEDIAN_INPUT} bp" "${MEDIAN_FILTERED} bp" >> "$OUTPUT_TABLE"
    echo "" >> "$OUTPUT_TABLE"
fi

# ============================================================================
# 6. BLACKLIST FILTERING
# ============================================================================
if [ -n "$BLACKLIST" ] && [ -f "$BLACKLIST" ]; then
    echo "┌────────────────────────────────────────────────────────────────────┐" >> "$OUTPUT_TABLE"
    echo "│  6. BLACKLIST FILTERING                                            │" >> "$OUTPUT_TABLE"
    echo "└────────────────────────────────────────────────────────────────────┘" >> "$OUTPUT_TABLE"
    echo "" >> "$OUTPUT_TABLE"
    
    READS_BL_INPUT=$(bedtools intersect -a "$INPUT_BAM" -b "$BLACKLIST" -u | samtools view -c -)
    READS_BL_FILTERED=$(bedtools intersect -a "$FILTERED_BAM" -b "$BLACKLIST" -u | samtools view -c -)
    
    BL_PCT_INPUT=$(calc_pct $READS_BL_INPUT $TOTAL_INPUT)
    BL_PCT_FILTERED=$(calc_pct $READS_BL_FILTERED $TOTAL_FILTERED)
    
    printf "%-40s %15s %15s %12s\n" "Metric" "Input BAM" "Filtered BAM" "Status" >> "$OUTPUT_TABLE"
    printf "%-40s %15s %15s %12s\n" "------" "---------" "------------" "------" >> "$OUTPUT_TABLE"
    
    if [ "$READS_BL_FILTERED" -eq 0 ]; then
        BL_STATUS="✅ PERFECT"
    elif (( $(echo "$BL_PCT_FILTERED < 0.1" | bc -l) )); then
        BL_STATUS="✅ GOOD"
    else
        BL_STATUS="⚠️  CHECK"
    fi
    
    printf "%-40s %15s %15s %12s\n" "Reads in blacklist regions" "$(format_number $READS_BL_INPUT) (${BL_PCT_INPUT}%)" "$(format_number $READS_BL_FILTERED) (${BL_PCT_FILTERED}%)" "$BL_STATUS" >> "$OUTPUT_TABLE"
    
    NUM_BL_REGIONS=$(wc -l < "$BLACKLIST")
    printf "%-40s %15s\n" "Total blacklist regions:" "$(format_number $NUM_BL_REGIONS)" >> "$OUTPUT_TABLE"
    
    echo "" >> "$OUTPUT_TABLE"
fi

# ============================================================================
# 7. CHROMOSOME COVERAGE (Top 5)
# ============================================================================
echo "┌────────────────────────────────────────────────────────────────────┐" >> "$OUTPUT_TABLE"
echo "│  7. CHROMOSOME COVERAGE (Top 5)                                    │" >> "$OUTPUT_TABLE"
echo "└────────────────────────────────────────────────────────────────────┘" >> "$OUTPUT_TABLE"
echo "" >> "$OUTPUT_TABLE"

printf "%-20s %18s %18s\n" "Chromosome" "Input BAM" "Filtered BAM" >> "$OUTPUT_TABLE"
printf "%-20s %18s %18s\n" "----------" "---------" "------------" >> "$OUTPUT_TABLE"

samtools idxstats "$INPUT_BAM" | sort -k3 -rn | head -5 | while read chr len count rest; do
    if [ "$chr" != "*" ]; then
        FILT_COUNT=$(samtools idxstats "$FILTERED_BAM" | awk -v chr="$chr" '$1==chr {print $3}')
        printf "%-20s %18s %18s\n" "$chr" "$(format_number $count)" "$(format_number $FILT_COUNT)" >> "$OUTPUT_TABLE"
    fi
done

echo "" >> "$OUTPUT_TABLE"

# ============================================================================
# 8. FINAL VERDICT
# ============================================================================
echo "┌────────────────────────────────────────────────────────────────────┐" >> "$OUTPUT_TABLE"
echo "│  8. FINAL VERDICT                                                  │" >> "$OUTPUT_TABLE"
echo "└────────────────────────────────────────────────────────────────────┘" >> "$OUTPUT_TABLE"
echo "" >> "$OUTPUT_TABLE"

PASS_COUNT=0
WARN_COUNT=0
FAIL_COUNT=0

# Check 1: Retention rate
if (( $(echo "$RETENTION_PCT >= 60" | bc -l) )); then
    echo "✅ Retention rate: ${RETENTION_PCT}% (PASS)" >> "$OUTPUT_TABLE"
    ((PASS_COUNT++))
elif (( $(echo "$RETENTION_PCT >= 40" | bc -l) )); then
    echo "⚠️  Retention rate: ${RETENTION_PCT}% (BORDERLINE)" >> "$OUTPUT_TABLE"
    ((WARN_COUNT++))
else
    echo "❌ Retention rate: ${RETENTION_PCT}% (LOW)" >> "$OUTPUT_TABLE"
    ((FAIL_COUNT++))
fi

# Check 2: Secondary alignments
if [ "$SECONDARY_FILTERED" -eq 0 ]; then
    echo "✅ Secondary alignments removed (PASS)" >> "$OUTPUT_TABLE"
    ((PASS_COUNT++))
else
    echo "❌ Secondary alignments present: $SECONDARY_FILTERED (FAIL)" >> "$OUTPUT_TABLE"
    ((FAIL_COUNT++))
fi

# Check 3: Supplementary alignments
if [ "$SUPP_FILTERED" -eq 0 ]; then
    echo "✅ Supplementary alignments removed (PASS)" >> "$OUTPUT_TABLE"
    ((PASS_COUNT++))
else
    echo "❌ Supplementary alignments present: $SUPP_FILTERED (FAIL)" >> "$OUTPUT_TABLE"
    ((FAIL_COUNT++))
fi

# Check 4: Multi-mappers (if keep_multi_map=false expected)
if [ "$MAPQ0_FILTERED" -eq 0 ]; then
    echo "✅ Multi-mappers (MAPQ=0) removed (PASS)" >> "$OUTPUT_TABLE"
    ((PASS_COUNT++))
else
    MAPQ0_PCT=$(calc_pct $MAPQ0_FILTERED $TOTAL_FILTERED)
    echo "⚠️  Multi-mappers present: $MAPQ0_FILTERED (${MAPQ0_PCT}%) - keep_multi_map=true?" >> "$OUTPUT_TABLE"
    ((WARN_COUNT++))
fi

# Check 5: Insert size (PE only)
if [ "$PAIRED_INPUT" -gt 0 ]; then
    if [ "$INSERT_FILT_500_PLUS" -eq 0 ]; then
        echo "✅ Fragments >500bp removed (PASS)" >> "$OUTPUT_TABLE"
        ((PASS_COUNT++))
    else
        echo "❌ Fragments >500bp present: $INSERT_FILT_500_PLUS (FAIL)" >> "$OUTPUT_TABLE"
        ((FAIL_COUNT++))
    fi
fi

# Check 6: Blacklist (if provided)
if [ -n "$BLACKLIST" ] && [ -f "$BLACKLIST" ]; then
    if [ "$READS_BL_FILTERED" -eq 0 ]; then
        echo "✅ Blacklist regions completely filtered (PERFECT)" >> "$OUTPUT_TABLE"
        ((PASS_COUNT++))
    elif (( $(echo "$BL_PCT_FILTERED < 0.1" | bc -l) )); then
        echo "✅ Blacklist filtering effective: ${BL_PCT_FILTERED}% residual (GOOD)" >> "$OUTPUT_TABLE"
        ((PASS_COUNT++))
    else
        echo "⚠️  Blacklist filtering partial: ${BL_PCT_FILTERED}% residual (CHECK)" >> "$OUTPUT_TABLE"
        ((WARN_COUNT++))
    fi
fi

echo "" >> "$OUTPUT_TABLE"
echo "────────────────────────────────────────────────────────────────────" >> "$OUTPUT_TABLE"

if [ "$FAIL_COUNT" -eq 0 ] && [ "$WARN_COUNT" -eq 0 ]; then
    echo "🎉 OVERALL VERDICT: ✅ EXCELLENT - All filters applied correctly" >> "$OUTPUT_TABLE"
elif [ "$FAIL_COUNT" -eq 0 ]; then
    echo "OVERALL VERDICT: ✅ PASS with $WARN_COUNT warnings" >> "$OUTPUT_TABLE"
else
    echo "OVERALL VERDICT: ❌ FAILED - $FAIL_COUNT critical issues detected" >> "$OUTPUT_TABLE"
fi

echo "────────────────────────────────────────────────────────────────────" >> "$OUTPUT_TABLE"
echo "" >> "$OUTPUT_TABLE"

echo "Summary: ${PASS_COUNT} checks passed, ${WARN_COUNT} warnings, ${FAIL_COUNT} failures" >> "$OUTPUT_TABLE"
echo "" >> "$OUTPUT_TABLE"
echo "======================================================================" >> "$OUTPUT_TABLE"

# Output su stdout e file
cat "$OUTPUT_TABLE"
echo ""
echo "📊 Summary table saved to: $OUTPUT_TABLE"
