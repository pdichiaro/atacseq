#!/bin/bash
# Script per verificare l'efficacia dei filtri post-alignment

SAMPLE=$1
OUTDIR=${2:-results}

echo "=== VERIFICA FILTRI POST-ALIGNMENT ==="
echo "Sample: $SAMPLE"
echo ""

# Percorsi dei file BAM
INPUT_BAM="${OUTDIR}/bowtie2/mergedLibrary/mdup/${SAMPLE}.mLb.mkD.sorted.bam"
FILTERED_BAM="${OUTDIR}/bowtie2/mergedLibrary/Final_BAM/${SAMPLE}.mLb.clN.sorted.bam"

# 1. FLAGSTAT COMPARISON
echo "### 1. CONFRONTO FLAGSTAT (Prima vs Dopo Filtri) ###"
echo ""

echo "--- INPUT BAM (dopo mark duplicates) ---"
samtools flagstat "$INPUT_BAM"
echo ""

echo "--- FILTERED BAM (dopo tutti i filtri) ---"
samtools flagstat "$FILTERED_BAM"
echo ""

# 2. CALCOLO PERCENTUALI DI RITENZIONE
echo "### 2. PERCENTUALI DI RITENZIONE ###"
echo ""

TOTAL_INPUT=$(samtools view -c "$INPUT_BAM")
TOTAL_FILTERED=$(samtools view -c "$FILTERED_BAM")
RETENTION=$(awk "BEGIN {printf \"%.2f\", ($TOTAL_FILTERED / $TOTAL_INPUT) * 100}")

echo "Reads totali INPUT:    $TOTAL_INPUT"
echo "Reads totali FILTERED: $TOTAL_FILTERED"
echo "Retention rate:        ${RETENTION}%"
echo ""

# 3. DISTRIBUZIONE MAPQ
echo "### 3. DISTRIBUZIONE MAPQ ###"
echo ""

echo "--- INPUT BAM ---"
samtools view "$INPUT_BAM" | \
    awk '{print $5}' | \
    sort -n | \
    uniq -c | \
    sort -rn | \
    head -10
echo ""

echo "--- FILTERED BAM ---"
samtools view "$FILTERED_BAM" | \
    awk '{print $5}' | \
    sort -n | \
    uniq -c | \
    sort -rn | \
    head -10
echo ""

# 4. DISTRIBUZIONE INSERT SIZE (PE only)
echo "### 4. DISTRIBUZIONE INSERT SIZE (PE) ###"
echo ""

echo "--- INPUT BAM ---"
samtools view -f 0x002 "$INPUT_BAM" | \
    awk '{if($9>0) print $9}' | \
    awk '{
        if($1<=100) bins[1]++
        else if($1<=200) bins[2]++
        else if($1<=300) bins[3]++
        else if($1<=400) bins[4]++
        else if($1<=500) bins[5]++
        else if($1<=600) bins[6]++
        else if($1<=800) bins[7]++
        else if($1<=1000) bins[8]++
        else bins[9]++
    }
    END {
        print "0-100bp:   ", bins[1]+0
        print "100-200bp: ", bins[2]+0
        print "200-300bp: ", bins[3]+0
        print "300-400bp: ", bins[4]+0
        print "400-500bp: ", bins[5]+0
        print "500-600bp: ", bins[6]+0
        print "600-800bp: ", bins[7]+0
        print "800-1000bp:", bins[8]+0
        print ">1000bp:   ", bins[9]+0
    }'
echo ""

echo "--- FILTERED BAM ---"
samtools view -f 0x002 "$FILTERED_BAM" | \
    awk '{if($9>0) print $9}' | \
    awk '{
        if($1<=100) bins[1]++
        else if($1<=200) bins[2]++
        else if($1<=300) bins[3]++
        else if($1<=400) bins[4]++
        else if($1<=500) bins[5]++
        else if($1<=600) bins[6]++
        else if($1<=800) bins[7]++
        else if($1<=1000) bins[8]++
        else bins[9]++
    }
    END {
        print "0-100bp:   ", bins[1]+0
        print "100-200bp: ", bins[2]+0
        print "200-300bp: ", bins[3]+0
        print "300-400bp: ", bins[4]+0
        print "400-500bp: ", bins[5]+0
        print "500-600bp: ", bins[6]+0
        print "600-800bp: ", bins[7]+0
        print "800-1000bp:", bins[8]+0
        print ">1000bp:   ", bins[9]+0
    }'
echo ""

# 5. COVERAGE PER CROMOSOMA
echo "### 5. COVERAGE PER CROMOSOMA (Top 10) ###"
echo ""

echo "--- INPUT BAM ---"
samtools idxstats "$INPUT_BAM" | \
    sort -k3 -rn | \
    head -10 | \
    awk '{printf "%-15s %10d reads\n", $1, $3}'
echo ""

echo "--- FILTERED BAM ---"
samtools idxstats "$FILTERED_BAM" | \
    sort -k3 -rn | \
    head -10 | \
    awk '{printf "%-15s %10d reads\n", $1, $3}'
echo ""

# 6. VERIFICHE SPECIFICHE FLAGS
echo "### 6. VERIFICHE FLAGS SPECIFICI ###"
echo ""

echo "Secondary alignments:"
echo "  INPUT:    $(samtools view -c -f 0x100 $INPUT_BAM)"
echo "  FILTERED: $(samtools view -c -f 0x100 $FILTERED_BAM)"
echo ""

echo "Supplementary alignments:"
echo "  INPUT:    $(samtools view -c -f 0x800 $INPUT_BAM)"
echo "  FILTERED: $(samtools view -c -f 0x800 $FILTERED_BAM)"
echo ""

echo "Unmapped reads:"
echo "  INPUT:    $(samtools view -c -f 0x004 $INPUT_BAM)"
echo "  FILTERED: $(samtools view -c -f 0x004 $FILTERED_BAM)"
echo ""

echo "Duplicates:"
echo "  INPUT:    $(samtools view -c -f 0x400 $INPUT_BAM)"
echo "  FILTERED: $(samtools view -c -f 0x400 $FILTERED_BAM)"
echo ""

echo "Multi-mappers (MAPQ=0):"
echo "  INPUT:    $(samtools view $INPUT_BAM | awk '$5==0' | wc -l)"
echo "  FILTERED: $(samtools view $FILTERED_BAM | awk '$5==0' | wc -l)"
echo ""

echo "=== FINE VERIFICA ==="
