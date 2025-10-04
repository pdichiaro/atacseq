#!/usr/bin/env python3
"""
Visualizza statistiche di filtering post-alignment
"""

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import sys
from pathlib import Path

def parse_flagstat(flagstat_file):
    """Estrae metriche da samtools flagstat output"""
    metrics = {}
    with open(flagstat_file) as f:
        for line in f:
            if 'in total' in line:
                metrics['total'] = int(line.split()[0])
            elif 'primary' in line:
                metrics['primary'] = int(line.split()[0])
            elif 'mapped (' in line and 'primary' not in line:
                metrics['mapped'] = int(line.split()[0])
            elif 'properly paired' in line:
                metrics['properly_paired'] = int(line.split()[0])
            elif 'singletons' in line:
                metrics['singletons'] = int(line.split()[0])
            elif 'duplicates' in line:
                metrics['duplicates'] = int(line.split()[0])
    return metrics

def plot_filtering_comparison(input_flagstat, filtered_flagstat, output_dir):
    """Genera grafici comparativi"""
    
    # Parse flagstat files
    input_metrics = parse_flagstat(input_flagstat)
    filtered_metrics = parse_flagstat(filtered_flagstat)
    
    # Crea figure
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle('Filtering Statistics Comparison', fontsize=16, fontweight='bold')
    
    # 1. BAR PLOT: Total reads comparison
    ax1 = axes[0, 0]
    categories = ['Total', 'Primary', 'Mapped', 'Proper Pairs']
    input_values = [
        input_metrics.get('total', 0),
        input_metrics.get('primary', 0),
        input_metrics.get('mapped', 0),
        input_metrics.get('properly_paired', 0)
    ]
    filtered_values = [
        filtered_metrics.get('total', 0),
        filtered_metrics.get('primary', 0),
        filtered_metrics.get('mapped', 0),
        filtered_metrics.get('properly_paired', 0)
    ]
    
    x = range(len(categories))
    width = 0.35
    ax1.bar([i - width/2 for i in x], input_values, width, label='Input BAM', alpha=0.8)
    ax1.bar([i + width/2 for i in x], filtered_values, width, label='Filtered BAM', alpha=0.8)
    ax1.set_ylabel('Number of Reads')
    ax1.set_title('Read Counts by Category')
    ax1.set_xticks(x)
    ax1.set_xticklabels(categories, rotation=45, ha='right')
    ax1.legend()
    ax1.grid(axis='y', alpha=0.3)
    
    # 2. PIE CHART: Retention rate
    ax2 = axes[0, 1]
    total_input = input_metrics.get('total', 1)
    total_filtered = filtered_metrics.get('total', 0)
    removed = total_input - total_filtered
    retention_pct = (total_filtered / total_input) * 100 if total_input > 0 else 0
    
    sizes = [total_filtered, removed]
    labels = [f'Retained\n{retention_pct:.1f}%', f'Removed\n{100-retention_pct:.1f}%']
    colors = ['#2ecc71', '#e74c3c']
    ax2.pie(sizes, labels=labels, colors=colors, autopct='%1.1f%%', startangle=90)
    ax2.set_title(f'Overall Retention Rate\n({total_filtered:,} / {total_input:,} reads)')
    
    # 3. BAR PLOT: Percentage comparison
    ax3 = axes[1, 0]
    input_total = input_metrics.get('total', 1)
    filtered_total = filtered_metrics.get('total', 1)
    
    pct_categories = ['Primary', 'Mapped', 'Proper Pairs']
    input_pcts = [
        (input_metrics.get('primary', 0) / input_total) * 100,
        (input_metrics.get('mapped', 0) / input_total) * 100,
        (input_metrics.get('properly_paired', 0) / input_total) * 100
    ]
    filtered_pcts = [
        (filtered_metrics.get('primary', 0) / filtered_total) * 100,
        (filtered_metrics.get('mapped', 0) / filtered_total) * 100,
        (filtered_metrics.get('properly_paired', 0) / filtered_total) * 100
    ]
    
    x = range(len(pct_categories))
    ax3.bar([i - width/2 for i in x], input_pcts, width, label='Input BAM', alpha=0.8)
    ax3.bar([i + width/2 for i in x], filtered_pcts, width, label='Filtered BAM', alpha=0.8)
    ax3.set_ylabel('Percentage (%)')
    ax3.set_title('Quality Metrics (%)')
    ax3.set_xticks(x)
    ax3.set_xticklabels(pct_categories, rotation=45, ha='right')
    ax3.set_ylim([0, 105])
    ax3.legend()
    ax3.grid(axis='y', alpha=0.3)
    
    # 4. TABLE: Summary statistics
    ax4 = axes[1, 1]
    ax4.axis('tight')
    ax4.axis('off')
    
    summary_data = [
        ['Metric', 'Input BAM', 'Filtered BAM', 'Change'],
        ['Total reads', f"{input_metrics.get('total', 0):,}", 
         f"{filtered_metrics.get('total', 0):,}", 
         f"{retention_pct:.1f}%"],
        ['Primary', f"{input_metrics.get('primary', 0):,}", 
         f"{filtered_metrics.get('primary', 0):,}", 
         f"{(filtered_metrics.get('primary', 0)/input_metrics.get('primary', 1)*100):.1f}%"],
        ['Proper pairs', f"{input_metrics.get('properly_paired', 0):,}", 
         f"{filtered_metrics.get('properly_paired', 0):,}", 
         f"{(filtered_metrics.get('properly_paired', 0)/input_metrics.get('properly_paired', 1)*100):.1f}%"],
        ['Singletons', f"{input_metrics.get('singletons', 0):,}", 
         f"{filtered_metrics.get('singletons', 0):,}", 
         '0.0%' if filtered_metrics.get('singletons', 0) == 0 else f"{(filtered_metrics.get('singletons', 0)/input_metrics.get('singletons', 1)*100):.1f}%"]
    ]
    
    table = ax4.table(cellText=summary_data, loc='center', cellLoc='right')
    table.auto_set_font_size(False)
    table.set_fontsize(9)
    table.scale(1, 2)
    
    # Stile header
    for i in range(4):
        table[(0, i)].set_facecolor('#3498db')
        table[(0, i)].set_text_props(weight='bold', color='white')
    
    ax4.set_title('Summary Statistics', fontweight='bold', pad=20)
    
    plt.tight_layout()
    
    # Salva
    output_file = Path(output_dir) / 'filtering_comparison.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Plot saved to: {output_file}")
    
    return fig

def plot_insert_size_distribution(insert_size_file, output_dir):
    """Plot insert size distribution"""
    
    # Leggi dati
    df = pd.read_csv(insert_size_file, sep='\t')
    
    fig, ax = plt.subplots(figsize=(10, 6))
    
    ax.hist(df['insert_size'], bins=50, alpha=0.7, edgecolor='black')
    ax.axvline(500, color='red', linestyle='--', linewidth=2, label='Filter threshold (500bp)')
    ax.set_xlabel('Insert Size (bp)')
    ax.set_ylabel('Frequency')
    ax.set_title('Insert Size Distribution')
    ax.legend()
    ax.grid(axis='y', alpha=0.3)
    
    plt.tight_layout()
    
    output_file = Path(output_dir) / 'insert_size_distribution.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Insert size plot saved to: {output_file}")

if __name__ == '__main__':
    if len(sys.argv) < 4:
        print("Usage: plot_filtering_stats.py <input.flagstat> <filtered.flagstat> <output_dir>")
        sys.exit(1)
    
    input_flagstat = sys.argv[1]
    filtered_flagstat = sys.argv[2]
    output_dir = sys.argv[3]
    
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    
    plot_filtering_comparison(input_flagstat, filtered_flagstat, output_dir)
    print("✅ Plots generated successfully!")
