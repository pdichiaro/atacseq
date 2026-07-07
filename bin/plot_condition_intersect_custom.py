#!/usr/bin/env python3

"""
Generate intersection file for condition-level UpSet plots.
This script creates a file compatible with plot_peak_intersect.r showing
overlaps between different conditions for the same antibody.
"""

import sys
import argparse
from collections import defaultdict

def main():
    parser = argparse.ArgumentParser(description='Generate condition intersection data for UpSet plot')
    parser.add_argument('merged_file', help='Merged peaks file with condition tags')
    parser.add_argument('condition_names', help='Comma-separated list of condition names')
    parser.add_argument('output_file', help='Output intersect file')
    args = parser.parse_args()
    
    conditions = args.condition_names.split(',')
    # Sort conditions by length (descending) to match longer names first
    # This prevents "WT" from matching "WT_TREATMENT_Interval_1" before "WT_TREATMENT" can
    conditions_sorted = sorted(conditions, key=len, reverse=True)
    condition_set = set(conditions)
    
    # Dictionary to count peak intersections
    # Key: frozenset of conditions, Value: count
    intersect_counts = defaultdict(int)
    
    # Debug output
    print(f"DEBUG: Processing {len(conditions)} conditions: {conditions}", file=sys.stderr)
    print(f"DEBUG: Sorted for matching: {conditions_sorted}", file=sys.stderr)
    
    # Track which conditions we actually find
    found_conditions = set()
    
    # Read merged peaks file
    line_num = 0
    with open(args.merged_file, 'r') as f:
        for line in f:
            line_num += 1
            fields = line.strip().split('\t')
            if len(fields) < 6:
                continue
                
            # Column 5 (index 5) contains comma-separated peak names
            # Format: CONDITION_Interval_N,CONDITION2_Interval_M,...
            peak_names = fields[5].split(',')
            
            # Debug first few lines
            if line_num <= 3:
                print(f"DEBUG Line {line_num}: peak_names = {peak_names[:3]}...", file=sys.stderr)
            
            # Extract condition names from peak names
            peak_conditions = set()
            for peak_name in peak_names:
                # Peak name format: CONDITION_Interval_N or CONDITION_...
                # Extract the first part before _Interval or before _peaks
                parts = peak_name.split('_')
                if len(parts) > 0:
                    # Try to match against known conditions (longest first to avoid partial matches)
                    matched = False
                    for cond in conditions_sorted:
                        if peak_name.startswith(cond + '_'):
                            peak_conditions.add(cond)
                            matched = True
                            break
                    
                    # Debug unmatched peak names
                    if not matched and line_num <= 10:
                        print(f"DEBUG: Could not match peak_name='{peak_name}' against conditions", file=sys.stderr)
            
            # Only count if we found valid conditions
            if peak_conditions:
                intersect_counts[frozenset(peak_conditions)] += 1
                # Track which conditions we found
                found_conditions.update(peak_conditions)
    
    # Write output in UpSetR compatible format
    # Format: condition1&condition2&condition3  count
    with open(args.output_file, 'w') as f:
        # Sort by count (descending), then by number of conditions, then alphabetically
        sorted_intersects = sorted(
            intersect_counts.items(),
            key=lambda x: (-x[1], -len(x[0]), '&'.join(sorted(x[0])))
        )
        
        for cond_set, count in sorted_intersects:
            cond_list = sorted(cond_set)
            f.write(f"{'&'.join(cond_list)}\t{count}\n")
    
    print(f"Generated intersection file with {len(intersect_counts)} unique combinations")
    print(f"Expected conditions: {sorted(conditions)}")
    print(f"Found conditions in data: {sorted(found_conditions)}")
    
    # Check for missing conditions
    missing_conditions = set(conditions) - found_conditions
    if missing_conditions:
        print(f"WARNING: Missing conditions in data: {sorted(missing_conditions)}", file=sys.stderr)
    
    for cond_set, count in sorted_intersects[:5]:  # Show top 5
        print(f"  {'&'.join(sorted(cond_set))}: {count} peaks")

if __name__ == '__main__':
    main()
