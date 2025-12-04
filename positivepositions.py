import pysam
import re
from collections import Counter
import pandas as pd 
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib
matplotlib.use('Agg')

from src import calculate_positional_proportions, plot_proportions_by_position 

def analyze_positive_positions(bam_file, file_label, positions_to_analyze, min_mapq):
    """
    Analyzes nucleotide proportions for specified positive positions, 
    using position-specific coverage as the denominator.
    """
    
    positive_positions = [p for p in positions_to_analyze if p > 0]
    if not positive_positions:
        return pd.DataFrame()

    # Counter for observed bases (numerator)
    positional_base_counts = {
        pos: Counter({'A': 0, 'C': 0, 'G': 0, 'T': 0})
        for pos in positive_positions
    }
    
    # Counter for total reads covering that position (denominator)
    total_coverage_at_pos = Counter({pos: 0 for pos in positive_positions})
    
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        
        for read in bam.fetch():
            
            # 1. FILTER READS
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            
            if read.mapping_quality < min_mapq:
                continue

            aligned_pairs = read.get_aligned_pairs(matches_only=False, with_seq=False)
            read_seq = read.seq 

            # 2. ANALYZE EACH DESIRED POSITION
            for pos in positive_positions:
                target_ref_idx = pos - 1
                for r_idx, ref_idx in aligned_pairs:
                    if ref_idx == target_ref_idx:
                        # Found the column in the alignment corresponding to target_ref_idx

                        if r_idx is not None:
                            # Base was OBSERVED (Matched/Mismatched/Substituted)
                            base = read_seq[r_idx]
                            positional_base_counts[pos][base] += 1
                       
                        total_coverage_at_pos[pos] += 1 

    # 3. Convert to DataFrame structure and CALCULATE PROPORTIONS
    data_for_df = []
    
    for pos in positive_positions:
        base_counts = positional_base_counts[pos]
        total_coverage = total_coverage_at_pos[pos] # Use position-specific denominator
        
        if total_coverage == 0:
            continue
            
        for base, count in base_counts.items():
            if count > 0:
                proportion = (count / total_coverage) * 100 # Corrected division
                data_for_df.append({
                    'Position': pos,
                    'File': file_label,
                    'Nucleotide': base,
                    'Count': count, 
                    'Proportion': proportion,
                    'Total_Coverage': total_coverage # Helpful for sanity check
                })

    if data_for_df:
        return pd.DataFrame(data_for_df)
    else:
        return pd.DataFrame()
    

POSITIONS_TO_ANALYZE = [1,2,3,4,5]
BAM_FILES = [f"./rsvwithcap/barcode0{i}.bam" for i in [1,2,3,4,5]]
BAM_FILES = [f"./rsvwithgg/FBE82653_pass_barcode0{i+1}_03a1f7c4_9c5a4d8a_0.bam" for i in range(5)]
MIN_OVERHANG_LENGTH = 1 
MIN_MAPQ = 50 

if __name__ == "__main__":
    
    all_dfs = [] 
    
    # Split the list into positive and negative sets for dedicated functions
    pos_positions = [p for p in POSITIONS_TO_ANALYZE if p > 0]
    
    for file_label in BAM_FILES:
        print(f"Analyzing {file_label}...")

        pos_df = analyze_positive_positions(
            file_label,
            file_label, 
            pos_positions, 
            MIN_MAPQ # Re-using the MAPQ filter
        )
        all_dfs.append(pos_df)
        
    if all_dfs:
        final_df = pd.concat(all_dfs, ignore_index=True)
        print(final_df)
        # Ensure the final plotting step sorts the positions correctly
        plot_proportions_by_position(final_df) 
    else:
        print("Analysis failed for all files. No plots generated.")