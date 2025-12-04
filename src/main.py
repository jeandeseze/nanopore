import pysam
import re
from collections import Counter
import pandas as pd 
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib
matplotlib.use('Agg')

# --- Configuration ---

# The minimum length of the soft-clipped region you will consider (e.g., 1 base)
MIN_OVERHANG_LENGTH = 1 
# The minimum mapping quality (MAPQ) to keep a 'well enough' aligned track.
# A common threshold for good alignments is 20 or 30.
MIN_MAPQ = 60 
# ---------------------

def analyze_5prime_overhang(bam_file, min_overhang_len, min_mapq):
    """
    Analyzes a BAM file for 5' soft-clipped sequences (the 'overhang').
    """
    # Use 'samtools view -c' equivalent to get total mapped reads for proportion calculation
    total_aligned_reads = 0
    overhang_sequences = Counter()
    
    # Pre-compiled regex to find the length of the leading Soft-Clip (e.g., '10S')
    # This only checks if the CIGAR *starts* with a soft-clip.
    cigar_regex = re.compile(r"^(\d+)S") 
    
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        
        # 1. & 2. FILTERING AND EXTRACTION
        for read in bam.fetch():
            
            # Skip unmapped, secondary, or supplementary alignments
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue

            total_aligned_reads += 1
            
            # Filter for alignment quality (your 'tracks where alignment could be done well enough')
            if read.mapping_quality < min_mapq:
                continue
            
            # Check the CIGAR string for a leading Soft-Clip ('S' operator)
            match = cigar_regex.match(read.cigarstring)
            
            if match:
                # Get the length of the soft-clip (the number captured in the regex)
                clip_length = int(match.group(1)) 
                
                if clip_length >= min_overhang_len:
                    
                    # 3. READ THE NUCLEOTIDES THAT ARE BEFORE
                    # The read.seq attribute holds the full read sequence.
                    # The soft-clipped bases are at the beginning (5' end) of the read.
                    # The sequence is read.seq[0 : clip_length]
                    overhang_seq = read.seq[:clip_length]
                    
                    # Tally the sequence
                    overhang_sequences[overhang_seq] += 1

    # --- RESULTS REPORTING ---
    # print(f"--- Analysis Summary for {BAM_FILE} ---")
    # print(f"Total Aligned Primary Reads: {total_aligned_reads}")
    # print(f"Reads with {min_overhang_len}+ nt 5'-Overhang: {sum(overhang_sequences.values())}")
    
    # print("\n--- Overhang Sequence Counts & Proportions ---")
    
    # Calculate total reads with an overhang for the denominator
    total_overhangs = sum(overhang_sequences.values())

    if total_aligned_reads > 0:
        for seq, count in overhang_sequences.most_common():
            proportion = (count / total_aligned_reads) * 100
            # print(f"Sequence: {seq:<10} | Count: {count:<8} | Proportion of Total Mapped Reads: {proportion:.3f}%")
    else:
        print("No reads passed the alignment filters.")
    
    return overhang_sequences

def calculate_positional_proportions(overhang_data, file_label):
    """
    Calculates nucleotide proportions for the last three positions of the overhang
    and returns the data as a Pandas DataFrame.
    """
    positional_counts = {
        -1: Counter({'A': 0, 'C': 0, 'G': 0, 'T': 0}),
        -2: Counter({'A': 0, 'C': 0, 'G': 0, 'T': 0}),
        -3: Counter({'A': 0, 'C': 0, 'G': 0, 'T': 0}),
        -4: Counter({'A': 0, 'C': 0, 'G': 0, 'T': 0})
    }
    
    # 1. Tally the bases for each position (same as before)
    for seq, count in overhang_data.items():
        L = len(seq)
        if L >= 1: positional_counts[-1][seq[-1]] += count
        if L >= 2: positional_counts[-2][seq[-2]] += count
        if L >= 3: positional_counts[-3][seq[-3]] += count
        if L >= 4: positional_counts[-4][seq[-4]] += count

    # 2. Convert to DataFrame structure
    data_for_df = []
    
    for pos, base_counts in positional_counts.items():
        total_bases_at_pos = sum(base_counts.values())
        
        if total_bases_at_pos == 0:
            continue
            
        for base, count in base_counts.items():
            proportion = (count / total_bases_at_pos) * 100
            data_for_df.append({
                'Position': pos,
                'File': file_label, # ADDED FILE LABEL HERE
                'Nucleotide': base,
                'Count': count,
                'Proportion': proportion
            })

    # Create the DataFrame
    if data_for_df:
        df = pd.DataFrame(data_for_df)
        return df
    else:
        return pd.DataFrame()

def plot_proportions_by_position(full_df):
    """
    Generates a separate stacked bar plot for each position, comparing all files.
    """
    if full_df.empty:
        print("Cannot plot: Global DataFrame is empty.")
        return

    # Identify all unique positions found in the data (e.g., -1, -2, -3)
    unique_positions = sorted(full_df['Position'].unique(), reverse=True)
    
    # 1. Loop through each position to create a separate plot
    for pos in unique_positions:
        # Filter the master DataFrame for the current position
        pos_df = full_df[full_df['Position'] == pos]
        
        # Pivot the data for stacking: Indices=File, Columns=Nucleotide, Values=Proportion
        plot_df = pos_df.pivot(index='File', columns='Nucleotide', values='Proportion').fillna(0)
        
        # Reorder columns for standard display (A, C, G, T)
        order = ['A', 'C', 'G', 'T']
        plot_df = plot_df[[col for col in order if col in plot_df.columns]]
        
        # --- Plotting ---
        plt.figure(figsize=(10, 5))
        
        plot_df.plot(
            kind='bar', 
            stacked=True, 
            color=sns.color_palette("Set1", n_colors=4),
            ax=plt.gca()
        )
        
        plt.title(f'Nucleotide Proportion at Position {pos} Across All Files')
        plt.xlabel('BAM File / Sample')
        plt.ylabel('Proportion (%)')
        plt.legend(title='Base', bbox_to_anchor=(1.05, 1), loc='upper left')
        # plt.xticks(rotation=45, ha='right') # Rotate file names for readability
        plt.tight_layout()
        
        filename = f'position_{pos}_cumulative_barplot.png'
        plt.savefig(filename)
        print(f"Plot for Position {pos} saved to {filename}")

if __name__ == "__main__":
    # BBAM_FILES = [f"./rsvwithcap/barcode0{i+1}.bam" for i in range(5)]
    BAM_FILES = [f"./rsvwithgg/FBE82653_pass_barcode0{i+1}_03a1f7c4_9c5a4d8a_0.bam" for i in range(5)]
    BAM_FILES = [f"./rsvwithgonly/FBE82653_pass_barcode0{i+1}_03a1f7c4_9c5a4d8a_0.bam" for i in range(5)]
    
    # List to collect all resulting DataFrames
    all_dfs = [] 
    
    # Loop through all BAM files
    for file_label in BAM_FILES:
        print(f"Analyzing {file_label}...")
        
        # 1. Analyze for overhang sequences
        overhang_results = analyze_5prime_overhang(file_label, MIN_OVERHANG_LENGTH, MIN_MAPQ)
        
        if overhang_results:
            # 2. Calculate proportions and capture the DataFrame
            df = calculate_positional_proportions(overhang_results, file_label)
            all_dfs.append(df)
        else:
            print(f"Warning: No valid overhang sequences found for {file_label}.")

    # 3. Concatenate all DataFrames into one master DataFrame
    if all_dfs:
        final_df = pd.concat(all_dfs, ignore_index=True)
        print("\n--- Consolidated DataFrame Head ---")
        print(final_df.head())
        
        # 4. Generate the plots
        plot_proportions_by_position(final_df)
        final_df.to_csv('filename.csv', index=False)
    else:
        print("Analysis failed for all files. No plots generated.")

