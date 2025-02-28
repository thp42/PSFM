import os
import csv
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from collections import Counter
from typing import List, Tuple

# Define a list of standard amino acids
STANDARD_AMINO_ACIDS = [
    'R', 'H', 'K', 'D', 'E', 'N', 'Q', 'S', 'T', 'C',
    'G', 'P', 'A', 'V', 'I', 'L', 'M', 'F', 'W', 'Y'
]

def parse_fasta(filepath: str) -> List[str]:
    """
    Parse a FASTA file and return a list of sequences.
    """
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"File not found: {filepath}")
    
    sequences = []
    with open(filepath, 'r') as file:
        sequence = []
        for line in file:
            line = line.strip()
            if line.startswith('>'):
                # Header line indicates a new sequence
                if sequence:
                    sequences.append(''.join(sequence))
                    sequence = []
            else:
                sequence.append(line)
        # Append the last sequence if present
        if sequence:
            sequences.append(''.join(sequence))
    return sequences

def calculate_aa_frequencies(sequences: List[str], 
                             amino_acids: List[str] = STANDARD_AMINO_ACIDS) -> Tuple[np.ndarray, List[str]]:
    """
    Calculate the frequency of each amino acid at each position across all sequences.
    """
    if not sequences:
        raise ValueError("No sequences provided for frequency calculation.")
    
    seq_length = len(sequences[0])
    if any(len(seq) != seq_length for seq in sequences):
        raise ValueError("Not all sequences are the same length.")
    
    matrix = np.zeros((len(amino_acids), seq_length))
    num_sequences = len(sequences)
    
    for i in range(seq_length):
        column = [seq[i] for seq in sequences]
        counts = Counter(column)
        for j, aa in enumerate(amino_acids):
            matrix[j, i] = counts.get(aa, 0) / num_sequences
            
    return matrix, amino_acids

def filter_matrix(matrix: np.ndarray, threshold: float = 0.00) -> np.ndarray:
    """
    Filter the frequency matrix, setting values below threshold to 0.
    """
    return np.where(matrix > threshold, matrix, 0)

def plot_heatmap(matrix: np.ndarray, 
                 amino_acids: List[str], 
                 title: str = 'Amino Acid Frequency Heatmap', 
                 highlight_sequence: str = None) -> None:
    """
    Plot a heatmap of amino acid frequencies with a publication-quality style.
    """
    # Create a colorblind-friendly colormap
    cmap = LinearSegmentedColormap.from_list('white_red', ['#ffffff', '#DB6B6A'])

    # Normalize frequencies between 0 and 1
    norm = Normalize(vmin=0, vmax=1)

    # Set up figure and axis with minimal style
    fig, ax = plt.subplots(figsize=(11, 8))

    # Increase DPI and ensure consistent font
    plt.rcParams['figure.dpi'] = 300
    plt.rcParams['font.family'] = 'arial'  # Consistent serif font
    plt.rcParams['font.size'] = 18

    # Plot the heatmap
    cax = ax.imshow(matrix, 
                    cmap=cmap, 
                    norm=norm, 
                    interpolation='nearest', 
                    aspect='auto')

    # Add horizontal colorbar below the graph
    cbar = fig.colorbar(cax, ax=ax, orientation='horizontal', fraction=0.05, pad=0.1)
    cbar.set_label("√Frequency", fontsize=18, family='arial', weight='bold')
    cbar.outline.set_linewidth(2)



    # Highlight specific amino acids in highlight_sequence if provided
    if highlight_sequence:
        for i, aa in enumerate(highlight_sequence):
            if aa in amino_acids:
                aa_index = amino_acids.index(aa)
                rect = plt.Rectangle((i - 0.5, aa_index - 0.5), 
                                     1, 1, fill=False, edgecolor='green', lw=3.0)  # Increased thickness
                ax.add_patch(rect)

    # Set ticks and labels
    num_positions = matrix.shape[1]
    x_labels = np.arange(-6, -6 + num_positions)
    ax.set_xticks(np.arange(num_positions))
    ax.set_xticklabels(x_labels, rotation=0, fontsize=22, family='arial')

    y_labels = amino_acids
    ax.set_yticks(np.arange(len(amino_acids)))
    ax.set_yticklabels(y_labels, fontsize=22, family='arial')

    # Adjust tick size for both x and y axes
    ax.tick_params(axis='x', labelsize=22, length=6, width=2)  # Adjust label size and tick appearance
    ax.tick_params(axis='y', labelsize=22, length=6, width=2)


    # Make the border of the graph thicker
    for spine in ax.spines.values():
        spine.set_linewidth(2)  # Set thickness to 2.5


    # Example positions to mark
    positions_to_mark = [1, 4, 8, 9]
    for pos in positions_to_mark:
        rect = plt.Rectangle((pos + 6 - 0.5, -0.5),
                             1, len(amino_acids),
                             fill=False, edgecolor='black', lw=1.5)  # Increased thickness
        ax.add_patch(rect)

    # Add axis labels and title with slightly larger fonts
    ax.set_xlabel('Position', fontsize=24, family='arial', fontweight='bold')
    ax.set_ylabel('Amino Acid', fontsize=24, family='arial', fontweight='bold')
    ax.set_title(title, fontsize=16, pad=20, family='arial')

    # Tight layout for better spacing
    plt.tight_layout()

    # Save as EPS file
    plt.savefig(r"G:\YOURPATH\output_heatmap.eps", format='eps', bbox_inches='tight')

    plt.show()



def power_transform(matrix: np.ndarray, power: float = 0.5) -> np.ndarray:
    """
    Apply a power-law transformation to flatten steepness.
    
    :param matrix: 2D numpy array of frequencies.
    :param power: Power factor to flatten steepness (default: square root, 0.5).
    :return: Power-transformed matrix.
    """
    return matrix ** power





def export_to_csv(matrix: np.ndarray, amino_acids: List[str], full_path: str) -> None:
    """
    Export the frequency matrix to a CSV file.
    """
    with open(full_path, mode='w', newline='') as file:
        writer = csv.writer(file)
        header = ['Amino Acid'] + [f'Position {i+1}' for i in range(matrix.shape[1])]
        writer.writerow(header)
        
        for aa, row in zip(amino_acids, matrix):
            writer.writerow([aa] + list(row))
    
    print(f"Frequency data exported to CSV file: {full_path}")


def main(family_file_paths: List[str]) -> None:
    """
    Main function to parse multiple alignment files, calculate average frequencies, 
    apply log transformation, plot them, and save the results to CSV.
    """
    frequency_matrices = []
    for file_path in family_file_paths:
        sequences = parse_fasta(file_path)
        if not sequences:
            print(f"No sequences found in file: {file_path}")
            continue
        
        frequency_matrix, _ = calculate_aa_frequencies(sequences)
        filtered_matrix = filter_matrix(frequency_matrix)
        frequency_matrices.append(filtered_matrix)

    if not frequency_matrices:
        print("No valid matrices to process.")
        return

    # Calculate average frequency across all matrices
    stacked_matrices = np.stack(frequency_matrices, axis=0)
    average_matrix = np.mean(stacked_matrices, axis=0)

    # Apply smoothed natural log transformation
#    log_transformed_matrix = power_transform(average_matrix, power=0.5)
    log_transformed_matrix = average_matrix
    # Highlight sequence for plotting
    highlight_seq = ''
    title = 'Log-Transformed Average Frequency Across Families - Protein Shroom3_Helix2'

    # Plot and save the figure
    plot_heatmap(log_transformed_matrix, STANDARD_AMINO_ACIDS, title=title, highlight_sequence=highlight_seq)

    # Export the log-transformed matrix to CSV
    output_csv_path = r"G:\YOURPATH\log_transformed_averaged_frequency_matrix.csv"
    export_to_csv(log_transformed_matrix, STANDARD_AMINO_ACIDS, output_csv_path)


if __name__ == "__main__":
    family_file_paths = [
        r'G:\YOURPATH\Protein1.fasta',
        r'G:\YOURPATH\Protein2.fasta',
    ]
    main(family_file_paths)
