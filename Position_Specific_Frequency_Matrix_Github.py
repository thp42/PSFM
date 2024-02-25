import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, Normalize
from collections import Counter
import csv

standard_amino_acids = ['R', 'H', 'K', 'S', 'T', 'N', 'Q', 'D', 'E', 'C', 'G', 'P', 'A', 'V', 'I', 'L', 'M', 'F', 'W', 'Y', ]   

#Parse Sequence from a fasta file or MSA as .aln
def parse_fasta(filepath):
    sequences = []
    with open(filepath, 'r') as file:
        sequence = ''
        for line in file:
            if line.startswith('>'):
                if sequence:
                    sequences.append(sequence)
                    sequence = ''
            else:
                sequence += line.strip()
        if sequence:
            sequences.append(sequence)
    return sequences

# Calculate the frequency of each amino acid at each position
def calculate_aa_frequencies(sequences):
    matrix = np.zeros((len(standard_amino_acids), len(sequences[0])))
    
    for i in range(len(sequences[0])):
        column = [seq[i] for seq in sequences]
        counts = Counter(column)
        for j, aa in enumerate(standard_amino_acids):
            matrix[j, i] = counts.get(aa, 0) / len(sequences)
    
    return matrix, standard_amino_acids

# Function to plot the heatmap
def plot_heatmap(matrix, all_amino_acids, title='Amino Acid Frequency Heatmap', highlight_sequence=None):
    masked_matrix = np.ma.masked_where(matrix == 0, matrix)
    cmap = plt.cm.Reds  
    cmaplist = [cmap(i) for i in range(cmap.N)]
    custom_cmap = ListedColormap(cmaplist)
    norm = Normalize(-0.25, vmax=np.max(matrix))
    fig, ax = plt.subplots()
    cax = ax.imshow(masked_matrix, cmap=custom_cmap, norm=norm, interpolation='nearest', aspect='auto')
    fig.colorbar(cax, ax=ax, extend='min')

    if highlight_sequence:
        for i, aa in enumerate(highlight_sequence):
            if aa in all_amino_acids:
                aa_index = all_amino_acids.index(aa)
                # Draw rectangles for highlighting, adjusted for 0.5 offset in imshow
                rect = plt.Rectangle((i-0.5, aa_index-0.5), 1, 1, fill=False, edgecolor='green', lw=2)
                ax.add_patch(rect)
    
    ax.set_xticks(np.arange(matrix.shape[1]))
    ax.set_yticks(np.arange(matrix.shape[0]))
    ax.set_xticklabels(np.arange(1, matrix.shape[1] + 1), rotation=90)
    ax.set_yticklabels(all_amino_acids)
    plt.xlabel('Position')
    plt.ylabel('Amino Acid')
    plt.title(title)
    plt.show()

# Function to export matrix to CSV
def export_to_csv(matrix, all_amino_acids, full_path):
    with open(full_path, mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['Amino Acid'] + [f'Position {i+1}' for i in range(matrix.shape[1])])
        for aa, row in zip(all_amino_acids, matrix):
            writer.writerow([aa] + list(row))
    print(f"Frequency data exported to CSV file: {full_path}")

def filter_matrix(matrix, threshold=0.002):
    filtered_matrix = np.where(matrix > threshold, matrix, 0)
    return filtered_matrix


def main(family_file_paths):
    frequency_matrices = []

    for file_path in family_file_paths:
        sequences = parse_fasta(file_path)
        if sequences: 
            frequency_matrix, _ = calculate_aa_frequencies(sequences)
            filtered_matrix = filter_matrix(frequency_matrix)
            frequency_matrices.append(filtered_matrix)

    # Check if there are matrices to process
    if not frequency_matrices:
        print("No valid matrices to process.")
        return

    # Calculate the average frequency across all matrices
    stacked_matrices = np.stack(frequency_matrices, axis=0)
    average_matrix = np.mean(stacked_matrices, axis=0)

    # Plot the aggregated matrix
    title = 'Average Frequency Across Families the 7 Proteins - USP 54'
    highlight_seq = ""  # Input the highlighting sequence
    plot_heatmap(average_matrix, standard_amino_acids, title=title, highlight_sequence=highlight_seq)


if __name__ == "__main__":
    family_file_paths = [
        #List of File Paths .aln
    ]
    main(family_file_paths)