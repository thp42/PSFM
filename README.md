# PSFM - Position Specific Frequency Matrix and Heatmap Generator

This Python script calculates the frequency of each amino acid at each position from multiple sequence alignments (MSAs, as `.aln`) and generates a heatmap visualization. It is designed to aid in the analysis of protein families by highlighting conservation and variability across sequences.

## Features

- Parses sequences from `.aln` files.
- Calculates amino acid frequencies at each position in the sequence alignment.
- Weights each family equally in the calculation of the average amino acid frequency, ensuring a balanced representation in the heatmap.
- Generates heatmaps to visualize the frequency of each amino acid.
- Implements a cut-off for underrepresentation of amino acids.
- Exports frequency data to a CSV file for further analysis.
- Allows for highlighting specific amino acids in the heatmap.

## Usage

To use the script, prepare a list of file paths to your `.aln` files containing the multiple sequence alignments. Modify the `family_file_paths` list in the `main` function accordingly. You can also add a highlighting sequence, to get an idea of how an individual sequence fits into your conserved sequence.

The heatmap will be displayed for the average frequency across all provided sequences. Note that each sequence family is weighted equally in the average calculation, ensuring that each family contributes identically to the final visualization, regardless of the number of sequences in each family.

