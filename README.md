# PSFM - Position Specific Frequency Matrix and Heatmap Generator

This Python script calculates the frequency of each amino acid at each position from multiple sequence alignments (MSAs, as `.fasta`) and generates a heatmap visualization. It is designed to aid in the analysis of protein families by highlighting conservation and variability across sequences.

---

## System Requirements

<details>

---

### Software Requirements

#### **OS Requirements**
This package should work for most operating systems, as it is a standalone python script. The package has been tested on the following systems:

  - Linux: Ubuntu 20.04

#### **Dependencies**
There is no ```environment.yml``` file provided. Please ensure the following Python packages are installed in your environment. The typical installation time should take some minutes. 
```  
python v.3.11.5
numpy v.1.26.2
matplotlib v.3.9.2
```

---

### Hardware Requirements


#### **Recommended System** and **Database Storage Requirements**

There are ***no specific hardware** requirements for running this script. A standard laptop or desktop system is sufficient for typical use cases.

</details>

---

## Installation & Environment

This is a standalone Python script — no compilation or packaging is required. Simply clone the repository and make sure you have the required Python packages installed.

---


## Function

<details>

- Parses sequences from `.fasta` files.
- Calculates amino acid frequencies at each position in the sequence alignment.
- Weights each family equally in the calculation of the average amino acid frequency, ensuring a balanced representation in the heatmap.
- Generates heatmaps to visualize the frequency of each amino acid.
- Implements a cut-off for underrepresentation of amino acids.
- Exports frequency data to a CSV file for further analysis.
- Allows for highlighting specific amino acids in the heatmap.

</details>

---

## Instructions

<details>

To use the script, prepare a list of file paths to your `.fasta` files containing the multiple sequence alignments. For this you can take MSAs generated either by jackhmmr (see also https://github.com/thp42/SLiMFold) or searched by BLAST (https://blast.ncbi.nlm.nih.gov/Blast.cgi). 
Modify the ```plt.savefig```, ```output_csv_path``` and ```family_file_paths``` list in the python script accordingly. You can also add a highlighting sequence, to get an idea of how an individual sequence fits into your conserved sequence.

The heatmap will be displayed for the average frequency across all provided sequences. Note that each sequence family is weighted equally in the average calculation, ensuring that each family contributes identically to the final visualization, regardless of the number of sequences in each family.

</details>

---


## Citation

If you use this script in published research, please cite:
- Manuscript
  - **Discovery of a new evolutionarily conserved short linear F-actin binding motif**: Themistoklis Paraschiakos, Biao Yuan, Kostiantyn Sopelniak, Michael Bucher, Lisa Simon, Ksenija Zonjic, Dominic Eggers, Franziska Selle, Jing Li, Stefan Linder, Thomas C. Marlovits, Sabine Windhorst. bioRxiv 2025.04.16.649135; [doi: https://doi.org/10.1101/2025.04.16.649135](https://www.biorxiv.org/content/10.1101/2025.04.16.649135v1)

---
