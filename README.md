# RNA-seq Differential Expression Analysis

This repository contains a Python script for performing differential expression analysis on RNA-seq data. The steps include data preprocessing, PCA visualization, TPM normalization, and generating a volcano plot for differentially expressed genes.

## Requirements

To run this script, you'll need the following Python libraries:
- `pandas`
- `numpy`
- `matplotlib`
- `seaborn`
- `scipy`
- `sklearn`

You can install the required libraries using:

## Files
- `rna_seq_analysis.py`: The main script for RNA-seq analysis.
- `metadata.csv`: Metadata file containing sample conditions (Control vs. Treatment).
- `counts.tsv`: RNA-seq count data (replace with your own data).

## How to Run

1. Clone the repository:
    ```
    git clone https://github.com/ibgsehriban/RNA-seq-analysis.git
    ```
2. Navigate to the directory:
    ```
    cd RNA-seq-analysis
    ```
3. Run the script:
    ```
    python rna_seq_analysis.py
    ```

## License
MIT License.
