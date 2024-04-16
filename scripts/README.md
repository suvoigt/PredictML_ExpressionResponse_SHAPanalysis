# Scripts Documentation

This document provides detailed descriptions and usage instructions for each script in the `scripts` folder. These scripts are designed to perform various tasks related to data processing, analysis, etc., as part of the project "Identifying Factors Involved in Temperature-Sensitive Expression of Polycomb Group-Regulated Genes Using Machine Learning".

## Table of Contents

- [wgs_trim_map_snps_cns.sh: Whole Genome Sequencing Pipeline](#wgs_trim_map_snps_cnssh-whole-genome-sequencing-pipeline)
- [merge_fastq.sh: Merge FASTQ Files](#merge_fastqsh-merge-fastq-files)
- [rnaseq_map_count.sh: RNA Sequencing Mapping and Counting](#rnaseq_map_countsh-rna-sequencing-mapping-and-counting)
- [get_fasta_per_gene.sh: Extract and Align FASTA Sequences by Gene](#get_fasta_per_genesh-extract-and-align-fasta-sequences-by-gene)
- [get_raw_motif_scores_pwmenrich.R: Motif Analysis for Sequence Groups](#get_raw_motif_scores_pwmenrichr-motif-analysis-for-sequence-groups)
- [get_expr_edger.R: Gene Expression Analysis with edgeR](#get_expr_edgerr-gene-expression-analysis-with-edger)
- [get_features.py: Build Feature Dataframe](#get_featurespy-build-feature-dataframe)
- [train.py: Train and Evaluate Machine Learning Models](#trainpy-train-and-evaluate-machine-learning-models)
- [predict.py: Generate Predictions and Evaluate Model](#predictpy-generate-predictions-and-evaluate-model)


## wgs_trim_map_snps_cns.sh: Whole Genome Sequencing Pipeline

### Purpose

This script performs a comprehensive pipeline for whole genome sequencing data, including trimming of sequencing reads, mapping to a reference genome, SNP calling, and consensus sequence generation.

### Usage

```bash
./wgs_trim_map_snps_cns.sh -i <path/to/fastq_ref/folder> -r <path/to/reference.fasta> [-t <number_of_threads>] [-h <help>]
```

- `-i <path/to/fastq_ref/folder>`: Directory containing FASTQ files for processing.
- `-r <path/to/reference.fasta>`: Path to the reference genome in FASTA format.
- `-t <number_of_threads>`: (Optional) Number of threads to use for processing. Default is 4.
- `-h`: Displays help information and exits.

### Dependencies

- **Trim Galore**: Used for quality and adapter trimming. [GPL v3 License](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/).
- **Samtools**: Utilized for manipulating alignments in the SAM format. [MIT License](http://www.htslib.org/).
- **BWA**: Software package for mapping low-divergent sequences against a large reference genome. [GPL v3 License](http://bio-bwa.sourceforge.net/).
- **GATK (Genome Analysis Toolkit)**: Provides tools for variant discovery in high-throughput sequencing data. [BSD 3-Clause License](https://gatk.broadinstitute.org/hc/en-us).
- **Bcftools**: Program for variant calling and manipulating VCFs and BCFs. [MIT License](http://www.htslib.org/).

### Notes

- Ensure each FASTQ input file has a unique name followed by an underscore, e.g., `sample1_R1.fastq.gz` and `sample1_R2.fastq.gz`.
- The default number of threads is set to 

### Author

- Susanne Voigt


## merge_fastq.sh: Merge FASTQ Files

### Purpose

This script combines multiple FASTQ files for the same sample into single R1 and R2 FASTQ files, respectively. It's designed to handle cases where sequencing outputs for a single sample are split across multiple files, allowing for streamlined downstream analysis.

### Usage

```bash
./merge_fastq.sh -i <path/to/info.csv> -d <path/to/folder_with_fastq_files>
```

- `-i <path/to/info.csv>`: Path to a CSV file where each row corresponds to a sample. The CSV should include: sample name for combined output files, and paths to R1 and R2 FASTQ files to be merged.
- `-d <path/to/folder_with_fastq_files>`: Directory containing FASTQ files to be processed.

The `info.csv` should be formatted with the following columns:
- **Column 1**: Sample name for the combined output files.
- **Column 2**: First R1 FASTQ file to merge.
- **Column 3**: First R2 FASTQ file to merge.
- **Column 4**: (Optional) Second R1 FASTQ file to merge.
- **Column 5**: (Optional) Second R2 FASTQ file to merge.

### Output Files

- For each sample, the script generates two output files:
  - `sample_name_R1.fastq.gz`: Combined R1 reads for the sample.
  - `sample_name_R2.fastq.gz`: Combined R2 reads for the sample.

### Dependencies

- **pigz**: A parallel implementation of gzip for modern multi-processor, multi-core machines. [Zlib License](https://zlib.net/zlib_license.html).

### Notes

- Ensure that the `info.csv` file correctly maps sample names to their corresponding FASTQ files.
- The script verifies that the combined R1 and R2 files for each sample have the same number of reads. A warning is issued for any discrepancies.
- Original FASTQ files are removed after successful verification of combined file integrity.

### Author

- Susanne Voigt


## rnaseq_map_count.sh: RNA Sequencing Mapping and Counting

### Purpose

This script automates the mapping of RNA sequencing reads to a reference genome and the counting of gene expression levels. It utilizes STAR for mapping and leverages subread for additional counting capabilities.

### Usage

```bash
./rnaseq_map_count.sh -i <path/to/fastq_ref/folder> -r <path/to/reference.fasta> -g <path/to/reference.gtf> [-t <number_of_threads>] [-h <help>]
```

- `-i <path/to/fastq_ref/folder>`: Directory containing FASTQ files for processing.
- `-r <path/to/reference.fasta>`: Path to the reference genome in FASTA format.
- `-g <path/to/reference.gtf>`: Path to the reference genome annotation file in GTF format.
- `-t <number_of_threads>`: (Optional) Number of threads to use for processing. Default is 4.
- `-h`: Displays help information and exits.

### Dependencies

- **STAR**: Used for efficient RNA sequencing read mapping. [GPLv3 License](https://github.com/alexdobin/STAR).
- **Subread (featureCounts)**: Utilized for counting gene features from mapped reads. [GPLv3 License](http://subread.sourceforge.net/).

### Notes

- Ensure each FASTQ input file has a unique name followed by an underscore, e.g., `sample1_R1.fastq.gz` and `sample1_R2.fastq.gz`.
- The script provides detailed gene counts for downstream expression analysis.
- The default number of threads is set to 4 unless specified otherwise by the `-t` option.

### Author

- Susanne Voigt


## get_fasta_per_gene.sh: Extract and Align FASTA Sequences by Gene

### Purpose

This script extracts FASTA sequences for specified genomic regions from one or more input FASTA files and then aligns these sequences. It's particularly useful for analyzing specific genes or regions across different sequences or species. The script uses a BED file to define the segments to be extracted, allowing for flexible, targeted sequence analysis.

### Usage

```bash
./get_fasta_per_gene.sh -i <path/to/input1.fas> [-i <path/to/input2.fas> ...] -r <path/to/regions.bed> [-h <help>]
```

- `-i <path/to/input.fas>`: One or more paths to input FASTA files. Multiple files are accepted by repeating the `-i` option.
- `-r <path/to/regions.bed>`: Path to a BED file defining the segments to extract. The BED file should have tab-separated columns with col1: chrom, col2: start, col3: end, col4: name.
- `-h`: Displays help information and exits.

### Dependencies

- **Bedtools**: Used for extracting sequences from FASTA files based on coordinates specified in a BED file. [MIT License](https://bedtools.readthedocs.io/en/latest/content/license.html).
- **MAFFT**: Utilized for aligning sequences extracted from the FASTA files. [MAFFT License](https://mafft.cbrc.jp/alignment/software/license.html) (BSD-like, very permissive).

### Notes

- The script generates a new directory `fasta_per_gene` within the same directory as the first input file, where it stores the extracted and aligned sequences for each region defined in the BED file.
- Each output file is named according to the gene and coordinates specified in the BED file, facilitating easy identification and subsequent analysis.

### Author

- Susanne Voigt


## get_raw_motif_scores_pwmenrich.R: Motif Analysis for Sequence Groups

### Purpose

This script performs motif analysis for transcription factor binding sites on provided sequence groups using PWMEnrich. It calculates raw motif scores for one or two groups of sequences, allowing for detailed examination of motif prevalence and strength across different sequence datasets.

### Usage

```bash
Rscript get_raw_motif_scores_pwmenrich.R --input <path/to/input.csv> [--output <path/to/output_raw_scores_pwmenrich.csv>] [--n_group1 <number_of_sequences_of_group1>] [--cores <number_of_cores>]
```

- `--input`: Path and name of the input CSV file containing a dataframe (col1: gene, col2: path/to/fasta).
- `--output`: (Optional) Path and name of the output file. Default is `"output_raw_scores_pwmenrich.csv"`.
- `--n_group1`: (Optional) Number of sequences of group 1. If not provided, motif scores are calculated across all sequences in the fasta.
- `--cores`: (Optional) Number of cores to utilize for the analysis. Default is 1.

### Dependencies

- **optparse**: Provides a mechanism for command-line parsing akin to what is available in Python with argparse. License: [GPL-3](https://www.gnu.org/licenses/gpl-3.0.en.html). More details: [CRAN - optparse](https://CRAN.R-project.org/package=optparse).
- **PWMEnrich**: A suite of tools for motif enrichment analyses in PWMs (Position Weight Matrices). License: [GPL-3](https://www.gnu.org/licenses/gpl-3.0.en.html). More information: [Bioconductor - PWMEnrich](https://bioconductor.org/packages/release/bioc/html/PWMEnrich.html).
- **PWMEnrich.Dmelanogaster.background**: Provides background datasets for PWMEnrich package specifically tailored for D. melanogaster. License: [Artistic-2.0](https://opensource.org/licenses/Artistic-2.0). More information: [Bioconductor - PWMEnrich.Dmelanogaster.background](https://bioconductor.org/packages/release/data/experiment/html/PWMEnrich.Dmelanogaster.background.html).


### Notes

- Ensure that `R` and the required libraries (`optparse`, `PWMEnrich`, `PWMEnrich.Dmelanogaster.background`) are installed and available.
- The script requires a CSV input specifying genes and paths to corresponding fasta files. Each fasta file should align with the structure anticipated by the analysis, potentially delineating groups of sequences for differential motif analysis.

### Author

- Susanne Voigt


## get_expr_edger.R: Gene Expression Analysis with edgeR

### Purpose

This script conducts differential gene expression analysis utilizing edgeR, processing gene count data alongside sample information to identify genes that are differentially expressed between groups. It leverages the power of contrasts to make specific group comparisons as defined by the user.

### Usage

```bash
Rscript get_expr_edger.R --input <path/to/gene_counts.csv> --samples <path/to/sample_info.csv> --output <path/to/output_directory> --annotation <number_of_annotation_columns> --contrast 'group2 - group1, group3 - group2'
```

- `--input`: Path and name of the input CSV file containing gene annotations followed by gene counts per sample.
- `--samples`: Path and name of the CSV file containing sample information (col1: replicates, col2: group).
- `--output`: Path for output files. Default is the current working directory.
- `--annotation`: Number of columns before the columns with gene counts in the input.csv that contain gene annotations.
- `--contrast`: Comma-separated list of contrasts to perform (e.g., "group2 - group1, group3 - group2"). This option is required for the analysis to specify the group comparisons.

### Dependencies

- **edgeR**: Used for differential expression analysis of RNA-seq expression profiles. License: [GPL-3](https://www.gnu.org/licenses/gpl-3.0.en.html). More information and citation: [Bioconductor - edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html).
- **optparse**: Provides a mechanism for command-line parsing akin to what is available in Python with argparse. License: [GPL-3](https://www.gnu.org/licenses/gpl-3.0.en.html). More details: [CRAN - optparse](https://CRAN.R-project.org/package=optparse).

### Notes

- It is essential that `R` and the necessary libraries (`optparse`, `edgeR`) are installed and accessible.
- The script necessitates a specific format for the input files: one for gene counts and another for sample information.

### Author

- Susanne Voigt


## get_features.py: Build Feature Dataframe

### Purpose

`get_features.py` creates a feature dataframe from transcription factor binding site (TFBS) motif scores for gene regions across different groups or populations. It's designed to support gene expression analysis with or without additional gene expression data. The primary output is a motif-features-only dataframe, suitable for analysis or machine learning prediction. If provided, gene expression data enriches the dataframe with expression features (logFCs), a binary target vector ('TSE', indicating significant overexpression in the positive direction (1) or none (0)), and a baseline predictor dataset consisting of gene activity (i.e., mean logCPM values per gene).

### Usage

```bash
python get_features.py --motif_scores <path/to/motif_scores.csv> --motif_info <path/to/motif_info.csv> --annotation <path/to/annotation.csv> --n_groups <number_of_groups> [--expression <path/to/expression_group1.csv> <path/to/expression_group2.csv> ...] --output <path/to/features.csv> [--output2 <path/to/gene_activity.csv>]
```

### Input Files and Format

#### `--motif_scores`
- **Description**: Required. CSV file with raw motif scores.
- **Format**: Output from `get_raw_motif_scores_pwmenrich.R`, listing raw motif scores for sequences analyzed.

#### `--motif_info`
- **Description**: Required. CSV file detailing information on non-redundant motifs.
- **Format**: Includes columns for TF name, motif ID, consensus sequence, data source, and combined TF and motif identifier.

#### `--annotation`
- **Description**: Required. CSV file providing gene annotations.
- **Format**: Contains gene ID, chromosome location, start and end positions, and gene name.

#### `--n_groups`
- **Description**: Required. Number of separate groups/lines/populations.
- **Format**: Integer indicating the distinct groups included in the analysis.

#### `--expression` (Optional)
- **Description**: Path(s) to CSV file(s) with gene expression data (edgeR output).
- **Format**: Each file should correspond to a group indicated by `--n_groups`, containing genome-wide expression results.

### Output Files

#### `--output`
- **Description**: Path for the output features dataframe CSV. Defaults to "features.csv".
- **Details**: The output includes a comprehensive dataframe, ready for further analysis or predictive modeling with machine learning techniques.

#### `--output2` (Optional)
- **Description**: Path for an output CSV with baseline predictor data (mean gene activity). Defaults to "gene_activity.csv".
- **Details**: Provides baseline gene activity levels for use as a baseline predictor dataset.

### Dependencies 

- **Python**: The core programming language used. Python is licensed under the [Python Software Foundation License](https://docs.python.org/3/license.html), a GPL-compatible but non-copyleft license.
- **pandas**: A powerful data manipulation and analysis library. pandas is open source and licensed under the [BSD 3-clause "New" or "Revised" License](https://pandas.pydata.org/pandas-docs/stable/getting_started/overview.html#license).
- **NumPy**: A fundamental package for scientific computing with Python. NumPy is also licensed under the [BSD 3-Clause "New" or "Revised" License](https://numpy.org/license.html).
- **argparse**: A standard Python module for parsing command-line options and arguments, included with Python. Since it's part of the Python Standard Library, its use is governed by Python's license.


### Notes

- The accuracy of the input files' format is crucial for the script to function as intended. Proper preparation of these files ensures seamless creation of the feature dataframe, facilitating robust gene expression analysis and modeling.

### Author

- Susanne Voigt


## train.py: Train and Evaluate Machine Learning Models

### Purpose

`train.py` facilitates the training and evaluation of machine learning models using a specified dataset. It selects specific features for training, loads an untrained model (.joblib file), trains this model with the dataset, and evaluates its performance using precision-recall and ROC metrics. If a baseline predictor dataset is provided, the script also compares the model's performance against this baseline.

### Usage

```bash
python train.py --input <path/to/dataset.csv> --features <path/to/features.csv> --target <target_column_name> --model <path/to/untrained_model.joblib> [--output <path/to/trained_model.joblib>] [--base <path/to/baseline_predictor.csv>]
```

### Input Files and Format

- `--input`: CSV file containing the dataset including features and target vector (get_features.py output).
- `--features`: CSV file listing the features to be selected. Features should be listed in the first column.
- `--target`: The name of the target vector column in the input dataset.
- `--model`: Path to the saved untrained model file (.joblib).

### Output Files

- `--output` (Optional): Path where the trained model will be saved. Defaults to "trained_model.joblib" in the current working directory.
- `--base` (Optional): CSV file of the baseline predictor dataset, which must also include the target vector alongside the baseline feature.

### Dependencies

- **Python**: The core programming language used for the script. Python's standard library, including modules like `argparse`, is covered under the [Python Software Foundation License Version 2](https://www.python.org/download/releases/2.7/license/).
- **pandas**: An open-source, BSD-licensed library providing high-performance, easy-to-use data structures and data analysis tools for Python. License: [BSD 3-Clause License](https://pandas.pydata.org/pandas-docs/stable/getting_started/overview.html#license).
- **joblib**: A set of tools to provide lightweight pipelining in Python. Mainly used for saving and loading Python objects that make use of NumPy data structures, efficiently. License: [BSD 3-Clause License](https://joblib.readthedocs.io/en/latest/license.html).
- **scikit-learn**: A simple and efficient tool for data mining and data analysis. Built on NumPy, SciPy, and matplotlib, it's open-source and commercially usable - BSD licensed. License: [BSD 3-Clause License](https://scikit-learn.org/stable/about.html#citing-scikit-learn).

### Notes

This script is useful for preparing datasets, training machine learning models, and evaluating their performance in a systematic manner. It supports the use of custom features and comparison against baseline predictors.

### Author

- Susanne Voigt


## predict.py: Generate Predictions and Evaluate Model

### Purpose

`predict.py` uses a trained machine learning model to generate predictions on a provided dataset. It can select specific features for prediction, load a trained model, and optionally evaluate predictions if a target vector is available. The script calculates the area under the Precision-Recall Curve (auPR) and the area under the Receiver Operating Characteristic Curve (auROC) for evaluation.

### Usage

```bash
python predict.py --input <path/to/dataset.csv> --features <path/to/features.csv> --model <path/to/trained_model.joblib> [--output <path/to/predictions.csv>] [--target <target_column_name>]
```

### Input Files and Format

- `--input`: CSV file containing the dataset including features and optional target vector (get_features.py output).
- `--features`: CSV file listing the features to be selected. Features should be listed in the first column.
- `--model`: Path to the saved trained model file (.joblib).

### Output Files

- `--output` (Optional): Path where the predictions will be saved. Defaults to "data_pred.csv" in the current working directory.

### Dependencies

- **Python**: The core programming language used. Python is licensed under the [Python Software Foundation License](https://docs.python.org/3/license.html), a GPL-compatible but non-copyleft license.
- **pandas**: A powerful data manipulation and analysis library. pandas is open source and licensed under the [BSD 3-clause "New" or "Revised" License](https://pandas.pydata.org/pandas-docs/stable/getting_started/overview.html#license).
- **joblib**: A set of tools to provide lightweight pipelining in Python. Mainly used for saving and loading Python objects that make use of NumPy data structures, efficiently. License: [BSD 3-Clause License](https://joblib.readthedocs.io/en/latest/license.html).
- **scikit-learn**: A simple and efficient tool for data mining and data analysis. Built on NumPy, SciPy, and matplotlib, it's open-source and commercially usable - BSD licensed. License: [BSD 3-Clause License](https://scikit-learn.org/stable/about.html#citing-scikit-learn).

### Notes

This script is instrumental in applying trained models to new datasets for generating predictions. It supports model evaluation through precision-recall and ROC metrics if a target vector is provided, facilitating an understanding of model performance on unseen data.

### Author

- Susanne Voigt
