# hStouffer: the Enhanced Meta-Analysis Method for the Comprehensive Analysis of Large-Scale RNA-seq Data

## Background
Advances in high-throughput sequencing have led to an exponential increase in publicly available RNA-seq datasets, enhancing opportunities for large-scale transcriptomic analyses. Compared to small-scale studies, large-scale RNA-seq meta-analyses offer greater robustness, reduce technical noise, and improve detection of subtle or novel gene expression patterns. However, while traditional meta-analysis methods—particularly p-value-based approaches like Stouffer's—have proven useful for integrating microarray and RNA-seq studies, they were typically designed for small datasets and face limitations when applied to large-scale data. To address this, we developed the hybrid Stouffer (hStouffer) method, an improved approach optimized for the integration of large RNA-seq datasets, enabling more reliable identification of differentially expressed genes across diverse studies.

<br>

## Prerequisites

To run this project, make sure you have the following Python packages installed (or newer):

| Package  | Version (>=) | Description                                                                 |
|----------|--------------|-----------------------------------------------------------------------------|
| `numpy`  | >=2.0.0       | Core library for numerical computing and array operations                  |
| `pandas` | >=2.2.2       | Powerful data structures for data analysis and manipulation                 |
| `joblib` | >=1.4.2       | Lightweight library for parallel processing and function caching            |
| `tqdm`   | >=4.66.4      | Fast, extensible progress bar for loops and CLI                            |
| `scipy`  | >=1.14.0      | Scientific computing library including statistics, optimization, and more   |

<br>
Prerequisities can be simply install using requirements.txt

<pre lang="markdown"> pip install -r requirements.txt </pre>
<br>

## Usage
| Argument | Required | Default | Description |
|----------|----------|---------|-------------|
| `-d`     | ✔        | –       | Directory containing DESeq2 result files |
| `-p`     | ✘        | `False` | p-value threshold (e.g., 3 = 1e-3) |
| `-o`     | ✔        | –       | Output file name |
| `-m`     | ✘        | `False` | Group size |
| `-r`     | ✘        | `1000`  | Number of repetitions |
| `-n`     | ✘        | `False` | Apply scaling (`True` / `False`) |
| `-t`     | ✘        | `"stouffer"` | Method for p-value combination:<br>`fisher`, `pearson`, `tippett`, `stouffer`, `mudholkar_george`, `median`, `percentile_70` |
| `-c`     | ✘        | `32`    | Number of CPU cores |
| `-w`     | ✘        | `True`  | Apply p-value capping (`True` / `False`) |
| `-l`     | ✘        | `False` | Max or min cutoff threshold |
<br>
*Datasets must be prepared in the format of DESeq2 output results <br><br>

Simple example :
<pre lang="markdown"> hStouffer.py -d dataset_directory -o meta_dataset_fat </pre>

More detail example :
<pre lang="markdown"> hStouffer.py -d dataset_directory -o meta_dataset_max -r 10000 -c 64 -l max </pre>

<br>


## Contact
Daehee Kim : rlarl0240@knu.ac.kr <br>
Jun-yeong Lee : junyeong@knu.ac.kr <br>
Laboratory of Genome Architecture and Regulation, School of Life Sciences, College of Natural Sciences, Kyungpook National University
