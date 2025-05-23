# hStouffer

## Background
Advances in high-throughput sequencing have led to an exponential increase in publicly available RNA-seq datasets, enhancing opportunities for large-scale transcriptomic analyses. Compared to small-scale studies, large-scale RNA-seq meta-analyses offer greater robustness, reduce technical noise, and improve detection of subtle or novel gene expression patterns. However, while traditional meta-analysis methods—particularly p-value-based approaches like Stouffer's—have proven useful for integrating microarray and RNA-seq studies, they were typically designed for small datasets and face limitations when applied to large-scale data. To address this, we developed the hybrid Stouffer (hStouffer) method, an improved approach optimized for the integration of large RNA-seq datasets, enabling more reliable identification of differentially expressed genes across diverse studies.

## Prerequisites

To run this project, make sure you have the following Python packages installed (or newer):

| Package  | Version (>=) | Description                                                                 |
|----------|--------------|-----------------------------------------------------------------------------|
| `numpy`  | >=2.0.0       | Core library for numerical computing and array operations                  |
| `pandas` | >=2.2.2       | Powerful data structures for data analysis and manipulation                 |
| `joblib` | >=1.4.2       | Lightweight library for parallel processing and function caching            |
| `tqdm`   | >=4.66.4      | Fast, extensible progress bar for loops and CLI                            |
| `scipy`  | >=1.14.0      | Scientific computing library including statistics, optimization, and more   |
