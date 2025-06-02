# AIRR-Pipeline

This project demonstrates a minimal immune repertoire analysis pipeline using Change-O formatted B-cell receptor (BCR) data. It performs clone size analysis and calculates key diversity metrics like **Shannon** and **Simpson’s reciprocal index**, with visualizations per sample. Designed as a lightweight, reproducible example to support immune repertoire data analysis in R.

---

## Features

- Clone size distribution analysis
- Shannon diversity index (overall & per sample)
- Simpson’s reciprocal diversity index (per sample)
- Modular R script
- Output plots and data tables
- Easy to extend or Dockerize

---

## Input

The pipeline uses a single Change-O formatted file (`example_db.tsv`) containing annotated BCR sequences. This file is sourced from the `alakazam` R package example dataset.

| Column | Description |
|--------|-------------|
| `SAMPLE` | Sample ID |
| `CLONE` | Clone assignment ID |
| `SEQUENCE_IMGT` | Full BCR sequence |
| `DUPCOUNT` | Number of reads for this sequence |
| ... | Other immune receptor annotations |

---

## 🚀 Quick Start

### 1. Clone the repo

```bash
git clone https://github.com/PSBioinfo/AIRR-Pipeline.git
cd AIRR-Pipeline
```
### 2. Install R Dependencies

```r
install.packages(c("dplyr", "ggplot2", "alakazam"))
```
### 3. Run the Analysis
```bash
Rscript scripts/clone_analysis.R
```

## Output

| File                           | Description                 |
| ------------------------------ | --------------------------- |
| `output/clone_histogram.png`   | Histogram of clone sizes    |
| `output/shannon_index.txt`     | Overall Shannon index       |
| `output/shannon_by_sample.tsv` | Per-sample Shannon indices  |
| `output/shannon_by_sample.png` | Bar plot of Shannon indices |
| `output/simpson_by_sample.tsv` | Per-sample Simpson indices  |
| `output/simpson_by_sample.png` | Bar plot of Simpson indices |

## Project Structure

AIRR-ToyPipeline/
├── data/
│   └── example_db.tsv          # Input BCR data
├── scripts/
│   └── clone_analysis.R        # Main analysis script
├── output/                     # All result files and plots
├── README.md                   # This file!

## Background

This project uses immune repertoire sequencing (AIRR-seq) data to analyze clonal diversity -- a key feature in immunological studies, cancer, and vaccine response. It mimics basic functionality from the Immcantation framework using base R packages.

## License

This project is open source and freely available under the MIT License.

