# Deep-EccScan

**Deep-EccScan** is a deep learning-based tool for detecting **eccDNA** (extrachromosomal DNA) from short-read sequencing data. This repository contains the evaluation scripts, benchmark datasets, and R notebooks used to compare Deep-EccScan against four existing eccDNA detection tools: **CircleMap**, **Circle_finder**, **ecc_finder**, and **ECCsplorer**.

## Overview

eccDNAs are circular DNA elements that exist independently of chromosomes and play important roles in cancer, gene amplification, and disease. Accurate detection of eccDNAs from long-read sequencing data remains challenging. Deep-EccScan leverages deep learning to improve detection sensitivity, precision, and computational efficiency compared to existing methods.

This project provides:
- A Python script for benchmarking detection performance using simulated ground truth
- R notebooks for generating publication-quality figures comparing tool performance
- Simulated positive/negative BED datasets for evaluation
- Aggregated metrics and timing data across multiple datasets and sequencing depths

## Repository Structure

```
├── Script/
│   └── comparison_bed.py        # Python benchmarking script
├── NoteBook/
│   ├── Figure3.ipynb            # Performance comparison plots (F1, Accuracy, Precision, etc.)
│   ├── Figure4B.ipynb           # Time and memory consumption comparison
│   ├── Figure5.ipynb            # Venn diagram comparison across 4 tools on real samples
│   └── data/
│       ├── Merged_Metric.tsv    # Aggregated detection metrics (F1, Accuracy, Precision, etc.)
│       └── Merged_Log.txt       # Runtime and memory usage logs
└── SimulatedData/
    ├── dataset.pos.bed         # Simulated positive (true eccDNA) regions
    └── dataset.neg.bed         # Simulated negative (non-eccDNA) regions
```

## Comparison Tools

| Tool | Type | Description |
|------|------|-------------|
| **Deep-EccScan** | Deep Learning | Novel eccDNA detector using deep learning (this project) |
| **CircleMap** | Split-read | Circular DNA mapper for long-read sequencing |
| **Circle_finder** | Algorithmic | Circular element finder for long-read data |
| **ecc_finder** | Algorithmic | Dedicated eccDNA detection from long reads |
| **ECCsplorer** | Algorithmic | eccDNA explorer for long-read sequencing |

## Evaluation Metrics

Performance is evaluated on simulated datasets at multiple sequencing depths (5X–50X) across 7 datasets, plus real biological samples. Key metrics include:

| Metric | Description |
|--------|-------------|
| **F1 Score** | Harmonic mean of precision and recall |
| **Accuracy** | Overall classification accuracy |
| **Precision** | Positive predictive value |
| **Sensitivity (Recall)** | True positive rate |
| **Specificity** | True negative rate |
| **Base Pair Difference** | Average junction breakpoint deviation (lower is better) |
| **Redundancy** | Ratio of detected to true positives (lower is better) |
| **Runtime** | Total processing time in minutes |
| **Peak Memory** | Maximum memory usage in GB |

## Scripts

### `Script/comparison_bed.py`

Evaluates a detection tool's output (BED format) against simulated ground truth using MUMmer for sequence-level alignment.

**Usage:**

```bash
python comparison_bed.py \
    -1 <positive_truth.bed> \
    -2 <negative_truth.bed> \
    -3 <detection_results.bed> \
    -o <output_directory> \
    [-t threads] [-c similarity_threshold]
```

**Arguments:**
- `-1, --bed1`: Positive ground truth BED file
- `-2, --bed2`: Negative ground truth BED file
- `-3, --bed3`: Tool detection results BED file
- `-o, --output`: Output directory for results
- `-t, --threads`: Number of parallel threads (default: 64)
- `-c, --threshold`: Similarity threshold percentage (default: 50)

**Outputs:**
- Matched coordinate files from MUMmer alignment
- Per-tool detection metrics (TP, FP, TN, FN, F1, etc.)
- Unmatched positive/negative BED files
- `Metrics.tsv` with all evaluation metrics

### R Notebooks

| Notebook | Description |
|----------|-------------|
| `Figure3.ipynb` | Scatter/line plots comparing F1 score, accuracy, precision, sensitivity, specificity, base pair difference, and redundancy across tools and sequencing depths |
| `Figure4B.ipynb` | Bar charts comparing runtime (minutes) and peak memory (GB) across tools and sample types |
| `Figure5.ipynb` | Venn diagrams showing overlap of detected eccDNAs between Deep-EccScan, CircleMap, Circle_finder, and ecc_finder across real samples (Hela, Mouse Muscle, etc.) |

## Data

### Simulated Data

`SimulatedData/` contains BED files of simulated eccDNA regions used as ground truth for benchmarking. These represent known positive (eccDNA) and negative (non-eccDNA) regions.

### Evaluation Results

- **`Merged_Metric.tsv`**: Detection performance metrics for all tools across 7 simulated datasets at depths from 5X to 50X, plus real samples. Contains TP, FP, TN, FN, Accuracy, Precision, Sensitivity, Specificity, F1 Score, False Positive/Negative Rates, Base Pair Difference, and Redundancy.
- **`Merged_Log.txt`**: Computational resource usage (runtime in minutes, peak memory in GB) for all tools across simulated datasets and real samples.

## Dependencies

### Python
- Python 3.x
- pandas
- MUMmer4 (nucmer, show-coords)
- bedtools (getfasta)

### R
- ggplot2
- ggpubr
- patchwork
- ggsci
- ggvenn
- dplyr / tidyverse
- readr
- RColorBrewer
- ggprism
- scales
- ggbeeswarm

## Citation

If you use Deep-EccScan in your research, please cite the corresponding publication.

## License

The project is licensed under the GNU General Public License.
