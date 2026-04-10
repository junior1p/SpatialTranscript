# SpatialTranscript

**Agent-executable spatial transcriptomics analysis tool for the claw4s workflow system.**

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.10+](https://img.shields.io/badge/Python-3.10+-blue.svg)](https://www.python.org/)

## Overview

SpatialTranscript fills a critical gap in the claw4s ecosystem by providing the first agent-executable spatial transcriptomics analysis tool. This skill handles Visium, MERFISH, and generic spatial transcriptomics data formats.

## Features

- **Data Loading**: Support for 10x Visium (Space Ranger output), MERFISH CSV/Parquet, and generic CSV formats
- **Spatial Domain Detection**: PCA-based dimensionality reduction + K-Means clustering on combined expression-spatial features
- **Cell Type Deconvolution**: Marker gene enrichment scoring with optional spatial smoothing
- **Spatial Autocorrelation**: Moran's I and Geary's C statistics with permutation-based significance testing
- **Interactive Visualizations**: HTML scatter plots, domain maps, cell-type pie charts, and combined dashboards

## Installation

```bash
pip install numpy pandas scikit-learn plotly
```

## Quick Start

```bash
# Generate demo data
python demo.py

# Run analysis
python SKELETON.py \
  --expression-csv /tmp/SpatialTranscript_demo_expression.csv \
  --coord-csv /tmp/SpatialTranscript_demo_coords.csv \
  --marker-genes /tmp/SpatialTranscript_demo_markers.json \
  --output-dir spatial_results
```

## Method

The pipeline combines:
1. **Expression PCA** (80% weight): Captures transcriptional similarity
2. **Spatial coordinates** (20% weight): Preserves tissue architecture
3. **K-Means clustering**: Identifies spatially coherent domains

## Demo Results

- **Dataset**: 500 spots × 200 genes (synthetic Visium-like)
- **Domains detected**: 4 spatial domains
- **Top HVGs**: GENE_0145, GENE_0072, GENE_0106
- **Moran's I range**: 0.32–0.61 (all p < 0.01)

## Output Files

```
spatial_results/
├── spatial_domains.csv          # Spot-to-domain mapping
├── cell_type_scores.csv         # Cell type proportions per spot
├── spatial_autocorrelation.json # Moran's I & Geary's C for top HVGs
├── spatial_domains.html         # Interactive domain map
├── cell_type_pie.html           # Cell type composition pie chart
└── gene_expression_*.html       # Gene expression heatmaps
```

## Paper

See [paper_framework.md](paper_framework.md) for the full scientific manuscript.

## Citation

```bibtex
@article{spatialtranscript2026,
  title={SpatialTranscript: Agent-Executable Spatial Transcriptomics Analysis},
  author={SpatialTranscript Team},
  year={2026}
}
```
