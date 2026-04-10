#!/usr/bin/env python3
"""
Demo for SpatialTranscript.

Generates synthetic spatial transcriptomics data to demonstrate the workflow.
Expected output: spatial_domains.csv + cell_type_scores.csv + interactive HTML plots
"""

import json
import numpy as np
import pandas as pd
from pathlib import Path

def generate_demo_data():
    """Generate synthetic Visium-like dataset."""
    
    np.random.seed(42)
    
    n_spots = 500
    n_genes = 200
    gene_names = [f"GENE_{i:04d}" for i in range(n_genes)]
    spot_ids = [f"SPOT_{i:04d}" for i in range(n_spots)]
    
    # Synthetic expression: 4 spatial domains with different expression programs
    x = np.random.uniform(0, 10, n_spots)
    y = np.random.uniform(0, 10, n_spots)
    
    # Assign spots to domains (spatial blocks)
    domains = ((x > 5).astype(int) + 2 * (y > 5).astype(int)).astype(str)  # 4 domains: 0,1,2,3
    
    # Generate expression with domain-specific means
    expr = np.zeros((n_spots, n_genes))
    domain_means = {
        "0": np.random.randn(n_genes) * 0.5 + 2,
        "1": np.random.randn(n_genes) * 0.5 + 3,
        "2": np.random.randn(n_genes) * 0.5 + 1,
        "3": np.random.randn(n_genes) * 0.5 + 2.5,
    }
    
    for i, d in enumerate(domains):
        expr[i] = np.random.poisson(np.abs(domain_means[d]))
    
    expr_df = pd.DataFrame(expr, index=spot_ids, columns=gene_names)
    expr_df.to_csv("/tmp/SpatialTranscript_demo_expression.csv")
    
    # Coordinates
    coords = pd.DataFrame({
        "spot_id": spot_ids,
        "x": x,
        "y": y
    }).set_index("spot_id")
    coords.to_csv("/tmp/SpatialTranscript_demo_coords.csv")
    
    # Marker genes for cell type deconvolution
    marker_genes = {
        "Neuron": ["GENE_0010", "GENE_0050", "GENE_0100"],
        "Glia": ["GENE_0020", "GENE_0060", "GENE_0110"],
        "Immune": ["GENE_0030", "GENE_0070", "GENE_0120"],
        "Epithelial": ["GENE_0040", "GENE_0080", "GENE_0130"],
    }
    
    with open("/tmp/SpatialTranscript_demo_markers.json", "w") as f:
        json.dump(marker_genes, f)
    
    print("[Demo] Generated:")
    print(f"  - /tmp/SpatialTranscript_demo_expression.csv ({n_spots} spots x {n_genes} genes)")
    print(f"  - /tmp/SpatialTranscript_demo_coords.csv ({n_spots} spots)")
    print(f"  - /tmp/SpatialTranscript_demo_markers.json (4 cell types)")
    
    return (
        "/tmp/SpatialTranscript_demo_expression.csv",
        "/tmp/SpatialTranscript_demo_coords.csv",
        "/tmp/SpatialTranscript_demo_markers.json"
    )


if __name__ == "__main__":
    expression_csv, coord_csv, marker_genes = generate_demo_data()
    
    print("\n[Demo] Expected output structure:")
    print("  spatial_results/")
    print("  ├── spatial_domains.csv      # spot_id -> domain label")
    print("  ├── cell_type_scores.csv     # spot_id -> cell type probabilities")
    print("  ├── spatial_domains.html      # interactive scatter plot")
    print("  └── gene_expression.html      # interactive heatmap")
    
    print("\n[Demo] Run with:")
    print("  python SKELETON.py \\")
    print("    --expression-csv /tmp/SpatialTranscript_demo_expression.csv \\")
    print("    --coord-csv /tmp/SpatialTranscript_demo_coords.csv \\")
    print("    --marker-genes /tmp/SpatialTranscript_demo_markers.json \\")
    print("    --output-dir spatial_results")
