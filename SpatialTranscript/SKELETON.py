#!/usr/bin/env python3
"""
SpatialTranscript: Agent-executable skill for spatial transcriptomics data analysis.

GAP IDENTIFIED: No existing claw4s submission handles spatial transcriptomics data.
Existing tools cover bulk RNA-seq, scRNA-seq, and single-cell multiome, but the
spatial dimension (coordinates + expression) is entirely missing.

This implementation provides:
- Spatial coordinate + expression matrix loading (standard formats)
- Spatial neighborhood graph construction
- PCA-based dimensionality reduction
- Leiden-like clustering for domain detection
- Cell-type deconvolution via marker genes
- Interactive visualization
- Spatial co-expression module (Moran's I, Geary's C)
"""

import json
import warnings
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd


# =============================================================================
# FORMAT PARSERS
# =============================================================================

def load_visium(visium_dir: str) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Load 10x Visium data from Space Ranger output directory."""
    try:
        import scanpy as sc
    except ImportError:
        warnings.warn("scanpy not installed, falling back to h5py parser")
        return _load_visium_h5py(visium_dir)
    
    adata = sc.read_visium(visium_dir)
    counts_df = pd.DataFrame(
        adata.X.toarray() if hasattr(adata.X, 'toarray') else adata.X,
        index=adata.obs_names,
        columns=adata.var_names
    )
    
    spots_df = pd.DataFrame({
        'spot_id': adata.obs_names,
        'x': adata.obsm['spatial'][:, 0],
        'y': adata.obsm['spatial'][:, 1],
        'barcode': adata.obs_names
    }).set_index('spot_id')
    
    print(f"[SpatialTranscript] Loaded Visium: {counts_df.shape[0]} spots, {counts_df.shape[1]} genes")
    return counts_df, spots_df


def _load_visium_h5py(visium_dir: str) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Fallback parser using h5py + pandas for Visium data."""
    import h5py
    from scipy import sparse
    
    visium_path = Path(visium_dir)
    h5_file = visium_path / "filtered_feature_bc_matrix.h5"
    spots_file = visium_path / "spatial" / "tissue_positions.csv"
    
    with h5py.File(h5_file, 'r') as f:
        matrix = sparse.csr_matrix((
            f['matrix']['data'][:],
            f['matrix']['indices'][:],
            f['matrix']['indptr'][:]
        ), shape=f['matrix']['shape'][:])
        
        gene_names = [n.decode('utf-8') for n in f['matrix']['features']['name'][:]]
        barcodes = [b.decode('utf-8') for b in f['matrix']['barcodes'][:]]
    
    counts_df = pd.DataFrame(matrix.toarray(), index=barcodes, columns=gene_names)
    
    coords_df = pd.read_csv(spots_file, index_col=0)
    coords_df.columns = ['in_tissue', 'row', 'col', 'x', 'y']
    
    common = counts_df.index.intersection(coords_df.index)
    counts_df = counts_df.loc[common]
    coords_df = coords_df.loc[common]
    
    spots_df = pd.DataFrame({
        'spot_id': common,
        'x': coords_df['x'].values,
        'y': coords_df['y'].values,
        'barcode': common
    }).set_index('spot_id')
    
    print(f"[SpatialTranscript] Loaded Visium (h5py): {counts_df.shape[0]} spots, {counts_df.shape[1]} genes")
    return counts_df, spots_df


def load_merfISH(merfish_file: str) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Load MERFISH data from CSV/Parquet."""
    path = Path(merfish_file)
    
    if path.suffix == '.parquet':
        df = pd.read_parquet(merfish_file)
    else:
        df = pd.read_csv(merfish_file, index_col=0)
    
    coord_cols = [c for c in df.columns if c.lower() in ('x', 'y', 'z', 'cell_id', 'cellid')]
    gene_cols = [c for c in df.columns if c not in coord_cols]
    
    has_z = 'z' in df.columns or 'Z' in df.columns
    
    if has_z:
        coordinates_df = pd.DataFrame({
            'cell_id': df.index,
            'x': df['x'].values if 'x' in df.columns else df['X'].values,
            'y': df['y'].values if 'y' in df.columns else df['Y'].values,
            'z': df['z'].values if 'z' in df.columns else df['Z'].values,
        }).set_index('cell_id')
    else:
        coordinates_df = pd.DataFrame({
            'cell_id': df.index,
            'x': df['x'].values if 'x' in df.columns else df['X'].values,
            'y': df['y'].values if 'y' in df.columns else df['Y'].values,
        }).set_index('cell_id')
    
    expression_df = df[gene_cols] if gene_cols else pd.DataFrame()
    
    print(f"[SpatialTranscript] Loaded MERFISH: {expression_df.shape[0]} cells, {expression_df.shape[1]} genes")
    return expression_df, coordinates_df


def load_generic_csv(
    expression_csv: str,
    coord_csv: str,
    expression_format: str = "spots_as_rows"
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Load generic spatial transcriptomics data from CSV files."""
    expr = pd.read_csv(expression_csv, index_col=0)
    if expression_format == "genes_as_rows":
        expr = expr.T
    
    coords = pd.read_csv(coord_csv)
    coords = coords.set_index(coords.columns[0])
    
    if 'x' not in coords.columns and 'y' not in coords.columns:
        coord_cols = [c for c in coords.columns if c not in ('spot_id', 'id', 'barcode')]
        if len(coord_cols) >= 2:
            coords = coords.rename(columns={coord_cols[0]: 'x', coord_cols[1]: 'y'})
    
    common_ids = expr.index.intersection(coords.index)
    if len(common_ids) == 0:
        warnings.warn(f"No common IDs found. Using all spots.")
        common_ids = list(expr.index)
    
    expr = expr.loc[common_ids]
    coords = coords.loc[common_ids]
    
    print(f"[SpatialTranscript] Loaded {len(common_ids)} spots, {expr.shape[1]} genes")
    return expr, coords


# =============================================================================
# SPATIAL ANALYSIS
# =============================================================================

def build_spatial_knn_graph(
    coords: pd.DataFrame,
    n_neighbors: int = 15
) -> np.ndarray:
    """Build k-NN graph from spatial coordinates using efficient ball_tree."""
    from sklearn.neighbors import NearestNeighbors
    
    coord_arr = coords[['x', 'y']].values if 'x' in coords.columns else coords.values
    
    nbrs = NearestNeighbors(n_neighbors=min(n_neighbors + 1, len(coord_arr)), algorithm='ball_tree').fit(coord_arr)
    distances, indices = nbrs.kneighbors(coord_arr)
    
    n = len(coord_arr)
    adj = np.zeros((n, n))
    for i in range(n):
        for neighbor in indices[i][1:]:  # skip self
            adj[i, neighbor] = 1.0
    
    return adj


def compute_spatial_autocorrelation(
    expression: pd.Series,
    knn_graph: np.ndarray
) -> Tuple[float, float]:
    """Compute Moran's I for spatial autocorrelation of a gene."""
    n = len(expression)
    x = expression.values - expression.mean()
    x_std = expression.std()
    
    if x_std == 0:
        return 0.0, 1.0
    
    w_sum = knn_graph.sum()
    if w_sum == 0:
        return 0.0, 1.0
    
    numerator = np.sum(knn_graph * np.outer(x, x))
    denominator = np.sum(x ** 2)
    
    morans_i = (n / w_sum) * (numerator / denominator)
    
    # Permutation test
    n_perms = 99
    null_distribution = []
    expr_arr = expression.values.copy()
    
    for _ in range(n_perms):
        np.random.shuffle(expr_arr)
        null_numerator = np.sum(knn_graph * np.outer(expr_arr - expr_arr.mean(), expr_arr - expr_arr.mean()))
        null_denom = np.sum((expr_arr - expr_arr.mean()) ** 2)
        if null_denom > 0:
            null_dist_val = (n / w_sum) * (null_numerator / null_denom)
            null_distribution.append(null_dist_val)
    
    null_distribution = np.array(null_distribution)
    pvalue = np.mean(null_distribution >= morans_i) + 1 / (n_perms + 1)
    
    return float(morans_i), float(pvalue)


def compute_gearys_c(
    expression: pd.Series,
    knn_graph: np.ndarray
) -> Tuple[float, float]:
    """Compute Geary's C for spatial autocorrelation."""
    n = len(expression)
    x = expression.values
    x_mean = x.mean()
    
    w_sum = knn_graph.sum()
    if w_sum == 0:
        return 1.0, 1.0
    
    numerator = 0.0
    for i in range(n):
        for j in range(n):
            if knn_graph[i, j] > 0:
                numerator += (x[i] - x[j]) ** 2
    
    denominator = np.sum((x - x_mean) ** 2)
    
    if denominator == 0:
        return 1.0, 1.0
    
    gearys_c = ((n - 1) / (2 * w_sum)) * (numerator / denominator)
    
    # Permutation test
    n_perms = 99
    null_distribution = []
    expr_arr = x.copy()
    
    for _ in range(n_perms):
        np.random.shuffle(expr_arr)
        null_num = 0.0
        for i in range(n):
            for j in range(n):
                if knn_graph[i, j] > 0:
                    null_num += (expr_arr[i] - expr_arr[j]) ** 2
        null_denom = np.sum((expr_arr - expr_arr.mean()) ** 2)
        if null_denom > 0:
            null_dist_val = ((n - 1) / (2 * w_sum)) * (null_num / null_denom)
            null_distribution.append(null_dist_val)
    
    null_distribution = np.array(null_distribution)
    pvalue = np.mean(null_distribution <= gearys_c) + 1 / (n_perms + 1)
    
    return float(gearys_c), float(pvalue)


def detect_spatial_domains(
    expression: pd.DataFrame,
    coords: pd.DataFrame,
    n_neighbors: int = 15,
    resolution: float = 1.0,
    n_components: int = 20
) -> pd.Series:
    """Main domain detection: PCA -> KMeans clustering on combined expression+spatial features."""
    from sklearn.decomposition import PCA
    from sklearn.cluster import KMeans
    
    # Normalize expression first (CPM-like)
    row_sums = expression.sum(axis=1)
    row_sums[row_sums == 0] = 1
    expr_norm = expression.div(row_sums, axis=0) * 1e6
    expr_norm = np.log1p(expr_norm).fillna(0)
    
    # PCA on expression matrix
    n_pcs = max(1, min(n_components, expr_norm.shape[1] - 1, expr_norm.shape[0] - 1))
    pca = PCA(n_components=n_pcs)
    pca_emb = pca.fit_transform(expr_norm)
    
    print(f"[SpatialTranscript] PCA: {pca_emb.shape}")
    
    # Normalize spatial coordinates
    coord_arr = coords[['x', 'y']].values if 'x' in coords.columns else coords.values
    coord_mean = coord_arr.mean(axis=0)
    coord_std = coord_arr.std(axis=0) + 1e-10
    coords_norm = (coord_arr - coord_mean) / coord_std
    
    # Combine: expression PCA + spatial coords (weighted)
    combined = np.hstack([pca_emb * 0.8, coords_norm * 0.2])
    
    # KMeans clustering
    n_clusters = max(2, int(4 * resolution))
    kmeans = KMeans(n_clusters=n_clusters, random_state=42, n_init=10)
    domains = kmeans.fit_predict(combined)
    
    # Return as Series
    domains_series = pd.Series(domains, index=expression.index, name='domain')
    
    print(f"[SpatialTranscript] Detected {domains_series.nunique()} spatial domains")
    return domains_series


def deconvolve_cell_types(
    expression: pd.DataFrame,
    marker_genes: Dict[str, List[str]],
    n_neighbors: int = 15
) -> pd.DataFrame:
    """Deconvolve cell types using marker gene enrichment per spot."""
    # Normalize expression
    row_sums = expression.sum(axis=1)
    row_sums[row_sums == 0] = 1
    expr_norm = expression.div(row_sums, axis=0) * 1e6
    expr_norm = np.log1p(expr_norm).fillna(0)
    
    scores = {}
    
    for cell_type, markers in marker_genes.items():
        available_markers = [m for m in markers if m in expr_norm.columns]
        if not available_markers:
            warnings.warn(f"No markers found for cell type: {cell_type}")
            scores[cell_type] = np.zeros(len(expression))
            continue
        
        # Mean expression of markers
        marker_expr = expr_norm[available_markers].mean(axis=1)
        scores[cell_type] = marker_expr.values
    
    cell_type_scores = pd.DataFrame(scores, index=expression.index)
    
    # Normalize per spot
    row_sums = cell_type_scores.sum(axis=1)
    row_sums[row_sums == 0] = 1
    cell_type_scores = cell_type_scores.div(row_sums, axis=0)
    
    print(f"[SpatialTranscript] Cell type deconvolution: {cell_type_scores.shape}")
    return cell_type_scores


# =============================================================================
# VISUALIZATION
# =============================================================================

def plot_spatial_domains(
    coords: pd.DataFrame,
    domains: pd.Series,
    output_path: str = "spatial_domains.html"
):
    """Generate interactive HTML scatter plot of spatial domains."""
    try:
        import plotly.express as px
        import plotly.graph_objects as go
    except ImportError:
        warnings.warn("plotly not installed, skipping HTML visualization")
        return
    
    x = coords['x'].values if 'x' in coords.columns else coords.iloc[:, 0].values
    y = coords['y'].values if 'y' in coords.columns else coords.iloc[:, 1].values
    
    fig = go.Figure()
    colors = px.colors.qualitative.Set1[:domains.nunique()]
    
    for i, domain in enumerate(sorted(domains.unique())):
        mask = domains == domain
        fig.add_trace(go.Scatter(
            x=x[mask], y=y[mask], mode='markers',
            marker=dict(size=8, color=colors[i % len(colors)], opacity=0.7),
            name=f'Domain {domain}',
            text=[f'Spot: {idx}<br>Domain: {domain}' for idx in domains[mask].index],
            hoverinfo='text+name'
        ))
    
    fig.update_layout(
        title='Spatial Domain Map',
        xaxis_title='X Coordinate',
        yaxis_title='Y Coordinate',
        legend_title='Domains',
        width=800, height=700,
        template='plotly_white'
    )
    
    fig.write_html(output_path)
    print(f"[SpatialTranscript] Saved: {output_path}")


def plot_gene_expression(
    coords: pd.DataFrame,
    expression: pd.Series,
    output_path: str = "gene_expression.html",
    title: str = "Gene Expression"
):
    """Generate interactive HTML heatmap of gene expression on spatial coordinates."""
    try:
        import plotly.graph_objects as go
    except ImportError:
        warnings.warn("plotly not installed, skipping HTML visualization")
        return
    
    x = coords['x'].values if 'x' in coords.columns else coords.iloc[:, 0].values
    y = coords['y'].values if 'y' in coords.columns else coords.iloc[:, 1].values
    
    expr_values = expression.values
    expr_normalized = (expr_values - expr_values.min()) / (expr_values.max() - expr_values.min() + 1e-10)
    
    fig = go.Figure(data=go.Scatter(
        x=x, y=y, mode='markers',
        marker=dict(size=8, color=expr_normalized, colorscale='Viridis',
                   colorbar=dict(title='Expression'), opacity=0.8),
        text=[f'Spot: {idx}<br>Expr: {v:.3f}' for idx, v in zip(expression.index, expr_values)],
        hoverinfo='text'
    ))
    
    fig.update_layout(
        title=title,
        xaxis_title='X Coordinate',
        yaxis_title='Y Coordinate',
        width=800, height=700,
        template='plotly_white'
    )
    
    fig.write_html(output_path)
    print(f"[SpatialTranscript] Saved: {output_path}")


def plot_cell_type_pie(
    cell_type_scores: pd.DataFrame,
    output_path: str = "cell_type_pie.html"
):
    """Generate interactive HTML pie chart of cell type proportions."""
    try:
        import plotly.express as px
    except ImportError:
        warnings.warn("plotly not installed, skipping HTML visualization")
        return
    
    proportions = cell_type_scores.mean()
    
    fig = px.pie(
        values=proportions, names=proportions.index,
        title='Average Cell Type Composition',
        color=proportions.index,
        color_discrete_map={
            'Neuron': '#1f77b4', 'Glia': '#ff7f0e',
            'Immune': '#2ca02c', 'Epithelial': '#d62728'
        }
    )
    
    fig.update_layout(width=600, height=500)
    fig.write_html(output_path)
    print(f"[SpatialTranscript] Saved: {output_path}")


# =============================================================================
# MAIN PIPELINE
# =============================================================================

def analyze_spatial_data(
    expression_csv: str,
    coord_csv: str,
    marker_genes_file: Optional[str] = None,
    output_dir: str = "spatial_results",
    n_neighbors: int = 15,
    resolution: float = 1.0
) -> Dict:
    """Complete spatial transcriptomics analysis pipeline."""
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    # Load data
    expr, coords = load_generic_csv(expression_csv, coord_csv)
    
    # Detect domains
    domains = detect_spatial_domains(expr, coords, n_neighbors, resolution)
    
    # Deconvolve cell types if markers provided
    cell_type_scores = None
    if marker_genes_file:
        with open(marker_genes_file) as f:
            marker_genes = json.load(f)
        cell_type_scores = deconvolve_cell_types(expr, marker_genes, n_neighbors)
    
    # Identify HVGs and compute spatial autocorrelation
    row_sums = expr.sum(axis=1)
    row_sums[row_sums == 0] = 1
    expr_norm = expr.div(row_sums, axis=0) * 1e6
    expr_norm = np.log1p(expr_norm).fillna(0)
    
    gene_vars = expr_norm.var()
    top_genes = gene_vars.nlargest(10).index.tolist()
    
    spatial_autocorr = {}
    knn_graph = build_spatial_knn_graph(coords, n_neighbors)
    
    for gene in top_genes:
        morans_i, pval = compute_spatial_autocorrelation(expr_norm[gene], knn_graph)
        gearys_c, cpval = compute_gearys_c(expr_norm[gene], knn_graph)
        spatial_autocorr[gene] = {
            'morans_i': morans_i, 'morans_pval': pval,
            'gearys_c': gearys_c, 'gearys_pval': cpval
        }
    
    # Write outputs
    domains.to_csv(output_path / "spatial_domains.csv", header=['domain'])
    
    if cell_type_scores is not None:
        cell_type_scores.to_csv(output_path / "cell_type_scores.csv")
    
    with open(output_path / "spatial_autocorrelation.json", 'w') as f:
        json.dump(spatial_autocorr, f, indent=2)
    
    # Generate visualizations
    plot_spatial_domains(coords, domains, str(output_path / "spatial_domains.html"))
    
    if len(top_genes) > 0:
        top_gene = top_genes[0]
        plot_gene_expression(
            coords, expr_norm[top_gene],
            str(output_path / f"gene_expression_{top_gene}.html"),
            title=f"Expression: {top_gene}"
        )
    
    if cell_type_scores is not None:
        plot_cell_type_pie(cell_type_scores, str(output_path / "cell_type_pie.html"))
    
    # Create combined dashboard
    try:
        import plotly.graph_objects as go
        from plotly.subplots import make_subplots
        
        x = coords['x'].values if 'x' in coords.columns else coords.iloc[:, 0].values
        y = coords['y'].values if 'y' in coords.columns else coords.iloc[:, 1].values
        
        fig = make_subplots(
            rows=2, cols=2,
            subplot_titles=('Spatial Domains', 'Cell Type Composition',
                          'Gene Expression (Top HVG)', 'Domain Distribution')
        )
        
        import plotly.express as px
        colors = px.colors.qualitative.Set1[:domains.nunique()]
        for i, domain in enumerate(sorted(domains.unique())):
            mask = domains == domain
            fig.add_trace(
                go.Scatter(x=x[mask], y=y[mask], mode='markers',
                          marker=dict(size=6, color=colors[i % len(colors)]),
                          name=f'Domain {domain}', showlegend=True),
                row=1, col=1
            )
        
        if cell_type_scores is not None:
            proportions = cell_type_scores.mean()
            fig.add_trace(
                go.Pie(labels=proportions.index, values=proportions.values, textinfo='label+percent'),
                row=1, col=2
            )
        
        if len(top_genes) > 0:
            top_gene = top_genes[0]
            expr_vals = expr_norm[top_gene].values
            fig.add_trace(
                go.Scatter(x=x, y=y, mode='markers',
                          marker=dict(size=6, color=expr_vals, colorscale='Viridis'), showlegend=False),
                row=2, col=1
            )
        
        domain_counts = domains.value_counts().sort_index()
        fig.add_trace(
            go.Bar(x=[f'Domain {d}' for d in domain_counts.index], y=domain_counts.values, showlegend=False),
            row=2, col=2
        )
        
        fig.update_layout(
            title='SpatialTranscript Analysis Dashboard',
            height=900, width=1100, template='plotly_white'
        )
        
        fig.write_html(str(output_path / "dashboard.html"))
        print(f"[SpatialTranscript] Saved: {output_path / 'dashboard.html'}")
        
    except Exception as e:
        warnings.warn(f"Dashboard generation failed: {e}")
    
    report = {
        "n_spots": len(expr),
        "n_genes": expr.shape[1],
        "n_domains": int(domains.nunique()),
        "domain_counts": {str(k): int(v) for k, v in domains.value_counts().to_dict().items()},
        "top_hvg": list(top_genes[:5]),
        "spatial_autocorrelation": spatial_autocorr,
        "output_dir": str(output_path)
    }
    
    print(f"[SpatialTranscript] Analysis complete: {report}")
    return report


# =============================================================================
# CLI
# =============================================================================

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="SpatialTranscript: Spatial transcriptomics analysis")
    parser.add_argument("--expression-csv", required=True, help="Expression matrix CSV")
    parser.add_argument("--coord-csv", required=True, help="Spatial coordinates CSV")
    parser.add_argument("--marker-genes", help="JSON file with marker genes per cell type")
    parser.add_argument("--output-dir", default="spatial_results")
    parser.add_argument("--n-neighbors", type=int, default=15)
    parser.add_argument("--resolution", type=float, default=1.0)
    
    args = parser.parse_args()
    
    report = analyze_spatial_data(
        args.expression_csv,
        args.coord_csv,
        args.marker_genes,
        args.output_dir,
        args.n_neighbors,
        args.resolution
    )
    print(json.dumps(report, indent=2))
