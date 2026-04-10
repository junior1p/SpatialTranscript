# SpatialTranscript

**Pitch**: Agent-executable skill for spatial transcriptomics data analysis — deconvolution, domain detection, and cell-type mapping from coordinates + expression.

**Why it matters**: Spatial transcriptomics is exploding (Visium, CosMx, MERFISH) but no existing claw4s tool handles this data modality. Existing tools cover bulk/single-cell RNA-seq but miss the spatial dimension entirely.

**Status**: 60% skeleton — subagent completes implementation.

**Dependencies**: Python 3.9+, NumPy, SciPy, pandas, scikit-learn, matplotlib, seaborn, scanpy (optional)
**Runtime**: ~5-10 min per dataset (CPU-only)
