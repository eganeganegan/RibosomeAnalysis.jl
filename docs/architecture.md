# Architecture

The top-level `RibosomeAnalysis` module includes focused source files:

- `structures`: stable atom/residue/chain records and selection helpers.
- `io`: BioStructures-backed parsing and reproducible result serialization.
- `binding`: ligand metadata and KD-tree site detection.
- `interactions`: atomic, hydrogen-bond, ionic, aromatic, and metal geometry.
- `geometry`: site descriptors without unsupported pocket-volume claims.
- `alignment`: sequence correspondence and Kabsch rigid-body fitting.
- `dynamics`: unit-aware trajectory observables, PCA, and clustering.
- `statistics`: independent and correlated-sample uncertainty methods.
- `visualization`: CairoMakie static scientific plots.
- `comparison`: paired-site and annotated-state summaries.

`AtomRecord.position` is a fixed-size vector; coordinate extraction produces a
contiguous `3 × N` matrix. Binding-site and contact queries build KD trees and
avoid full ribosome-by-ligand distance matrices. Trajectories are represented in
memory for deterministic core algorithms; external readers should stream or
slice before construction for very large simulations.

The MDAnalysis bridge is a Julia package extension. Static workflows do not load
PythonCall, and failure to install the optional dependency produces a readable
error rather than a partial analysis.
