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
- `visualization`: optional CairoMakie static scientific plots.
- `comparison`: correspondence-aware paired-site contact changes and
  annotated-state summaries.

`AtomRecord.position` is a fixed-size vector; coordinate extraction produces a
contiguous `3 × N` matrix. Binding-site and contact queries build KD trees and
avoid full ribosome-by-ligand distance matrices. Trajectories are represented in
memory for deterministic core algorithms; external readers should stream or
slice before construction for very large simulations.

The MDAnalysis bridge and CairoMakie plotting backend are Julia package
extensions. Static workflows do not load PythonCall or CairoMakie, and a
missing optional dependency produces a readable error rather than a partial
analysis.

Polymer residue changes are compared only through sequence-derived
correspondence. Non-polymer site components such as ions or solvent are paired
only when their complete residue identifiers and classified molecule types
agree, and that weaker correspondence is labeled separately in output tables.
