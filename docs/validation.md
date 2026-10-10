# Validation status

Passing unit tests are evidence about implementation behavior, not a claim of
complete scientific validation.

## Unit-tested functionality

- BioStructures parsing of a controlled PDB fixture; preservation of atom,
  residue, chain, model, occupancy, B-factor, element, and hetero fields.
- Exact synthetic KD-tree contact inclusion and monotonic cutoff behavior.
- Kabsch recovery of a known rotation/translation and near-zero fitted RMSD.
- Centroid, radius of gyration, extent, contact density, and pairwise counts for
  analytic coordinate sets.
- Explicit-hydrogen 180° D–H–A geometry.
- Nanometer-to-Angstrom conversion (factor 10).
- RMSD invariance to rigid rotation/translation, RMSF non-negativity, contact
  persistence, distance series, PCA shapes, and deterministic clustering.
- Autocorrelation normalization, effective sample-size bounds, bootstrap and
  moving-block-bootstrap behavior, and replica-level summaries.
- Immutable result-directory behavior and machine-readable exports.

## Cross-validated functionality

Structure parsing is delegated to BioStructures.jl 4.9.0, whose current API is
used directly. Atom counts are checked independently through the BioStructures
API in the public-data validation suite.

MDAnalysis 2.10.0 independently reproduced RMSD, aligned RMSF, and distance
values for a deterministic analytic trajectory. The maximum observed difference
from the Julia implementation was below `5e-8` Å. The checked-in reference also
uses MDAnalysis to analyze the official PDB-format RCSB 1J7T entry. At a 5 Å
cutoff, MDAnalysis and this package both report:

- 356 paromomycin–environment atom-pair contacts;
- 88 contacting non-ligand atoms;
- 18 contacting residues; and
- a minimum distance of approximately 2.51381 Å.

An end-to-end run on the compressed mmCIF RCSB 4V7U entry loaded 284,525 atoms
in 104 author chains, selected ERY at author chain BA/residue 3136, and produced
the standard metadata, summary, contact table, and visualization selections.
At 5 Å it reported 283 atom-pair contacts, 96 contacting atoms, and 12 contacting
residues. The input SHA-256 checksums and exact acquisition URLs are recorded in
`validation/public_structures.toml`; the data are downloaded separately rather
than bundled.

## Functionality requiring further validation

- A curated, versioned chemical-component database broader than the initial
  candidate mapping.
- Hydrogen-bond assignments against a protonation/topology-aware reference.
- Ionic/aromatic/metal candidates against quantum-chemical or curated
  interaction annotations.
- MDAnalysis XTC/TRR/DCD streaming against GROMACS tools across unit-cell and
  periodic-boundary conventions. The in-memory trajectory calculations are
  cross-validated, but production-format periodic-boundary behavior is not.
- Sequence correspondence across divergent species and paralogous chains.
- Solvent-accessibility and pocket-occupancy algorithms (not implemented).
- Time-lagged dimensionality reduction (not implemented).
- Interactive molecular rendering (exports are provided instead).
- Independent reproduction of publication figures from a real public ribosome
  dataset. Real-structure tables and selections are covered by the public-data
  suite; figure comparison remains qualitative.

## Experimental or exploratory methods

Cartesian PCA, deterministic Lloyd clustering, aromatic ring orientation, and
metal-mediated geometry are exploratory. They require domain-specific feature
selection, convergence checks, sensitivity analysis, and external validation.

## Interpretation boundaries

A geometric contact indicates proximity. A hydrogen-bond geometry indicates
compatibility with user-supplied topology. Interaction persistence indicates
cutoff occupancy over sampled frames. Structural stability describes limited
coordinate observables. None alone establishes affinity, binding free energy,
causality, clinical efficacy, or resistance mechanism.
