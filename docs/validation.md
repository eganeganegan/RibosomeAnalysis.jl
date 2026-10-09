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
used directly. Kabsch and RMSD were checked against analytically generated
coordinates independent of structure parsing. No production ribosome dataset
is bundled, so no claim is made that every PDBx/mmCIF edge case has been tested.

## Functionality requiring further validation

- A curated, versioned chemical-component database broader than the initial
  candidate mapping.
- Hydrogen-bond assignments against a protonation/topology-aware reference.
- Ionic/aromatic/metal candidates against quantum-chemical or curated
  interaction annotations.
- MDAnalysis XTC/TRR/DCD integration against GROMACS tools across unit-cell and
  periodic-boundary conventions.
- Sequence correspondence across divergent species and paralogous chains.
- Solvent-accessibility and pocket-occupancy algorithms (not implemented).
- Time-lagged dimensionality reduction (not implemented).
- Interactive molecular rendering (exports are provided instead).
- Independent reproduction of figures and tables from a real public ribosome
  dataset.

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
