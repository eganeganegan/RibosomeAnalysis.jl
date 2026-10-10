# Ribosome structural analysis in Julia

[![CI](https://github.com/eganeganegan/RibosomeAnalysis.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/eganeganegan/RibosomeAnalysis.jl/actions/workflows/ci.yml)
[![Documentation](https://github.com/eganeganegan/RibosomeAnalysis.jl/actions/workflows/documentation.yml/badge.svg)](https://eganeganegan.github.io/RibosomeAnalysis.jl/stable/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

This repository is a Julia package and command-line toolkit for reproducible
analysis of ribosomal structures, antibiotic binding sites, structural
correspondence, geometric interactions, and molecular-dynamics trajectories.
It is intended for bacterial 30S/50S work, but its core geometry is not tied to
one organism or ribosomal state.

The package deliberately keeps distinct concepts separate. A distance cutoff
defines a geometric contact, not a hydrogen bond, favorable interaction,
binding affinity, binding free energy, or antibiotic efficacy. Functional-state
labels come from user metadata; they are not inferred from coordinates.

## Install

Use Julia 1.10 or newer:

```bash
git clone https://github.com/eganeganegan/RibosomeAnalysis.jl.git
cd RibosomeAnalysis.jl
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. -e 'using Pkg; Pkg.test()'
```

Load the API:

```julia
using RibosomeAnalysis
structure = load_structure("ribosome.cif")
ligand = select_ligand(structure; residue_name="ACTUAL_CCD_CODE",
                       chain_id="actual_chain", residue_number=123)
site = identify_binding_site(structure, ligand; cutoff=5.0)
contacts = atomic_contacts(structure, ligand; cutoff=4.0)
```

Use mmCIF for ribosomes when possible. Both author and label chain/residue
identifiers are retained. Explicit ligand identifiers should be taken from the
input file; no residue code is assumed to be universal.

## Command-line workflows

```bash
julia --project=. scripts/analyze_structure.jl \
  --structure structure.cif --output results/structure-run

julia --project=. scripts/analyze_binding_site.jl \
  --structure structure.cif --ligand ACTUAL_CCD_CODE \
  --ligand-chain L --ligand-residue 123 --cutoff 5.0 \
  --output results/site-run

julia --project=. scripts/compare_structures.jl \
  --reference state1.cif --target state2.cif \
  --output results/comparison-run

julia --project=. scripts/analyze_trajectory.jl \
  --topology system.pdb --trajectory simulation.xtc \
  --selection "all" --output results/trajectory-run
```

Every workflow creates `metadata.json`, `summary.json`, and `tables/`,
`figures/`, `structures/`, and `logs/` directories. Existing output directories
are not silently overwritten.

## Implemented scope

- PDB and mmCIF parsing through BioStructures.jl, including multiple models,
  hetero atoms, alternate locations, occupancies, B factors, and mmCIF
  author/label identifiers.
- RNA, protein, chain, residue, ligand, and coordinate selections.
- Explicit and metadata-assisted antibiotic identification with user CSV maps.
- KD-tree binding sites and atom contacts at arbitrary Angstrom cutoffs.
- Explicit-topology hydrogen-bond geometry, explicit-group ionic/aromatic
  geometry, and candidate metal-mediated geometry.
- Sequence-aware residue matching, explicit chain maps, Kabsch fitting, RMSD,
  local geometry descriptors, and state-comparison tables.
- In-memory trajectory RMSD/RMSF, contact persistence, distances, PCA,
  deterministic clustering, and representative-frame selection.
- Optional XTC/TRR/DCD/PDB-trajectory streaming through PythonCall and
  MDAnalysis.
- Autocorrelation, effective sample size, deterministic bootstrap and moving
  block-bootstrap intervals, and replica-level summaries.
- CairoMakie PNG/SVG/PDF plots and PyMOL/ChimeraX selection exports.

See the [hosted documentation](https://eganeganegan.github.io/RibosomeAnalysis.jl/stable/),
[methodology](docs/src/methodology.md), [validation](docs/src/validation.md), and
[known limitations](docs/src/validation.md#functionality-requiring-further-validation)
before interpreting results.

## Validation

The validation suite includes independently generated MDAnalysis reference
values and optional integration checks against public RCSB PDB entries 1J7T
and 4V7U. Public structures are downloaded separately and verified by SHA-256,
so core tests remain offline and lightweight:

```bash
julia --project=. validation/download_public_structures.jl validation/data
julia --project=. validation/validate_public_structures.jl validation/data
```

The full 4V7U workflow exercises a 284,525-atom antibiotic-bound bacterial
ribosome. See [validation/README.md](validation/README.md) for provenance,
independent MDAnalysis reproduction, and expected results.

## Repository map

`src/` contains the public API; `scripts/` contains reproducible CLI workflows;
`examples/` contains runnable controlled examples; `test/fixtures/` contains a
synthetic, non-biological structure; `benchmarks/` contains executable
benchmarks; and `docs/` contains scientific and operational documentation.

## Citation and provenance

No biological result or benchmark is bundled as a scientific claim. Analyses
record absolute input paths, SHA-256 checksums, Julia/package versions,
parameters, timestamps, selection definitions, and random seeds. Users remain
responsible for citing primary structures, simulations, BioStructures.jl,
MDAnalysis (if used), and the relevant experimental literature.
