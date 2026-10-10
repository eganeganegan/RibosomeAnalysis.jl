# RibosomeAnalysis.jl

RibosomeAnalysis.jl is a Julia package and command-line toolkit for
reproducible analysis of ribosomal structures, antibiotic binding sites,
structural correspondence, geometric interactions, and molecular-dynamics
trajectories.

The software distinguishes geometric observations from chemical or biological
interpretation. A distance cutoff identifies a geometric contact; it does not
establish a hydrogen bond, favorable interaction, binding affinity, binding
free energy, or antibiotic efficacy. Functional-state labels are supplied by
the user and are not inferred from coordinates.

## Where to begin

- Follow [Installation](installation.md) to set up the package.
- See [Usage](usage.md) for API and command-line examples.
- Read [Methodology](methodology.md) for definitions and assumptions.
- Review [Validation](validation.md) before interpreting results.
- Consult the [API reference](reference.md) for the exported interface.

## Validated scope

The offline test suite covers numerical kernels, structure parsing, selections,
binding sites, contacts, comparisons, trajectory analysis, result provenance,
and error handling. Separately runnable public-data checks use checksum-pinned
RCSB entries 1J7T and 4V7U. Independent MDAnalysis 2.10.0 reference values are
checked for RMSD, RMSF, distances, and the 1J7T contact set.

See [Validation](validation.md#functionality-requiring-further-validation) for
functionality that still requires production-scale cross-validation.
