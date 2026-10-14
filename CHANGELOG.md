# Changelog

## 0.2.0 - 2026-10-13

- Reject unknown, duplicated, missing, and malformed command-line options with
  concise user-facing errors and exit status 2.
- Exercise complete structure, binding-site, and comparison CLI workflows in
  the offline test suite.
- Expand regression coverage for optional extensions, unit-cell input forms,
  modified RNA classification, metal-mediated candidates, and reproducible
  output metadata.
- Make residue correspondence chemically conservative by validating explicit
  chain maps, separating RNA from protein identities, and requiring exact names
  for unrecognized residues.
- Compare binding-site contacts through sequence-derived polymer
  correspondence, with separately labeled exact-identifier matching for ions,
  solvent, and other non-polymers.
- Add per-residue RMSF, inter-group minimum-distance series, aligned group
  centroid displacement, stricter trajectory validation, and zero-variance PCA
  handling.
- Add distance-distribution, cluster-population, binding-site-comparison, and
  state-metric figures through the optional CairoMakie extension.
- Expand state-comparison tables, add explicit CLI chain maps, and include a
  controlled functional-state comparison example.

## 0.1.0 - 2026-10-09

- Initial Julia package and command-line workflows for ribosome structure,
  binding-site, interaction, structural-comparison, and trajectory analysis.
- Offline unit, integration, regression, precision, and error-handling tests.
- Independent MDAnalysis RMSD, RMSF, distance, and real-structure contact
  references.
- Per-frame orthorhombic and triclinic unit-cell preservation, minimum-image
  distances, and periodic contact persistence.
- Callback-based one-frame-at-a-time XTC/TRR/DCD processing, cross-validated
  against MDAnalysis 2.10.0 and GROMACS.
- Checksum-pinned public-data validation using RCSB PDB entries 1J7T and 4V7U.
- Reproducible JSON/CSV outputs, scientific plots, and PyMOL/ChimeraX exports.
