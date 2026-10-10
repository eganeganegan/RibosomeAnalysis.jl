# Changelog

## Unreleased

- Preserve per-frame orthorhombic and triclinic unit-cell vectors when reading
  trajectories through MDAnalysis.
- Add minimum-image distances and periodic contact persistence with explicit
  opt-in behavior and missing-cell validation.
- Add callback-based one-frame-at-a-time XTC/TRR/DCD processing.
- Cross-validate deterministic XTC, TRR, and DCD fixtures against MDAnalysis
  2.10.0 and periodic distances against GROMACS.

## 0.1.0 - 2026-10-09

- Initial Julia package and command-line workflows for ribosome structure,
  binding-site, interaction, structural-comparison, and trajectory analysis.
- Offline unit, integration, regression, precision, and error-handling tests.
- Independent MDAnalysis RMSD, RMSF, distance, and real-structure contact
  references.
- Checksum-pinned public-data validation using RCSB PDB entries 1J7T and 4V7U.
- Reproducible JSON/CSV outputs, scientific plots, and PyMOL/ChimeraX exports.
