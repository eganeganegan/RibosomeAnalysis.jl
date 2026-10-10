# API reference

The exported API is intentionally small and direct.

## Structures and selections

- `load_structure(path)`: load PDB/mmCIF, including all models and disorder.
- `structure_summary`, `incomplete_information`: composition and data warnings.
- `select_atoms`, `select_chains`, `select_residues`: explicit selection helpers.
- `rna_atoms`, `protein_atoms`, `ligand_atoms`, `coordinates`: common accessors.
- `residue_records`, `chain_records`: stable hierarchy summaries.

## Ligands, sites, and interactions

- `identify_ligands`, `select_ligand`: metadata-assisted or explicit ligands.
- `load_ligand_mapping`: CSV extension point.
- `identify_binding_site`, `binding_site_summary`: KD-tree proximity sites.
- `atomic_contacts`: heavy-atom contact records.
- `explicit_hydrogen_bonds`: topology- and hydrogen-dependent geometry.
- `ionic_proximities`, `aromatic_interactions`: explicit chemical-group geometry.
- `metal_mediated_interactions`: candidate metal bridging geometry.
- `export_selection`: PyMOL or ChimeraX residue selection.

## Alignment and geometry

- `match_residues`: chain-controlled global sequence correspondence.
- `kabsch`, `apply_transform`, `rmsd`: numerical rigid-body operations.
- `align_structures`: correspondence, fit, coverage, and unmatched-residue report.
- `centroid`, `radius_of_gyration`, `spatial_extent`, `pairwise_distances`,
  `local_contact_density`, `geometry_descriptors`: geometric observables.
- `compare_binding_sites`: global/local site RMSD, ligand displacement, contacts.

## Dynamics and statistics

- `TrajectoryData`: unit-normalized `3 × atoms × frames` representation with
  optional `3 × 3 × frames` unit-cell vectors.
- `TrajectoryFrame`, `foreach_trajectory_frame`: one-frame-at-a-time optional
  MDAnalysis processing.
- `load_trajectory`: optional in-memory MDAnalysis reader.
- `trajectory_rmsd`, `trajectory_rmsf`, `contact_persistence`,
  `distance_series`: trajectory observables.
- `minimum_image_displacement`, `periodic_distance`: explicit orthorhombic and
  triclinic periodic geometry.
- `trajectory_pca`, `cluster_frames`, `representative_frames`: exploratory
  ensemble reduction.
- `autocorrelation`, `effective_sample_size`, `bootstrap_ci`,
  `block_bootstrap_ci`, `summarize_samples`, `compare_replicates`: uncertainty.

## Figures and results

- `plot_rmsd`, `plot_rmsf`, `plot_contact_heatmap`, `plot_pca`: CairoMakie
  PNG/SVG/PDF output.
- `prepare_results`, `write_metadata`, `write_summary`, `write_table`:
  non-overwriting reproducible output.

Use Julia's help mode (`?function_name`) for signatures and method docstrings.
