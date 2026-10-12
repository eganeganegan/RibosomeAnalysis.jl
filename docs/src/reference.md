# API reference

The exported API is grouped below. Each entry is sourced from the same
docstring shown by Julia's help mode.

## Structures and selections

```@docs
RibosomeAnalysis
AtomRecord
ResidueKey
ResidueRecord
ChainRecord
StructureMetadata
StructureData
load_structure
structure_summary
incomplete_information
select_atoms
select_chains
select_residues
rna_atoms
protein_atoms
ligand_atoms
coordinates
residue_records
chain_records
```

## Ligands, sites, and interactions

```@docs
Ligand
BindingSite
ContactRecord
InteractionRecord
antibiotic_metadata
load_ligand_mapping
identify_ligands
select_ligand
identify_binding_site
binding_site_summary
export_selection
atomic_contacts
explicit_hydrogen_bonds
ionic_proximities
aromatic_interactions
metal_mediated_interactions
```

## Alignment and geometry

```@docs
AlignmentResult
match_residues
kabsch
apply_transform
rmsd
align_structures
centroid
radius_of_gyration
spatial_extent
pairwise_distances
local_contact_density
geometry_descriptors
compare_binding_sites
```

## Dynamics and statistics

```@docs
TrajectoryData
TrajectoryFrame
load_trajectory
foreach_trajectory_frame
trajectory_rmsd
trajectory_rmsf
contact_persistence
distance_series
minimum_image_displacement
periodic_distance
trajectory_pca
cluster_frames
representative_frames
autocorrelation
effective_sample_size
bootstrap_ci
block_bootstrap_ci
summarize_samples
compare_replicates
```

## Figures and results

The plotting functions require the optional CairoMakie extension.

```@docs
plot_rmsd
plot_rmsf
plot_contact_heatmap
plot_pca
prepare_results
write_metadata
write_summary
write_table
```

## State comparisons

```@docs
StateMetadata
ComparisonResult
compare_states
```
