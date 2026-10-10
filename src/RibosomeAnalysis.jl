module RibosomeAnalysis

using LinearAlgebra
using Statistics
using SparseArrays
using Random
using SHA
using Dates
using Printf
using StaticArrays
using NearestNeighbors
using StatsBase
using DataFrames
using CSV
using JSON3
using BioStructures

include("structures/types.jl")
include("io/structures.jl")
include("structures/selections.jl")
include("binding/ligands.jl")
include("binding/sites.jl")
include("interactions/interactions.jl")
include("geometry/descriptors.jl")
include("alignment/alignment.jl")
include("dynamics/trajectories.jl")
include("statistics/ensemble.jl")
include("visualization/plots.jl")
include("io/results.jl")
include("comparison/comparison.jl")

export AtomRecord,
    ResidueKey,
    ResidueRecord,
    ChainRecord,
    StructureMetadata,
    StructureData,
    Ligand,
    BindingSite,
    ContactRecord,
    InteractionRecord,
    AlignmentResult,
    TrajectoryData,
    StateMetadata,
    ComparisonResult
export load_structure,
    structure_summary,
    incomplete_information,
    select_atoms,
    select_chains,
    select_residues,
    rna_atoms,
    protein_atoms,
    ligand_atoms,
    coordinates,
    residue_records,
    chain_records
export antibiotic_metadata, load_ligand_mapping, identify_ligands, select_ligand
export identify_binding_site, binding_site_summary, export_selection
export atomic_contacts,
    explicit_hydrogen_bonds,
    ionic_proximities,
    aromatic_interactions,
    metal_mediated_interactions
export centroid,
    radius_of_gyration,
    spatial_extent,
    pairwise_distances,
    local_contact_density,
    geometry_descriptors
export kabsch, apply_transform, rmsd, match_residues, align_structures
export trajectory_rmsd,
    trajectory_rmsf,
    contact_persistence,
    distance_series,
    trajectory_pca,
    cluster_frames,
    representative_frames,
    load_trajectory
export autocorrelation,
    effective_sample_size,
    bootstrap_ci,
    block_bootstrap_ci,
    summarize_samples,
    compare_replicates
export plot_rmsd, plot_rmsf, plot_contact_heatmap, plot_pca
export prepare_results, write_metadata, write_table, write_summary
export compare_binding_sites, compare_states

end
