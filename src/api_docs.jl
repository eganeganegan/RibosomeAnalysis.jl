@doc """
    ResidueKey

Stable residue identity containing model, chain, residue number, insertion code,
and residue name. ResidueKey values are sortable and suitable as dictionary
keys.
""" ResidueKey

@doc """
    ResidueRecord

A residue identity, its indices into StructureData.atoms, and its classified
molecule type (:rna, :protein, :ligand, or :other).
""" ResidueRecord

@doc """
    ChainRecord

A model/chain identity and the indices of its residues in residue_records.
""" ChainRecord

@doc """
    StructureMetadata

Source, format, parser, identifier, and free-form annotations recorded while a
structure is loaded.
""" StructureMetadata

@doc """
    StructureData

A collection of AtomRecords together with provenance in StructureMetadata.
Coordinates are stored in Angstrom.
""" StructureData

@doc """
    Ligand

An explicitly identified ligand residue and the atom indices that belong to it.
The identification method and any antibiotic annotation remain attached to the
record.
""" Ligand

@doc """
    BindingSite

Residues and atoms found within a stated cutoff of a Ligand, including the
minimum ligand distance for each residue.
""" BindingSite

@doc """
    ContactRecord

One structure-atom/ligand-atom contact with residue identity and distance in
Angstrom.
""" ContactRecord

@doc """
    InteractionRecord

A geometry-based interaction candidate with its participants, measured distance
or angle, confidence label, and interpretive limitation.
""" InteractionRecord

@doc """
    AlignmentResult

Rigid-body transform and correspondence diagnostics returned by
align_structures, including RMSD, coverage, matches, and unmatched residues.
""" AlignmentResult

@doc """
    StateMetadata

User-supplied biological and experimental annotations for a structure state.
These labels are never inferred from coordinates.
""" StateMetadata

@doc """
    ComparisonResult

A labeled pairwise state comparison containing the alignment, scalar metrics,
changed-contact table, and explicit limitations.
""" ComparisonResult

@doc """
    structure_summary(structure)

Return a dictionary summarizing atom, residue, chain, model, and molecule-type
counts for structure.
""" structure_summary

@doc """
    incomplete_information(structure)

Return warnings describing missing or ambiguous structural metadata that may
limit downstream interpretation.
""" incomplete_information

@doc """
    select_atoms(structure, predicate)

Return atom indices for which predicate(atom) is true.
""" select_atoms

@doc """
    select_chains(structure, chain_ids; model_id=nothing)

Return atom indices belonging to the requested chain identifiers, optionally
restricted to one model.
""" select_chains

@doc """
    select_residues(structure; kwargs...)

Return atom indices matching explicit residue fields such as chain, residue
number, insertion code, residue name, or model.
""" select_residues

@doc """
    rna_atoms(structure)

Return the indices of atoms classified as RNA.
""" rna_atoms

@doc """
    protein_atoms(structure)

Return the indices of atoms classified as protein.
""" protein_atoms

@doc """
    ligand_atoms(structure)

Return the indices of atoms classified as non-polymer ligands.
""" ligand_atoms

@doc """
    coordinates(atoms)
    coordinates(structure, indices)

Return coordinates as a contiguous 3 × N matrix in Angstrom.
""" coordinates

@doc """
    residue_records(structure)

Return stable residue records assembled from the atoms in structure.
""" residue_records

@doc """
    chain_records(structure)

Return model-aware chain records referencing the structure's residue records.
""" chain_records

@doc """
    antibiotic_metadata()

Return a copy of the package's conservative residue-name metadata. The mapping
is an annotation aid, not universal ligand identification.
""" antibiotic_metadata

@doc """
    select_ligand(structure; kwargs...)

Select exactly one ligand using explicit identifying fields. An error is raised
when the criteria match zero or multiple residues.
""" select_ligand

@doc """
    identify_binding_site(structure, ligand; cutoff=5.0,
                          include_hydrogen=false, include_ligand=false)

Find residues with atoms within cutoff Angstrom of ligand using a KD tree.
By default hydrogens and the ligand residue itself are excluded.
""" identify_binding_site

@doc """
    binding_site_summary(structure, site)

Return counts and minimum-distance diagnostics for a BindingSite.
""" binding_site_summary

@doc """
    atomic_contacts(structure, ligand; cutoff=4.0, heavy_atoms_only=true)

Return atom-pair contacts within cutoff Angstrom, sorted by distance. A contact
is a geometric proximity and does not imply binding energy.
""" atomic_contacts

@doc """
    metal_mediated_interactions(structure, ligand; metal_cutoff=3.0,
                                ribosome_cutoff=3.0)

Report distance-based candidate bridges between recognized metal elements, the
ligand, and ribosomal N/O/S atoms. Coordination chemistry is not inferred.
""" metal_mediated_interactions

@doc """
    centroid(xyz)

Return the centroid of a 3 × N coordinate matrix.
""" centroid

@doc """
    radius_of_gyration(xyz; masses=nothing)

Return the unweighted or mass-weighted radius of gyration in the coordinate
matrix's length unit.
""" radius_of_gyration

@doc """
    spatial_extent(xyz)

Return the axis-aligned coordinate extent along each of the three dimensions.
""" spatial_extent

@doc """
    pairwise_distances(xyz)

Return the symmetric Euclidean distance matrix for columns of xyz.
""" pairwise_distances

@doc """
    local_contact_density(xyz; cutoff=5.0)

Return the number of neighboring points within cutoff of each coordinate,
excluding self contacts.
""" local_contact_density

@doc """
    geometry_descriptors(structure, site)

Return transparent geometric descriptors for a binding site, including size,
centroid, extent, radius of gyration, and local contact density.
""" geometry_descriptors

@doc """
    kabsch(reference, mobile)

Compute the proper rotation and translation that least-squares superimpose the
3 × N mobile coordinates on the reference coordinates.
""" kabsch

@doc """
    apply_transform(xyz, rotation, translation)

Apply a rigid-body rotation and translation to a 3 × N coordinate matrix.
""" apply_transform

@doc """
    rmsd(reference, target; superimpose=false)

Compute coordinate RMSD, optionally performing a Kabsch superposition first.
""" rmsd

@doc """
    align_structures(reference, target; kwargs...)

Match polymer residues, fit corresponding atoms, and return an AlignmentResult
with coverage and unmatched-residue diagnostics.
""" align_structures

@doc """
    trajectory_rmsd(trajectory, reference; kwargs...)

Compute RMSD for each trajectory frame, with optional alignment on a specified
atom subset.
""" trajectory_rmsd

@doc """
    trajectory_rmsf(trajectory; kwargs...)

Compute per-atom root-mean-square fluctuations after optional frame alignment.
""" trajectory_rmsf

@doc """
    residue_rmsf(structure, trajectory; kwargs...)

Aggregate aligned per-atom RMSF by residue using the root mean square of the
atom RMSF values. Trajectory atom identifiers must map to atom indices in the
supplied structure.
""" residue_rmsf

@doc """
    contact_persistence(trajectory, group_a, group_b; kwargs...)

Return the fraction of frames in which each selected atom pair lies within the
specified cutoff, optionally using periodic minimum-image distances.
""" contact_persistence

@doc """
    residue_contact_persistence(structure, trajectory, ligand_indices; kwargs...)

Return a residue-level geometric contact fingerprint over trajectory frames.
Indices are trajectory-local and residue identities are recovered through the
trajectory-to-structure atom mapping. Contact fractions do not measure binding
affinity or interaction energy.
""" residue_contact_persistence

@doc """
    distance_series(trajectory, atom_i, atom_j; periodic=false)

Return the interatomic distance for every frame, optionally using each frame's
periodic unit cell.
""" distance_series

@doc """
    minimum_distance_series(trajectory, group_a, group_b; periodic=false)

Return the minimum inter-group atom distance for each trajectory frame,
optionally using periodic minimum-image distances.
""" minimum_distance_series

@doc """
    group_centroid_displacement(trajectory, group_indices; kwargs...)

Align each frame on the requested atoms and return displacement of the selected
group centroid from its reference-frame position. This geometric displacement
does not establish ligand unbinding or binding affinity.
""" group_centroid_displacement

@doc """
    minimum_image_displacement(displacement, cell)

Map a Cartesian displacement to the nearest periodic image for a validated
three-dimensional cell matrix.
""" minimum_image_displacement

@doc """
    periodic_distance(a, b, cell)

Return the minimum-image distance between two Cartesian positions.
""" periodic_distance

@doc """
    trajectory_pca(trajectory; kwargs...)

Perform PCA on flattened trajectory coordinates and return means, components,
scores, variances, and explained-variance ratios.
""" trajectory_pca

@doc """
    representative_frames(features, clustering)

Return the frame nearest each cluster center in feature space.
""" representative_frames

@doc """
    autocorrelation(values; max_lag=nothing)

Estimate the normalized autocorrelation sequence through the requested lag.
""" autocorrelation

@doc """
    effective_sample_size(values; max_lag=nothing)

Estimate effective sample size using the initial positive sequence of
autocorrelations.
""" effective_sample_size

@doc """
    bootstrap_ci(values; statistic=mean, confidence=0.95, n_resamples=2000,
                 rng=MersenneTwister(1))

Compute a deterministic-by-default percentile bootstrap confidence interval for
independent observations.
""" bootstrap_ci

@doc """
    block_bootstrap_ci(values; block_length, statistic=mean, confidence=0.95,
                       n_resamples=2000, rng=MersenneTwister(1))

Compute a moving-block bootstrap interval for correlated observations.
""" block_bootstrap_ci

@doc """
    summarize_samples(values; correlated=false, kwargs...)

Return descriptive statistics, effective sample size, and an appropriate
bootstrap interval for a sample series.
""" summarize_samples

@doc """
    compare_replicates(replicates; statistic=mean)

Summarize each independent replicate and the between-replicate distribution of
the chosen statistic.
""" compare_replicates

@doc """
    summarize_state_replicas(observations; state=:state, replica=:replica,
                             value=:value, time=nothing, kwargs...)

Reduce each correlated trajectory replica to one estimate, calculate
within-replica moving-block intervals, and summarize independent replica
estimates with equal replica weighting.
When `time` names a column, observations are sorted chronologically within each
replica before block resampling; otherwise input row order is preserved.
""" summarize_state_replicas

@doc """
    pairwise_state_differences(replica_summary; kwargs...)

Calculate target-minus-reference differences between state-level replica
estimates. Intervals are omitted unless both states have at least two
independent replicas.
""" pairwise_state_differences

@doc """
    compare_state_ensembles(observations; kwargs...)

Run replica-aware state summarization and pairwise comparison for a per-frame
metric table. The result includes per-replica, per-state, and difference tables
plus explicit statistical assumptions.
""" compare_state_ensembles

@doc """
    prepare_results(base; overwrite=false)

Create a non-overwriting result directory with tables, figures, structures, and
logs subdirectories.
""" prepare_results

@doc """
    write_metadata(output_dir; kwargs...)

Write metadata.json containing inputs, checksums, parameters, environment
versions, timestamp, and optional annotations.
""" write_metadata

@doc """
    write_summary(output_dir, summary)

Write summary.json without replacing an existing file.
""" write_summary

@doc """
    write_table(output_dir, filename, table)

Write a CSV table beneath the result directory's tables folder without
replacing an existing file.
""" write_table

@doc """
    compare_binding_sites(reference, target, reference_ligand, target_ligand; kwargs...)

Align two structures and compare their binding-site geometry, ligand placement,
and atom-contact sets.
""" compare_binding_sites

@doc """
    compare_states(items; cutoff=5.0)

Summarize binding-site geometry and supplied metadata for each
`(structure, ligand, state_metadata)` tuple in a DataFrame. State labels and
biological annotations are supplied by the caller, not inferred.
""" compare_states
