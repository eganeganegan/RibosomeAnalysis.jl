# Usage

## Structure selection

Selections use explicit Julia arguments rather than an undocumented expression
language:

```julia
select_chains(structure, ["A", "B"]; model_id=1)
select_residues(structure; chain_id="A", numbers=1490:1500, names=["A", "G"])
rna_atoms(structure)
protein_atoms(structure)
ligand_atoms(structure)
select_atoms(structure, atom -> atom.element == "MG")
```

Residue identity includes model, chain, residue number, insertion code, and
residue name. Alternate locations remain separate atoms. This prevents silent
chain merging or renumbering.

## Ligands and sites

Use explicit selection when possible:

```julia
ligand = select_ligand(structure; residue_name="PAR", chain_id="L",
                       residue_number=100, model_id=1)
sites = Dict(r => identify_binding_site(structure, ligand; cutoff=r)
             for r in (3.0, 4.0, 5.0, 8.0))
```

`PAR` is only an example from the controlled fixture. Read the actual component
identifier from the structure. Metadata-assisted identification returns only
hetero residues in its mapping and never classifies every nonstandard residue
as an antibiotic. Supply mappings as CSV with `residue_name,name,class`.

## Alignment

```julia
fit = align_structures(reference, target;
                       chain_map=Dict("author_chain_A" => "author_chain_B"))
```

Residues are globally sequence-aligned inside paired chains; only chemically
compatible matches are retained. Explicit chain maps must be one-to-one and
cannot pair RNA with protein. Without a map, identical chain IDs are paired
only when their polymer types agree. Unrecognized or modified residues match
only when their residue names agree. Protein N/CA/C and RNA P/C4'/C3'/O3' atoms
are fitted by default. Use `atom_names` or explicit residue pairs for a local
binding-site fit. Review unmatched residues and coverage before using the RMSD.

The structure-comparison CLI accepts explicit chain pairs as a comma-separated
`reference:target` map:

```bash
julia --project=. scripts/compare_structures.jl \
  --reference state1.cif --target state2.cif \
  --chain-map A:B,C:D --output results/comparison-run
```

`compare_binding_sites` uses the same sequence-derived polymer correspondence
for its contact-change table. Non-polymers are paired only by complete matching
residue identifiers and are labeled `exact_identifier` rather than `sequence`.

## Trajectory arrays

Coordinates have shape `3 × atoms × frames`. Internal length units are Angstrom
and times are picoseconds:

```julia
trajectory = TrajectoryData(array; times_ps=times, input_unit=:nanometer)
rmsd_values = trajectory_rmsd(trajectory; align_indices=backbone,
                              measure_indices=binding_site)
rmsf_values = trajectory_rmsf(trajectory; align_indices=backbone)
persistence = contact_persistence(trajectory, ligand_indices, site_indices;
                                  cutoff=4.0)
minimum_distances = minimum_distance_series(trajectory, ligand_indices,
                                            site_indices)
ligand_displacement = group_centroid_displacement(
    trajectory, ligand_indices; align_indices=backbone)
per_residue = residue_rmsf(structure, trajectory; align_indices=backbone)
```

Nanometers are multiplied by exactly 10. Align before interpreting RMSF unless
whole-body motion is intentionally part of the observable.

`group_centroid_displacement` is measured after fitting the selected alignment
atoms. It is a coordinate observable, not evidence of ligand dissociation or a
binding free-energy estimate. `residue_rmsf` groups atom RMSF values by the
trajectory-to-structure atom-index mapping and reports their root mean square.

When a loaded trajectory supplies unit-cell vectors, periodic observables are
explicitly enabled rather than silently assumed:

```julia
distances = distance_series(trajectory, ligand_atom, site_atom; periodic=true)
persistence = contact_persistence(trajectory, ligand_indices, site_indices;
                                  cutoff=4.0, periodic=true)
```

Both calls fail if `periodic=true` and unit-cell data are unavailable. Cell
columns are lattice vectors in Angstrom. The minimum-image implementation
checks neighboring lattice translations and supports orthorhombic and reduced
triclinic MD cells. It does not reconstruct molecules split across boundaries;
perform topology-aware unwrapping before RMSD/RMSF when required.

For trajectories too large to hold in memory, process independent frame copies
with a callback:

```julia
foreach_trajectory_frame("system.tpr", "trajectory.xtc";
                         selection="nucleic or resname ERY") do frame
    # frame.coordinates is 3 × selected_atoms in Angstrom
    # frame.unit_cell is either a 3 × 3 matrix or nothing
end
```

## Command-line validation

The command-line workflows reject unknown and duplicate options before an
analysis starts. Required options, numeric ligand residue identifiers, and
distance cutoffs are validated with concise error messages. Invalid user input
exits with status code 2 and does not emit an internal Julia stack trace.

## Outputs

`prepare_results` refuses an existing directory unless `overwrite=true` is
explicit. CLI workflows do not expose overwrite, encouraging immutable runs.
Tables are CSV; metadata and summaries are JSON; figures have physical units.

## State comparisons and figures

State labels, organism, subunit, experimental method, resolution, antibiotic,
and translation state are user annotations. They are never inferred from
coordinates:

```julia
states = compare_states([(structure, ligand, state_metadata)]; cutoff=5.0)
plot_state_metric(states, :radius_of_gyration_angstrom, "state-rg.svg")
plot_binding_site_comparison(comparison, "site-comparison.pdf")
plot_distance_distribution(minimum_distances, "distances.png")
plot_cluster_populations(clustering.labels, "clusters.svg")
```

The plotting calls require CairoMakie. Missing state metrics are omitted rather
than imputed.
