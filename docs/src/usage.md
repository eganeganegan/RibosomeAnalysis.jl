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

Residues are globally sequence-aligned inside explicitly paired chains; only
matching residue symbols are retained. Protein N/CA/C and RNA P/C4'/C3'/O3'
atoms are fitted by default. Use `atom_names` or explicit residue pairs for a
local binding-site fit. Review unmatched residues and coverage before using the
RMSD.

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
```

Nanometers are multiplied by exactly 10. Align before interpreting RMSF unless
whole-body motion is intentionally part of the observable.

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

## Outputs

`prepare_results` refuses an existing directory unless `overwrite=true` is
explicit. CLI workflows do not expose overwrite, encouraging immutable runs.
Tables are CSV; metadata and summaries are JSON; figures have physical units.
