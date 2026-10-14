using RibosomeAnalysis

path = joinpath(@__DIR__, "..", "test", "fixtures", "synthetic.pdb")
reference = load_structure(path; identifier = "reference")
target = load_structure(path; identifier = "target")
reference_ligand = select_ligand(reference; residue_name = "PAR")
target_ligand = select_ligand(target; residue_name = "PAR")

reference_state = StateMetadata(
    "reference state",
    "controlled fixture",
    "synthetic assembly",
    "synthetic coordinates",
    missing,
    "fixture ligand",
    "user annotation A",
    Dict("provenance" => "test/fixtures/synthetic.pdb"),
)
target_state = StateMetadata(
    "target state",
    "controlled fixture",
    "synthetic assembly",
    "synthetic coordinates",
    missing,
    "fixture ligand",
    "user annotation B",
    Dict("provenance" => "test/fixtures/synthetic.pdb"),
)

comparison = compare_binding_sites(reference, target, reference_ligand, target_ligand)
states = compare_states([
    (reference, reference_ligand, reference_state),
    (target, target_ligand, target_state),
])

display(comparison)
display(states)
