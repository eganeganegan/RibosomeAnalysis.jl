using RibosomeAnalysis

path = joinpath(@__DIR__, "..", "test", "fixtures", "synthetic.pdb")
reference = load_structure(path; identifier = "reference")
target = load_structure(path; identifier = "target")
reference_ligand = select_ligand(reference; residue_name = "PAR")
target_ligand = select_ligand(target; residue_name = "PAR")

display(align_structures(reference, target))
display(compare_binding_sites(reference, target, reference_ligand, target_ligand))
