using RibosomeAnalysis

structure = load_structure(joinpath(@__DIR__, "..", "test", "fixtures", "synthetic.pdb"))
display(structure_summary(structure))
display(identify_ligands(structure))
