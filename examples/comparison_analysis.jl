using RibosomeAnalysis

path = joinpath(@__DIR__, "..", "test", "fixtures", "synthetic.pdb")
reference = load_structure(path; identifier="reference")
target = load_structure(path; identifier="target")
display(align_structures(reference, target))
