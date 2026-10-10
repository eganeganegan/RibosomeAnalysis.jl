using RibosomeAnalysis

structure = load_structure(joinpath(@__DIR__, "..", "test", "fixtures", "synthetic.pdb"))
ligand =
    select_ligand(structure; residue_name = "PAR", chain_id = "L", residue_number = 100)
site = identify_binding_site(structure, ligand; cutoff = 5.0)
display(binding_site_summary(structure, site))
display(atomic_contacts(structure, ligand; cutoff = 4.0))
