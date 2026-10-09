@testset "Ligands and binding sites" begin
    structure = load_structure(joinpath(@__DIR__, "..", "fixtures", "synthetic.pdb"))
    candidates = identify_ligands(structure)
    @test length(candidates) == 1
    @test candidates[1].residue_name == "PAR"
    ligand = select_ligand(structure; residue_name="PAR", chain_id="L", residue_number=100)
    @test ligand.identification == :explicit
    site3 = identify_binding_site(structure, ligand; cutoff=3.0)
    site5 = identify_binding_site(structure, ligand; cutoff=5.0)
    @test length(site3.residue_keys) <= length(site5.residue_keys)
    @test any(k -> k.chain_id == "A", site5.residue_keys)
    @test !any(k -> k.chain_id == "L", site5.residue_keys)
    @test binding_site_summary(structure, site5)["rna_residues"] == 1
    @test_throws ArgumentError identify_binding_site(structure, ligand; cutoff=-1)
    mktempdir() do directory
        mapping_path = joinpath(directory, "ligands.csv")
        write(mapping_path, "residue_name,name,class\nPAR,paromomycin,aminoglycoside\n")
        @test load_ligand_mapping(mapping_path)["PAR"]["class"] == "aminoglycoside"
    end
end
