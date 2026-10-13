@testset "Interaction geometry" begin
    structure = load_structure(joinpath(@__DIR__, "..", "fixtures", "synthetic.pdb"))
    ligand = select_ligand(structure; residue_name = "PAR")
    contacts = atomic_contacts(structure, ligand; cutoff = 4.0)
    @test !isempty(contacts)
    @test all(c.distance_angstrom <= 4.0 for c in contacts)
    ionic = ionic_proximities(structure, [12], [5]; cutoff = 5.0)
    @test length(ionic) == 1
    @test ionic[1].confidence == :candidate
    metals = metal_mediated_interactions(
        structure,
        ligand;
        metal_cutoff = 3.0,
        ribosome_cutoff = 6.0,
    )
    @test !isempty(metals)
    @test all(i.interaction_type == :metal_mediated_candidate for i in metals)
    @test all(i.confidence == :candidate for i in metals)
    @test isempty(explicit_hydrogen_bonds(structure, ligand, Tuple{Int,Int}[]))
end

@testset "Controlled hydrogen-bond angle" begin
    atoms = [
        AtomRecord(
            1,
            "N",
            "LIG",
            "L",
            "L",
            "L",
            1,
            1,
            1,
            ' ',
            ' ',
            1,
            1,
            0,
            "N",
            true,
            [0.0, 0, 0],
        ),
        AtomRecord(
            2,
            "H",
            "LIG",
            "L",
            "L",
            "L",
            1,
            1,
            1,
            ' ',
            ' ',
            1,
            1,
            0,
            "H",
            true,
            [1.0, 0, 0],
        ),
        AtomRecord(
            3,
            "O",
            "A",
            "A",
            "A",
            "A",
            2,
            2,
            2,
            ' ',
            ' ',
            1,
            1,
            0,
            "O",
            false,
            [2.8, 0, 0],
        ),
    ]
    structure =
        StructureData(atoms, StructureMetadata("hbond", "", :synthetic, "test", Dict()))
    ligand = Ligand(
        "LIG:L:1:model1",
        "LIG",
        "L",
        1,
        ' ',
        1,
        [1, 2],
        missing,
        :explicit,
        String[],
    )
    bonds = [(1, 2)]
    hbonds = explicit_hydrogen_bonds(structure, ligand, bonds)
    @test length(hbonds) == 1
    @test hbonds[1].angle_degrees ≈ 180.0
end
