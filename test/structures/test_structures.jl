@testset "Structure loading and selection" begin
    fixture = joinpath(@__DIR__, "..", "fixtures", "synthetic.pdb")
    structure = load_structure(fixture)
    @test length(structure.atoms) == 14
    @test structure.metadata.format == :pdb
    @test length(rna_atoms(structure)) == 5
    @test length(protein_atoms(structure)) == 4
    @test length(ligand_atoms(structure)) == 3
    @test length(select_chains(structure, "A")) == 5
    @test length(select_residues(structure; chain_id = "L", numbers = [100])) == 3
    summary = structure_summary(structure)
    @test summary["models"] == [1]
    @test summary["residues"] == 5
    @test_throws ArgumentError load_structure("missing.pdb")
    mktempdir() do directory
        unknown = joinpath(directory, "structure.xyz")
        write(unknown, "not a supported structure\n")
        @test_throws ArgumentError load_structure(unknown)

        compressed_pdb = joinpath(directory, "structure.pdb.gz")
        write(compressed_pdb, "not actually compressed\n")
        @test_throws ArgumentError load_structure(compressed_pdb)
    end
end

@testset "mmCIF and compressed mmCIF" begin
    fixture = joinpath(@__DIR__, "..", "fixtures", "synthetic.pdb")
    raw = BioStructures.read(fixture, BioStructures.PDBFormat)
    mktempdir() do directory
        cif_path = joinpath(directory, "synthetic.cif")
        compressed_path = joinpath(directory, "synthetic.cif.gz")
        BioStructures.writemmcif(cif_path, raw)
        BioStructures.writemmcif(compressed_path, raw; gzip = true)
        cif = load_structure(cif_path)
        compressed = load_structure(compressed_path)
        @test cif.metadata.format == :mmcif
        @test length(cif.atoms) == 14
        @test cif.atoms[1].author_chain_id == "A"
        @test length(compressed.atoms) == 14
        @test compressed.metadata.annotations["compressed"]
    end
end

@testset "Modified RNA classification" begin
    atoms = [
        AtomRecord(
            serial,
            atom_name,
            "PSU",
            "R",
            "R",
            "R",
            1,
            1,
            1,
            ' ',
            ' ',
            1,
            1.0,
            0.0,
            element,
            false,
            [Float64(serial), 0.0, 0.0],
        ) for
        (serial, atom_name, element) in ((1, "C1'", "C"), (2, "P", "P"), (3, "O2", "O"))
    ]
    structure = StructureData(
        atoms,
        StructureMetadata("modified-rna", "", :synthetic, "test", Dict()),
    )
    @test only(residue_records(structure)).molecule_type == :rna
end
