@testset "Kabsch validation" begin
    reference = [0.0 1 0 0; 0 0 1 0; 0 0 0 1]
    rotation_true = [0.0 -1 0; 1 0 0; 0 0 1]
    translation_true = [4.0, -2.0, 3.0]
    mobile = rotation_true * reference .+ translation_true
    rotation, translation, value = kabsch(reference, mobile)
    @test value < 1e-12
    @test apply_transform(mobile, rotation, translation) ≈ reference atol = 1e-12
    @test rmsd(reference, mobile; superimpose = true) < 1e-12
end

@testset "Sequence alignment gaps" begin
    @test RibosomeAnalysis._needleman_wunsch(['A', 'C'], ['A']) == [(1, 1)]
    @test RibosomeAnalysis._needleman_wunsch(['A'], ['A', 'C']) == [(1, 1)]
end

function synthetic_polymer(chain_id, residue_name, atom_names)
    atoms = [
        AtomRecord(
            serial,
            atom_name,
            residue_name,
            chain_id,
            chain_id,
            chain_id,
            1,
            1,
            1,
            ' ',
            ' ',
            1,
            1.0,
            0.0,
            first(atom_name),
            false,
            [Float64(serial), 0.0, 0.0],
        ) for (serial, atom_name) in enumerate(atom_names)
    ]
    StructureData(
        atoms,
        StructureMetadata("synthetic-polymer", "", :synthetic, "test", Dict()),
    )
end

@testset "Chemically conservative residue matching" begin
    path = joinpath(@__DIR__, "..", "fixtures", "synthetic.pdb")
    structure = load_structure(path)

    @test_throws ArgumentError match_residues(
        structure,
        structure;
        chain_map = Dict("missing" => "A"),
    )
    @test_throws ArgumentError match_residues(
        structure,
        structure;
        chain_map = Dict("A" => "missing"),
    )
    @test_throws ArgumentError match_residues(
        structure,
        structure;
        chain_map = Dict("A" => "B"),
    )

    pseudouridine = synthetic_polymer("R", "PSU", ["C1'", "P", "O2"])
    second_reference = synthetic_polymer("Q", "PSU", ["C1'", "P", "O2"])
    same_modification = synthetic_polymer("S", "psu", ["C1'", "P", "O2"])
    second_target = synthetic_polymer("T", "PSU", ["C1'", "P", "O2"])
    different_modification = synthetic_polymer("S", "M2G", ["C1'", "P", "O2"])
    two_reference_chains = StructureData(
        vcat(pseudouridine.atoms, second_reference.atoms),
        pseudouridine.metadata,
    )
    two_target_chains = StructureData(
        vcat(same_modification.atoms, second_target.atoms),
        same_modification.metadata,
    )
    @test_throws ArgumentError match_residues(
        two_reference_chains,
        same_modification;
        chain_map = Dict("R" => "S", "Q" => "S"),
    )
    @test_throws ArgumentError match_residues(
        pseudouridine,
        two_target_chains;
        chain_map = Dict{Any,Any}(:R => "S", "R" => "T"),
    )
    same_pairs, same_unmatched_reference, same_unmatched_target =
        match_residues(pseudouridine, same_modification; chain_map = Dict("R" => "S"))
    @test length(same_pairs) == 1
    @test isempty(same_unmatched_reference)
    @test isempty(same_unmatched_target)

    different_pairs, different_unmatched_reference, different_unmatched_target =
        match_residues(pseudouridine, different_modification; chain_map = Dict("R" => "S"))
    @test isempty(different_pairs)
    @test length(different_unmatched_reference) == 1
    @test length(different_unmatched_target) == 1

    alanine = synthetic_polymer("A", "ALA", ["N", "CA", "C"])
    adenosine = synthetic_polymer("A", "A", ["C1'", "P", "O2"])
    automatic_pairs, _, _ = match_residues(alanine, adenosine)
    @test isempty(automatic_pairs)
    @test_throws ArgumentError match_residues(
        alanine,
        adenosine;
        chain_map = Dict("A" => "A"),
    )
end

@testset "Sequence-aware structure alignment" begin
    path = joinpath(@__DIR__, "..", "fixtures", "synthetic.pdb")
    reference = load_structure(path)
    target = load_structure(path)
    result = align_structures(reference, target)
    @test result.rmsd_angstrom < 1e-12
    @test result.aligned_atoms == 7
    @test result.correspondence_coverage == 1.0
end
