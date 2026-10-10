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

@testset "Sequence-aware structure alignment" begin
    path = joinpath(@__DIR__, "..", "fixtures", "synthetic.pdb")
    reference = load_structure(path)
    target = load_structure(path)
    result = align_structures(reference, target)
    @test result.rmsd_angstrom < 1e-12
    @test result.aligned_atoms == 7
    @test result.correspondence_coverage == 1.0
end
