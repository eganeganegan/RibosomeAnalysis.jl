using JSON3

@testset "Independent MDAnalysis trajectory reference" begin
    reference_path =
        joinpath(@__DIR__, "..", "..", "validation", "reference", "mdanalysis-2.10.0.json")
    reference_data = JSON3.read(read(reference_path, String))
    expected = reference_data.trajectory

    coordinates = [0.0 1 0 0; 0 0 1 0; 0 0 0 1]
    rotation = [0.0 -1 0; 1 0 0; 0 0 1]
    translated = rotation * coordinates .+ [2.0, 3.0, 4.0]
    deformed = copy(translated)
    deformed[:, 4] .+= [0.1, 0.0, 0.0]
    trajectory = TrajectoryData(cat(coordinates, translated, deformed; dims = 3))

    @test trajectory_rmsd(trajectory) ≈ Float64.(expected.rmsd_angstrom) atol = 5e-8
    @test trajectory_rmsf(trajectory) ≈ Float64.(expected.rmsf_angstrom) atol = 5e-8
    @test distance_series(trajectory, 1, 2) ≈ Float64.(expected.distance_1_2_angstrom) atol =
        5e-8
end
