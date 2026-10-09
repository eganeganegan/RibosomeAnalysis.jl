@testset "Trajectory analysis" begin
    reference = [0.0 1 0 0; 0 0 1 0; 0 0 0 1]
    rotation = [0.0 -1 0; 1 0 0; 0 0 1]
    translated = rotation * reference .+ [2.0, 3.0, 4.0]
    deformed = copy(translated); deformed[:, 4] .+= [0.1, 0.0, 0.0]
    trajectory = TrajectoryData(cat(reference, translated, deformed; dims=3); times_ps=[0, 5, 10])
    values = trajectory_rmsd(trajectory)
    @test values[1] < 1e-12
    @test values[2] < 1e-12
    @test values[3] > 0
    fluctuations = trajectory_rmsf(trajectory)
    @test length(fluctuations) == 4
    @test all(>=(0), fluctuations)
    distances = distance_series(trajectory, 1, 2)
    @test length(distances) == 3
    persistence = contact_persistence(trajectory, [1], [2]; cutoff=1.1)
    @test persistence[1,1] == 1.0
    pca = trajectory_pca(trajectory)
    @test size(pca.scores, 2) == 3
    clustering = cluster_frames(pca.scores[1:2, :], 2)
    @test length(clustering.labels) == 3
    @test length(representative_frames(pca.scores[1:2, :], clustering)) == 2
end

@testset "Unit conversion" begin
    trajectory = TrajectoryData(ones(3, 3, 2); input_unit=:nanometer)
    @test all(==(10.0), trajectory.coordinates)
    @test trajectory.unit_conversion_to_angstrom == 10.0
end
