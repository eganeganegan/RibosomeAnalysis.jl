@testset "Trajectory analysis" begin
    reference = [0.0 1 0 0; 0 0 1 0; 0 0 0 1]
    rotation = [0.0 -1 0; 1 0 0; 0 0 1]
    translated = rotation * reference .+ [2.0, 3.0, 4.0]
    deformed = copy(translated)
    deformed[:, 4] .+= [0.1, 0.0, 0.0]
    trajectory = TrajectoryData(
        cat(reference, translated, deformed; dims = 3);
        times_ps = [0, 5, 10],
    )
    values = trajectory_rmsd(trajectory)
    @test values[1] < 1e-12
    @test values[2] < 1e-12
    @test values[3] > 0
    fluctuations = trajectory_rmsf(trajectory)
    @test length(fluctuations) == 4
    @test all(>=(0), fluctuations)
    distances = distance_series(trajectory, 1, 2)
    @test length(distances) == 3
    persistence = contact_persistence(trajectory, [1], [2]; cutoff = 1.1)
    @test persistence[1, 1] == 1.0
    pca = trajectory_pca(trajectory)
    @test size(pca.scores, 2) == 3
    clustering = cluster_frames(pca.scores[1:2, :], 2)
    @test length(clustering.labels) == 3
    @test length(representative_frames(pca.scores[1:2, :], clustering)) == 2
end

@testset "Unit conversion" begin
    trajectory = TrajectoryData(ones(3, 3, 2); input_unit = :nanometer)
    @test all(==(10.0), trajectory.coordinates)
    @test trajectory.unit_conversion_to_angstrom == 10.0
end

@testset "Periodic trajectory geometry" begin
    coordinates = zeros(3, 2, 3)
    coordinates[:, 1, :] .= [0.5, 0.5, 0.5]
    coordinates[:, 2, 1] .= [9.7, 0.5, 0.5]
    coordinates[:, 2, 2] .= [9.4, 0.5, 0.5]
    coordinates[:, 2, 3] .= [9.8, 0.5, 0.5]
    dimensions = repeat([10.0, 10.0, 10.0, 90.0, 90.0, 90.0], 1, 3)
    trajectory = TrajectoryData(coordinates; unit_cells = dimensions)

    @test size(trajectory.unit_cells) == (3, 3, 3)
    @test distance_series(trajectory, 1, 2) ≈ [9.2, 8.9, 9.3]
    @test distance_series(trajectory, 1, 2; periodic = true) ≈ [0.8, 1.1, 0.7]
    @test contact_persistence(trajectory, [1], [2]; cutoff = 1.1)[1, 1] == 0.0
    @test contact_persistence(trajectory, [1], [2]; cutoff = 1.1, periodic = true)[1, 1] ==
          1.0

    triclinic = [10.0 2.0 0.0; 0.0 9.0 0.0; 0.0 0.0 8.0]
    point_a = triclinic * [0.05, 0.05, 0.5]
    point_b = triclinic * [0.95, 0.95, 0.5]
    expected = norm(triclinic * [0.1, 0.1, 0.0])
    @test periodic_distance(point_a, point_b, triclinic) ≈ expected atol = 1e-12
    @test norm(minimum_image_displacement(point_a - point_b, triclinic)) ≈ expected atol =
        1e-12

    nanometer_trajectory = TrajectoryData(
        zeros(3, 1, 2);
        input_unit = :nanometer,
        unit_cells = [1.0, 2.0, 3.0, 90.0, 90.0, 90.0],
    )
    @test nanometer_trajectory.unit_cells[:, :, 1] ≈ Diagonal([10.0, 20.0, 30.0])

    without_cell = TrajectoryData(zeros(3, 2, 1))
    @test_throws ArgumentError distance_series(without_cell, 1, 2; periodic = true)
    @test_throws ArgumentError contact_persistence(without_cell, [1], [2]; periodic = true)
end

@testset "Unit-cell validation" begin
    @test_throws DimensionMismatch TrajectoryData(zeros(3, 1, 2); unit_cells = ones(6, 3))
    @test_throws ArgumentError TrajectoryData(
        zeros(3, 1, 1);
        unit_cells = [10.0, 10.0, 0.0, 90.0, 90.0, 90.0],
    )
    @test_throws ArgumentError periodic_distance(zeros(3), ones(3), zeros(3, 3))
end
