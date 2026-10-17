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
    @test minimum_distance_series(trajectory, [1, 2], [3, 4]) ≈ ones(3) atol = 1e-12
    separated = TrajectoryData(reshape([0.0, 0, 0, 4, 0, 0], 3, 2, 1))
    @test minimum_distance_series(separated, [1], [2]) == [4.0]
    displacement = group_centroid_displacement(trajectory, [4]; align_indices = [1, 2, 3])
    @test displacement[1:2] ≈ zeros(2) atol = 1e-12
    @test displacement[3] > 0
    pca = trajectory_pca(trajectory)
    @test size(pca.scores, 2) == 3
    clustering = cluster_frames(pca.scores[1:2, :], 2)
    @test length(clustering.labels) == 3
    @test length(representative_frames(pca.scores[1:2, :], clustering)) == 2

    identical = TrajectoryData(repeat(reshape(reference, 3, 4, 1), 1, 1, 2))
    @test all(iszero, trajectory_pca(identical).explained_fraction)

    structure = load_structure(joinpath(@__DIR__, "..", "fixtures", "synthetic.pdb"))
    structure_trajectory = TrajectoryData(
        repeat(reshape(coordinates(structure), 3, length(structure.atoms), 1), 1, 1, 2),
    )
    residue_fluctuations = residue_rmsf(structure, structure_trajectory)
    @test nrow(residue_fluctuations) == length(residue_records(structure))
    @test all(iszero, residue_fluctuations.rmsf_angstrom)
end

@testset "Residue contact persistence" begin
    path = joinpath(@__DIR__, "..", "fixtures", "synthetic.pdb")
    structure = load_structure(path)
    ligand =
        select_ligand(structure; residue_name = "PAR", chain_id = "L", residue_number = 100)
    reference = Matrix(coordinates(structure))
    displaced = copy(reference)
    displaced[:, ligand.atom_indices] .+= 40.0
    trajectory = TrajectoryData(
        cat(reference, displaced; dims = 3);
        times_ps = [0.0, 10.0],
        atom_indices = collect(eachindex(structure.atoms)),
    )
    fingerprint = residue_contact_persistence(
        structure,
        trajectory,
        ligand.atom_indices;
        cutoff = 5.0,
    )
    @test nrow(fingerprint) == length(residue_records(structure)) - 1
    @test all(==(2), fingerprint.sampled_frames)
    @test all(0 .<= fingerprint.contact_fraction .<= 1)
    @test any(==(0.5), fingerprint.contact_fraction)
    @test !("PAR" in fingerprint.residue_name)
    @test_throws ArgumentError residue_contact_persistence(
        structure,
        trajectory,
        ligand.atom_indices;
        partner_indices = ligand.atom_indices,
    )
    @test_throws ArgumentError residue_contact_persistence(
        structure,
        trajectory,
        ligand.atom_indices;
        periodic = true,
    )
end

@testset "Unit conversion" begin
    trajectory = TrajectoryData(ones(3, 3, 2); input_unit = :nanometer)
    @test all(==(10.0), trajectory.coordinates)
    @test trajectory.unit_conversion_to_angstrom == 10.0
    @test_throws ArgumentError TrajectoryData(zeros(3, 0, 1))
    @test_throws ArgumentError TrajectoryData(zeros(3, 1, 0))
    @test_throws ArgumentError TrajectoryData(fill(NaN, 3, 1, 1))
    @test_throws ArgumentError TrajectoryData(zeros(3, 1, 2); times_ps = [1.0, 0.0])
    @test_throws ArgumentError TrajectoryData(zeros(3, 2, 1); atom_indices = [1, 1])
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
    @test minimum_distance_series(trajectory, [1], [2]; periodic = true) ≈ [0.8, 1.1, 0.7]

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
    @test_throws ArgumentError minimum_distance_series(
        without_cell,
        [1],
        [2];
        periodic = true,
    )
    @test_throws ArgumentError contact_persistence(without_cell, Int[], [2])
    @test_throws ArgumentError trajectory_rmsd(without_cell; align_indices = [1, 2])
end

@testset "Unit-cell validation" begin
    matrix_cell = Diagonal([10.0, 11.0, 12.0])
    repeated = TrajectoryData(zeros(3, 1, 2); unit_cells = matrix_cell)
    @test repeated.unit_cells[:, :, 1] == matrix_cell
    @test repeated.unit_cells[:, :, 2] == matrix_cell

    frame_cells = cat(matrix_cell, Diagonal([20.0, 21.0, 22.0]); dims = 3)
    per_frame = TrajectoryData(zeros(3, 1, 2); unit_cells = frame_cells)
    @test per_frame.unit_cells == frame_cells

    @test_throws DimensionMismatch TrajectoryData(zeros(3, 1, 2); unit_cells = ones(6, 3))
    @test_throws DimensionMismatch TrajectoryData(
        zeros(3, 1, 2);
        unit_cells = ones(3, 3, 3),
    )
    @test_throws DimensionMismatch TrajectoryData(zeros(3, 1, 1); unit_cells = ones(2, 2))
    @test_throws ArgumentError TrajectoryData(
        zeros(3, 1, 1);
        unit_cells = [10.0, 10.0, 0.0, 90.0, 90.0, 90.0],
    )
    @test_throws ArgumentError TrajectoryData(
        zeros(3, 1, 1);
        unit_cells = [10.0 NaN 0.0; 0.0 10.0 0.0; 0.0 0.0 10.0],
    )
    @test_throws ArgumentError periodic_distance(zeros(3), ones(3), zeros(3, 3))
end
