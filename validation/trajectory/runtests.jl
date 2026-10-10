using JSON3
using PythonCall
using RibosomeAnalysis
using Test

data_directory = length(ARGS) >= 1 ? abspath(ARGS[1]) : joinpath(@__DIR__, "data")
gromacs = if length(ARGS) >= 2
    ARGS[2]
else
    executable = Sys.which("gmx")
    isnothing(executable) ? Sys.which("gmx_mpi") : executable
end
isdir(data_directory) || error("trajectory validation data do not exist: $data_directory")
isnothing(gromacs) && error("gmx or gmx_mpi is required for independent validation")

reference = JSON3.read(read(joinpath(data_directory, "reference.json"), String))
frames = length(reference.coordinates_angstrom)
atoms = length(reference.coordinates_angstrom[1])
expected_coordinates = Array{Float64}(undef, 3, atoms, frames)
for frame = 1:frames, atom = 1:atoms, dimension = 1:3
    expected_coordinates[dimension, atom, frame] =
        reference.coordinates_angstrom[frame][atom][dimension]
end
expected_times = Float64.(reference.times_ps)
expected_periodic = Float64.(reference.periodic_distance_1_2_angstrom)
topology = joinpath(data_directory, "topology.pdb")

@testset "Validation tool provenance" begin
    @test string(reference.mdanalysis_version) == "2.10.0"
    @test pyconvert(String, pyimport("MDAnalysis").__version__) == "2.10.0"
    gromacs_version = read(`$(gromacs) --version`, String)
    version_match = match(r"GROMACS version:\s+(\S+)", gromacs_version)
    @test !isnothing(version_match)
    println("GROMACS validation version: ", version_match.captures[1])
end

@testset "MDAnalysis production trajectory formats" begin
    tolerances = Dict("xtc" => 1.1e-2, "trr" => 2e-5, "dcd" => 2e-5)
    for format in ("xtc", "trr", "dcd")
        path = joinpath(data_directory, "trajectory.$format")
        trajectory = load_trajectory(topology, path)
        @test size(trajectory.coordinates) == (3, atoms, frames)
        @test trajectory.atom_indices == collect(1:atoms)
        @test trajectory.times_ps ≈ expected_times atol = 2e-4
        @test trajectory.unit_conversion_to_angstrom == 1.0
        @test !isnothing(trajectory.unit_cells)
        @test maximum(abs, trajectory.coordinates - expected_coordinates) <=
              tolerances[format]
        @test distance_series(trajectory, 1, 2; periodic = true) ≈ expected_periodic atol =
            1.2e-2
    end
end

@testset "Frame slicing, atom selection, and streaming" begin
    path = joinpath(data_directory, "trajectory.xtc")
    sliced = load_trajectory(
        topology,
        path;
        selection = "index 0 1",
        start = 2,
        stop = 5,
        step = 2,
    )
    @test size(sliced.coordinates) == (3, 2, 2)
    @test sliced.atom_indices == [1, 2]
    @test sliced.times_ps ≈ expected_times[[2, 4]] atol = 2e-4

    visited = Int[]
    streamed_distances = Float64[]
    foreach_trajectory_frame(topology, path; start = 2, stop = 5, step = 2) do frame
        push!(visited, frame.frame_index)
        @test size(frame.coordinates) == (3, atoms)
        @test !isnothing(frame.unit_cell)
        push!(
            streamed_distances,
            periodic_distance(
                frame.coordinates[:, 1],
                frame.coordinates[:, 2],
                frame.unit_cell,
            ),
        )
    end
    @test visited == [2, 4]
    @test streamed_distances ≈ expected_periodic[[2, 4]] atol = 1.2e-2

    @test_throws ArgumentError load_trajectory(
        topology,
        path;
        selection = "name DOES_NOT_EXIST",
    )
    @test_throws ArgumentError load_trajectory(topology, path; start = 0)
    @test_throws ArgumentError load_trajectory(topology, path; step = 0)
    @test_throws ArgumentError load_trajectory(
        joinpath(data_directory, "wrong-topology.pdb"),
        path,
    )
end

function read_xvg_distances(path)
    values = Float64[]
    for line in eachline(path)
        stripped = strip(line)
        (isempty(stripped) || startswith(stripped, ('#', '@'))) && continue
        columns = split(stripped)
        push!(values, parse(Float64, columns[2]))
    end
    values
end

@testset "Independent GROMACS periodic-distance reference" begin
    output = joinpath(data_directory, "gromacs-distance.xvg")
    command = `$(gromacs) distance -f $(joinpath(data_directory, "trajectory.xtc")) -s $topology -select "atomnr 1 plus atomnr 2" -oall $output -xvg none -pbc yes`
    run(pipeline(command; stdout = devnull, stderr = devnull))
    gromacs_distances_angstrom = 10 .* read_xvg_distances(output)
    julia_trajectory = load_trajectory(topology, joinpath(data_directory, "trajectory.xtc"))
    julia_distances = distance_series(julia_trajectory, 1, 2; periodic = true)
    @test gromacs_distances_angstrom ≈ expected_periodic atol = 1.2e-2
    @test julia_distances ≈ gromacs_distances_angstrom atol = 1.2e-2
end
