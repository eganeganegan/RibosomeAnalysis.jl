module RibosomeAnalysisPythonCallExt

using RibosomeAnalysis
using PythonCall
using LinearAlgebra

function _universe(topology::AbstractString, trajectory::AbstractString)
    isfile(topology) || throw(ArgumentError("topology file does not exist: $topology"))
    isfile(trajectory) ||
        throw(ArgumentError("trajectory file does not exist: $trajectory"))
    mda = try
        pyimport("MDAnalysis")
    catch error
        throw(
            ArgumentError(
                "MDAnalysis is not importable in PythonCall's environment: $(sprint(showerror, error))",
            ),
        )
    end
    try
        mda.Universe(topology, trajectory)
    catch error
        throw(
            ArgumentError(
                "MDAnalysis could not combine topology $(abspath(topology)) with trajectory $(abspath(trajectory)): $(sprint(showerror, error))",
            ),
        )
    end
end

function _frame_indices(total::Int, start::Integer, stop, step::Integer)
    start >= 1 || throw(ArgumentError("start must be at least 1"))
    step >= 1 || throw(ArgumentError("step must be positive"))
    final_frame = isnothing(stop) ? total : min(Int(stop), total)
    final_frame >= start || throw(ArgumentError("frame slice is empty"))
    collect(Int(start):Int(step):final_frame)
end

function _unit_cell(time_step)
    vectors = time_step.triclinic_dimensions
    pyis(vectors, pybuiltins.None) && return nothing
    row_vectors = pyconvert(Matrix{Float64}, vectors)
    size(row_vectors) == (3, 3) ||
        throw(DimensionMismatch("MDAnalysis returned a unit cell that is not 3 × 3"))
    all(isfinite, row_vectors) && abs(det(row_vectors)) > eps(Float64) || return nothing
    Matrix(row_vectors')
end

function _selection(universe, selection::AbstractString)
    atoms = universe.select_atoms(selection)
    pyconvert(Int, length(atoms)) > 0 ||
        throw(ArgumentError("MDAnalysis selection matched no atoms: $selection"))
    atoms
end


function load_mdanalysis(
    topology::AbstractString,
    trajectory::AbstractString;
    selection::AbstractString = "all",
    start::Integer = 1,
    stop = nothing,
    step::Integer = 1,
)
    universe = _universe(topology, trajectory)
    atoms = _selection(universe, selection)
    total = pyconvert(Int, length(universe.trajectory))
    frames = _frame_indices(total, start, stop, step)
    atom_count = pyconvert(Int, length(atoms))
    xyz = Array{Float64}(undef, 3, atom_count, length(frames))
    times = Vector{Float64}(undef, length(frames))
    cells = Vector{Union{Nothing,Matrix{Float64}}}(undef, length(frames))
    for (out_frame, one_based_frame) in enumerate(frames)
        time_step = universe.trajectory[one_based_frame-1]
        positions = pyconvert(Matrix{Float64}, atoms.positions)
        xyz[:, :, out_frame] = positions'
        times[out_frame] = pyconvert(Float64, time_step.time)
        cells[out_frame] = _unit_cell(time_step)
    end
    has_cells = .!isnothing.(cells)
    all(has_cells) ||
        all(.!has_cells) ||
        throw(ArgumentError("unit-cell data are present for only some selected frames"))
    unit_cells =
        all(has_cells) ? cat((cell for cell in cells if !isnothing(cell))...; dims = 3) :
        nothing
    atom_ids = pyconvert(Vector{Int}, atoms.indices) .+ 1
    RibosomeAnalysis.TrajectoryData(
        xyz;
        times_ps = times,
        atom_indices = atom_ids,
        source = abspath(trajectory),
        input_unit = :angstrom,
        unit_cells = unit_cells,
    )
end

function foreach_mdanalysis_frame(
    callback::Function,
    topology::AbstractString,
    trajectory::AbstractString;
    selection::AbstractString = "all",
    start::Integer = 1,
    stop = nothing,
    step::Integer = 1,
)
    universe = _universe(topology, trajectory)
    atoms = _selection(universe, selection)
    total = pyconvert(Int, length(universe.trajectory))
    frames = _frame_indices(total, start, stop, step)
    atom_ids = pyconvert(Vector{Int}, atoms.indices) .+ 1
    for one_based_frame in frames
        time_step = universe.trajectory[one_based_frame-1]
        positions = pyconvert(Matrix{Float64}, atoms.positions)
        frame = RibosomeAnalysis.TrajectoryFrame(
            Matrix(positions'),
            pyconvert(Float64, time_step.time),
            copy(atom_ids),
            one_based_frame,
            _unit_cell(time_step),
        )
        callback(frame)
    end
    nothing
end

end
