module RibosomeAnalysisPythonCallExt

using RibosomeAnalysis
using PythonCall

function load_mdanalysis(topology::AbstractString, trajectory::AbstractString;
                         selection::AbstractString="all",
                         start::Integer=1, stop=nothing, step::Integer=1)
    isfile(topology) || throw(ArgumentError("topology file does not exist: $topology"))
    isfile(trajectory) || throw(ArgumentError("trajectory file does not exist: $trajectory"))
    mda = try
        pyimport("MDAnalysis")
    catch error
        throw(ArgumentError("MDAnalysis is not importable in PythonCall's environment: $(sprint(showerror, error))"))
    end
    universe = mda.Universe(topology, trajectory)
    atoms = universe.select_atoms(selection)
    total = pyconvert(Int, length(universe.trajectory))
    final_frame = isnothing(stop) ? total : min(Int(stop), total)
    frames = collect(Int(start):Int(step):final_frame)
    isempty(frames) && throw(ArgumentError("frame slice is empty"))
    atom_count = pyconvert(Int, length(atoms))
    xyz = Array{Float64}(undef, 3, atom_count, length(frames))
    times = Vector{Float64}(undef, length(frames))
    for (out_frame, one_based_frame) in enumerate(frames)
        time_step = universe.trajectory[one_based_frame - 1]
        positions = pyconvert(Matrix{Float64}, atoms.positions)
        xyz[:, :, out_frame] = positions'
        times[out_frame] = pyconvert(Float64, time_step.time)
    end
    atom_ids = pyconvert(Vector{Int}, atoms.indices) .+ 1
    RibosomeAnalysis.TrajectoryData(xyz; times_ps=times, atom_indices=atom_ids,
                                    source=abspath(trajectory), input_unit=:angstrom)
end

end
