function TrajectoryData(coordinates::AbstractArray{<:Real,3}; times_ps=nothing,
                        atom_indices=collect(1:size(coordinates, 2)), source="in-memory",
                        input_unit::Symbol=:angstrom)
    size(coordinates, 1) == 3 || throw(DimensionMismatch("trajectory must have shape 3 × atoms × frames"))
    factor = input_unit === :angstrom ? 1.0 : input_unit === :nanometer ? 10.0 :
             throw(ArgumentError("input_unit must be :angstrom or :nanometer"))
    nframes = size(coordinates, 3)
    times = isnothing(times_ps) ? Float64.(0:nframes-1) : Float64.(times_ps)
    length(times) == nframes || throw(DimensionMismatch("one time value is required per frame"))
    length(atom_indices) == size(coordinates, 2) || throw(DimensionMismatch("one atom index is required per trajectory atom"))
    TrajectoryData(Float64.(coordinates) .* factor, times, Int.(atom_indices), string(source), factor)
end

"""Extension point for trajectory readers. Install PythonCall and MDAnalysis for XTC/TRR/DCD."""
function load_trajectory(topology::AbstractString, trajectory::AbstractString; kwargs...)
    extension = Base.get_extension(@__MODULE__, :RibosomeAnalysisPythonCallExt)
    isnothing(extension) && throw(ArgumentError("trajectory integration is unavailable; install PythonCall and MDAnalysis, then restart Julia"))
    extension.load_mdanalysis(topology, trajectory; kwargs...)
end

function _fit_frame(reference, frame, align_indices)
    rotation, translation, _ = kabsch(reference[:, align_indices], frame[:, align_indices])
    apply_transform(frame, rotation, translation)
end

function trajectory_rmsd(trajectory::TrajectoryData;
                         reference_frame::Int=1,
                         align_indices=axes(trajectory.coordinates, 2),
                         measure_indices=align_indices)
    checkbounds(trajectory.coordinates, :, :, reference_frame)
    reference = trajectory.coordinates[:, :, reference_frame]
    [rmsd(reference[:, measure_indices],
          _fit_frame(reference, trajectory.coordinates[:, :, frame], align_indices)[:, measure_indices])
     for frame in axes(trajectory.coordinates, 3)]
end

function trajectory_rmsf(trajectory::TrajectoryData;
                         reference_frame::Int=1,
                         align_indices=axes(trajectory.coordinates, 2),
                         measure_indices=axes(trajectory.coordinates, 2))
    reference = trajectory.coordinates[:, :, reference_frame]
    fitted = Array{Float64}(undef, 3, length(measure_indices), size(trajectory.coordinates, 3))
    for frame in axes(trajectory.coordinates, 3)
        fitted[:, :, frame] = _fit_frame(reference, trajectory.coordinates[:, :, frame], align_indices)[:, measure_indices]
    end
    average = mean(fitted; dims=3)
    [sqrt(mean(sum(abs2, fitted[:, i, frame] - average[:, i, 1]) for frame in axes(fitted, 3)))
     for i in axes(fitted, 2)]
end

function contact_persistence(trajectory::TrajectoryData, group_a, group_b; cutoff::Real=4.0)
    a, b = collect(group_a), collect(group_b)
    counts = zeros(Int, length(a), length(b))
    for frame in axes(trajectory.coordinates, 3)
        tree = KDTree(trajectory.coordinates[:, b, frame])
        for (local_i, i) in enumerate(a)
            for local_j in inrange(tree, trajectory.coordinates[:, i, frame], Float64(cutoff))
                counts[local_i, local_j] += 1
            end
        end
    end
    counts ./ size(trajectory.coordinates, 3)
end

function distance_series(trajectory::TrajectoryData, atom_i::Int, atom_j::Int)
    [norm(trajectory.coordinates[:, atom_i, frame] - trajectory.coordinates[:, atom_j, frame])
     for frame in axes(trajectory.coordinates, 3)]
end

function trajectory_pca(trajectory::TrajectoryData; atom_indices=axes(trajectory.coordinates, 2),
                        align::Bool=true, reference_frame::Int=1)
    frames = size(trajectory.coordinates, 3)
    matrix = Matrix{Float64}(undef, 3 * length(atom_indices), frames)
    reference = trajectory.coordinates[:, :, reference_frame]
    for frame in 1:frames
        xyz = align ? _fit_frame(reference, trajectory.coordinates[:, :, frame], atom_indices) : trajectory.coordinates[:, :, frame]
        matrix[:, frame] = vec(xyz[:, atom_indices])
    end
    means = mean(matrix; dims=2)
    centered = matrix .- means
    decomposition = svd(centered / sqrt(max(frames - 1, 1)); full=false)
    scores = decomposition.U' * centered
    variance = decomposition.S .^ 2
    (scores=scores, components=decomposition.U, variance=variance,
     explained_fraction=variance ./ sum(variance), mean_coordinates=vec(means))
end

"""Deterministic Lloyd clustering of columns in a feature matrix."""
function cluster_frames(features::AbstractMatrix, k::Integer; max_iterations::Integer=100)
    n = size(features, 2)
    1 <= k <= n || throw(ArgumentError("k must be between 1 and the number of frames"))
    seeds = unique(round.(Int, range(1, n; length=k)))
    centers = copy(features[:, seeds])
    labels = ones(Int, n)
    for _ in 1:max_iterations
        new_labels = [argmin([sum(abs2, features[:, i] - centers[:, j]) for j in 1:k]) for i in 1:n]
        new_labels == labels && break
        labels = new_labels
        for j in 1:k
            members = findall(==(j), labels)
            !isempty(members) && (centers[:, j] = vec(mean(features[:, members]; dims=2)))
        end
    end
    (labels=labels, centers=centers)
end

function representative_frames(features::AbstractMatrix, clustering)
    [argmin([sum(abs2, features[:, i] - clustering.centers[:, j]) for i in axes(features, 2)])
     for j in axes(clustering.centers, 2)]
end
