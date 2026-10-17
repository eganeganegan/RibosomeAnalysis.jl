function TrajectoryData(
    coordinates::AbstractArray{<:Real,3};
    times_ps = nothing,
    atom_indices = collect(1:size(coordinates, 2)),
    source = "in-memory",
    input_unit::Symbol = :angstrom,
    unit_cells = nothing,
)
    size(coordinates, 1) == 3 ||
        throw(DimensionMismatch("trajectory must have shape 3 × atoms × frames"))
    factor =
        input_unit === :angstrom ? 1.0 :
        input_unit === :nanometer ? 10.0 :
        throw(ArgumentError("input_unit must be :angstrom or :nanometer"))
    size(coordinates, 2) > 0 || throw(ArgumentError("trajectory contains no atoms"))
    nframes = size(coordinates, 3)
    nframes > 0 || throw(ArgumentError("trajectory contains no frames"))
    times = isnothing(times_ps) ? Float64.(0:(nframes-1)) : Float64.(times_ps)
    length(times) == nframes ||
        throw(DimensionMismatch("one time value is required per frame"))
    all(isfinite, times) || throw(ArgumentError("trajectory times must be finite"))
    issorted(times) || throw(ArgumentError("trajectory times must be nondecreasing"))
    length(atom_indices) == size(coordinates, 2) ||
        throw(DimensionMismatch("one atom index is required per trajectory atom"))
    indices = Int.(atom_indices)
    all(>(0), indices) || throw(ArgumentError("trajectory atom indices must be positive"))
    allunique(indices) || throw(ArgumentError("trajectory atom indices must be unique"))
    converted_coordinates = Float64.(coordinates) .* factor
    all(isfinite, converted_coordinates) ||
        throw(ArgumentError("trajectory coordinates must be finite"))
    cells = _normalize_unit_cells(unit_cells, nframes, factor)
    TrajectoryData(converted_coordinates, times, indices, string(source), factor, cells)
end

function _cell_from_lengths_angles(dimensions::AbstractVector{<:Real})
    length(dimensions) == 6 || throw(
        DimensionMismatch("unit-cell dimensions must contain a, b, c, alpha, beta, gamma"),
    )
    a, b, c, alpha, beta, gamma = Float64.(dimensions)
    all(>(0), (a, b, c)) || throw(ArgumentError("unit-cell lengths must be positive"))
    all(angle -> 0 < angle < 180, (alpha, beta, gamma)) ||
        throw(ArgumentError("unit-cell angles must lie strictly between 0 and 180 degrees"))
    alpha_r, beta_r, gamma_r = deg2rad.((alpha, beta, gamma))
    sin_gamma = sin(gamma_r)
    abs(sin_gamma) > sqrt(eps(Float64)) ||
        throw(ArgumentError("unit-cell gamma angle produces a singular cell"))
    c_x = c * cos(beta_r)
    c_y = c * (cos(alpha_r) - cos(beta_r) * cos(gamma_r)) / sin_gamma
    c_z_squared = c^2 - c_x^2 - c_y^2
    c_z_squared > 0 || throw(ArgumentError("unit-cell dimensions produce zero volume"))
    [a b*cos(gamma_r) c_x; 0.0 b*sin_gamma c_y; 0.0 0.0 sqrt(c_z_squared)]
end

function _validate_cell(cell::AbstractMatrix{<:Real})
    size(cell) == (3, 3) || throw(DimensionMismatch("a unit cell must be a 3 × 3 matrix"))
    all(isfinite, cell) || throw(ArgumentError("unit-cell vectors must be finite"))
    volume = abs(det(Matrix{Float64}(cell)))
    scale = maximum(abs, cell)^3
    volume > max(eps(Float64), scale * 1e-12) ||
        throw(ArgumentError("unit-cell vectors must span a nonzero volume"))
    nothing
end

function _normalize_unit_cells(unit_cells, nframes::Int, factor::Float64)
    isnothing(unit_cells) && return nothing
    cells = if ndims(unit_cells) == 2 && size(unit_cells) == (3, 3)
        repeat(reshape(Float64.(unit_cells), 3, 3, 1), 1, 1, nframes)
    elseif ndims(unit_cells) == 3 && size(unit_cells, 1) == 3 && size(unit_cells, 2) == 3
        size(unit_cells, 3) == nframes ||
            throw(DimensionMismatch("one 3 × 3 unit cell is required per frame"))
        Float64.(unit_cells)
    elseif ndims(unit_cells) == 1 && length(unit_cells) == 6
        repeat(reshape(_cell_from_lengths_angles(unit_cells), 3, 3, 1), 1, 1, nframes)
    elseif ndims(unit_cells) == 2 && size(unit_cells, 1) == 6
        size(unit_cells, 2) == nframes ||
            throw(DimensionMismatch("one six-value unit cell is required per frame"))
        cat(
            (_cell_from_lengths_angles(unit_cells[:, frame]) for frame = 1:nframes)...;
            dims = 3,
        )
    else
        throw(
            DimensionMismatch(
                "unit_cells must be 3 × 3, 3 × 3 × frames, length 6, or 6 × frames",
            ),
        )
    end
    cells .*= factor
    for frame = 1:nframes
        _validate_cell(@view cells[:, :, frame])
    end
    cells
end

"""Extension point for trajectory readers. Install PythonCall and MDAnalysis for XTC/TRR/DCD."""
function load_trajectory(topology::AbstractString, trajectory::AbstractString; kwargs...)
    extension = Base.get_extension(@__MODULE__, :RibosomeAnalysisPythonCallExt)
    isnothing(extension) && throw(
        ArgumentError(
            "trajectory integration is unavailable; install PythonCall and MDAnalysis, then restart Julia",
        ),
    )
    extension.load_mdanalysis(topology, trajectory; kwargs...)
end

"""
Visit trajectory frames one at a time through the optional MDAnalysis reader.

The callback receives a `TrajectoryFrame`. Its coordinate and cell arrays are
independent copies and remain valid after the callback returns.
"""
function foreach_trajectory_frame(
    callback::Function,
    topology::AbstractString,
    trajectory::AbstractString;
    kwargs...,
)
    extension = Base.get_extension(@__MODULE__, :RibosomeAnalysisPythonCallExt)
    isnothing(extension) && throw(
        ArgumentError(
            "trajectory integration is unavailable; install PythonCall and MDAnalysis, then restart Julia",
        ),
    )
    extension.foreach_mdanalysis_frame(callback, topology, trajectory; kwargs...)
end

"""
Return the shortest periodic displacement for a displacement vector and cell.

Cell columns are lattice vectors in the same unit as the displacement. The
27 neighboring lattice translations around the fractional-coordinate estimate
are searched, which handles orthorhombic and reduced triclinic MD cells.
"""
function minimum_image_displacement(
    displacement::AbstractVector{<:Real},
    unit_cell::AbstractMatrix{<:Real},
)
    length(displacement) == 3 ||
        throw(DimensionMismatch("a displacement must have length 3"))
    _validate_cell(unit_cell)
    cell = Matrix{Float64}(unit_cell)
    delta = Float64.(displacement)
    fractional = cell \ delta
    nearest = round.(Int, fractional)
    best = copy(delta)
    best_squared = sum(abs2, best)
    for i = -1:1, j = -1:1, k = -1:1
        candidate = delta - cell * (nearest .+ (i, j, k))
        candidate_squared = sum(abs2, candidate)
        if candidate_squared < best_squared
            best = candidate
            best_squared = candidate_squared
        end
    end
    best
end

function periodic_distance(
    point_a::AbstractVector{<:Real},
    point_b::AbstractVector{<:Real},
    unit_cell::AbstractMatrix{<:Real},
)
    length(point_a) == 3 && length(point_b) == 3 ||
        throw(DimensionMismatch("periodic-distance points must have length 3"))
    norm(minimum_image_displacement(point_a - point_b, unit_cell))
end

function _fit_frame(reference, frame, align_indices)
    rotation, translation, _ = kabsch(reference[:, align_indices], frame[:, align_indices])
    apply_transform(frame, rotation, translation)
end

function _trajectory_indices(
    trajectory::TrajectoryData,
    indices,
    label::AbstractString;
    minimum_count::Int = 1,
)
    selected = Int.(collect(indices))
    length(selected) >= minimum_count || throw(
        ArgumentError("$label requires at least $minimum_count trajectory atom indices"),
    )
    allunique(selected) || throw(ArgumentError("$label contains duplicate atom indices"))
    atom_count = size(trajectory.coordinates, 2)
    all(i -> 1 <= i <= atom_count, selected) ||
        throw(ArgumentError("$label must contain indices between 1 and $atom_count"))
    return selected
end

function trajectory_rmsd(
    trajectory::TrajectoryData;
    reference_frame::Int = 1,
    align_indices = axes(trajectory.coordinates, 2),
    measure_indices = align_indices,
)
    checkbounds(trajectory.coordinates, :, :, reference_frame)
    alignment_atoms =
        _trajectory_indices(trajectory, align_indices, "align_indices"; minimum_count = 3)
    measured_atoms = _trajectory_indices(trajectory, measure_indices, "measure_indices")
    reference = trajectory.coordinates[:, :, reference_frame]
    [
        rmsd(
            reference[:, measured_atoms],
            _fit_frame(reference, trajectory.coordinates[:, :, frame], alignment_atoms)[
                :,
                measured_atoms,
            ],
        ) for frame in axes(trajectory.coordinates, 3)
    ]
end

function trajectory_rmsf(
    trajectory::TrajectoryData;
    reference_frame::Int = 1,
    align_indices = axes(trajectory.coordinates, 2),
    measure_indices = axes(trajectory.coordinates, 2),
)
    checkbounds(trajectory.coordinates, :, :, reference_frame)
    alignment_atoms =
        _trajectory_indices(trajectory, align_indices, "align_indices"; minimum_count = 3)
    measured_atoms = _trajectory_indices(trajectory, measure_indices, "measure_indices")
    reference = trajectory.coordinates[:, :, reference_frame]
    fitted =
        Array{Float64}(undef, 3, length(measured_atoms), size(trajectory.coordinates, 3))
    for frame in axes(trajectory.coordinates, 3)
        fitted[:, :, frame] =
            _fit_frame(reference, trajectory.coordinates[:, :, frame], alignment_atoms)[
                :,
                measured_atoms,
            ]
    end
    average = mean(fitted; dims = 3)
    [
        sqrt(
            mean(
                sum(abs2, fitted[:, i, frame] - average[:, i, 1]) for
                frame in axes(fitted, 3)
            ),
        ) for i in axes(fitted, 2)
    ]
end

function contact_persistence(
    trajectory::TrajectoryData,
    group_a,
    group_b;
    cutoff::Real = 4.0,
    periodic::Bool = false,
)
    a = _trajectory_indices(trajectory, group_a, "group_a")
    b = _trajectory_indices(trajectory, group_b, "group_b")
    cutoff >= 0 || throw(ArgumentError("cutoff must be nonnegative"))
    periodic &&
        isnothing(trajectory.unit_cells) &&
        throw(ArgumentError("periodic contact persistence requires unit-cell data"))
    counts = zeros(Int, length(a), length(b))
    for frame in axes(trajectory.coordinates, 3)
        if periodic
            cell = @view trajectory.unit_cells[:, :, frame]
            for (local_i, i) in enumerate(a), (local_j, j) in enumerate(b)
                periodic_distance(
                    @view(trajectory.coordinates[:, i, frame]),
                    @view(trajectory.coordinates[:, j, frame]),
                    cell,
                ) <= cutoff && (counts[local_i, local_j] += 1)
            end
        else
            tree = KDTree(trajectory.coordinates[:, b, frame])
            for (local_i, i) in enumerate(a)
                for local_j in
                    inrange(tree, trajectory.coordinates[:, i, frame], Float64(cutoff))
                    counts[local_i, local_j] += 1
                end
            end
        end
    end
    counts ./ size(trajectory.coordinates, 3)
end

"""
Calculate a residue-level geometric contact fingerprint over a trajectory.

`ligand_indices` and `partner_indices` refer to trajectory-local atom indices.
The trajectory's `atom_indices` mapping is used to recover residue identities
from `structure`. A residue is present in a frame when any selected ligand atom
lies within `cutoff` Angstrom of any selected atom in that residue.
"""
function residue_contact_persistence(
    structure::StructureData,
    trajectory::TrajectoryData,
    ligand_indices;
    partner_indices = nothing,
    cutoff::Real = 4.0,
    periodic::Bool = false,
)
    cutoff >= 0 || throw(ArgumentError("cutoff must be nonnegative"))
    ligand = _trajectory_indices(trajectory, ligand_indices, "ligand_indices")
    partners = _trajectory_indices(
        trajectory,
        isnothing(partner_indices) ?
        setdiff(collect(axes(trajectory.coordinates, 2)), ligand) : partner_indices,
        "partner_indices",
    )
    isempty(intersect(Set(ligand), Set(partners))) ||
        throw(ArgumentError("ligand_indices and partner_indices must not overlap"))
    periodic &&
        isnothing(trajectory.unit_cells) &&
        throw(ArgumentError("periodic residue contact persistence requires unit-cell data"))
    all(i -> 1 <= i <= length(structure.atoms), trajectory.atom_indices) ||
        throw(ArgumentError("trajectory atom indices do not map to the supplied structure"))

    molecule_type_by_key =
        Dict(residue.key => residue.molecule_type for residue in residue_records(structure))
    grouped = Dict{ResidueKey,Vector{Int}}()
    for local_index in partners
        topology_index = trajectory.atom_indices[local_index]
        key = residue_key(structure.atoms[topology_index])
        push!(get!(grouped, key, Int[]), local_index)
    end
    isempty(grouped) && throw(ArgumentError("partner_indices contain no residues"))

    residue_keys = sort!(collect(keys(grouped)))
    frame_count = size(trajectory.coordinates, 3)
    contact_counts = Dict(key => 0 for key in residue_keys)
    minimum_distances = Dict(key => Inf for key in residue_keys)
    distance_sums = Dict(key => 0.0 for key in residue_keys)
    for frame in axes(trajectory.coordinates, 3)
        ligand_tree = periodic ? nothing : KDTree(trajectory.coordinates[:, ligand, frame])
        for key in residue_keys
            residue_atoms = grouped[key]
            distance = if periodic
                cell = @view trajectory.unit_cells[:, :, frame]
                minimum(
                    periodic_distance(
                        @view(trajectory.coordinates[:, ligand_atom, frame]),
                        @view(trajectory.coordinates[:, residue_atom, frame]),
                        cell,
                    ) for ligand_atom in ligand for residue_atom in residue_atoms
                )
            else
                minimum(
                    last(
                        knn(ligand_tree, trajectory.coordinates[:, residue_atom, frame], 1),
                    )[1] for residue_atom in residue_atoms
                )
            end
            minimum_distances[key] = min(minimum_distances[key], distance)
            distance_sums[key] += distance
            distance <= cutoff && (contact_counts[key] += 1)
        end
    end

    rows = DataFrame(
        model_id = Int[],
        chain_id = String[],
        residue_number = Int[],
        insertion_code = Char[],
        residue_name = String[],
        molecule_type = Symbol[],
        atom_count = Int[],
        contact_frames = Int[],
        sampled_frames = Int[],
        contact_fraction = Float64[],
        minimum_distance_angstrom = Float64[],
        mean_minimum_distance_angstrom = Float64[],
    )
    for key in residue_keys
        residue_atoms = grouped[key]
        contacts = contact_counts[key]
        push!(
            rows,
            (
                key.model_id,
                key.chain_id,
                key.number,
                key.insertion_code,
                key.name,
                molecule_type_by_key[key],
                length(residue_atoms),
                contacts,
                frame_count,
                contacts / frame_count,
                minimum_distances[key],
                distance_sums[key] / frame_count,
            ),
        )
    end
    return rows
end

function distance_series(
    trajectory::TrajectoryData,
    atom_i::Int,
    atom_j::Int;
    periodic::Bool = false,
)
    _trajectory_indices(trajectory, (atom_i,), "atom_i")
    _trajectory_indices(trajectory, (atom_j,), "atom_j")
    periodic &&
        isnothing(trajectory.unit_cells) &&
        throw(ArgumentError("periodic distances require unit-cell data"))
    [
        periodic ?
        periodic_distance(
            @view(trajectory.coordinates[:, atom_i, frame]),
            @view(trajectory.coordinates[:, atom_j, frame]),
            @view(trajectory.unit_cells[:, :, frame]),
        ) :
        norm(
            @view(trajectory.coordinates[:, atom_i, frame]) -
            @view(trajectory.coordinates[:, atom_j, frame]),
        ) for frame in axes(trajectory.coordinates, 3)
    ]
end

function minimum_distance_series(
    trajectory::TrajectoryData,
    group_a,
    group_b;
    periodic::Bool = false,
)
    a = _trajectory_indices(trajectory, group_a, "group_a")
    b = _trajectory_indices(trajectory, group_b, "group_b")
    periodic &&
        isnothing(trajectory.unit_cells) &&
        throw(ArgumentError("periodic distances require unit-cell data"))
    values = Vector{Float64}(undef, size(trajectory.coordinates, 3))
    for frame in axes(trajectory.coordinates, 3)
        if periodic
            cell = @view trajectory.unit_cells[:, :, frame]
            values[frame] = minimum(
                periodic_distance(
                    @view(trajectory.coordinates[:, i, frame]),
                    @view(trajectory.coordinates[:, j, frame]),
                    cell,
                ) for i in a for j in b
            )
        else
            tree = KDTree(trajectory.coordinates[:, b, frame])
            values[frame] = minimum(
                last(knn(tree, trajectory.coordinates[:, i, frame], 1))[1] for i in a
            )
        end
    end
    return values
end

function group_centroid_displacement(
    trajectory::TrajectoryData,
    group_indices;
    align_indices = axes(trajectory.coordinates, 2),
    reference_frame::Int = 1,
)
    checkbounds(trajectory.coordinates, :, :, reference_frame)
    group = _trajectory_indices(trajectory, group_indices, "group_indices")
    alignment_atoms =
        _trajectory_indices(trajectory, align_indices, "align_indices"; minimum_count = 3)
    reference = trajectory.coordinates[:, :, reference_frame]
    reference_centroid = centroid(reference[:, group])
    [
        norm(
            centroid(
                _fit_frame(reference, trajectory.coordinates[:, :, frame], alignment_atoms)[
                    :,
                    group,
                ],
            ) - reference_centroid,
        ) for frame in axes(trajectory.coordinates, 3)
    ]
end

function residue_rmsf(
    structure::StructureData,
    trajectory::TrajectoryData;
    reference_frame::Int = 1,
    align_indices = axes(trajectory.coordinates, 2),
    measure_indices = axes(trajectory.coordinates, 2),
)
    measured_atoms = _trajectory_indices(trajectory, measure_indices, "measure_indices")
    topology_indices = trajectory.atom_indices[measured_atoms]
    all(i -> 1 <= i <= length(structure.atoms), topology_indices) ||
        throw(ArgumentError("trajectory atom indices do not map to the supplied structure"))
    atom_values = trajectory_rmsf(
        trajectory;
        reference_frame = reference_frame,
        align_indices = align_indices,
        measure_indices = measured_atoms,
    )
    grouped_values = Dict{ResidueKey,Vector{Float64}}()
    molecule_types = Dict{ResidueKey,Symbol}()
    type_by_key =
        Dict(residue.key => residue.molecule_type for residue in residue_records(structure))
    for (topology_index, value) in zip(topology_indices, atom_values)
        key = residue_key(structure.atoms[topology_index])
        push!(get!(grouped_values, key, Float64[]), value)
        molecule_types[key] = type_by_key[key]
    end
    rows = DataFrame(
        model_id = Int[],
        chain_id = String[],
        residue_number = Int[],
        insertion_code = Char[],
        residue_name = String[],
        molecule_type = Symbol[],
        atom_count = Int[],
        rmsf_angstrom = Float64[],
    )
    for key in sort!(collect(keys(grouped_values)))
        values = grouped_values[key]
        push!(
            rows,
            (
                key.model_id,
                key.chain_id,
                key.number,
                key.insertion_code,
                key.name,
                molecule_types[key],
                length(values),
                sqrt(mean(abs2, values)),
            ),
        )
    end
    return rows
end

function trajectory_pca(
    trajectory::TrajectoryData;
    atom_indices = axes(trajectory.coordinates, 2),
    align::Bool = true,
    reference_frame::Int = 1,
)
    frames = size(trajectory.coordinates, 3)
    selected_atoms = _trajectory_indices(trajectory, atom_indices, "atom_indices")
    align &&
        length(selected_atoms) < 3 &&
        throw(ArgumentError("aligned PCA requires at least 3 atom indices"))
    matrix = Matrix{Float64}(undef, 3 * length(selected_atoms), frames)
    reference = trajectory.coordinates[:, :, reference_frame]
    for frame = 1:frames
        xyz =
            align ?
            _fit_frame(reference, trajectory.coordinates[:, :, frame], selected_atoms) :
            trajectory.coordinates[:, :, frame]
        matrix[:, frame] = vec(xyz[:, selected_atoms])
    end
    means = mean(matrix; dims = 2)
    centered = matrix .- means
    decomposition = svd(centered / sqrt(max(frames - 1, 1)); full = false)
    scores = decomposition.U' * centered
    variance = decomposition.S .^ 2
    total_variance = sum(variance)
    (
        scores = scores,
        components = decomposition.U,
        variance = variance,
        explained_fraction = iszero(total_variance) ? zeros(length(variance)) :
                             variance ./ total_variance,
        mean_coordinates = vec(means),
    )
end

"""Deterministic Lloyd clustering of columns in a feature matrix."""
function cluster_frames(features::AbstractMatrix, k::Integer; max_iterations::Integer = 100)
    size(features, 1) > 0 || throw(ArgumentError("features must contain at least one row"))
    all(isfinite, features) || throw(ArgumentError("features must be finite"))
    max_iterations > 0 || throw(ArgumentError("max_iterations must be positive"))
    n = size(features, 2)
    1 <= k <= n || throw(ArgumentError("k must be between 1 and the number of frames"))
    seeds = unique(round.(Int, range(1, n; length = k)))
    centers = copy(features[:, seeds])
    labels = ones(Int, n)
    for _ = 1:max_iterations
        new_labels =
            [argmin([sum(abs2, features[:, i] - centers[:, j]) for j = 1:k]) for i = 1:n]
        new_labels == labels && break
        labels = new_labels
        for j = 1:k
            members = findall(==(j), labels)
            !isempty(members) && (centers[:, j] = vec(mean(features[:, members]; dims = 2)))
        end
    end
    (labels = labels, centers = centers)
end

function representative_frames(features::AbstractMatrix, clustering)
    [
        argmin([
            sum(abs2, features[:, i] - clustering.centers[:, j]) for i in axes(features, 2)
        ]) for j in axes(clustering.centers, 2)
    ]
end
