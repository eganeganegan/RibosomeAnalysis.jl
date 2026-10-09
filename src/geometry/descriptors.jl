centroid(xyz::AbstractMatrix) = vec(mean(xyz; dims=2))

function radius_of_gyration(xyz::AbstractMatrix; masses::Union{Nothing,AbstractVector}=nothing)
    size(xyz, 1) == 3 || throw(DimensionMismatch("coordinates must have shape 3 × N"))
    n = size(xyz, 2)
    n > 0 || throw(ArgumentError("coordinates are empty"))
    weights = isnothing(masses) ? fill(1 / n, n) : Float64.(masses) ./ sum(masses)
    length(weights) == n || throw(DimensionMismatch("one mass is required per coordinate"))
    center = xyz * weights
    sqrt(sum(weights[i] * sum(abs2, view(xyz, :, i) .- center) for i in 1:n))
end

function spatial_extent(xyz::AbstractMatrix)
    size(xyz, 2) > 0 || throw(ArgumentError("coordinates are empty"))
    vec(maximum(xyz; dims=2) - minimum(xyz; dims=2))
end

function pairwise_distances(xyz::AbstractMatrix)
    n = size(xyz, 2)
    result = Vector{Float64}(undef, n * (n - 1) ÷ 2)
    k = 1
    for j in 2:n, i in 1:j-1
        result[k] = norm(view(xyz, :, i) - view(xyz, :, j))
        k += 1
    end
    result
end

function local_contact_density(xyz::AbstractMatrix; cutoff::Real=5.0)
    size(xyz, 2) == 0 && return 0.0
    tree = KDTree(xyz)
    contacts = sum(length(inrange(tree, view(xyz, :, i), Float64(cutoff))) - 1 for i in axes(xyz, 2)) ÷ 2
    contacts / size(xyz, 2)
end

function geometry_descriptors(structure::StructureData, site::BindingSite)
    isempty(site.residue_atom_indices) && throw(ArgumentError("binding site is empty at the selected cutoff"))
    indices = unique(vcat(site.residue_atom_indices...))
    xyz = coordinates(structure, indices)
    distances = pairwise_distances(xyz)
    Dict{String,Any}(
        "centroid_angstrom" => centroid(xyz),
        "radius_of_gyration_angstrom" => radius_of_gyration(xyz),
        "extent_angstrom" => spatial_extent(xyz),
        "pairwise_distance_mean_angstrom" => isempty(distances) ? 0.0 : mean(distances),
        "pairwise_distance_median_angstrom" => isempty(distances) ? 0.0 : median(distances),
        "local_contact_density_per_atom" => local_contact_density(xyz; cutoff=5.0),
        "atom_count" => length(indices),
        "limitations" => "geometric descriptors are not pocket volumes or binding-affinity estimates",
    )
end
