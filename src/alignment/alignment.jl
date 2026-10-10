const RESIDUE_SYMBOL = Dict(
    "ALA" => 'A',
    "ARG" => 'R',
    "ASN" => 'N',
    "ASP" => 'D',
    "CYS" => 'C',
    "GLN" => 'Q',
    "GLU" => 'E',
    "GLY" => 'G',
    "HIS" => 'H',
    "ILE" => 'I',
    "LEU" => 'L',
    "LYS" => 'K',
    "MET" => 'M',
    "PHE" => 'F',
    "PRO" => 'P',
    "SER" => 'S',
    "THR" => 'T',
    "TRP" => 'W',
    "TYR" => 'Y',
    "VAL" => 'V',
    "A" => 'A',
    "C" => 'C',
    "G" => 'G',
    "U" => 'U',
    "I" => 'I',
)

function kabsch(reference::AbstractMatrix, mobile::AbstractMatrix)
    size(reference) == size(mobile) ||
        throw(DimensionMismatch("coordinate matrices must have equal shape"))
    size(reference, 1) == 3 || throw(DimensionMismatch("coordinates must have shape 3 × N"))
    size(reference, 2) >= 3 ||
        throw(ArgumentError("at least three corresponding atoms are required"))
    reference_center = vec(mean(reference; dims = 2))
    mobile_center = vec(mean(mobile; dims = 2))
    x = mobile .- mobile_center
    y = reference .- reference_center
    decomposition = svd(x * y')
    correction = Diagonal([1.0, 1.0, sign(det(decomposition.V * decomposition.U'))])
    rotation = decomposition.V * correction * decomposition.U'
    translation = reference_center - rotation * mobile_center
    fitted = rotation * mobile .+ translation
    value = sqrt(sum(abs2, fitted - reference) / size(reference, 2))
    rotation, translation, value
end

apply_transform(
    xyz::AbstractMatrix,
    rotation::AbstractMatrix,
    translation::AbstractVector,
) = rotation * xyz .+ translation

function rmsd(reference::AbstractMatrix, target::AbstractMatrix; superimpose::Bool = false)
    size(reference) == size(target) ||
        throw(DimensionMismatch("coordinate matrices must have equal shape"))
    if superimpose
        _, _, value = kabsch(reference, target)
        value
    else
        sqrt(sum(abs2, reference - target) / size(reference, 2))
    end
end

function _needleman_wunsch(
    a::Vector{Char},
    b::Vector{Char};
    match = 2,
    mismatch = -1,
    gap = -2,
)
    n, m = length(a), length(b)
    score = zeros(Int, n + 1, m + 1)
    trace = zeros(Int8, n + 1, m + 1) # 1 diagonal, 2 up, 3 left
    for i = 1:n
        score[i+1, 1] = i * gap
        trace[i+1, 1] = 2
    end
    for j = 1:m
        score[1, j+1] = j * gap
        trace[1, j+1] = 3
    end
    for i = 1:n, j = 1:m
        diagonal = score[i, j] + (a[i] == b[j] ? match : mismatch)
        up = score[i, j+1] + gap
        left = score[i+1, j] + gap
        best = max(diagonal, up, left)
        score[i+1, j+1] = best
        trace[i+1, j+1] = best == diagonal ? 1 : (best == up ? 2 : 3)
    end
    pairs = Tuple{Int,Int}[]
    i, j = n, m
    while i > 0 || j > 0
        direction = trace[i+1, j+1]
        if direction == 1
            push!(pairs, (i, j))
            i -= 1
            j -= 1
        elseif direction == 2
            i -= 1
        else
            j -= 1
        end
    end
    reverse!(pairs)
end

function _polymer_residues(structure::StructureData, chain_id::String)
    [
        r for r in residue_records(structure) if
        r.key.chain_id == chain_id && r.molecule_type in (:rna, :protein)
    ]
end

"""
Sequence-aware residue matching. Chain equivalence is explicit: identical chain
IDs are used only when no `chain_map` is supplied.
"""
function match_residues(
    reference::StructureData,
    target::StructureData;
    chain_map::Union{Nothing,AbstractDict} = nothing,
)
    ref_chains = sort!(
        unique(
            r.key.chain_id for
            r in residue_records(reference) if r.molecule_type in (:rna, :protein)
        ),
    )
    mapping =
        isnothing(chain_map) ? Dict(c => c for c in ref_chains) :
        Dict(string(k) => string(v) for (k, v) in chain_map)
    matched = Tuple{ResidueRecord,ResidueRecord}[]
    used_reference = Set{ResidueKey}()
    used_target = Set{ResidueKey}()
    for ref_chain in ref_chains
        haskey(mapping, ref_chain) || continue
        target_chain = mapping[ref_chain]
        a = _polymer_residues(reference, ref_chain)
        b = _polymer_residues(target, target_chain)
        isempty(a) && continue
        isempty(b) && continue
        seq_a = [get(RESIDUE_SYMBOL, uppercase(r.key.name), 'X') for r in a]
        seq_b = [get(RESIDUE_SYMBOL, uppercase(r.key.name), 'X') for r in b]
        for (i, j) in _needleman_wunsch(seq_a, seq_b)
            seq_a[i] == seq_b[j] || continue
            push!(matched, (a[i], b[j]))
            push!(used_reference, a[i].key)
            push!(used_target, b[j].key)
        end
    end
    all_ref = Set(
        r.key for r in residue_records(reference) if r.molecule_type in (:rna, :protein)
    )
    all_target =
        Set(r.key for r in residue_records(target) if r.molecule_type in (:rna, :protein))
    matched,
    sort!(collect(setdiff(all_ref, used_reference))),
    sort!(collect(setdiff(all_target, used_target)))
end

function _atom_by_name(structure::StructureData, residue::ResidueRecord)
    result = Dict{String,Int}()
    for i in residue.atom_indices
        name = structure.atoms[i].name
        if !haskey(result, name) ||
           structure.atoms[i].occupancy > structure.atoms[result[name]].occupancy
            result[name] = i
        end
    end
    result
end

function align_structures(
    reference::StructureData,
    target::StructureData;
    chain_map::Union{Nothing,AbstractDict} = nothing,
    atom_names::Union{Nothing,AbstractVector} = nothing,
    residue_pairs = nothing,
)
    pairs, unmatched_ref, unmatched_target =
        isnothing(residue_pairs) ?
        match_residues(reference, target; chain_map = chain_map) :
        (collect(residue_pairs), ResidueKey[], ResidueKey[])
    ref_indices, target_indices = Int[], Int[]
    for (ref_residue, target_residue) in pairs
        ref_atoms = _atom_by_name(reference, ref_residue)
        target_atoms = _atom_by_name(target, target_residue)
        defaults =
            ref_residue.molecule_type === :rna ? ["P", "C4'", "C3'", "O3'"] :
            ["N", "CA", "C"]
        names_to_use = isnothing(atom_names) ? defaults : string.(atom_names)
        for name in names_to_use
            haskey(ref_atoms, name) && haskey(target_atoms, name) || continue
            push!(ref_indices, ref_atoms[name])
            push!(target_indices, target_atoms[name])
        end
    end
    length(ref_indices) >= 3 || throw(
        ArgumentError(
            "fewer than three corresponding atoms were found; provide a chain or residue correspondence map",
        ),
    )
    rotation, translation, value =
        kabsch(coordinates(reference, ref_indices), coordinates(target, target_indices))
    reference_polymer_count =
        count(r -> r.molecule_type in (:rna, :protein), residue_records(reference))
    target_polymer_count =
        count(r -> r.molecule_type in (:rna, :protein), residue_records(target))
    denominator = max(reference_polymer_count, target_polymer_count)
    coverage = denominator == 0 ? 0.0 : length(pairs) / denominator
    AlignmentResult(
        rotation,
        translation,
        value,
        length(ref_indices),
        coverage,
        collect(zip(ref_indices, target_indices)),
        unmatched_ref,
        unmatched_target,
    )
end
