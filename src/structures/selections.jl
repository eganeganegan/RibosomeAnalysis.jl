const STANDARD_RNA = Set(["A", "C", "G", "U", "I", "DA", "DC", "DG", "DT"])
const STANDARD_PROTEIN = Set([
    "ALA",
    "ARG",
    "ASN",
    "ASP",
    "CYS",
    "GLN",
    "GLU",
    "GLY",
    "HIS",
    "ILE",
    "LEU",
    "LYS",
    "MET",
    "PHE",
    "PRO",
    "SER",
    "THR",
    "TRP",
    "TYR",
    "VAL",
    "ASX",
    "GLX",
    "SEC",
    "PYL",
])
const WATER_NAMES = Set(["HOH", "WAT", "H2O", "DOD"])
const ION_NAMES = Set(["MG", "MN", "ZN", "CA", "NA", "K", "CL", "FE", "CO", "CU", "NI"])

coordinates(atoms::AbstractVector{AtomRecord}) =
    isempty(atoms) ? zeros(3, 0) : reduce(hcat, (a.position for a in atoms))
coordinates(
    structure::StructureData,
    indices::AbstractVector{<:Integer} = eachindex(structure.atoms),
) = coordinates(structure.atoms[indices])

select_atoms(structure::StructureData, predicate::Function) =
    findall(predicate, structure.atoms)

function select_chains(
    structure::StructureData,
    chain_ids;
    model_id::Union{Int,Nothing} = nothing,
)
    ids = Set(string.(chain_ids isa AbstractString ? [chain_ids] : chain_ids))
    findall(
        a -> a.chain_id in ids && (isnothing(model_id) || a.model_id == model_id),
        structure.atoms,
    )
end

function select_residues(
    structure::StructureData;
    chain_id::Union{String,Nothing} = nothing,
    numbers = nothing,
    names = nothing,
    model_id::Union{Int,Nothing} = nothing,
)
    number_set = isnothing(numbers) ? nothing : Set(Int.(numbers))
    name_set = isnothing(names) ? nothing : Set(uppercase.(string.(names)))
    findall(structure.atoms) do a
        (isnothing(chain_id) || a.chain_id == chain_id) &&
            (isnothing(number_set) || a.residue_number in number_set) &&
            (isnothing(name_set) || uppercase(a.residue_name) in name_set) &&
            (isnothing(model_id) || a.model_id == model_id)
    end
end

function _molecule_type(atoms::Vector{AtomRecord}, indices::Vector{Int})
    name = uppercase(atoms[first(indices)].residue_name)
    name in WATER_NAMES && return :water
    name in ION_NAMES && return :ion
    name in STANDARD_RNA && return :rna
    name in STANDARD_PROTEIN && return :protein
    atom_names = Set(atoms[i].name for i in indices)
    has_sugar =
        any(in(atom_names), ("C1'", "C2'", "C3'", "C4'", "C1*", "C2*", "C3*", "C4*"))
    has_sugar &&
        any(n -> startswith(n, "P") || startswith(n, "O"), atom_names) &&
        return :rna
    all(in(atom_names), ("N", "CA", "C")) && return :protein
    atoms[first(indices)].hetero ? :ligand : :other
end

function residue_records(structure::StructureData)
    groups = Dict{ResidueKey,Vector{Int}}()
    for (i, atom) in enumerate(structure.atoms)
        push!(get!(groups, residue_key(atom), Int[]), i)
    end
    keys_sorted = sort!(collect(keys(groups)))
    [
        ResidueRecord(key, groups[key], _molecule_type(structure.atoms, groups[key])) for
        key in keys_sorted
    ]
end

function chain_records(structure::StructureData)
    residues = residue_records(structure)
    groups = Dict{Tuple{Int,String},Vector{Int}}()
    for (i, residue) in enumerate(residues)
        push!(get!(groups, (residue.key.model_id, residue.key.chain_id), Int[]), i)
    end
    [ChainRecord(k[1], k[2], v) for (k, v) in sort!(collect(groups); by = first)]
end

function _indices_of_type(structure::StructureData, kind::Symbol)
    reduce(
        vcat,
        (r.atom_indices for r in residue_records(structure) if r.molecule_type === kind);
        init = Int[],
    )
end

rna_atoms(structure::StructureData) = _indices_of_type(structure, :rna)
protein_atoms(structure::StructureData) = _indices_of_type(structure, :protein)
ligand_atoms(structure::StructureData) = _indices_of_type(structure, :ligand)
