"""
Curated candidate residue names. These mappings support identification only;
they do not imply equivalent binding modes within an antibiotic class.
"""
const DEFAULT_ANTIBIOTIC_METADATA = Dict(
    "PAR" => Dict("name" => "paromomycin", "class" => "aminoglycoside"),
    "NEO" => Dict("name" => "neomycin", "class" => "aminoglycoside"),
    "KAN" => Dict("name" => "kanamycin", "class" => "aminoglycoside"),
    "STR" => Dict("name" => "streptomycin", "class" => "aminoglycoside"),
    "ERY" => Dict("name" => "erythromycin", "class" => "macrolide"),
    "AZM" => Dict("name" => "azithromycin", "class" => "macrolide"),
    "TET" => Dict("name" => "tetracycline", "class" => "tetracycline"),
    "DOX" => Dict("name" => "doxycycline", "class" => "tetracycline"),
    "LZD" => Dict("name" => "linezolid", "class" => "oxazolidinone"),
    "CLI" => Dict("name" => "clindamycin", "class" => "lincosamide"),
    "CAM" => Dict("name" => "chloramphenicol", "class" => "chloramphenicol"),
    "RET" => Dict("name" => "retapamulin", "class" => "pleuromutilin"),
)

antibiotic_metadata() = deepcopy(DEFAULT_ANTIBIOTIC_METADATA)

"""Load a CSV ligand mapping with required columns `residue_name`, `name`, `class`."""
function load_ligand_mapping(path::AbstractString)
    isfile(path) || throw(ArgumentError("ligand mapping file does not exist: $path"))
    table = CSV.read(path, DataFrame)
    required = [:residue_name, :name, :class]
    available = Set(Symbol.(names(table)))
    all(in(available), required) ||
        throw(ArgumentError("mapping must contain residue_name,name,class columns"))
    Dict(
        uppercase(string(row.residue_name)) =>
            Dict("name" => string(row.name), "class" => string(row.class)) for
        row in eachrow(table)
    )
end

function _make_ligand(residue::ResidueRecord, mapping, method::Symbol)
    info = get(mapping, uppercase(residue.key.name), nothing)
    class = isnothing(info) ? missing : string(info["class"])
    identifier = string(
        residue.key.name,
        ":",
        residue.key.chain_id,
        ":",
        residue.key.number,
        residue.key.insertion_code == ' ' ? "" : residue.key.insertion_code,
        ":model",
        residue.key.model_id,
    )
    notes =
        isnothing(info) ? ["explicitly selected; no antibiotic class asserted"] :
        [
            "metadata mapping identifies a candidate; identity should be checked against structure metadata",
        ]
    Ligand(
        identifier,
        residue.key.name,
        residue.key.chain_id,
        residue.key.number,
        residue.key.insertion_code,
        residue.key.model_id,
        copy(residue.atom_indices),
        class,
        method,
        notes,
    )
end

"""
    identify_ligands(structure; residue_name=nothing, chain_id=nothing,
                     residue_number=nothing, mapping=antibiotic_metadata())

With `residue_name`, performs explicit selection. Without it, returns only
hetero residues present in the supplied metadata mapping; arbitrary nonstandard
residues are never labeled as antibiotics.
"""
function identify_ligands(
    structure::StructureData;
    residue_name::Union{String,Nothing} = nothing,
    chain_id::Union{String,Nothing} = nothing,
    residue_number::Union{Int,Nothing} = nothing,
    model_id::Union{Int,Nothing} = nothing,
    mapping = antibiotic_metadata(),
)
    result = Ligand[]
    for residue in residue_records(structure)
        atom = structure.atoms[first(residue.atom_indices)]
        atom.hetero || continue
        explicit = !isnothing(residue_name)
        explicit && uppercase(residue.key.name) != uppercase(residue_name) && continue
        !isnothing(chain_id) && residue.key.chain_id != chain_id && continue
        !isnothing(residue_number) && residue.key.number != residue_number && continue
        !isnothing(model_id) && residue.key.model_id != model_id && continue
        if !explicit && !haskey(mapping, uppercase(residue.key.name))
            continue
        end
        push!(
            result,
            _make_ligand(residue, mapping, explicit ? :explicit : :metadata_assisted),
        )
    end
    result
end

function select_ligand(structure::StructureData; kwargs...)
    ligands = identify_ligands(structure; kwargs...)
    isempty(ligands) && throw(ArgumentError("no ligand matched the supplied identifiers"))
    length(ligands) == 1 || throw(
        ArgumentError(
            "ligand selection is ambiguous ($(length(ligands)) matches); specify chain, residue number, and model",
        ),
    )
    only(ligands)
end
