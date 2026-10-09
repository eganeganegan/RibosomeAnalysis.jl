"""Stable atom representation. Coordinates are stored in Angstrom."""
struct AtomRecord
    serial::Int
    name::String
    residue_name::String
    chain_id::String
    author_chain_id::String
    label_chain_id::String
    residue_number::Int
    author_residue_number::Union{Int,Missing}
    label_residue_number::Union{Int,Missing}
    insertion_code::Char
    alternate_location::Char
    model_id::Int
    occupancy::Float64
    bfactor::Float64
    element::String
    hetero::Bool
    position::SVector{3,Float64}
end

function AtomRecord(serial, name, residue_name, chain_id, author_chain_id, label_chain_id,
                    residue_number, author_residue_number, label_residue_number,
                    insertion_code, alternate_location, model_id, occupancy, bfactor,
                    element, hetero, position::AbstractVector)
    length(position) == 3 || throw(DimensionMismatch("atom position must have length 3"))
    AtomRecord(Int(serial), string(name), string(residue_name), string(chain_id),
               string(author_chain_id), string(label_chain_id), Int(residue_number),
               author_residue_number, label_residue_number, Char(insertion_code),
               Char(alternate_location), Int(model_id), Float64(occupancy), Float64(bfactor),
               string(element), Bool(hetero), SVector{3,Float64}(position))
end

struct ResidueKey
    model_id::Int
    chain_id::String
    number::Int
    insertion_code::Char
    name::String
end

Base.isless(a::ResidueKey, b::ResidueKey) =
    (a.model_id, a.chain_id, a.number, a.insertion_code, a.name) <
    (b.model_id, b.chain_id, b.number, b.insertion_code, b.name)

struct ResidueRecord
    key::ResidueKey
    atom_indices::Vector{Int}
    molecule_type::Symbol
end

struct ChainRecord
    model_id::Int
    id::String
    residue_indices::Vector{Int}
end

struct StructureMetadata
    identifier::String
    source_path::String
    format::Symbol
    parser::String
    annotations::Dict{String,Any}
end

struct StructureData
    atoms::Vector{AtomRecord}
    metadata::StructureMetadata
end

struct Ligand
    identifier::String
    residue_name::String
    chain_id::String
    residue_number::Int
    insertion_code::Char
    model_id::Int
    atom_indices::Vector{Int}
    antibiotic_class::Union{String,Missing}
    identification::Symbol
    notes::Vector{String}
end

struct BindingSite
    ligand::Ligand
    cutoff_angstrom::Float64
    residue_keys::Vector{ResidueKey}
    residue_atom_indices::Vector{Vector{Int}}
    contacting_atom_indices::Vector{Int}
    minimum_distances::Dict{ResidueKey,Float64}
end

struct ContactRecord
    structure_atom_index::Int
    ligand_atom_index::Int
    residue::ResidueKey
    distance_angstrom::Float64
end

struct InteractionRecord
    interaction_type::Symbol
    participant_indices::Vector{Int}
    participant_residues::Vector{ResidueKey}
    distance_angstrom::Union{Float64,Missing}
    angle_degrees::Union{Float64,Missing}
    confidence::Symbol
    limitation::String
end

struct AlignmentResult
    rotation::Matrix{Float64}
    translation::Vector{Float64}
    rmsd_angstrom::Float64
    aligned_atoms::Int
    correspondence_coverage::Float64
    matched_indices::Vector{Tuple{Int,Int}}
    unmatched_reference::Vector{ResidueKey}
    unmatched_target::Vector{ResidueKey}
end

"""In-memory trajectory. `coordinates` has shape 3 × atoms × frames in Angstrom."""
struct TrajectoryData
    coordinates::Array{Float64,3}
    times_ps::Vector{Float64}
    atom_indices::Vector{Int}
    source::String
    unit_conversion_to_angstrom::Float64
end

struct StateMetadata
    label::String
    organism::Union{String,Missing}
    subunit::Union{String,Missing}
    experimental_method::Union{String,Missing}
    resolution_angstrom::Union{Float64,Missing}
    antibiotic::Union{String,Missing}
    translation_state::Union{String,Missing}
    annotations::Dict{String,Any}
end

struct ComparisonResult
    reference_label::String
    target_label::String
    alignment::AlignmentResult
    metrics::Dict{String,Any}
    changed_contacts::DataFrame
    limitations::Vector{String}
end

residue_key(a::AtomRecord) = ResidueKey(a.model_id, a.chain_id, a.residue_number,
                                        a.insertion_code, a.residue_name)
