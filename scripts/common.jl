using RibosomeAnalysis
using DataFrames

function parse_cli(arguments; allowed_options = nothing)
    options = Dict{String,String}()
    i = 1
    while i <= length(arguments)
        startswith(arguments[i], "--") && length(arguments[i]) > 2 || throw(
            ArgumentError("unexpected argument: $(arguments[i]); expected --option value"),
        )
        key = arguments[i][3:end]
        i == length(arguments) && throw(ArgumentError("missing value for --$key"))
        haskey(options, key) &&
            throw(ArgumentError("option --$key was supplied more than once"))
        options[key] = arguments[i+1]
        i += 2
    end
    if !isnothing(allowed_options)
        allowed = Set(string.(allowed_options))
        unknown = sort!(collect(setdiff(keys(options), allowed)))
        isempty(unknown) || throw(
            ArgumentError(
                "unknown option(s): " * join(("--$key" for key in unknown), ", "),
            ),
        )
    end
    options
end

function require_option(options, key)
    haskey(options, key) || throw(ArgumentError("required option --$key was not supplied"))
    options[key]
end

function parse_option(options, key, ::Type{T}; default = nothing) where {T<:Real}
    haskey(options, key) || return default
    value = tryparse(T, options[key])
    isnothing(value) && throw(
        ArgumentError("option --$key must be a valid $(nameof(T)); got '$(options[key])'"),
    )
    value
end

function parse_chain_map(value::AbstractString)
    isempty(strip(value)) && throw(ArgumentError("option --chain-map cannot be empty"))
    mapping = Dict{String,String}()
    for entry in split(value, ',')
        parts = split(entry, ':'; limit = 2)
        length(parts) == 2 || throw(
            ArgumentError(
                "option --chain-map must use reference:target pairs separated by commas",
            ),
        )
        reference_chain, target_chain = strip.(parts)
        isempty(reference_chain) &&
            throw(ArgumentError("option --chain-map contains an empty reference chain"))
        isempty(target_chain) &&
            throw(ArgumentError("option --chain-map contains an empty target chain"))
        haskey(mapping, reference_chain) && throw(
            ArgumentError(
                "option --chain-map contains reference chain '$reference_chain' more than once",
            ),
        )
        mapping[reference_chain] = target_chain
    end
    return mapping
end

function run_cli(main_function, arguments = ARGS)
    try
        main_function(arguments)
        0
    catch error
        if error isa ArgumentError || error isa DimensionMismatch
            println(stderr, "Error: ", sprint(showerror, error))
            return 2
        end
        rethrow()
    end
end

function contact_dataframe(structure, contacts)
    DataFrame(
        chain = [c.residue.chain_id for c in contacts],
        residue_number = [c.residue.number for c in contacts],
        residue_name = [c.residue.name for c in contacts],
        structure_atom = [structure.atoms[c.structure_atom_index].name for c in contacts],
        ligand_atom = [structure.atoms[c.ligand_atom_index].name for c in contacts],
        distance_angstrom = [c.distance_angstrom for c in contacts],
    )
end
