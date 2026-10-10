using RibosomeAnalysis
using DataFrames

function parse_cli(arguments)
    options = Dict{String,String}()
    i = 1
    while i <= length(arguments)
        startswith(arguments[i], "--") ||
            throw(ArgumentError("unexpected argument: $(arguments[i])"))
        key = arguments[i][3:end]
        i == length(arguments) && throw(ArgumentError("missing value for --$key"))
        options[key] = arguments[i+1]
        i += 2
    end
    options
end

function require_option(options, key)
    haskey(options, key) || throw(ArgumentError("required option --$key was not supplied"))
    options[key]
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
