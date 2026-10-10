function identify_binding_site(
    structure::StructureData,
    ligand::Ligand;
    cutoff::Real = 5.0,
    include_hydrogen::Bool = false,
    include_ligand::Bool = false,
)
    cutoff > 0 || throw(ArgumentError("cutoff must be positive"))
    ligand_set = Set(ligand.atom_indices)
    ligand_coords = coordinates(structure, ligand.atom_indices)
    tree = KDTree(ligand_coords)
    residues = residue_records(structure)
    keys = ResidueKey[]
    residue_atoms = Vector{Int}[]
    contacting_atoms = Int[]
    minimum_distances = Dict{ResidueKey,Float64}()
    for residue in residues
        (!include_ligand && any(in(ligand_set), residue.atom_indices)) && continue
        candidates =
            include_hydrogen ? residue.atom_indices :
            [
                i for
                i in residue.atom_indices if uppercase(structure.atoms[i].element) != "H"
            ]
        isempty(candidates) && continue
        local_contacts = Int[]
        min_distance = Inf
        for i in candidates
            idxs = inrange(tree, structure.atoms[i].position, Float64(cutoff))
            if !isempty(idxs)
                push!(local_contacts, i)
                min_distance = min(
                    min_distance,
                    minimum(
                        norm(
                            structure.atoms[i].position -
                            structure.atoms[ligand.atom_indices[j]].position,
                        ) for j in idxs
                    ),
                )
            end
        end
        if !isempty(local_contacts)
            push!(keys, residue.key)
            push!(residue_atoms, copy(residue.atom_indices))
            append!(contacting_atoms, local_contacts)
            minimum_distances[residue.key] = min_distance
        end
    end
    BindingSite(
        ligand,
        Float64(cutoff),
        keys,
        residue_atoms,
        sort!(unique(contacting_atoms)),
        minimum_distances,
    )
end

function binding_site_summary(structure::StructureData, site::BindingSite)
    types = Dict(r.key => r.molecule_type for r in residue_records(structure))
    chains = countmap(k.chain_id for k in site.residue_keys)
    Dict{String,Any}(
        "ligand" => site.ligand.identifier,
        "cutoff_angstrom" => site.cutoff_angstrom,
        "contacting_residues" => length(site.residue_keys),
        "contacting_atoms" => length(site.contacting_atom_indices),
        "rna_residues" => count(k -> get(types, k, :other) === :rna, site.residue_keys),
        "protein_residues" =>
            count(k -> get(types, k, :other) === :protein, site.residue_keys),
        "chain_contributions" => Dict(string(k) => v for (k, v) in chains),
        "minimum_distance_angstrom" =>
            isempty(site.minimum_distances) ? nothing :
            minimum(values(site.minimum_distances)),
    )
end

"""Export a site as a PyMOL or ChimeraX residue-selection expression."""
function export_selection(site::BindingSite, path::AbstractString; format::Symbol = :pymol)
    isempty(site.residue_keys) &&
        throw(ArgumentError("cannot export an empty binding site"))
    expressions = if format === :pymol
        [
            "(chain $(k.chain_id) and resi $(k.number)$(k.insertion_code == ' ' ? "" : k.insertion_code))"
            for k in site.residue_keys
        ]
    elseif format === :chimerax
        [
            "/$(k.chain_id):$(k.number)$(k.insertion_code == ' ' ? "" : k.insertion_code)"
            for k in site.residue_keys
        ]
    else
        throw(ArgumentError("format must be :pymol or :chimerax"))
    end
    ispath(path) && throw(ArgumentError("selection output already exists: $path"))
    open(path, "w") do io
        if format === :pymol
            println(io, "select binding_site, ", join(expressions, " or "))
        else
            println(io, "select ", join(expressions, " "))
        end
    end
    path
end
