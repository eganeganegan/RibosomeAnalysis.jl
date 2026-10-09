function atomic_contacts(structure::StructureData, ligand::Ligand; cutoff::Real=4.0,
                         heavy_atoms_only::Bool=true)
    cutoff > 0 || throw(ArgumentError("cutoff must be positive"))
    ligand_indices = heavy_atoms_only ?
        [i for i in ligand.atom_indices if uppercase(structure.atoms[i].element) != "H"] : ligand.atom_indices
    isempty(ligand_indices) && return ContactRecord[]
    tree = KDTree(coordinates(structure, ligand_indices))
    ligand_set = Set(ligand.atom_indices)
    contacts = ContactRecord[]
    for (i, atom) in enumerate(structure.atoms)
        i in ligand_set && continue
        heavy_atoms_only && uppercase(atom.element) == "H" && continue
        for local_j in inrange(tree, atom.position, Float64(cutoff))
            j = ligand_indices[local_j]
            push!(contacts, ContactRecord(i, j, residue_key(atom),
                                          norm(atom.position - structure.atoms[j].position)))
        end
    end
    sort!(contacts; by=c -> (c.distance_angstrom, c.structure_atom_index, c.ligand_atom_index))
end

_angle_degrees(a, b, c) = begin
    u, v = a - b, c - b
    denom = norm(u) * norm(v)
    denom == 0 && return NaN
    acosd(clamp(dot(u, v) / denom, -1.0, 1.0))
end

"""
Find hydrogen bonds only from explicit covalent bonds and explicit hydrogens.
`bonds` contains atom-index pairs. No bond or charge assignment is inferred.
"""
function explicit_hydrogen_bonds(structure::StructureData, ligand::Ligand,
                                 bonds::AbstractVector{<:Tuple};
                                 donor_acceptor_cutoff::Real=3.5,
                                 minimum_dha_angle::Real=120.0)
    adjacency = Dict{Int,Vector{Int}}()
    for (i, j) in bonds
        push!(get!(adjacency, i, Int[]), j)
        push!(get!(adjacency, j, Int[]), i)
    end
    ligand_set = Set(ligand.atom_indices)
    acceptors = Set(i for (i, a) in enumerate(structure.atoms) if uppercase(a.element) in ("N", "O", "S"))
    interactions = InteractionRecord[]
    for (donor, neighbors) in adjacency
        uppercase(structure.atoms[donor].element) in ("N", "O", "S") || continue
        for h in neighbors
            uppercase(structure.atoms[h].element) == "H" || continue
            for acceptor in acceptors
                acceptor == donor && continue
                ((donor in ligand_set) ⊻ (acceptor in ligand_set)) || continue
                d = norm(structure.atoms[donor].position - structure.atoms[acceptor].position)
                d <= donor_acceptor_cutoff || continue
                angle = _angle_degrees(structure.atoms[donor].position,
                                       structure.atoms[h].position,
                                       structure.atoms[acceptor].position)
                angle >= minimum_dha_angle || continue
                push!(interactions, InteractionRecord(:hydrogen_bond, [donor, h, acceptor],
                    unique(residue_key(structure.atoms[i]) for i in (donor, acceptor)),
                    d, angle, :geometry_with_explicit_hydrogen,
                    "geometric assignment; protonation and energetics are not established"))
            end
        end
    end
    interactions
end

"""Report proximity between explicitly supplied positive and negative atom groups."""
function ionic_proximities(structure::StructureData, positive_indices, negative_indices;
                           cutoff::Real=4.0)
    interactions = InteractionRecord[]
    for i in positive_indices, j in negative_indices
        d = norm(structure.atoms[i].position - structure.atoms[j].position)
        d <= cutoff || continue
        push!(interactions, InteractionRecord(:ionic_proximity, [i, j],
            [residue_key(structure.atoms[i]), residue_key(structure.atoms[j])],
            d, missing, :candidate,
            "charge groups were user-supplied; distance alone does not establish an energetic interaction"))
    end
    interactions
end

function _ring_geometry(structure::StructureData, indices)
    length(indices) >= 3 || throw(ArgumentError("an aromatic ring requires at least three atom indices"))
    xyz = coordinates(structure, collect(indices))
    center = vec(mean(xyz; dims=2))
    centered = xyz .- center
    decomposition = svd(centered')
    normal = decomposition.V[:, end]
    center, normal
end

"""Analyze explicitly defined aromatic ring atom groups; aromaticity is not inferred."""
function aromatic_interactions(structure::StructureData, rings_a, rings_b;
                               centroid_cutoff::Real=6.0)
    out = InteractionRecord[]
    for ra in rings_a, rb in rings_b
        ca, na = _ring_geometry(structure, ra)
        cb, nb = _ring_geometry(structure, rb)
        d = norm(ca - cb)
        d <= centroid_cutoff || continue
        angle = acosd(clamp(abs(dot(na, nb)), 0.0, 1.0))
        participants = vcat(collect(ra), collect(rb))
        push!(out, InteractionRecord(:aromatic_geometry, participants,
            unique(residue_key(structure.atoms[i]) for i in participants), d, angle,
            :candidate, "ring groups were user-supplied; geometry does not establish interaction energy"))
    end
    out
end

function metal_mediated_interactions(structure::StructureData, ligand::Ligand;
                                     metal_cutoff::Real=3.0,
                                     ribosome_cutoff::Real=3.0)
    metals = findall(a -> uppercase(a.element) in ("MG", "MN", "ZN", "CA", "FE", "CO", "CU", "NI"), structure.atoms)
    ligand_set = Set(ligand.atom_indices)
    donors = findall(a -> uppercase(a.element) in ("N", "O", "S"), structure.atoms)
    out = InteractionRecord[]
    for m in metals
        near_ligand = [i for i in ligand.atom_indices if norm(structure.atoms[m].position - structure.atoms[i].position) <= metal_cutoff]
        isempty(near_ligand) && continue
        near_site = [i for i in donors if !(i in ligand_set) && norm(structure.atoms[m].position - structure.atoms[i].position) <= ribosome_cutoff]
        for li in near_ligand, si in near_site
            d = norm(structure.atoms[m].position - structure.atoms[li].position)
            angle = _angle_degrees(structure.atoms[li].position, structure.atoms[m].position, structure.atoms[si].position)
            push!(out, InteractionRecord(:metal_mediated_candidate, [li, m, si],
                unique(residue_key(structure.atoms[i]) for i in (li, m, si)), d, angle,
                :candidate, "distance-based candidate only; coordination chemistry and bond order are not assigned"))
        end
    end
    out
end
