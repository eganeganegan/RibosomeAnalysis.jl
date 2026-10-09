function _contact_keys(contacts)
    Set((c.residue, c.structure_atom_index, c.ligand_atom_index) for c in contacts)
end

function compare_binding_sites(reference::StructureData, target::StructureData,
                               reference_ligand::Ligand, target_ligand::Ligand;
                               cutoff::Real=5.0, chain_map=nothing)
    alignment = align_structures(reference, target; chain_map=chain_map)
    reference_site = identify_binding_site(reference, reference_ligand; cutoff=cutoff)
    target_site = identify_binding_site(target, target_ligand; cutoff=cutoff)
    reference_contacts = atomic_contacts(reference, reference_ligand; cutoff=cutoff)
    target_contacts = atomic_contacts(target, target_ligand; cutoff=cutoff)
    reference_site_keys = Set(reference_site.residue_keys)
    target_site_keys = Set(target_site.residue_keys)
    site_atom_pairs = [(i, j) for (i, j) in alignment.matched_indices
                       if residue_key(reference.atoms[i]) in reference_site_keys &&
                          residue_key(target.atoms[j]) in target_site_keys]
    global_site_rmsd = missing
    local_site_rmsd = missing
    if !isempty(site_atom_pairs)
        ref_site_xyz = coordinates(reference, first.(site_atom_pairs))
        target_site_xyz = coordinates(target, last.(site_atom_pairs))
        globally_fitted = apply_transform(target_site_xyz, alignment.rotation, alignment.translation)
        global_site_rmsd = rmsd(ref_site_xyz, globally_fitted)
        length(site_atom_pairs) >= 3 && (local_site_rmsd = rmsd(ref_site_xyz, target_site_xyz; superimpose=true))
    end
    target_ligand_fitted = apply_transform(coordinates(target, target_ligand.atom_indices),
                                           alignment.rotation, alignment.translation)
    ligand_centroid_displacement = norm(centroid(coordinates(reference, reference_ligand.atom_indices)) -
                                        centroid(target_ligand_fitted))
    ref_residues = Set((k.chain_id, k.number, k.insertion_code, k.name) for k in reference_site.residue_keys)
    target_residues = Set((k.chain_id, k.number, k.insertion_code, k.name) for k in target_site.residue_keys)
    changed = DataFrame(
        change=vcat(fill("reference_only", length(setdiff(ref_residues, target_residues))),
                    fill("target_only", length(setdiff(target_residues, ref_residues)))),
        residue=string.(vcat(collect(setdiff(ref_residues, target_residues)),
                             collect(setdiff(target_residues, ref_residues)))),
    )
    metrics = Dict{String,Any}(
        "alignment_rmsd_angstrom" => alignment.rmsd_angstrom,
        "aligned_atoms" => alignment.aligned_atoms,
        "correspondence_coverage" => alignment.correspondence_coverage,
        "reference_binding_site_residues" => length(reference_site.residue_keys),
        "target_binding_site_residues" => length(target_site.residue_keys),
        "reference_atomic_contacts" => length(reference_contacts),
        "target_atomic_contacts" => length(target_contacts),
        "matched_binding_site_atoms" => length(site_atom_pairs),
        "binding_site_rmsd_after_global_fit_angstrom" => global_site_rmsd,
        "binding_site_local_rmsd_angstrom" => local_site_rmsd,
        "ligand_centroid_displacement_after_global_fit_angstrom" => ligand_centroid_displacement,
        "cutoff_angstrom" => Float64(cutoff),
    )
    limitations = ["contact differences use matched residue identifiers and require review when numbering or homology differs",
                   "geometric changes do not establish changes in binding free energy or antibiotic efficacy"]
    ComparisonResult(reference.metadata.identifier, target.metadata.identifier,
                     alignment, metrics, changed, limitations)
end

function compare_states(items::AbstractVector{<:Tuple})
    rows = DataFrame(label=String[], antibiotic=Union{Missing,String}[],
                     subunit=Union{Missing,String}[], site_residues=Int[],
                     radius_of_gyration_angstrom=Float64[])
    for (structure, ligand, state) in items
        site = identify_binding_site(structure, ligand)
        geometry = geometry_descriptors(structure, site)
        push!(rows, (state.label, state.antibiotic, state.subunit,
                     length(site.residue_keys), geometry["radius_of_gyration_angstrom"]))
    end
    rows
end
