function _contact_summary(contacts)
    counts = Dict{ResidueKey,Int}()
    minimum_distances = Dict{ResidueKey,Float64}()
    for contact in contacts
        counts[contact.residue] = get(counts, contact.residue, 0) + 1
        minimum_distances[contact.residue] =
            min(get(minimum_distances, contact.residue, Inf), contact.distance_angstrom)
    end
    counts, minimum_distances
end

function _comparison_residue_pairs(
    reference::StructureData,
    target::StructureData;
    chain_map = nothing,
)
    pairs, _, _ = match_residues(reference, target; chain_map = chain_map)
    target_by_key = Dict(residue.key => residue for residue in residue_records(target))
    for reference_residue in residue_records(reference)
        reference_residue.molecule_type in (:rna, :protein) && continue
        target_residue = get(target_by_key, reference_residue.key, nothing)
        isnothing(target_residue) && continue
        target_residue.molecule_type == reference_residue.molecule_type || continue
        push!(pairs, (reference_residue, target_residue))
    end
    return pairs
end

function _contact_change_table(
    residue_pairs,
    reference_site::BindingSite,
    target_site::BindingSite,
    reference_contacts,
    target_contacts,
)
    reference_site_keys = Set(reference_site.residue_keys)
    target_site_keys = Set(target_site.residue_keys)
    reference_counts, reference_distances = _contact_summary(reference_contacts)
    target_counts, target_distances = _contact_summary(target_contacts)
    rows = DataFrame(
        change = String[],
        correspondence = Union{Missing,String}[],
        reference_residue = Union{Missing,String}[],
        target_residue = Union{Missing,String}[],
        reference_contact_count = Union{Missing,Int}[],
        target_contact_count = Union{Missing,Int}[],
        contact_count_difference = Union{Missing,Int}[],
        reference_minimum_distance_angstrom = Union{Missing,Float64}[],
        target_minimum_distance_angstrom = Union{Missing,Float64}[],
        minimum_distance_difference_angstrom = Union{Missing,Float64}[],
    )
    mapped_reference = Set{ResidueKey}()
    mapped_target = Set{ResidueKey}()
    for (reference_residue, target_residue) in residue_pairs
        reference_key = reference_residue.key
        target_key = target_residue.key
        reference_present = reference_key in reference_site_keys
        target_present = target_key in target_site_keys
        reference_present || target_present || continue
        push!(mapped_reference, reference_key)
        push!(mapped_target, target_key)
        reference_count =
            reference_present ? get(reference_counts, reference_key, 0) : missing
        target_count = target_present ? get(target_counts, target_key, 0) : missing
        reference_distance =
            reference_present ? get(reference_distances, reference_key, Inf) : missing
        target_distance = target_present ? get(target_distances, target_key, Inf) : missing
        push!(
            rows,
            (
                reference_present && target_present ? "matched" :
                reference_present ? "reference_only" : "target_only",
                reference_residue.molecule_type in (:rna, :protein) ? "sequence" :
                "exact_identifier",
                reference_present ? string(reference_key) : missing,
                target_present ? string(target_key) : missing,
                reference_count,
                target_count,
                reference_present && target_present ? target_count - reference_count :
                missing,
                reference_distance,
                target_distance,
                reference_present && target_present ? target_distance - reference_distance :
                missing,
            ),
        )
    end
    for reference_key in sort!(collect(setdiff(reference_site_keys, mapped_reference)))
        push!(
            rows,
            (
                "reference_only",
                missing,
                string(reference_key),
                missing,
                get(reference_counts, reference_key, 0),
                missing,
                missing,
                get(reference_distances, reference_key, Inf),
                missing,
                missing,
            ),
        )
    end
    for target_key in sort!(collect(setdiff(target_site_keys, mapped_target)))
        push!(
            rows,
            (
                "target_only",
                missing,
                missing,
                string(target_key),
                missing,
                get(target_counts, target_key, 0),
                missing,
                missing,
                get(target_distances, target_key, Inf),
                missing,
            ),
        )
    end
    sort!(rows, [:change, :reference_residue, :target_residue])
end

function compare_binding_sites(
    reference::StructureData,
    target::StructureData,
    reference_ligand::Ligand,
    target_ligand::Ligand;
    cutoff::Real = 5.0,
    chain_map = nothing,
)
    alignment = align_structures(reference, target; chain_map = chain_map)
    reference_site = identify_binding_site(reference, reference_ligand; cutoff = cutoff)
    target_site = identify_binding_site(target, target_ligand; cutoff = cutoff)
    reference_contacts = atomic_contacts(reference, reference_ligand; cutoff = cutoff)
    target_contacts = atomic_contacts(target, target_ligand; cutoff = cutoff)
    residue_pairs = _comparison_residue_pairs(reference, target; chain_map = chain_map)
    reference_site_keys = Set(reference_site.residue_keys)
    target_site_keys = Set(target_site.residue_keys)
    site_atom_pairs = [
        (i, j) for (i, j) in alignment.matched_indices if
        residue_key(reference.atoms[i]) in reference_site_keys &&
        residue_key(target.atoms[j]) in target_site_keys
    ]
    global_site_rmsd = missing
    local_site_rmsd = missing
    if !isempty(site_atom_pairs)
        ref_site_xyz = coordinates(reference, first.(site_atom_pairs))
        target_site_xyz = coordinates(target, last.(site_atom_pairs))
        globally_fitted =
            apply_transform(target_site_xyz, alignment.rotation, alignment.translation)
        global_site_rmsd = rmsd(ref_site_xyz, globally_fitted)
        length(site_atom_pairs) >= 3 &&
            (local_site_rmsd = rmsd(ref_site_xyz, target_site_xyz; superimpose = true))
    end
    target_ligand_fitted = apply_transform(
        coordinates(target, target_ligand.atom_indices),
        alignment.rotation,
        alignment.translation,
    )
    ligand_centroid_displacement = norm(
        centroid(coordinates(reference, reference_ligand.atom_indices)) -
        centroid(target_ligand_fitted),
    )
    changed = _contact_change_table(
        residue_pairs,
        reference_site,
        target_site,
        reference_contacts,
        target_contacts,
    )
    matched_site_residues = count(==("matched"), changed.change)
    metrics = Dict{String,Any}(
        "alignment_rmsd_angstrom" => alignment.rmsd_angstrom,
        "aligned_atoms" => alignment.aligned_atoms,
        "correspondence_coverage" => alignment.correspondence_coverage,
        "reference_binding_site_residues" => length(reference_site.residue_keys),
        "target_binding_site_residues" => length(target_site.residue_keys),
        "reference_atomic_contacts" => length(reference_contacts),
        "target_atomic_contacts" => length(target_contacts),
        "matched_binding_site_residues" => matched_site_residues,
        "reference_only_binding_site_residues" =>
            count(==("reference_only"), changed.change),
        "target_only_binding_site_residues" => count(==("target_only"), changed.change),
        "matched_binding_site_atoms" => length(site_atom_pairs),
        "binding_site_rmsd_after_global_fit_angstrom" => global_site_rmsd,
        "binding_site_local_rmsd_angstrom" => local_site_rmsd,
        "ligand_centroid_displacement_after_global_fit_angstrom" =>
            ligand_centroid_displacement,
        "cutoff_angstrom" => Float64(cutoff),
    )
    limitations = [
        "polymer contact differences use sequence-derived correspondence; chain assignment and homology still require expert review",
        "non-polymer residues correspond only when their complete residue identifiers and classified types agree",
        "geometric changes do not establish changes in binding free energy or antibiotic efficacy",
    ]
    ComparisonResult(
        reference.metadata.identifier,
        target.metadata.identifier,
        alignment,
        metrics,
        changed,
        limitations,
    )
end

function compare_states(items::AbstractVector{<:Tuple}; cutoff::Real = 5.0)
    cutoff > 0 || throw(ArgumentError("cutoff must be positive"))
    rows = DataFrame(
        label = String[],
        organism = Union{Missing,String}[],
        antibiotic = Union{Missing,String}[],
        subunit = Union{Missing,String}[],
        experimental_method = Union{Missing,String}[],
        resolution_angstrom = Union{Missing,Float64}[],
        translation_state = Union{Missing,String}[],
        site_residues = Int[],
        contacting_atoms = Int[],
        radius_of_gyration_angstrom = Union{Missing,Float64}[],
        extent_x_angstrom = Union{Missing,Float64}[],
        extent_y_angstrom = Union{Missing,Float64}[],
        extent_z_angstrom = Union{Missing,Float64}[],
        local_contact_density_per_atom = Union{Missing,Float64}[],
        minimum_ligand_distance_angstrom = Union{Missing,Float64}[],
    )
    for (structure, ligand, state) in items
        site = identify_binding_site(structure, ligand; cutoff = cutoff)
        geometry =
            isempty(site.residue_keys) ? nothing : geometry_descriptors(structure, site)
        extent = isnothing(geometry) ? fill(missing, 3) : geometry["extent_angstrom"]
        push!(
            rows,
            (
                state.label,
                state.organism,
                state.antibiotic,
                state.subunit,
                state.experimental_method,
                state.resolution_angstrom,
                state.translation_state,
                length(site.residue_keys),
                length(site.contacting_atom_indices),
                isnothing(geometry) ? missing : geometry["radius_of_gyration_angstrom"],
                extent[1],
                extent[2],
                extent[3],
                isnothing(geometry) ? missing : geometry["local_contact_density_per_atom"],
                isempty(site.minimum_distances) ? missing :
                minimum(values(site.minimum_distances)),
            ),
        )
    end
    return rows
end
