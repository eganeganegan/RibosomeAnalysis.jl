const _PDB_SUFFIXES = (".pdb", ".ent")
const _CIF_SUFFIXES = (".cif", ".mmcif")

_clean_missing(x) = x in (".", "?", "", nothing) ? missing : x

function _parse_int_or_missing(x)
    y = _clean_missing(x)
    ismissing(y) && return missing
    try
        return parse(Int, string(y))
    catch
        return missing
    end
end

function _detect_structure_format(path::AbstractString)
    p = lowercase(path)
    p = endswith(p, ".gz") ? chop(p; tail = 3) : p
    any(endswith(p, s) for s in _PDB_SUFFIXES) && return :pdb
    any(endswith(p, s) for s in _CIF_SUFFIXES) && return :mmcif
    throw(ArgumentError("unsupported structure format for '$path'; expected PDB or mmCIF"))
end

function _cif_column(cif, name::String, n::Int; default = missing)
    try
        values = cif[name]
        return length(values) == n ? collect(values) : fill(default, n)
    catch
        return fill(default, n)
    end
end

"""Read author/label identifiers directly from mmCIF without replacing either namespace."""
function _mmcif_identifiers(path::AbstractString; gzip::Bool = false)
    cif = BioStructures.MMCIFDict(path; gzip = gzip)
    ids = try
        collect(cif["_atom_site.id"])
    catch
        Any[]
    end
    n = length(ids)
    n == 0 && return Dict{Int,NamedTuple}()
    auth_chain = _cif_column(cif, "_atom_site.auth_asym_id", n)
    label_chain = _cif_column(cif, "_atom_site.label_asym_id", n)
    auth_seq = _cif_column(cif, "_atom_site.auth_seq_id", n)
    label_seq = _cif_column(cif, "_atom_site.label_seq_id", n)
    out = Dict{Int,NamedTuple}()
    for i in eachindex(ids)
        serial = _parse_int_or_missing(ids[i])
        ismissing(serial) && continue
        out[serial] = (
            author_chain = ismissing(_clean_missing(auth_chain[i])) ? "" :
                           string(auth_chain[i]),
            label_chain = ismissing(_clean_missing(label_chain[i])) ? "" :
                          string(label_chain[i]),
            author_number = _parse_int_or_missing(auth_seq[i]),
            label_number = _parse_int_or_missing(label_seq[i]),
        )
    end
    out
end

"""
    load_structure(path; identifier=basename(path))

Load all models, hetero atoms, and alternate locations from PDB or mmCIF. mmCIF
is recommended for ribosomes. Coordinates are normalized to Angstrom by the
format parser. Author and label identifiers are retained where mmCIF provides them.
"""
function load_structure(path::AbstractString; identifier::AbstractString = basename(path))
    isfile(path) || throw(ArgumentError("structure file does not exist: $path"))
    format = _detect_structure_format(path)
    gzip = endswith(lowercase(path), ".gz")
    format === :pdb &&
        gzip &&
        throw(
            ArgumentError(
                "compressed PDB is not supported by BioStructures; decompress it or use compressed mmCIF",
            ),
        )
    format_type = format === :pdb ? BioStructures.PDBFormat : BioStructures.MMCIFFormat
    raw = if format === :pdb
        BioStructures.read(
            path,
            format_type;
            remove_disorder = false,
            read_std_atoms = true,
            read_het_atoms = true,
        )
    else
        BioStructures.read(
            path,
            format_type;
            remove_disorder = false,
            read_std_atoms = true,
            read_het_atoms = true,
            gzip = gzip,
        )
    end
    raw_atoms = BioStructures.collectatoms(raw; expand_disordered = true)
    cif_ids =
        format === :mmcif ? _mmcif_identifiers(path; gzip = gzip) : Dict{Int,NamedTuple}()
    atoms = Vector{AtomRecord}(undef, length(raw_atoms))
    for (i, atom) in enumerate(raw_atoms)
        serial = BioStructures.serial(atom)
        ids = get(
            cif_ids,
            serial,
            (
                author_chain = BioStructures.chainid(atom),
                label_chain = BioStructures.chainid(atom),
                author_number = BioStructures.resnumber(atom),
                label_number = missing,
            ),
        )
        chain =
            !isempty(ids.author_chain) ? ids.author_chain :
            (!isempty(ids.label_chain) ? ids.label_chain : BioStructures.chainid(atom))
        xyz = BioStructures.coords(atom)
        atoms[i] = AtomRecord(
            serial,
            BioStructures.atomname(atom),
            BioStructures.resname(atom),
            chain,
            ids.author_chain,
            ids.label_chain,
            BioStructures.resnumber(atom),
            ids.author_number,
            ids.label_number,
            BioStructures.inscode(atom),
            BioStructures.altlocid(atom),
            BioStructures.modelnumber(atom),
            BioStructures.occupancy(atom),
            BioStructures.tempfactor(atom),
            uppercase(BioStructures.element(atom)),
            BioStructures.ishetero(atom),
            SVector{3,Float64}(xyz),
        )
    end
    annotations = Dict{String,Any}(
        "compressed" => gzip,
        "author_and_label_ids_retained" => format === :mmcif,
    )
    meta = StructureMetadata(
        string(identifier),
        abspath(path),
        format,
        "BioStructures.jl",
        annotations,
    )
    StructureData(atoms, meta)
end

function structure_summary(structure::StructureData)
    residues = residue_records(structure)
    Dict{String,Any}(
        "identifier" => structure.metadata.identifier,
        "format" => String(structure.metadata.format),
        "atoms" => length(structure.atoms),
        "residues" => length(residues),
        "chains" => length(unique((a.model_id, a.chain_id) for a in structure.atoms)),
        "models" => sort(unique(a.model_id for a in structure.atoms)),
        "hetero_atoms" => count(a -> a.hetero, structure.atoms),
        "alternate_locations" =>
            count(a -> a.alternate_location ∉ (' ', '\0'), structure.atoms),
        "molecule_types" => countmap(String(r.molecule_type) for r in residues),
    )
end

function incomplete_information(structure::StructureData)
    warnings = String[]
    any(isempty(a.element) for a in structure.atoms) &&
        push!(warnings, "one or more atom elements are missing")
    any(a.occupancy == 0.0 for a in structure.atoms) &&
        push!(warnings, "one or more atoms have zero occupancy")
    if structure.metadata.format === :pdb
        push!(
            warnings,
            "PDB cannot represent some large-assembly identifiers; prefer source mmCIF",
        )
    end
    for residue in residue_records(structure)
        names = Set(structure.atoms[i].name for i in residue.atom_indices)
        if residue.molecule_type === :protein && !all(in(names), ("N", "CA", "C"))
            push!(
                warnings,
                "protein residue $(residue.key) lacks one or more N/CA/C backbone atoms",
            )
        elseif residue.molecule_type === :rna && !any(in(names), ("P", "C4'", "C4*"))
            push!(warnings, "RNA residue $(residue.key) lacks expected backbone atoms")
        end
    end
    unique(warnings)
end
