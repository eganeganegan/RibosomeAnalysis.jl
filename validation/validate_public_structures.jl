using BioStructures
using RibosomeAnalysis
using SHA
using TOML
using Test

function file_sha256(path)
    bytes2hex(open(SHA.sha256, path))
end

function verified_path(data_directory, manifest, key)
    record = manifest["structures"][key]
    path = joinpath(data_directory, record["filename"])
    isfile(path) ||
        error("missing validation input $path; run download_public_structures.jl first")
    file_sha256(path) == record["sha256"] || error("checksum mismatch for $path")
    path
end

function brute_force_contact_count(structure, ligand; cutoff = 5.0)
    ligand_set = Set(ligand.atom_indices)
    count = 0
    minimum_distance = Inf
    for (index, atom) in enumerate(structure.atoms)
        index in ligand_set && continue
        uppercase(atom.element) == "H" && continue
        for ligand_index in ligand.atom_indices
            ligand_atom = structure.atoms[ligand_index]
            uppercase(ligand_atom.element) == "H" && continue
            distance = sqrt(sum(abs2, atom.position - ligand_atom.position))
            minimum_distance = min(minimum_distance, distance)
            count += distance <= cutoff
        end
    end
    count, minimum_distance
end

function main(arguments = ARGS)
    data_directory =
        isempty(arguments) ? joinpath(@__DIR__, "data") : abspath(first(arguments))
    manifest = TOML.parsefile(joinpath(@__DIR__, "public_structures.toml"))
    one_j7t_path = verified_path(data_directory, manifest, "1J7T_cif")
    four_v7u_path = verified_path(data_directory, manifest, "4V7U_cif")

    @testset "RCSB 1J7T" begin
        structure = load_structure(one_j7t_path; identifier = "1J7T")
        raw = BioStructures.read(
            one_j7t_path,
            BioStructures.MMCIFFormat;
            remove_disorder = false,
            read_std_atoms = true,
            read_het_atoms = true,
            gzip = true,
        )
        @test length(BioStructures.collectatoms(raw; expand_disordered = true)) == 1036
        @test length(structure.atoms) == 1036
        ligand = select_ligand(
            structure;
            residue_name = "PAR",
            chain_id = "A",
            residue_number = 45,
        )
        site = identify_binding_site(structure, ligand; cutoff = 5.0)
        contacts = atomic_contacts(structure, ligand; cutoff = 5.0)
        brute_count, brute_minimum = brute_force_contact_count(structure, ligand)
        @test length(ligand.atom_indices) == 42
        @test length(contacts) == brute_count == 356
        @test length(site.contacting_atom_indices) == 88
        @test length(site.residue_keys) == 18
        @test minimum(contact.distance_angstrom for contact in contacts) ≈ brute_minimum atol =
            1e-12
        @test brute_minimum ≈ 2.513807072947321 atol = 1e-9
    end

    @testset "RCSB 4V7U" begin
        structure = load_structure(four_v7u_path; identifier = "4V7U")
        summary = structure_summary(structure)
        @test summary["atoms"] == 284525
        @test summary["residues"] == 22506
        @test summary["chains"] == 104
        ligand = select_ligand(
            structure;
            residue_name = "ERY",
            chain_id = "BA",
            residue_number = 3136,
        )
        site = identify_binding_site(structure, ligand; cutoff = 5.0)
        contacts = atomic_contacts(structure, ligand; cutoff = 5.0)
        @test length(ligand.atom_indices) == 51
        @test length(contacts) == 283
        @test length(site.contacting_atom_indices) == 96
        @test length(site.residue_keys) == 12
        @test minimum(contact.distance_angstrom for contact in contacts) ≈ 2.525122967302783 atol =
            1e-9
    end
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()
