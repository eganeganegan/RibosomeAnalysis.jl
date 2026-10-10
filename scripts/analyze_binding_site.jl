include(joinpath(@__DIR__, "common.jl"))

function main(arguments = ARGS)
    options = parse_cli(arguments)
    input = require_option(options, "structure")
    output = require_option(options, "output")
    ligand_name = require_option(options, "ligand")
    cutoff = parse(Float64, get(options, "cutoff", "5.0"))
    structure = load_structure(input)
    ligand = select_ligand(
        structure;
        residue_name = ligand_name,
        chain_id = get(options, "ligand-chain", nothing),
        residue_number = haskey(options, "ligand-residue") ?
                         parse(Int, options["ligand-residue"]) : nothing,
    )
    site = identify_binding_site(structure, ligand; cutoff = cutoff)
    contacts = atomic_contacts(structure, ligand; cutoff = cutoff)
    directory = prepare_results(output)
    summary =
        merge(binding_site_summary(structure, site), geometry_descriptors(structure, site))
    summary["scientific_limitation"] = "contacts are geometric; they do not establish affinity, free energy, or efficacy"
    write_summary(directory, summary)
    write_metadata(
        directory;
        inputs = [input],
        parameters = Dict(
            "workflow" => "binding_site",
            "ligand" => ligand.identifier,
            "cutoff_angstrom" => cutoff,
        ),
    )
    write_table(directory, "atomic_contacts.csv", contact_dataframe(structure, contacts))
    export_selection(
        site,
        joinpath(directory, "structures", "binding_site.pml");
        format = :pymol,
    )
    export_selection(
        site,
        joinpath(directory, "structures", "binding_site.cxc");
        format = :chimerax,
    )
    println("Wrote binding-site analysis to $directory")
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()
