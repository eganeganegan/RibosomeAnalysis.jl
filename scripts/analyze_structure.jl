include(joinpath(@__DIR__, "common.jl"))

function main(arguments=ARGS)
    options = parse_cli(arguments)
    input = require_option(options, "structure")
    output = require_option(options, "output")
    structure = load_structure(input)
    directory = prepare_results(output)
    summary = structure_summary(structure)
    summary["data_quality_warnings"] = incomplete_information(structure)
    summary["candidate_antibiotics"] = [l.identifier for l in identify_ligands(structure)]
    write_summary(directory, summary)
    write_metadata(directory; inputs=[input], parameters=Dict("workflow" => "structure_analysis"))
    println("Wrote structure analysis to $directory")
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()
