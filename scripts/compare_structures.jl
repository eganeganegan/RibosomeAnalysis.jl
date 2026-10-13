include(joinpath(@__DIR__, "common.jl"))

function main(arguments = ARGS)
    options = parse_cli(arguments; allowed_options = ("reference", "target", "output"))
    reference_path = require_option(options, "reference")
    target_path = require_option(options, "target")
    output = require_option(options, "output")
    reference = load_structure(reference_path; identifier = "reference")
    target = load_structure(target_path; identifier = "target")
    alignment = align_structures(reference, target)
    directory = prepare_results(output)
    summary = Dict(
        "rmsd_angstrom" => alignment.rmsd_angstrom,
        "aligned_atoms" => alignment.aligned_atoms,
        "correspondence_coverage" => alignment.correspondence_coverage,
        "unmatched_reference_residues" => string.(alignment.unmatched_reference),
        "unmatched_target_residues" => string.(alignment.unmatched_target),
        "rotation" => alignment.rotation,
        "translation_angstrom" => alignment.translation,
        "scientific_limitation" => "alignment is sequence-aware but chain equivalence and cross-species homology require expert review",
    )
    write_summary(directory, summary)
    write_metadata(
        directory;
        inputs = [reference_path, target_path],
        parameters = Dict("workflow" => "structure_comparison"),
    )
    println("Wrote structural comparison to $directory")
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && exit(run_cli(main))
