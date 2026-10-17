include(joinpath(@__DIR__, "common.jl"))
using CSV
using CairoMakie

function main(arguments = ARGS)
    options = parse_cli(
        arguments;
        allowed_options = (
            "input",
            "output",
            "state-column",
            "replica-column",
            "value-column",
            "time-column",
            "metric",
            "unit",
            "confidence",
            "iterations",
            "seed",
            "block-length",
        ),
    )
    input = require_option(options, "input")
    output = require_option(options, "output")
    isfile(input) || throw(ArgumentError("input file does not exist: $input"))
    state_column = Symbol(get(options, "state-column", "state"))
    replica_column = Symbol(get(options, "replica-column", "replica"))
    value_column = Symbol(get(options, "value-column", "value"))
    metric = get(options, "metric", string(value_column))
    unit = get(options, "unit", "unspecified")
    confidence = parse_option(options, "confidence", Float64; default = 0.95)
    iterations = parse_option(options, "iterations", Int; default = 2000)
    seed = parse_option(options, "seed", Int; default = 1)
    block_length = parse_option(options, "block-length", Int; default = nothing)

    observations = CSV.read(input, DataFrame)
    time_column = if haskey(options, "time-column")
        Symbol(options["time-column"])
    elseif :time_ps in propertynames(observations)
        :time_ps
    else
        nothing
    end
    analysis = compare_state_ensembles(
        observations;
        state = state_column,
        replica = replica_column,
        value = value_column,
        time = time_column,
        confidence = confidence,
        iterations = iterations,
        seed = seed,
        block_length = block_length,
    )
    directory = prepare_results(output)
    write_table(directory, "replica_summary.csv", analysis.replicas)
    write_table(directory, "state_summary.csv", analysis.states)
    write_table(directory, "pairwise_state_differences.csv", analysis.differences)
    plot_state_ensemble(
        analysis.states,
        joinpath(directory, "figures", "state_ensemble.pdf");
        ylabel = unit == "unspecified" ? metric : "$metric ($unit)",
        title = "Replica-aware comparison: $metric",
    )
    write_summary(
        directory,
        Dict(
            "metric" => metric,
            "unit" => unit,
            "states" => nrow(analysis.states),
            "replicas" => nrow(analysis.replicas),
            "observations" => nrow(observations),
            "confidence" => confidence,
            "assumptions" => analysis.assumptions,
        ),
    )
    write_metadata(
        directory;
        inputs = [input],
        parameters = Dict(
            "workflow" => "state_ensemble_comparison",
            "state_column" => string(state_column),
            "replica_column" => string(replica_column),
            "value_column" => string(value_column),
            "time_column" => isnothing(time_column) ? nothing : string(time_column),
            "metric" => metric,
            "unit" => unit,
            "confidence" => confidence,
            "iterations" => iterations,
            "seed" => seed,
            "block_length" => block_length,
        ),
    )
    println("Wrote state-ensemble comparison to $directory")
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && exit(run_cli(main))
