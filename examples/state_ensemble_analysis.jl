using RibosomeAnalysis
using DataFrames
using CairoMakie

output = isempty(ARGS) ? joinpath(pwd(), "results", "state-ensemble-example") : ARGS[1]

# Controlled deterministic values illustrate the workflow only. They are not
# measurements from a ribosome simulation and carry no biological conclusion.
observations = DataFrame(state = String[], replica = String[], time_ps = Float64[], value = Float64[])
for (state_index, state) in enumerate(("annotated state A", "annotated state B"))
    for replica in 1:3
        for frame in 0:39
            value =
                2.0 + 0.25 * (state_index - 1) + 0.04 * replica +
                0.12 * sin(2pi * frame / 20 + replica)
            push!(observations, (state, "replica $replica", 5.0 * frame, value))
        end
    end
end

analysis = compare_state_ensembles(
    observations;
    time = :time_ps,
    confidence = 0.95,
    iterations = 1000,
    seed = 17,
    block_length = 6,
)
directory = prepare_results(output)
write_table(directory, "controlled_observations.csv", observations)
write_table(directory, "replica_summary.csv", analysis.replicas)
write_table(directory, "state_summary.csv", analysis.states)
write_table(directory, "pairwise_state_differences.csv", analysis.differences)
plot_state_ensemble(
    analysis.states,
    joinpath(directory, "figures", "state_ensemble.png");
    ylabel = "Controlled distance-like metric (Å)",
)
write_summary(
    directory,
    Dict(
        "controlled_example" => true,
        "states" => nrow(analysis.states),
        "replicas" => nrow(analysis.replicas),
        "assumptions" => analysis.assumptions,
        "interpretation" => "Synthetic demonstration only; not evidence of affinity or efficacy.",
    ),
)
println("Wrote controlled state-ensemble example to $directory")
