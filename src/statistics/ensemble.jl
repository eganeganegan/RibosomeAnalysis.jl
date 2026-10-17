function autocorrelation(
    values::AbstractVector;
    maxlag::Integer = min(length(values) - 1, 1000),
)
    n = length(values)
    n > 1 || throw(ArgumentError("at least two observations are required"))
    0 <= maxlag < n || throw(ArgumentError("maxlag must be between 0 and n-1"))
    centered = Float64.(values) .- mean(values)
    variance_sum = sum(abs2, centered)
    variance_sum == 0 && return [lag == 0 ? 1.0 : 0.0 for lag = 0:maxlag]
    [
        dot(view(centered, 1:(n-lag)), view(centered, (1+lag):n)) / variance_sum for
        lag = 0:maxlag
    ]
end

function effective_sample_size(
    values::AbstractVector;
    maxlag::Integer = min(length(values) - 1, 1000),
)
    ac = autocorrelation(values; maxlag = maxlag)
    positive = Float64[]
    for value in ac[2:end]
        value <= 0 && break
        push!(positive, value)
    end
    clamp(length(values) / (1 + 2sum(positive)), 1.0, Float64(length(values)))
end

function _quantile_interval(samples, confidence)
    alpha = (1 - confidence) / 2
    (quantile(samples, alpha), quantile(samples, 1 - alpha))
end

function _validate_resampling_options(confidence::Real, iterations::Integer)
    0 < confidence < 1 || throw(ArgumentError("confidence must be between 0 and 1"))
    iterations > 0 || throw(ArgumentError("iterations must be positive"))
    nothing
end

function bootstrap_ci(
    values::AbstractVector;
    statistic = mean,
    confidence::Real = 0.95,
    iterations::Integer = 2000,
    seed::Integer = 1,
)
    isempty(values) && throw(ArgumentError("values are empty"))
    all(isfinite, values) || throw(ArgumentError("values must be finite"))
    _validate_resampling_options(confidence, iterations)
    rng = MersenneTwister(seed)
    draws = [statistic(rand(rng, values, length(values))) for _ = 1:iterations]
    estimate = statistic(values)
    lower, upper = _quantile_interval(draws, confidence)
    (
        estimate = estimate,
        lower = lower,
        upper = upper,
        confidence = Float64(confidence),
        iterations = iterations,
        seed = seed,
        assumption = "independent observations",
    )
end

function block_bootstrap_ci(
    values::AbstractVector;
    statistic = mean,
    block_length::Integer = max(1, round(Int, sqrt(length(values)))),
    confidence::Real = 0.95,
    iterations::Integer = 2000,
    seed::Integer = 1,
)
    n = length(values)
    n > 0 || throw(ArgumentError("values are empty"))
    all(isfinite, values) || throw(ArgumentError("values must be finite"))
    1 <= block_length <= n || throw(ArgumentError("block_length must be between 1 and n"))
    _validate_resampling_options(confidence, iterations)
    rng = MersenneTwister(seed)
    starts = 1:(n-block_length+1)
    draws = Vector{Float64}(undef, iterations)
    for iteration = 1:iterations
        sampled = eltype(values)[]
        while length(sampled) < n
            start = rand(rng, starts)
            append!(sampled, values[start:(start+block_length-1)])
        end
        draws[iteration] = statistic(view(sampled, 1:n))
    end
    lower, upper = _quantile_interval(draws, confidence)
    (
        estimate = statistic(values),
        lower = lower,
        upper = upper,
        confidence = Float64(confidence),
        iterations = iterations,
        seed = seed,
        block_length = block_length,
        assumption = "moving-block bootstrap preserves correlation within each block",
    )
end

function summarize_samples(
    values::AbstractVector;
    correlated::Bool = false,
    block_length::Union{Nothing,Int} = nothing,
    seed::Int = 1,
)
    isempty(values) && throw(ArgumentError("values are empty"))
    ci =
        correlated ?
        block_bootstrap_ci(
            values;
            block_length = something(
                block_length,
                max(1, round(Int, sqrt(length(values)))),
            ),
            seed = seed,
        ) : bootstrap_ci(values; seed = seed)
    Dict{String,Any}(
        "n" => length(values),
        "mean" => mean(values),
        "median" => median(values),
        "standard_deviation" => length(values) > 1 ? std(values) : missing,
        "quantiles" => Dict(
            "0.05" => quantile(values, 0.05),
            "0.25" => quantile(values, 0.25),
            "0.75" => quantile(values, 0.75),
            "0.95" => quantile(values, 0.95),
        ),
        "effective_sample_size" =>
            correlated && length(values) > 1 ? effective_sample_size(values) :
            length(values),
        "confidence_interval" => Dict(string(k) => v for (k, v) in pairs(ci)),
    )
end

"""Compare independent replicas using one summary value per replica."""
function compare_replicates(replicates::AbstractVector{<:AbstractVector}; statistic = mean)
    isempty(replicates) && throw(ArgumentError("at least one replica is required"))
    all(replica -> !isempty(replica), replicates) ||
        throw(ArgumentError("replicas must not be empty"))
    all(replica -> all(isfinite, replica), replicates) ||
        throw(ArgumentError("replica values must be finite"))
    estimates = statistic.(replicates)
    Dict{String,Any}(
        "replicas" => length(replicates),
        "frames_per_replica" => length.(replicates),
        "replica_estimates" => estimates,
        "between_replica_mean" => mean(estimates),
        "between_replica_standard_deviation" =>
            length(estimates) > 1 ? std(estimates) : missing,
        "assumption" => "replicas are independent; trajectory frames are not treated as independent replicas",
    )
end

function _require_columns(table::AbstractDataFrame, columns)
    available = Set(propertynames(table))
    missing_columns = [column for column in columns if column ∉ available]
    isempty(missing_columns) || throw(
        ArgumentError(
            "table is missing required column(s): " * join(string.(missing_columns), ", "),
        ),
    )
    nothing
end

function _state_replica_observations(
    observations::AbstractDataFrame,
    state_column::Symbol,
    replica_column::Symbol,
    value_column::Symbol,
    time_column::Union{Nothing,Symbol},
)
    required =
        isnothing(time_column) ? (state_column, replica_column, value_column) :
        (state_column, replica_column, value_column, time_column)
    _require_columns(observations, required)
    nrow(observations) > 0 || throw(ArgumentError("observations are empty"))
    for column in (state_column, replica_column, value_column)
        any(ismissing, observations[!, column]) &&
            throw(ArgumentError("column $(repr(string(column))) contains missing values"))
    end
    values = try
        Float64.(observations[!, value_column])
    catch
        throw(ArgumentError("column $(repr(string(value_column))) must be numeric"))
    end
    all(isfinite, values) ||
        throw(ArgumentError("column $(repr(string(value_column))) must be finite"))
    normalized = DataFrame(
        state = string.(observations[!, state_column]),
        replica = string.(observations[!, replica_column]),
        value = values,
        time = Float64.(axes(observations, 1)),
    )
    if !isnothing(time_column)
        any(ismissing, observations[!, time_column]) && throw(
            ArgumentError("column $(repr(string(time_column))) contains missing values"),
        )
        times = try
            Float64.(observations[!, time_column])
        catch
            throw(ArgumentError("column $(repr(string(time_column))) must be numeric"))
        end
        all(isfinite, times) ||
            throw(ArgumentError("column $(repr(string(time_column))) must be finite"))
        normalized.time = times
    end
    return normalized
end

"""
Summarize a per-frame metric without pooling independent replicas.

Each replica is reduced to one estimate before state-level statistics are
calculated. Within-replica intervals use a moving-block bootstrap. State-level
intervals resample the independent replica estimates and are reported only when
at least two replicas are available.
"""
function summarize_state_replicas(
    observations::AbstractDataFrame;
    state::Symbol = :state,
    replica::Symbol = :replica,
    value::Symbol = :value,
    time::Union{Nothing,Symbol} = nothing,
    statistic = mean,
    confidence::Real = 0.95,
    iterations::Integer = 2000,
    seed::Integer = 1,
    block_length::Union{Nothing,Integer} = nothing,
)
    _validate_resampling_options(confidence, iterations)
    isnothing(block_length) ||
        block_length > 0 ||
        throw(ArgumentError("block_length must be positive"))
    normalized = _state_replica_observations(observations, state, replica, value, time)
    sort!(normalized, [:state, :replica, :time])
    replica_rows = DataFrame(
        state = String[],
        replica = String[],
        frames = Int[],
        estimate = Float64[],
        effective_sample_size = Float64[],
        block_length = Int[],
        within_replica_lower = Union{Missing,Float64}[],
        within_replica_upper = Union{Missing,Float64}[],
    )
    for (group_index, group) in enumerate(groupby(normalized, [:state, :replica]))
        values = group.value
        !isnothing(block_length) &&
            block_length > length(values) &&
            throw(
                ArgumentError(
                    "block_length exceeds the $(length(values)) observations in state $(repr(group.state[1])) replica $(repr(group.replica[1]))",
                ),
            )
        used_block_length =
            isnothing(block_length) ? max(1, round(Int, sqrt(length(values)))) :
            Int(block_length)
        estimate = Float64(statistic(values))
        isfinite(estimate) || throw(ArgumentError("statistic must return a finite value"))
        lower, upper = if length(values) > 1
            interval = block_bootstrap_ci(
                values;
                statistic = statistic,
                block_length = used_block_length,
                confidence = confidence,
                iterations = iterations,
                seed = seed + group_index - 1,
            )
            (Float64(interval.lower), Float64(interval.upper))
        else
            (missing, missing)
        end
        push!(
            replica_rows,
            (
                group.state[1],
                group.replica[1],
                length(values),
                estimate,
                length(values) > 1 ? effective_sample_size(values) : 1.0,
                used_block_length,
                lower,
                upper,
            ),
        )
    end

    state_rows = DataFrame(
        state = String[],
        replicas = Int[],
        total_frames = Int[],
        mean_replica_estimate = Float64[],
        median_replica_estimate = Float64[],
        between_replica_standard_deviation = Union{Missing,Float64}[],
        lower = Union{Missing,Float64}[],
        upper = Union{Missing,Float64}[],
        confidence = Float64[],
    )
    for (state_index, group) in enumerate(groupby(replica_rows, :state))
        estimates = group.estimate
        lower, upper = if length(estimates) > 1
            interval = bootstrap_ci(
                estimates;
                statistic = mean,
                confidence = confidence,
                iterations = iterations,
                seed = seed + nrow(replica_rows) + state_index - 1,
            )
            (Float64(interval.lower), Float64(interval.upper))
        else
            (missing, missing)
        end
        push!(
            state_rows,
            (
                group.state[1],
                length(estimates),
                sum(group.frames),
                mean(estimates),
                median(estimates),
                length(estimates) > 1 ? std(estimates) : missing,
                lower,
                upper,
                Float64(confidence),
            ),
        )
    end
    return (replicas = replica_rows, states = state_rows)
end

"""
Calculate target-minus-reference differences from independent replica estimates.

Intervals are omitted when either state has fewer than two independent
replicas. No hypothesis test or significance claim is made.
"""
function pairwise_state_differences(
    replica_summary::AbstractDataFrame;
    confidence::Real = 0.95,
    iterations::Integer = 2000,
    seed::Integer = 1,
)
    _validate_resampling_options(confidence, iterations)
    _require_columns(replica_summary, (:state, :replica, :estimate))
    nrow(replica_summary) > 0 || throw(ArgumentError("replica summary is empty"))
    any(ismissing, replica_summary.estimate) &&
        throw(ArgumentError("replica estimates contain missing values"))
    estimates = Float64.(replica_summary.estimate)
    all(isfinite, estimates) || throw(ArgumentError("replica estimates must be finite"))
    normalized = DataFrame(
        state = string.(replica_summary.state),
        replica = string.(replica_summary.replica),
        estimate = estimates,
    )
    length(unique(collect(zip(normalized.state, normalized.replica)))) ==
    nrow(normalized) ||
        throw(ArgumentError("each state and replica pair must appear exactly once"))
    state_labels = sort!(unique(normalized.state))
    rows = DataFrame(
        reference_state = String[],
        target_state = String[],
        reference_replicas = Int[],
        target_replicas = Int[],
        difference = Float64[],
        lower = Union{Missing,Float64}[],
        upper = Union{Missing,Float64}[],
        confidence = Float64[],
    )
    pair_index = 0
    for reference_index = 1:(length(state_labels)-1)
        for target_index = (reference_index+1):length(state_labels)
            pair_index += 1
            reference_state = state_labels[reference_index]
            target_state = state_labels[target_index]
            reference_values = normalized.estimate[normalized.state.==reference_state]
            target_values = normalized.estimate[normalized.state.==target_state]
            difference = mean(target_values) - mean(reference_values)
            lower, upper = if length(reference_values) > 1 && length(target_values) > 1
                rng = MersenneTwister(seed + pair_index - 1)
                draws = Vector{Float64}(undef, iterations)
                for iteration = 1:iterations
                    draws[iteration] =
                        mean(rand(rng, target_values, length(target_values))) -
                        mean(rand(rng, reference_values, length(reference_values)))
                end
                _quantile_interval(draws, confidence)
            else
                (missing, missing)
            end
            push!(
                rows,
                (
                    reference_state,
                    target_state,
                    length(reference_values),
                    length(target_values),
                    difference,
                    lower,
                    upper,
                    Float64(confidence),
                ),
            )
        end
    end
    return rows
end

"""
Run the complete replica-aware state comparison workflow.

The returned named tuple contains per-replica summaries, per-state summaries,
pairwise target-minus-reference differences, and explicit statistical
assumptions.
"""
function compare_state_ensembles(observations::AbstractDataFrame; kwargs...)
    summaries = summarize_state_replicas(observations; kwargs...)
    selected = Dict(kwargs)
    differences = pairwise_state_differences(
        summaries.replicas;
        confidence = get(selected, :confidence, 0.95),
        iterations = get(selected, :iterations, 2000),
        seed = get(selected, :seed, 1),
    )
    return (
        replicas = summaries.replicas,
        states = summaries.states,
        differences = differences,
        assumptions = [
            "trajectory frames within a replica may be correlated and are summarized with a moving-block bootstrap",
            "independent replicas receive equal weight regardless of frame count",
            "state-level uncertainty resamples replica estimates and requires at least two replicas per state",
            "pairwise differences are target minus reference and do not establish biological or clinical significance",
        ],
    )
end
