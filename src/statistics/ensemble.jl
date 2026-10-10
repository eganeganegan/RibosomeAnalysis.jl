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
        dot(view(centered, 1:n-lag), view(centered, 1+lag:n)) / variance_sum for
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

function bootstrap_ci(
    values::AbstractVector;
    statistic = mean,
    confidence::Real = 0.95,
    iterations::Integer = 2000,
    seed::Integer = 1,
)
    isempty(values) && throw(ArgumentError("values are empty"))
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
    1 <= block_length <= n || throw(ArgumentError("block_length must be between 1 and n"))
    rng = MersenneTwister(seed)
    starts = 1:(n-block_length+1)
    draws = Vector{Float64}(undef, iterations)
    for iteration = 1:iterations
        sampled = eltype(values)[]
        while length(sampled) < n
            start = rand(rng, starts)
            append!(sampled, values[start:start+block_length-1])
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
