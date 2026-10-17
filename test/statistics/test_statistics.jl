@testset "Time-series statistics" begin
    values = collect(1.0:100.0)
    ac = autocorrelation(values; maxlag = 5)
    @test ac[1] ≈ 1.0
    @test 1 <= effective_sample_size(values) <= length(values)
    ci = bootstrap_ci(values; iterations = 200, seed = 42)
    @test ci.lower < mean(values) < ci.upper
    block_ci = block_bootstrap_ci(values; iterations = 200, block_length = 10, seed = 42)
    @test block_ci.lower < mean(values) < block_ci.upper
    summary = summarize_samples(values; correlated = true, block_length = 10)
    @test summary["n"] == 100
    replicas = compare_replicates([values, values .+ 1])
    @test replicas["replicas"] == 2
    @test_throws ArgumentError bootstrap_ci([1.0, Inf])
    @test_throws ArgumentError bootstrap_ci([1.0, 2.0]; confidence = 1.0)
    @test_throws ArgumentError block_bootstrap_ci([1.0, 2.0]; iterations = 0)
    @test_throws ArgumentError compare_replicates([[1.0], Float64[]])
end

@testset "Replica-aware state ensembles" begin
    observations = DataFrame(state = String[], replica = String[], value = Float64[])
    for state in ("A", "B"), replica = 1:3, frame = 1:12
        offset = state == "A" ? 0.0 : 2.0
        push!(observations, (state, "r$replica", offset + replica + 0.01 * frame))
    end
    analysis =
        compare_state_ensembles(observations; iterations = 200, seed = 41, block_length = 3)
    @test size(analysis.replicas) == (6, 8)
    @test size(analysis.states) == (2, 9)
    @test size(analysis.differences) == (1, 8)
    @test all(==(3), analysis.states.replicas)
    @test analysis.differences.difference[1] ≈ 2.0
    @test analysis.differences.lower[1] <=
          analysis.differences.difference[1] <=
          analysis.differences.upper[1]
    repeated =
        compare_state_ensembles(observations; iterations = 200, seed = 41, block_length = 3)
    @test isequal(analysis.states, repeated.states)
    @test isequal(analysis.differences, repeated.differences)

    unbalanced = DataFrame(
        state = [fill("A", 101); fill("B", 2)],
        replica = [fill("long", 100); "short"; "one"; "two"],
        value = [zeros(100); 10.0; 20.0; 22.0],
    )
    summary = summarize_state_replicas(unbalanced; iterations = 100, seed = 2)
    @test summary.states.mean_replica_estimate[summary.states.state.=="A"] == [5.0]
    @test summary.states.mean_replica_estimate[summary.states.state.=="B"] == [21.0]
    @test all(==(2), summary.states.replicas)

    timed = DataFrame(
        state = fill("A", 4),
        replica = fill("r1", 4),
        time_ps = [3.0, 1.0, 4.0, 2.0],
        value = [3.0, 1.0, 4.0, 2.0],
    )
    timed_summary =
        summarize_state_replicas(timed; time = :time_ps, iterations = 20, block_length = 2)
    @test timed_summary.replicas.estimate == [2.5]

    single = summarize_state_replicas(
        DataFrame(state = ["A"], replica = ["one"], value = [1.0]);
        iterations = 20,
    )
    @test ismissing(single.states.lower[1])
    @test ismissing(single.replicas.within_replica_lower[1])
    @test_throws ArgumentError summarize_state_replicas(
        DataFrame(state = ["A"], value = [1.0]),
    )
    @test_throws ArgumentError summarize_state_replicas(
        DataFrame(state = ["A"], replica = ["one"], value = [NaN]),
    )
    @test_throws ArgumentError summarize_state_replicas(
        DataFrame(state = ["A"], replica = ["one"], value = [1.0]);
        block_length = 2,
    )
    @test_throws ArgumentError pairwise_state_differences(
        DataFrame(state = ["A", "A"], replica = ["one", "one"], estimate = [1.0, 2.0]),
    )
end
