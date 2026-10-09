@testset "Time-series statistics" begin
    values = collect(1.0:100.0)
    ac = autocorrelation(values; maxlag=5)
    @test ac[1] ≈ 1.0
    @test 1 <= effective_sample_size(values) <= length(values)
    ci = bootstrap_ci(values; iterations=200, seed=42)
    @test ci.lower < mean(values) < ci.upper
    block_ci = block_bootstrap_ci(values; iterations=200, block_length=10, seed=42)
    @test block_ci.lower < mean(values) < block_ci.upper
    summary = summarize_samples(values; correlated=true, block_length=10)
    @test summary["n"] == 100
    replicas = compare_replicates([values, values .+ 1])
    @test replicas["replicas"] == 2
end
