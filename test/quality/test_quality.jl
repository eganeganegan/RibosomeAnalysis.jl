@testset "Package quality" begin
    @test isempty(Test.detect_unbound_args(RibosomeAnalysis))
    @test isempty(Test.detect_ambiguities(RibosomeAnalysis))
    @test all(isdefined(RibosomeAnalysis, name) for name in names(RibosomeAnalysis))
end
