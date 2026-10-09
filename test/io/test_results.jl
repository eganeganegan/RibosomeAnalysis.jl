@testset "Reproducible outputs" begin
    mktempdir() do root
        output = prepare_results(joinpath(root, "analysis"))
        @test isdir(joinpath(output, "tables"))
        write_summary(output, Dict("value" => 1))
        write_metadata(output; parameters=Dict("cutoff_angstrom" => 5.0), random_seed=1)
        write_table(output, "values.csv", DataFrame(x=[1, 2]))
        @test isfile(joinpath(output, "summary.json"))
        @test isfile(joinpath(output, "metadata.json"))
        @test isfile(joinpath(output, "tables", "values.csv"))
        @test_throws ArgumentError prepare_results(output)
    end
end
