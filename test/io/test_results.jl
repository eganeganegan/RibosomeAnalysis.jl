@testset "Reproducible outputs" begin
    mktempdir() do root
        output = prepare_results(joinpath(root, "analysis"))
        @test isdir(joinpath(output, "tables"))
        write_summary(output, Dict("value" => 1))
        metadata = write_metadata(
            output;
            parameters = Dict("cutoff_angstrom" => 5.0),
            random_seed = 1,
        )
        write_table(output, "values.csv", DataFrame(x = [1, 2]))
        @test isfile(joinpath(output, "summary.json"))
        @test isfile(joinpath(output, "metadata.json"))
        @test isfile(joinpath(output, "tables", "values.csv"))
        @test metadata["package_version"] == "0.1.0"
        @test metadata["dependency_versions"]["BioStructures"] ==
              string(pkgversion(BioStructures))
        @test length(metadata["dependency_versions"]) == 8
        @test_throws ArgumentError prepare_results(output)
    end
end
