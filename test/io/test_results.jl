@testset "Reproducible outputs" begin
    mktempdir() do root
        input = joinpath(root, "input.txt")
        write(input, "reproducible input\n")
        output = prepare_results(joinpath(root, "analysis"))
        @test isdir(joinpath(output, "tables"))
        @test write_summary(output, Dict("value" => 1)) == joinpath(output, "summary.json")
        metadata = write_metadata(
            output;
            inputs = [input, joinpath(root, "missing.txt")],
            parameters = Dict("cutoff_angstrom" => 5.0),
            random_seed = 1,
        )
        @test write_table(output, "values.csv", DataFrame(x = [1, 2])) ==
              joinpath(output, "tables", "values.csv")
        @test isfile(joinpath(output, "summary.json"))
        @test isfile(joinpath(output, "metadata.json"))
        @test isfile(joinpath(output, "tables", "values.csv"))
        @test metadata["inputs"][1]["sha256"] == bytes2hex(open(SHA.sha256, input))
        @test isnothing(metadata["inputs"][2]["sha256"])
        @test metadata["package_version"] == string(pkgversion(RibosomeAnalysis))
        @test metadata["dependency_versions"]["BioStructures"] ==
              string(pkgversion(BioStructures))
        @test Set(keys(metadata["dependency_versions"])) == Set([
            "BioStructures",
            "CSV",
            "DataFrames",
            "JSON3",
            "NearestNeighbors",
            "StaticArrays",
            "StatsBase",
        ],)
        @test_throws ArgumentError prepare_results(output)
        @test prepare_results(output; overwrite = true) == output
        @test_throws ArgumentError write_table(output, "values.tsv", DataFrame(x = [1]))
    end
end
