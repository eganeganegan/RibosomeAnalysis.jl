module AnalyzeStructureCLI
include(joinpath(@__DIR__, "..", "..", "scripts", "analyze_structure.jl"))
end

module AnalyzeBindingSiteCLI
include(joinpath(@__DIR__, "..", "..", "scripts", "analyze_binding_site.jl"))
end

module AnalyzeTrajectoryCLI
include(joinpath(@__DIR__, "..", "..", "scripts", "analyze_trajectory.jl"))
end

module CompareStructuresCLI
include(joinpath(@__DIR__, "..", "..", "scripts", "compare_structures.jl"))
end

module CompareStateEnsemblesCLI
include(joinpath(@__DIR__, "..", "..", "scripts", "compare_state_ensembles.jl"))
end

@testset "CLI argument validation" begin
    parse_cli = AnalyzeStructureCLI.parse_cli
    @test parse_cli(["--structure", "input.cif"]; allowed_options = ("structure",)) ==
          Dict("structure" => "input.cif")
    @test_throws ArgumentError parse_cli(["structure", "input.cif"])
    @test_throws ArgumentError parse_cli(["--structure"])
    @test_throws ArgumentError parse_cli([
        "--structure",
        "one.cif",
        "--structure",
        "two.cif",
    ],)
    @test_throws ArgumentError parse_cli(
        ["--structuer", "input.cif"];
        allowed_options = ("structure",),
    )
    @test_throws ArgumentError AnalyzeBindingSiteCLI.parse_option(
        Dict("cutoff" => "nearby"),
        "cutoff",
        Float64,
    )
    @test CompareStructuresCLI.parse_chain_map("A:B, C:D") == Dict("A" => "B", "C" => "D")
    @test_throws ArgumentError CompareStructuresCLI.parse_chain_map("")
    @test_throws ArgumentError CompareStructuresCLI.parse_chain_map("A")
    @test_throws ArgumentError CompareStructuresCLI.parse_chain_map("A:B,A:C")
end

@testset "CLI workflows" begin
    fixture = joinpath(@__DIR__, "..", "fixtures", "synthetic.pdb")
    mktempdir() do root
        structure_output = joinpath(root, "structure")
        @test isnothing(
            AnalyzeStructureCLI.main([
                "--structure",
                fixture,
                "--output",
                structure_output,
            ]),
        )
        @test all(
            isfile,
            (
                joinpath(structure_output, "metadata.json"),
                joinpath(structure_output, "summary.json"),
            ),
        )

        binding_output = joinpath(root, "binding")
        @test isnothing(
            AnalyzeBindingSiteCLI.main([
                "--structure",
                fixture,
                "--ligand",
                "PAR",
                "--ligand-chain",
                "L",
                "--ligand-residue",
                "100",
                "--cutoff",
                "5.0",
                "--output",
                binding_output,
            ]),
        )
        @test all(
            isfile,
            (
                joinpath(binding_output, "tables", "atomic_contacts.csv"),
                joinpath(binding_output, "structures", "binding_site.pml"),
                joinpath(binding_output, "structures", "binding_site.cxc"),
            ),
        )

        comparison_output = joinpath(root, "comparison")
        @test isnothing(
            CompareStructuresCLI.main([
                "--reference",
                fixture,
                "--target",
                fixture,
                "--chain-map",
                "A:A,B:B",
                "--output",
                comparison_output,
            ]),
        )
        @test isfile(joinpath(comparison_output, "summary.json"))

        ensemble_input = joinpath(root, "ensemble.csv")
        CSV.write(
            ensemble_input,
            DataFrame(
                state = repeat(["A", "B"], inner = 12),
                replica = repeat(repeat(["r1", "r2"], inner = 6), 2),
                value = vcat(collect(1.0:12.0), collect(3.0:14.0)),
            ),
        )
        ensemble_output = joinpath(root, "ensemble")
        @test isnothing(
            CompareStateEnsemblesCLI.main([
                "--input",
                ensemble_input,
                "--output",
                ensemble_output,
                "--metric",
                "controlled metric",
                "--unit",
                "angstrom",
                "--iterations",
                "50",
                "--seed",
                "7",
            ]),
        )
        @test all(
            isfile,
            (
                joinpath(ensemble_output, "tables", "replica_summary.csv"),
                joinpath(ensemble_output, "tables", "state_summary.csv"),
                joinpath(ensemble_output, "tables", "pairwise_state_differences.csv"),
                joinpath(ensemble_output, "figures", "state_ensemble.pdf"),
                joinpath(ensemble_output, "summary.json"),
                joinpath(ensemble_output, "metadata.json"),
            ),
        )

        @test_throws ArgumentError AnalyzeTrajectoryCLI.main([
            "--topology",
            fixture,
            "--trajectory",
            joinpath(root, "missing.xtc"),
            "--output",
            joinpath(root, "trajectory"),
        ],)
    end
end

@testset "CLI user-facing failures" begin
    mktemp() do _, error_output
        code = redirect_stderr(error_output) do
            AnalyzeStructureCLI.run_cli(AnalyzeStructureCLI.main, ["--structure", "input.pdb"])
        end
        flush(error_output)
        seekstart(error_output)
        message = read(error_output, String)
        @test code == 2
        @test occursin("Error: required option --output was not supplied", message)
        @test !occursin("Stacktrace", message)
    end
end

@testset "Optional extension failures" begin
    error = try
        RibosomeAnalysis.load_trajectory("topology.pdb", "trajectory.xtc")
        nothing
    catch caught
        caught
    end
    @test error isa ArgumentError
    @test occursin("install PythonCall and MDAnalysis", sprint(showerror, error))

    project = dirname(Base.active_project())
    expression = """
    using RibosomeAnalysis
    @assert isnothing(Base.get_extension(RibosomeAnalysis, :RibosomeAnalysisCairoMakieExt))
    for plotter in (
        plot_rmsd,
        plot_rmsf,
        plot_contact_heatmap,
        plot_pca,
        plot_distance_distribution,
        plot_cluster_populations,
        plot_binding_site_comparison,
        plot_state_metric,
        plot_state_ensemble,
    )
        failed_as_expected = false
        try
            plotter()
        catch error
            error isa ArgumentError || rethrow()
            occursin("optional CairoMakie extension", sprint(showerror, error)) || rethrow()
            failed_as_expected = true
        end
        failed_as_expected || exit(1)
    end
    exit(0)
    """
    command = `$(Base.julia_cmd()) --startup-file=no --project=$project -e $expression`
    @test success(command)
end
