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
                "--output",
                comparison_output,
            ]),
        )
        @test isfile(joinpath(comparison_output, "summary.json"))

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
    try
        plot_rmsd([0.0], [0.0], "unused.png")
    catch error
        error isa ArgumentError || rethrow()
        occursin("optional CairoMakie extension", sprint(showerror, error)) || rethrow()
        exit(0)
    end
    exit(1)
    """
    command = `$(Base.julia_cmd()) --startup-file=no --project=$project -e $expression`
    @test success(command)
end
