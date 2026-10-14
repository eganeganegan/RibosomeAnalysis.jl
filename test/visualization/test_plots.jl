@testset "CairoMakie plotting extension" begin
    @test !isnothing(Base.get_extension(RibosomeAnalysis, :RibosomeAnalysisCairoMakieExt))
    mktempdir() do directory
        rmsd_path = joinpath(directory, "rmsd.png")
        rmsf_path = joinpath(directory, "rmsf.svg")
        heatmap_path = joinpath(directory, "contacts.pdf")
        pca_path = joinpath(directory, "pca.png")
        distances_path = joinpath(directory, "distances.svg")
        clusters_path = joinpath(directory, "clusters.pdf")
        comparison_path = joinpath(directory, "comparison.png")
        states_path = joinpath(directory, "states.svg")

        @test plot_rmsd(0.0:1.0:2.0, [0.0, 0.2, 0.1], rmsd_path) == rmsd_path
        @test plot_rmsf(1:3, [0.1, 0.3, 0.2], rmsf_path) == rmsf_path
        @test plot_contact_heatmap([0.0 0.5; 1.0 0.25], heatmap_path) == heatmap_path
        @test plot_pca([0.0 1.0 2.0; 2.0 1.0 0.0], pca_path) == pca_path
        @test plot_distance_distribution([2.0, 2.5, 3.0], distances_path; bins = 2) ==
              distances_path
        @test plot_cluster_populations([1, 1, 2], clusters_path) == clusters_path
        comparison = (
            metrics = Dict(
                "reference_binding_site_residues" => 3,
                "target_binding_site_residues" => 4,
                "reference_atomic_contacts" => 8,
                "target_atomic_contacts" => 10,
            ),
        )
        @test plot_binding_site_comparison(comparison, comparison_path) == comparison_path
        states = DataFrame(
            label = ["state A", "state B"],
            radius_of_gyration_angstrom = [4.0, 4.5],
        )
        @test plot_state_metric(states, :radius_of_gyration_angstrom, states_path) ==
              states_path
        @test all(
            isfile,
            (
                rmsd_path,
                rmsf_path,
                heatmap_path,
                pca_path,
                distances_path,
                clusters_path,
                comparison_path,
                states_path,
            ),
        )
        @test_throws ArgumentError plot_rmsd(
            [0.0, 1.0],
            [0.0, 0.1],
            joinpath(directory, "rmsd.txt"),
        )
        @test_throws ArgumentError plot_distance_distribution(Float64[], distances_path)
        @test_throws ArgumentError plot_cluster_populations([0, 1], clusters_path)
        @test_throws ArgumentError plot_state_metric(states, :missing_metric, states_path)
    end
end
