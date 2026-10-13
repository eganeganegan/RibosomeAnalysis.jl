@testset "CairoMakie plotting extension" begin
    @test !isnothing(Base.get_extension(RibosomeAnalysis, :RibosomeAnalysisCairoMakieExt))
    mktempdir() do directory
        rmsd_path = joinpath(directory, "rmsd.png")
        rmsf_path = joinpath(directory, "rmsf.svg")
        heatmap_path = joinpath(directory, "contacts.pdf")
        pca_path = joinpath(directory, "pca.png")

        @test plot_rmsd(0.0:1.0:2.0, [0.0, 0.2, 0.1], rmsd_path) == rmsd_path
        @test plot_rmsf(1:3, [0.1, 0.3, 0.2], rmsf_path) == rmsf_path
        @test plot_contact_heatmap([0.0 0.5; 1.0 0.25], heatmap_path) == heatmap_path
        @test plot_pca([0.0 1.0 2.0; 2.0 1.0 0.0], pca_path) == pca_path
        @test all(isfile, (rmsd_path, rmsf_path, heatmap_path, pca_path))
        @test_throws ArgumentError plot_rmsd(
            [0.0, 1.0],
            [0.0, 0.1],
            joinpath(directory, "rmsd.txt"),
        )
    end
end
