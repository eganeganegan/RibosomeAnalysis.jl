@testset "Geometry descriptors" begin
    xyz = [1.0 -1 0 0; 0 0 1 -1; 0 0 0 0]
    @test centroid(xyz) ≈ zeros(3)
    @test radius_of_gyration(xyz) ≈ 1.0
    @test spatial_extent(xyz) ≈ [2.0, 2.0, 0.0]
    @test length(pairwise_distances(xyz)) == 6
    @test local_contact_density(xyz; cutoff=1.5) ≈ 1.0
end
