@testset "Binding-site comparison" begin
    path = joinpath(@__DIR__, "..", "fixtures", "synthetic.pdb")
    reference = load_structure(path; identifier = "reference")
    target = load_structure(path; identifier = "target")
    reference_ligand = select_ligand(reference; residue_name = "PAR")
    target_ligand = select_ligand(target; residue_name = "PAR")
    comparison = compare_binding_sites(reference, target, reference_ligand, target_ligand)
    @test comparison.metrics["binding_site_rmsd_after_global_fit_angstrom"] < 1e-12
    @test comparison.metrics["binding_site_local_rmsd_angstrom"] < 1e-12
    @test comparison.metrics["ligand_centroid_displacement_after_global_fit_angstrom"] <
          1e-12
end
