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
    @test comparison.metrics["matched_binding_site_residues"] ==
          comparison.metrics["reference_binding_site_residues"]
    @test all(==("matched"), comparison.changed_contacts.change)
    @test Set(comparison.changed_contacts.correspondence) ==
          Set(["sequence", "exact_identifier"])
    @test all(==(0), comparison.changed_contacts.contact_count_difference)
    @test all(
        x -> abs(x) < 1e-12,
        comparison.changed_contacts.minimum_distance_difference_angstrom,
    )

    state = StateMetadata(
        "reference state",
        "Escherichia coli",
        "70S",
        "X-ray diffraction",
        2.5,
        "paromomycin",
        "annotated state",
        Dict("source" => "controlled fixture"),
    )
    states = compare_states([(reference, reference_ligand, state)]; cutoff = 5.0)
    @test size(states) == (1, 15)
    @test states.organism == ["Escherichia coli"]
    @test states.site_residues == [
        length(
            identify_binding_site(reference, reference_ligand; cutoff = 5.0).residue_keys,
        ),
    ]
    @test states.minimum_ligand_distance_angstrom[1] > 0
    @test_throws ArgumentError compare_states(
        [(reference, reference_ligand, state)];
        cutoff = 0.0,
    )

    reference_mapped = ResidueKey(1, "A", 1, ' ', "A")
    target_mapped = ResidueKey(1, "A", 11, ' ', "A")
    reference_other = ResidueKey(1, "X", 2, ' ', "MG")
    target_other = ResidueKey(1, "Y", 3, ' ', "MG")
    residue_pairs = [
        (
            ResidueRecord(reference_mapped, [1], :rna),
            ResidueRecord(ResidueKey(1, "A", 10, ' ', "A"), [1], :rna),
        ),
        (
            ResidueRecord(ResidueKey(1, "A", 12, ' ', "A"), [1], :rna),
            ResidueRecord(target_mapped, [1], :rna),
        ),
    ]
    reference_site = BindingSite(
        reference_ligand,
        5.0,
        [reference_mapped, reference_other],
        [[1], [1]],
        [1],
        Dict(reference_mapped => 2.0, reference_other => 3.0),
    )
    target_site = BindingSite(
        target_ligand,
        5.0,
        [target_mapped, target_other],
        [[1], [1]],
        [1],
        Dict(target_mapped => 2.5, target_other => 3.5),
    )
    edge_changes = RibosomeAnalysis._contact_change_table(
        residue_pairs,
        reference_site,
        target_site,
        [
            ContactRecord(1, 1, reference_mapped, 2.0),
            ContactRecord(1, 1, reference_other, 3.0),
        ],
        [ContactRecord(1, 1, target_mapped, 2.5), ContactRecord(1, 1, target_other, 3.5)],
    )
    @test count(==("reference_only"), edge_changes.change) == 2
    @test count(==("target_only"), edge_changes.change) == 2
    @test count(ismissing, edge_changes.correspondence) == 2
end
