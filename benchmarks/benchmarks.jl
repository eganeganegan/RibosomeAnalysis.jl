using BenchmarkTools
using RibosomeAnalysis

fixture = joinpath(@__DIR__, "..", "test", "fixtures", "synthetic.pdb")
structure = load_structure(fixture)
ligand = select_ligand(structure; residue_name = "PAR")
reference = [0.0 1 0 0; 0 0 1 0; 0 0 0 1]
trajectory = TrajectoryData(repeat(reshape(reference, 3, 4, 1), 1, 1, 100))
structure_trajectory = TrajectoryData(
    repeat(reshape(coordinates(structure), 3, length(structure.atoms), 1), 1, 1, 100),
)

suite = BenchmarkGroup()
suite["structure_loading"] = @benchmarkable load_structure($fixture)
suite["binding_site"] =
    @benchmarkable identify_binding_site($structure, $ligand; cutoff = 5.0)
suite["atomic_contacts"] = @benchmarkable atomic_contacts($structure, $ligand; cutoff = 4.0)
suite["kabsch"] = @benchmarkable kabsch($reference, $reference)
suite["trajectory_rmsd_100_frames"] = @benchmarkable trajectory_rmsd($trajectory)
suite["minimum_distance_100_frames"] =
    @benchmarkable minimum_distance_series($structure_trajectory, 1:5, 10:12)
suite["residue_rmsf_100_frames"] =
    @benchmarkable residue_rmsf($structure, $structure_trajectory)
suite["residue_contact_persistence_100_frames"] =
    @benchmarkable residue_contact_persistence(
        $structure,
        $structure_trajectory,
        $(ligand.atom_indices);
        cutoff = 5.0,
    )
suite["binding_site_comparison"] =
    @benchmarkable compare_binding_sites($structure, $structure, $ligand, $ligand)

results = run(suite; seconds = 2)
show(stdout, MIME("text/plain"), results)
println()
