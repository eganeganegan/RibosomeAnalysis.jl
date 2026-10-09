using BenchmarkTools
using RibosomeAnalysis

fixture = joinpath(@__DIR__, "..", "test", "fixtures", "synthetic.pdb")
structure = load_structure(fixture)
ligand = select_ligand(structure; residue_name="PAR")
reference = [0.0 1 0 0; 0 0 1 0; 0 0 0 1]
trajectory = TrajectoryData(repeat(reshape(reference, 3, 4, 1), 1, 1, 100))

suite = BenchmarkGroup()
suite["structure_loading"] = @benchmarkable load_structure($fixture)
suite["binding_site"] = @benchmarkable identify_binding_site($structure, $ligand; cutoff=5.0)
suite["atomic_contacts"] = @benchmarkable atomic_contacts($structure, $ligand; cutoff=4.0)
suite["kabsch"] = @benchmarkable kabsch($reference, $reference)
suite["trajectory_rmsd_100_frames"] = @benchmarkable trajectory_rmsd($trajectory)

results = run(suite; seconds=2)
show(stdout, MIME("text/plain"), results)
println()
