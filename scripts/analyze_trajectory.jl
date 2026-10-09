include(joinpath(@__DIR__, "common.jl"))

function main(arguments=ARGS)
    options = parse_cli(arguments)
    topology = require_option(options, "topology")
    trajectory_path = require_option(options, "trajectory")
    output = require_option(options, "output")
    selection = get(options, "selection", "all")
    trajectory = load_trajectory(topology, trajectory_path; selection=selection)
    values = trajectory_rmsd(trajectory)
    fluctuations = trajectory_rmsf(trajectory)
    directory = prepare_results(output)
    write_table(directory, "rmsd.csv", DataFrame(time_ps=trajectory.times_ps, rmsd_angstrom=values))
    write_table(directory, "rmsf.csv", DataFrame(atom_index=trajectory.atom_indices, rmsf_angstrom=fluctuations))
    plot_rmsd(trajectory.times_ps, values, joinpath(directory, "figures", "rmsd.pdf"))
    plot_rmsf(trajectory.atom_indices, fluctuations, joinpath(directory, "figures", "rmsf.pdf"))
    write_summary(directory, Dict("frames" => length(trajectory.times_ps), "atoms" => length(trajectory.atom_indices),
                                  "rmsd" => summarize_samples(values; correlated=true),
                                  "unit_conversion_to_angstrom" => trajectory.unit_conversion_to_angstrom))
    write_metadata(directory; inputs=[topology, trajectory_path], parameters=Dict("workflow" => "trajectory_analysis", "selection" => selection))
    println("Wrote trajectory analysis to $directory")
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()
