using RibosomeAnalysis

reference = [0.0 1 0 0; 0 0 1 0; 0 0 0 1]
frames =
    cat(reference, reference .+ [0.1, 0.0, 0.0], reference .+ [0.2, 0.0, 0.0]; dims = 3)
trajectory = TrajectoryData(frames; times_ps = [0.0, 10.0, 20.0])
display(trajectory_rmsd(trajectory))
display(trajectory_rmsf(trajectory))
