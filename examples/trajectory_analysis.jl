using RibosomeAnalysis

reference = [0.0 1 0 0; 0 0 1 0; 0 0 0 1]
translated = reference .+ [0.1, 0.0, 0.0]
deformed = reference .+ [0.2, 0.0, 0.0]
deformed[:, 4] .+= [0.1, 0.0, 0.0]
frames = cat(reference, translated, deformed; dims = 3)
trajectory = TrajectoryData(frames; times_ps = [0.0, 10.0, 20.0])
display(trajectory_rmsd(trajectory))
display(trajectory_rmsf(trajectory))
display(minimum_distance_series(trajectory, [1, 2], [3, 4]))
display(group_centroid_displacement(trajectory, [4]; align_indices = [1, 2, 3]))
