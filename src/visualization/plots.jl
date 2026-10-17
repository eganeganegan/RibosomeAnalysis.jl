function _plotting_unavailable(name::Symbol)
    throw(
        ArgumentError(
            "$name requires the optional CairoMakie extension; install CairoMakie and load it with `using CairoMakie`",
        ),
    )
end

"""
    plot_rmsd(times_ps, values_angstrom, path; title="Trajectory RMSD")

Write an RMSD line plot to a PNG, SVG, or PDF file. Install and load
`CairoMakie` before calling this optional plotting function.
"""
plot_rmsd(args...; kwargs...) = _plotting_unavailable(:plot_rmsd)

"""
    plot_rmsf(residue_or_atom_ids, values_angstrom, path;
              title="Structural fluctuations")

Write an RMSF line plot to a PNG, SVG, or PDF file. Install and load
`CairoMakie` before calling this optional plotting function.
"""
plot_rmsf(args...; kwargs...) = _plotting_unavailable(:plot_rmsf)

"""
    plot_contact_heatmap(matrix, path; xlabels=nothing, ylabels=nothing,
                         title="Contact persistence")

Write a contact-persistence heatmap to a PNG, SVG, or PDF file. Values are
displayed on a fixed zero-to-one color scale. Install and load `CairoMakie`
before calling this optional plotting function.
"""
plot_contact_heatmap(args...; kwargs...) = _plotting_unavailable(:plot_contact_heatmap)

"""
    plot_pca(scores, path; components=(1, 2), title="Trajectory PCA")

Write a scatter plot of two PCA score rows to a PNG, SVG, or PDF file. Install
and load `CairoMakie` before calling this optional plotting function.
"""
plot_pca(args...; kwargs...) = _plotting_unavailable(:plot_pca)

"""
    plot_distance_distribution(values_angstrom, path; bins=30,
                               title="Distance distribution")

Write a distance histogram with Angstrom units to a PNG, SVG, or PDF file.
Install and load `CairoMakie` before calling this optional plotting function.
"""
plot_distance_distribution(args...; kwargs...) =
    _plotting_unavailable(:plot_distance_distribution)

"""
    plot_cluster_populations(labels, path; title="Cluster populations")

Write frame counts per cluster to a PNG, SVG, or PDF file. Install and load
`CairoMakie` before calling this optional plotting function.
"""
plot_cluster_populations(args...; kwargs...) =
    _plotting_unavailable(:plot_cluster_populations)

"""
    plot_binding_site_comparison(comparison, path;
                                 title="Binding-site comparison")

Plot reference and target binding-site residue and atomic-contact counts from a
`ComparisonResult`. Install and load `CairoMakie` before calling this optional
plotting function.
"""
plot_binding_site_comparison(args...; kwargs...) =
    _plotting_unavailable(:plot_binding_site_comparison)

"""
    plot_state_metric(states, metric, path; title=nothing, ylabel=nothing)

Plot one numeric column from a state-comparison table against its `label`
column. Rows with missing values are omitted. Install and load `CairoMakie`
before calling this optional plotting function.
"""
plot_state_metric(args...; kwargs...) = _plotting_unavailable(:plot_state_metric)

"""
    plot_state_ensemble(states, path; kwargs...)

Plot equally weighted state-level replica estimates with bootstrap confidence
intervals when available. Install and load `CairoMakie` before calling this
optional plotting function.
"""
plot_state_ensemble(args...; kwargs...) = _plotting_unavailable(:plot_state_ensemble)
