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
