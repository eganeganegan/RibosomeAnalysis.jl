module RibosomeAnalysisCairoMakieExt

using CairoMakie
import RibosomeAnalysis: plot_contact_heatmap, plot_pca, plot_rmsd, plot_rmsf

function _save_figure(figure, path)
    extension = lowercase(splitext(path)[2])
    extension in (".png", ".svg", ".pdf") ||
        throw(ArgumentError("figure path must end in .png, .svg, or .pdf"))
    mkpath(dirname(abspath(path)))
    CairoMakie.save(path, figure)
    path
end

function plot_rmsd(
    times_ps,
    values_angstrom,
    path::AbstractString;
    title = "Trajectory RMSD",
)
    figure = Figure(size = (900, 550))
    axis = Axis(figure[1, 1], xlabel = "Time (ps)", ylabel = "RMSD (Å)", title = title)
    lines!(axis, times_ps, values_angstrom, linewidth = 2, label = "RMSD")
    axislegend(axis)
    _save_figure(figure, path)
end

function plot_rmsf(
    residue_or_atom_ids,
    values_angstrom,
    path::AbstractString;
    title = "Structural fluctuations",
)
    figure = Figure(size = (900, 550))
    axis = Axis(
        figure[1, 1],
        xlabel = "Atom or residue index",
        ylabel = "RMSF (Å)",
        title = title,
    )
    lines!(axis, residue_or_atom_ids, values_angstrom, linewidth = 2, label = "RMSF")
    axislegend(axis)
    _save_figure(figure, path)
end

function plot_contact_heatmap(
    matrix,
    path::AbstractString;
    xlabels = nothing,
    ylabels = nothing,
    title = "Contact persistence",
)
    figure = Figure(size = (800, 650))
    axis =
        Axis(figure[1, 1], xlabel = "Group B atom", ylabel = "Group A atom", title = title)
    isnothing(xlabels) || (axis.xticks = (collect(axes(matrix, 2)), xlabels))
    isnothing(ylabels) || (axis.yticks = (collect(axes(matrix, 1)), ylabels))
    heat = heatmap!(axis, matrix; colorrange = (0, 1))
    Colorbar(figure[1, 2], heat, label = "Fraction of sampled frames")
    _save_figure(figure, path)
end

function plot_pca(
    scores,
    path::AbstractString;
    components = (1, 2),
    title = "Trajectory PCA",
)
    i, j = components
    figure = Figure(size = (750, 650))
    axis = Axis(figure[1, 1], xlabel = "PC$i score", ylabel = "PC$j score", title = title)
    scatter!(axis, scores[i, :], scores[j, :], color = axes(scores, 2), colormap = :viridis)
    _save_figure(figure, path)
end

end
