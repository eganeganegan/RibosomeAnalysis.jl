module RibosomeAnalysisCairoMakieExt

using CairoMakie
import RibosomeAnalysis:
    plot_binding_site_comparison,
    plot_cluster_populations,
    plot_contact_heatmap,
    plot_distance_distribution,
    plot_pca,
    plot_rmsd,
    plot_rmsf,
    plot_state_metric

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

function plot_distance_distribution(
    values_angstrom,
    path::AbstractString;
    bins = 30,
    title = "Distance distribution",
)
    values = Float64.(collect(values_angstrom))
    isempty(values) && throw(ArgumentError("distance values are empty"))
    all(isfinite, values) || throw(ArgumentError("distance values must be finite"))
    bins isa Integer && bins <= 0 && throw(ArgumentError("bins must be positive"))
    figure = Figure(size = (800, 550))
    axis = Axis(figure[1, 1], xlabel = "Distance (Å)", ylabel = "Count", title = title)
    hist!(axis, values; bins = bins, color = :steelblue)
    _save_figure(figure, path)
end

function plot_cluster_populations(
    labels,
    path::AbstractString;
    title = "Cluster populations",
)
    values = Int.(collect(labels))
    isempty(values) && throw(ArgumentError("cluster labels are empty"))
    all(>(0), values) || throw(ArgumentError("cluster labels must be positive integers"))
    clusters = sort!(unique(values))
    counts = [count(==(cluster), values) for cluster in clusters]
    figure = Figure(size = (800, 550))
    axis = Axis(
        figure[1, 1],
        xlabel = "Cluster",
        ylabel = "Frames",
        title = title,
        xticks = (clusters, string.(clusters)),
    )
    barplot!(axis, clusters, counts; color = :steelblue)
    _save_figure(figure, path)
end

function plot_binding_site_comparison(
    comparison,
    path::AbstractString;
    title = "Binding-site comparison",
)
    metrics = comparison.metrics
    categories = ["Site residues", "Atomic contacts"]
    reference_values = Float64.([
        metrics["reference_binding_site_residues"],
        metrics["reference_atomic_contacts"],
    ])
    target_values = Float64.([
        metrics["target_binding_site_residues"],
        metrics["target_atomic_contacts"],
    ])
    positions = collect(1:length(categories))
    figure = Figure(size = (800, 550))
    axis = Axis(
        figure[1, 1],
        xlabel = "Measure",
        ylabel = "Count",
        title = title,
        xticks = (positions, categories),
    )
    barplot!(axis, positions .- 0.2, reference_values; width = 0.38, label = "Reference")
    barplot!(axis, positions .+ 0.2, target_values; width = 0.38, label = "Target")
    axislegend(axis)
    _save_figure(figure, path)
end

function plot_state_metric(
    states,
    metric,
    path::AbstractString;
    title = nothing,
    ylabel = nothing,
)
    column = Symbol(metric)
    :label in propertynames(states) ||
        throw(ArgumentError("states must contain a label column"))
    column in propertynames(states) ||
        throw(ArgumentError("state metric column $(repr(string(column))) was not found"))
    keep = findall(!ismissing, states[!, column])
    isempty(keep) && throw(ArgumentError("state metric contains no non-missing values"))
    values = Float64.(states[keep, column])
    all(isfinite, values) || throw(ArgumentError("state metric values must be finite"))
    labels = string.(states[keep, :label])
    positions = collect(eachindex(values))
    metric_label = replace(string(column), '_' => ' ')
    figure = Figure(size = (max(800, 110 * length(values)), 550))
    axis = Axis(
        figure[1, 1],
        xlabel = "State",
        ylabel = isnothing(ylabel) ? metric_label : string(ylabel),
        title = isnothing(title) ? "State comparison: $metric_label" : string(title),
        xticks = (positions, labels),
    )
    barplot!(axis, positions, values; color = :steelblue)
    _save_figure(figure, path)
end

end
