#!/usr/bin/env julia

include(joinpath(@__DIR__, "..", "src", "RodentDMDAnalysis.jl"))

using .RodentDMDAnalysis
using CairoMakie
using CSV
using DataFrames
using DMDAnalysis
using LinearAlgebra
using ProgressBars
using Random
using Serialization
using Statistics

const ROOT = normpath(joinpath(@__DIR__, ".."))
const RESULTS_DIR = joinpath(ROOT, "results")
const FIG_DIR = joinpath(ROOT, "figures")
const GENOTYPE_GROUPS = ["Cntnap_neg", "Cntnap_het", "Cntnap_hom", "TSC_neg", "TSC_het", "TSC_hom"]
const CURVE_GROUPS = ["Cntnap_neg", "Cntnap_het", "Cntnap_hom", "TSC_neg", "TSC_het", "TSC_hom", "new_preds"]
const IO_GROUPS = CURVE_GROUPS
const GROUP_LABELS = Dict(
    "Cntnap_neg" => "Cntnap control",
    "Cntnap_het" => "Cntnap het",
    "Cntnap_hom" => "Cntnap hom/KO",
    "TSC_neg" => "TSC control",
    "TSC_het" => "TSC het",
    "TSC_hom" => "TSC hom/mut",
)
const GROUP_COLORS = Dict(
    "Cntnap_neg" => RGBf(0.00, 0.45, 0.70),
    "Cntnap_het" => RGBf(0.35, 0.70, 0.90),
    "Cntnap_hom" => RGBf(0.00, 0.62, 0.45),
    "TSC_neg" => RGBf(0.80, 0.47, 0.65),
    "TSC_het" => RGBf(0.90, 0.62, 0.00),
    "TSC_hom" => RGBf(0.84, 0.37, 0.00),
)

label_for(g) = get(GROUP_LABELS, String(g), String(g))
color_for(g) = get(GROUP_COLORS, String(g), RGBf(0.3, 0.3, 0.3))

function finite_mean(x)
    vals = Float64.(collect(skipmissing(x)))
    vals = vals[isfinite.(vals)]
    isempty(vals) ? NaN : mean(vals)
end

function io_matrix(row)
    n = Int(row.n_features)
    M = fill(NaN, n, n)
    for i in 1:n, j in 1:n
        col = Symbol("io_", i, "_from_", j)
        hasproperty(row, col) && (M[i, j] = Float64(row[col]))
    end
    mx = maximum(M)
    isfinite(mx) && mx > 0 || return M
    return M ./ mx
end

function blockmean(M, outs, ins)
    vals = vec(M[outs, ins])
    vals = vals[isfinite.(vals)]
    isempty(vals) ? NaN : mean(vals)
end

function limb_blocks(n_features)
    if n_features == 4
        return (
            fore=1:2,
            hind=3:4,
            left=[1, 4],
            right=[2, 3],
        )
    elseif n_features == 8
        return (
            fore=1:4,
            hind=5:8,
            left=[1, 2, 7, 8],
            right=[3, 4, 5, 6],
        )
    else
        error("Cannot assign fore/hind and left/right IO blocks for $n_features state features")
    end
end

function geometry_blocks(n_features)
    if n_features == 8
        return (
            radius=[1, 3, 5, 7],
            angle=[2, 4, 6, 8],
        )
    elseif n_features == 4
        return (
            radius=Int[],
            angle=collect(1:4),
        )
    else
        error("Cannot assign radius/angle IO blocks for $n_features state features")
    end
end

function subspace_coupling_labels(state_space::AbstractString)
    if state_space == "centroid_geometry"
        return (
            selectivity_title="Angle subspace selectivity",
            selectivity_ylabel="Angle-angle / angle-radius IO",
            directional_title="Angle/radius input-output",
            within_label="angle -> angle",
            cross_label="radius -> angle",
        )
    elseif state_space == "body_aligned_xy" || state_space == "centroid_xy" || state_space == "raw_xy"
        return (
            selectivity_title="Y subspace selectivity",
            selectivity_ylabel="Y-Y / Y-X IO",
            directional_title="X/Y input-output",
            within_label="y -> y",
            cross_label="x -> y",
        )
    elseif state_space == "body_ap"
        return (
            selectivity_title="AP self-selectivity",
            selectivity_ylabel="Self / cross-paw IO",
            directional_title="AP input-output",
            within_label="same paw",
            cross_label="other paw",
        )
    else
        return (
            selectivity_title="Subspace selectivity",
            selectivity_ylabel="within / cross IO",
            directional_title="Subspace input-output",
            within_label="within subspace",
            cross_label="cross subspace",
        )
    end
end

function io_summary_table(df)
    rows = NamedTuple[]
    for row in eachrow(df)
        M = io_matrix(row)
        n = size(M, 1)
        blocks = limb_blocks(n)
        geom = geometry_blocks(n)
        diagvals = diag(M)
        offvals = [M[i, j] for i in 1:n, j in 1:n if i != j]
        diag_gain = finite_mean(diagvals)
        off_gain = finite_mean(offvals)
        fore_hind = mean([
            blockmean(M, blocks.fore, blocks.hind),
            blockmean(M, blocks.hind, blocks.fore),
        ])
        left_right = mean([
            blockmean(M, blocks.left, blocks.right),
            blockmean(M, blocks.right, blocks.left),
        ])
        state_space = hasproperty(row, :state_space) ? String(row.state_space) : "unknown"
        if state_space == "body_ap"
            angle_angle = diag_gain
            radius_radius = NaN
            radius_to_angle = off_gain
            angle_to_radius = NaN
            angle_radius_cross = off_gain
        else
            angle_angle = isempty(geom.angle) ? NaN : blockmean(M, geom.angle, geom.angle)
            radius_radius = isempty(geom.radius) ? NaN : blockmean(M, geom.radius, geom.radius)
            radius_to_angle = isempty(geom.radius) ? NaN : blockmean(M, geom.angle, geom.radius)
            angle_to_radius = isempty(geom.radius) ? NaN : blockmean(M, geom.radius, geom.angle)
            angle_radius_cross = mean([radius_to_angle, angle_to_radius])
        end
        push!(rows, (
            group=row.group,
            subject=row.subject,
            speed_bin=row.speed_bin,
            mean_speed=row.mean_speed,
            model_index=row.model_index,
            diagonal_dominance=diag_gain / off_gain,
            self_gain=diag_gain,
            cross_gain=off_gain,
            fore_hind_coupling=fore_hind,
            left_right_coupling=left_right,
            angle_angle_coupling=angle_angle,
            radius_radius_coupling=radius_radius,
            radius_to_angle_coupling=radius_to_angle,
            angle_to_radius_coupling=angle_to_radius,
            angle_radius_cross_coupling=angle_radius_cross,
            angle_selectivity=angle_angle / angle_radius_cross,
        ))
    end
    out = DataFrame(rows)
    CSV.write(joinpath(RESULTS_DIR, "rodent_io_summary_metrics.csv"), out)
    return out
end

function hierarchical_bootstrap(group_df, metric; nboot=1000, seed=20260724)
    rng = MersenneTwister(seed)
    subjects = unique(group_df.subject)
    finite_all = Float64.(collect(skipmissing(group_df[!, metric])))
    finite_all = finite_all[isfinite.(finite_all)]
    isempty(finite_all) && return (mean=NaN, sem=NaN, lo=NaN, hi=NaN)
    vals = Float64[]
    for _ in 1:nboot
        sampled_subjects = rand(rng, subjects, length(subjects))
        sample_vals = Float64[]
        for subject in sampled_subjects
            sub = Float64.(collect(skipmissing(group_df[group_df.subject .== subject, metric])))
            sub = sub[isfinite.(sub)]
            isempty(sub) && continue
            push!(sample_vals, rand(rng, sub))
        end
        isempty(sample_vals) || push!(vals, mean(sample_vals))
    end
    isempty(vals) && return (mean=mean(finite_all), sem=NaN, lo=NaN, hi=NaN)
    return (mean=mean(finite_all), sem=std(vals), lo=quantile(vals, 0.025), hi=quantile(vals, 0.975))
end

function bootstrap_summary(io_df; nboot=1000)
    rows = NamedTuple[]
    for g in IO_GROUPS
        sub = io_df[io_df.group .== g, :]
        nrow(sub) == 0 && continue
        for metric in [:diagonal_dominance, :fore_hind_coupling, :left_right_coupling,
                       :cross_gain, :angle_angle_coupling, :radius_radius_coupling,
                       :radius_to_angle_coupling, :angle_to_radius_coupling,
                       :angle_radius_cross_coupling, :angle_selectivity]
            b = hierarchical_bootstrap(sub, metric; nboot=nboot)
            push!(rows, (
                group=g,
                metric=String(metric),
                n_models=nrow(sub),
                n_subjects=length(unique(sub.subject)),
                mean=b.mean,
                sem=b.sem,
                ci_low=b.lo,
                ci_high=b.hi,
            ))
        end
    end
    out = DataFrame(rows)
    CSV.write(joinpath(RESULTS_DIR, "rodent_io_summary_bootstrap.csv"), out)
    return out
end

function model_stride_frequency(df, i)
    f0 = hasproperty(df, :stride_frequency) ? Float64(df.stride_frequency[i]) : NaN
    isfinite(f0) && f0 > 0 && return f0
    vals = Float64.(df.stride_frequency)
    vals = vals[isfinite.(vals) .& (vals .> 0)]
    isempty(vals) && error("Cannot compute stride-normalized resolvent curves without finite stride frequencies")
    return median(vals)
end

function maybe_load_resolvent_curves(df; force=false,
                                     stride_multiples=collect(range(0.25, 5.0; length=90)))
    cache = joinpath(RESULTS_DIR, "rodent_resolvent_curves.csv")
    if !force && isfile(cache)
        curves = CSV.read(cache, DataFrame)
        if hasproperty(curves, :stride_multiple) &&
           length(unique(curves.stride_multiple)) == length(stride_multiples)
            return curves
        end
    end
    results = deserialize(joinpath(RESULTS_DIR, "rodent_dmd_results.jls"))
    length(results) == nrow(df) || error("DMD result count $(length(results)) does not match metrics rows $(nrow(df))")
    rows = NamedTuple[]
    println("Computing stride-normalized resolvent curves for $(nrow(df)) models; caching to $cache")
    for i in ProgressBar(1:nrow(df))
        result = results[i]
        f_stride = model_stride_frequency(df, i)
        freqs = stride_multiples .* f_stride
        gains = DMDAnalysis.resolvent_gain_at(result, freqs; fs=RodentDMDAnalysis.FS)
        for (m, f, g) in zip(stride_multiples, freqs, gains)
            push!(rows, (
                group=df.group[i],
                subject=df.subject[i],
                speed_bin=df.speed_bin[i],
                model_index=df.model_index[i],
                stride_frequency=f_stride,
                stride_multiple=m,
                frequency=f,
                log_gain=log10(max(g, eps(Float64))),
            ))
        end
    end
    curves = DataFrame(rows)
    CSV.write(cache, curves)
    return curves
end

function curve_subject_summary(curves)
    rows = NamedTuple[]
    for sub in groupby(curves, [:group, :subject, :stride_multiple])
        push!(rows, (
            group=sub.group[1],
            subject=sub.subject[1],
            stride_multiple=sub.stride_multiple[1],
            log_gain=mean(sub.log_gain),
        ))
    end
    return DataFrame(rows)
end

function curve_group_summary(curves)
    subj = curve_subject_summary(curves)
    rows = NamedTuple[]
    for sub in groupby(subj, [:group, :stride_multiple])
        vals = Float64.(sub.log_gain)
        push!(rows, (
            group=sub.group[1],
            stride_multiple=sub.stride_multiple[1],
            mean=mean(vals),
            sem=length(vals) > 1 ? std(vals) / sqrt(length(vals)) : 0.0,
            n_subjects=length(unique(sub.subject)),
        ))
    end
    return DataFrame(rows)
end

function resolvent_summary_table(curves)
    rows = NamedTuple[]
    for sub in groupby(curves, [:group, :subject, :speed_bin, :model_index])
        vals = Float64.(sub.log_gain)
        freqs = Float64.(sub.frequency)
        multiples = Float64.(sub.stride_multiple)
        peak_i = argmax(vals)
        band = (multiples .>= 0.75) .& (multiples .<= 2.0)
        push!(rows, (
            group=sub.group[1],
            subject=sub.subject[1],
            speed_bin=sub.speed_bin[1],
            model_index=sub.model_index[1],
            peak_log_gain=vals[peak_i],
            peak_frequency=freqs[peak_i],
            peak_stride_multiple=multiples[peak_i],
            mean_log_gain=mean(vals),
            band_0p75_2x_log_gain=any(band) ? mean(vals[band]) : NaN,
        ))
    end
    out = DataFrame(rows)
    CSV.write(joinpath(RESULTS_DIR, "rodent_resolvent_summary_metrics.csv"), out)
    return out
end

function bootstrap_resolvent_summary(res_df; nboot=1000)
    rows = NamedTuple[]
    for g in CURVE_GROUPS
        sub = res_df[res_df.group .== g, :]
        nrow(sub) == 0 && continue
        for metric in [:peak_log_gain, :peak_frequency, :peak_stride_multiple,
                       :mean_log_gain, :band_0p75_2x_log_gain]
            b = hierarchical_bootstrap(sub, metric; nboot=nboot)
            push!(rows, (
                group=g,
                metric=String(metric),
                n_models=nrow(sub),
                n_subjects=length(unique(sub.subject)),
                mean=b.mean,
                sem=b.sem,
                ci_low=b.lo,
                ci_high=b.hi,
            ))
        end
    end
    out = DataFrame(rows)
    CSV.write(joinpath(RESULTS_DIR, "rodent_resolvent_summary_bootstrap.csv"), out)
    return out
end

function panel_label!(fig, pos, label)
    Label(fig[pos..., TopLeft()], label; fontsize=24, font=:bold,
        padding=(0, 8, 8, 0), halign=:right, tellwidth=false, tellheight=false)
end

function grouped_bar_panel!(ax, boot, metric; ylabel)
    ax.ylabel = ylabel
    ax.xticks = (1:3, ["Cntnap", "TSC", "new_preds"])
    ax.xlabel = ""
    offsets = [-0.22, 0.0, 0.22]
    widths = 0.16
    for (line_i, line) in enumerate(["Cntnap", "TSC"])
        groups = ["$(line)_neg", "$(line)_het", "$(line)_hom"]
        for (j, g) in enumerate(groups)
            row = boot[(boot.group .== g) .& (boot.metric .== String(metric)), :]
            nrow(row) == 0 && continue
            x = line_i + offsets[j]
            y = row.mean[1]
            barplot!(ax, [x], [y]; width=widths, color=color_for(g))
            lines!(ax, [x, x], [y - row.sem[1], y + row.sem[1]]; color=:black, linewidth=1.4)
        end
    end
    row = boot[(boot.group .== "new_preds") .& (boot.metric .== String(metric)), :]
    if nrow(row) > 0
        y = row.mean[1]
        x = 3.0
        barplot!(ax, [x], [y]; width=0.18, color=color_for("new_preds"))
        lines!(ax, [x, x], [y - row.sem[1], y + row.sem[1]]; color=:black, linewidth=1.4)
    end
    xlims!(ax, 0.55, 3.45)
end

function directional_angle_radius_panel!(ax, boot; labels=subspace_coupling_labels("centroid_geometry"))
    ax.ylabel = "Normalized IO gain"
    ax.xticks = (1:3, ["Cntnap", "TSC", "new_preds"])
    ax.xlabel = ""
    offsets = [-0.22, 0.0, 0.22]
    widths = 0.10
    metrics = [
        (:angle_angle_coupling, -0.055, :solid),
        (:radius_to_angle_coupling, 0.055, :open),
    ]
    for (line_i, line) in enumerate(["Cntnap", "TSC"])
        groups = ["$(line)_neg", "$(line)_het", "$(line)_hom"]
        for (j, g) in enumerate(groups)
            base_x = line_i + offsets[j]
            for (metric, dx, style) in metrics
                row = boot[(boot.group .== g) .& (boot.metric .== String(metric)), :]
                nrow(row) == 0 && continue
                x = base_x + dx
                y = row.mean[1]
                c = color_for(g)
                if style == :solid
                    barplot!(ax, [x], [y]; width=widths, color=c,
                        strokecolor=c, strokewidth=0.0)
                else
                    barplot!(ax, [x], [y]; width=widths, color=(c, 0.18),
                        strokecolor=c, strokewidth=2.0)
                end
                lines!(ax, [x, x], [y - row.sem[1], y + row.sem[1]]; color=:black, linewidth=1.2)
            end
        end
    end
    for (metric, dx, style) in metrics
        row = boot[(boot.group .== "new_preds") .& (boot.metric .== String(metric)), :]
        nrow(row) == 0 && continue
        x = 3.0 + dx
        y = row.mean[1]
        c = color_for("new_preds")
        if style == :solid
            barplot!(ax, [x], [y]; width=widths, color=c,
                strokecolor=c, strokewidth=0.0)
        else
            barplot!(ax, [x], [y]; width=widths, color=(c, 0.18),
                strokecolor=c, strokewidth=2.0)
        end
        lines!(ax, [x, x], [y - row.sem[1], y + row.sem[1]]; color=:black, linewidth=1.2)
    end
    xlims!(ax, 0.55, 3.45)
    axislegend(ax,
        [PolyElement(color=:gray55, strokecolor=:gray55, strokewidth=0),
         PolyElement(color=(:gray55, 0.18), strokecolor=:gray35, strokewidth=2)],
        [labels.within_label, labels.cross_label];
        framevisible=false,
        position=:rt)
end

function resolvent_panel!(ax, curve_summary, groups; title, reference_group=nothing)
    ax.title = title
    ax.xlabel = "Forcing frequency (f / f_stride)"
    ax.ylabel = "log10 resolvent gain"
    if reference_group !== nothing
        ref = curve_summary[curve_summary.group .== reference_group, :]
        if nrow(ref) > 0
            sort!(ref, :stride_multiple)
            c = color_for(reference_group)
            band!(ax, ref.stride_multiple, ref.mean .- ref.sem, ref.mean .+ ref.sem;
                color=(c, 0.10))
            lines!(ax, ref.stride_multiple, ref.mean; color=c, linewidth=2.4,
                linestyle=:dash, label=label_for(reference_group))
        end
    end
    for g in groups
        sub = curve_summary[curve_summary.group .== g, :]
        sort!(sub, :stride_multiple)
        c = color_for(g)
        band!(ax, sub.stride_multiple, sub.mean .- sub.sem, sub.mean .+ sub.sem;
            color=(c, 0.18))
        lines!(ax, sub.stride_multiple, sub.mean; color=c, linewidth=2.6, label=label_for(g))
    end
    vlines!(ax, [1, 2, 3, 4]; color=(:black, 0.18), linewidth=1.0, linestyle=:dash)
end

function make_figure3(; nboot=1000, force_curves=false)
    RodentDMDAnalysis.pub_theme!()
    df = RodentDMDAnalysis.load_metrics()
    state_space = RodentDMDAnalysis.metric_state_space(df)
    subspace_labels = subspace_coupling_labels(state_space)
    io_df = io_summary_table(df)
    boot = bootstrap_summary(io_df; nboot=nboot)
    curves = maybe_load_resolvent_curves(df; force=force_curves)
    res_df = resolvent_summary_table(curves)
    bootstrap_resolvent_summary(res_df; nboot=nboot)
    curve_summary = curve_group_summary(curves)

    fig = Figure(size=(1180, 1280))
    axA = Axis(fig[1, 1], title="Self-dominance")
    panel_label!(fig, (1, 1), "A")
    grouped_bar_panel!(axA, boot, :diagonal_dominance; ylabel="Diagonal / off-diagonal IO")

    axB = Axis(fig[1, 2], title=subspace_labels.selectivity_title)
    panel_label!(fig, (1, 2), "B")
    grouped_bar_panel!(axB, boot, :angle_selectivity; ylabel=subspace_labels.selectivity_ylabel)

    axC = Axis(fig[2, 1], title=subspace_labels.directional_title)
    panel_label!(fig, (2, 1), "C")
    directional_angle_radius_panel!(axC, boot; labels=subspace_labels)

    axD = Axis(fig[2, 2], title="Fore-hind coupling")
    panel_label!(fig, (2, 2), "D")
    grouped_bar_panel!(axD, boot, :fore_hind_coupling; ylabel="Normalized IO gain")

    axE = Axis(fig[3, 1], title="Left-right coupling")
    panel_label!(fig, (3, 1), "E")
    grouped_bar_panel!(axE, boot, :left_right_coupling; ylabel="Normalized IO gain")

    axF = Axis(fig[3, 2], title="Distributed IO")
    panel_label!(fig, (3, 2), "F")
    grouped_bar_panel!(axF, boot, :cross_gain; ylabel="Mean off-diagonal IO")

    axG = Axis(fig[4, 1])
    panel_label!(fig, (4, 1), "G")
    resolvent_panel!(axG, curve_summary, ["Cntnap_neg", "Cntnap_het", "Cntnap_hom"];
        title="Cntnap resolvent gain", reference_group="new_preds")
    axislegend(axG; position=:rt, framevisible=false)

    axH = Axis(fig[4, 2])
    panel_label!(fig, (4, 2), "H")
    resolvent_panel!(axH, curve_summary, ["TSC_neg", "TSC_het", "TSC_hom"];
        title="TSC resolvent gain", reference_group="new_preds")
    axislegend(axH; position=:rt, framevisible=false)

    linkyaxes!(axG, axH)
    save(joinpath(FIG_DIR, "Figure 3.svg"), fig)
    println("Saved: ", joinpath(FIG_DIR, "Figure 3.svg"))
    return fig
end

function main(args=ARGS)
    force_curves = "--force-curves" in args
    nboot = begin
        i = findfirst(==("--nboot"), args)
        isnothing(i) ? 1000 : parse(Int, args[i + 1])
    end
    make_figure3(; nboot=nboot, force_curves=force_curves)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
