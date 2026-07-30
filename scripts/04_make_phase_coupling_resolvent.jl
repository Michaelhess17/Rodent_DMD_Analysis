#!/usr/bin/env julia

include(joinpath(@__DIR__, "..", "src", "RodentDMDAnalysis.jl"))

using .RodentDMDAnalysis
using CairoMakie
using CSV
using DataFrames
using LinearAlgebra
using ProgressBars
using Serialization
using Statistics

const ROOT = normpath(joinpath(@__DIR__, ".."))
const RESULTS_DIR = joinpath(ROOT, "results")
const FIG_DIR = joinpath(ROOT, "figures")
const PAWS = ["LF", "RF", "RH", "LH"]
const GENOTYPE_GROUPS = ["Cntnap_neg", "Cntnap_het", "Cntnap_hom", "TSC_neg", "TSC_het", "TSC_hom"]
const GROUP_LABELS = RodentDMDAnalysis.GROUP_LABELS
const GROUP_COLORS = RodentDMDAnalysis.GROUP_COLORS

label_for(g) = get(GROUP_LABELS, String(g), String(g))
color_for(g) = get(GROUP_COLORS, String(g), RGBf(0.3, 0.3, 0.3))
finite_mean(v) = (x = RodentDMDAnalysis.finite_vec(v); isempty(x) ? NaN : mean(x))
finite_sem(v) = (x = RodentDMDAnalysis.finite_vec(v); length(x) <= 1 ? 0.0 : std(x) / sqrt(length(x)))

function wrap_stride_fraction(x)
    y = mod(x + 0.5, 1.0) - 0.5
    return isapprox(y, -0.5; atol=1e-12) ? 0.5 : y
end

function state_space_from_cache(df, cfg)
    if hasproperty(df, :state_space) && nrow(df) > 0 && df.state_space[1] != "unknown"
        return String(df.state_space[1])
    end
    hasproperty(cfg, :state_space) && return String(cfg.state_space)
    return "unknown"
end

function limb_feature_groups(state_space::AbstractString, n_features::Int)
    if n_features == 4
        return [[1], [2], [3], [4]]
    elseif n_features == 8
        # All current 8D spaces are packed by paw: LF variables, RF variables, RH variables, LH variables.
        return [[1, 2], [3, 4], [5, 6], [7, 8]]
    end
    error("Cannot map $n_features features from state_space=$state_space onto four limbs")
end

function state_component_labels(state_space::AbstractString, n_features::Int)
    labels = RodentDMDAnalysis.state_labels(state_space)
    length(labels) == n_features && return labels
    return ["x$i" for i in 1:n_features]
end

function stride_frequency(row)
    if hasproperty(row, :stride_frequency)
        f0 = Float64(row.stride_frequency)
        isfinite(f0) && f0 > 0 && return f0
    end
    return 1.0
end

function phase_resolvent_matrix(result, frequency; fs=RodentDMDAnalysis.FS)
    A, U = result.A, result.U
    z = exp(2im * pi * frequency / fs)
    R = inv(z * I(size(A, 1)) - A)
    H = abs.(U * R * U')
    mx = maximum(H)
    isfinite(mx) && mx > 0 || return fill(NaN, size(H))
    return Float64.(H ./ mx)
end

hankel_index(feature::Int, lag::Int, n_features::Int) = feature + lag * n_features

function directed_limb_lag_profiles(H, result, limb_groups, phase_steps::Int)
    n_features = result.n_features
    n_delays = result.hankel.n_delays
    delay_interval = result.hankel.delay_interval
    deltas = collect(-(n_delays - 1):(n_delays - 1))
    rows = NamedTuple[]
    for jout in 1:4, jin in 1:4, delta in deltas
        vals = Float64[]
        for in_lag in 0:n_delays-1
            out_lag = in_lag + delta
            0 <= out_lag <= n_delays - 1 || continue
            for fout in limb_groups[jout], fin in limb_groups[jin]
                v = H[hankel_index(fout, out_lag, n_features), hankel_index(fin, in_lag, n_features)]
                isfinite(v) && push!(vals, v)
            end
        end
        push!(rows, (
            output_limb=PAWS[jout],
            input_limb=PAWS[jin],
            pair_label="$(PAWS[jout]) <- $(PAWS[jin])",
            lag_samples=delta * delay_interval,
            raw_lag_fraction=delta * delay_interval / phase_steps,
            lag_fraction=wrap_stride_fraction(delta * delay_interval / phase_steps),
            normalized_gain=isempty(vals) ? NaN : mean(vals),
            n_cells=length(vals),
        ))
    end
    return rows
end

function category_pairs(name::Symbol)
    name == :same && return [("LF", "LF"), ("RF", "RF"), ("RH", "RH"), ("LH", "LH")]
    name == :fore_homologs && return [("LF", "RF"), ("RF", "LF")]
    name == :hind_homologs && return [("LH", "RH"), ("RH", "LH")]
    name == :ipsilateral && return [("LF", "LH"), ("LH", "LF"), ("RF", "RH"), ("RH", "RF")]
    name == :diagonal && return [("LF", "RH"), ("RH", "LF"), ("RF", "LH"), ("LH", "RF")]
    name == :all_cross && return [(o, i) for o in PAWS for i in PAWS if o != i]
    error("Unknown coupling category $name")
end

function category_label(name::Symbol)
    labels = Dict(
        :same => "same limb",
        :fore_homologs => "fore homologs",
        :hind_homologs => "hind homologs",
        :ipsilateral => "ipsilateral fore-hind",
        :diagonal => "diagonal",
        :all_cross => "all cross-limb",
    )
    return labels[name]
end

function phase_vector_strength(lags, gains)
    x = Float64.(lags)
    w = Float64.(gains)
    ok = isfinite.(x) .& isfinite.(w) .& (w .> 0)
    count(ok) >= 2 || return NaN
    theta = 2pi .* x[ok]
    z = sum(w[ok] .* exp.(im .* theta)) / sum(w[ok])
    return abs(z)
end

function profile_metrics(profile_rows, row, frequency_multiple)
    df = DataFrame(profile_rows)
    pairs = Pair{Symbol,Any}[
        :group => row.group,
        :subject => row.subject,
        :speed_bin => row.speed_bin,
        :model_index => row.model_index,
        :mean_speed => row.mean_speed,
        :frequency_multiple => frequency_multiple,
        :forcing_frequency => frequency_multiple * stride_frequency(row),
    ]
    for cat in [:same, :fore_homologs, :hind_homologs, :ipsilateral, :diagonal, :all_cross]
        sub = DataFrame()
        cat_pairs = category_pairs(cat)
        mask = [(String(o), String(i)) in cat_pairs for (o, i) in zip(df.output_limb, df.input_limb)]
        sub = df[mask, :]
        bylag = combine(groupby(sub, :lag_fraction), :normalized_gain => finite_mean => :gain)
        vals = RodentDMDAnalysis.finite_vec(bylag.gain)
        if isempty(vals)
            append!(pairs, [
                Symbol(cat, "_mean_gain") => NaN,
                Symbol(cat, "_peak_gain") => NaN,
                Symbol(cat, "_peak_lag_fraction") => NaN,
                Symbol(cat, "_phase_selectivity") => NaN,
            ])
        else
            idx = argmax(bylag.gain)
            append!(pairs, [
                Symbol(cat, "_mean_gain") => mean(vals),
                Symbol(cat, "_peak_gain") => bylag.gain[idx],
                Symbol(cat, "_peak_lag_fraction") => bylag.lag_fraction[idx],
                Symbol(cat, "_phase_selectivity") => phase_vector_strength(bylag.lag_fraction, bylag.gain),
            ])
        end
    end
    same = get(Dict(pairs), :same_mean_gain, NaN)
    cross = get(Dict(pairs), :all_cross_mean_gain, NaN)
    push!(pairs, :self_to_cross_ratio => same / cross)
    return NamedTuple(pairs)
end

function subject_balanced_profiles(model_profiles)
    model_df = DataFrame(model_profiles)
    CSV.write(joinpath(RESULTS_DIR, "rodent_phase_coupling_model_profiles.csv"), model_df)

    subject_speed = combine(groupby(model_df, [:state_space, :n_features, :n_delays, :group, :subject, :speed_bin, :output_limb, :input_limb, :pair_label, :lag_samples, :lag_fraction]),
        :normalized_gain => finite_mean => :normalized_gain,
        :n_cells => sum => :n_cells)
    CSV.write(joinpath(RESULTS_DIR, "rodent_phase_coupling_subject_profiles.csv"), subject_speed)

    group_speed = combine(groupby(subject_speed, [:state_space, :n_features, :n_delays, :group, :speed_bin, :output_limb, :input_limb, :pair_label, :lag_fraction]),
        :normalized_gain => finite_mean => :mean,
        :normalized_gain => finite_sem => :sem,
        :subject => (x -> length(unique(x))) => :n_subjects)
    group_speed.level .= "group_speed"

    subject_all = combine(groupby(subject_speed, [:state_space, :n_features, :n_delays, :group, :subject, :output_limb, :input_limb, :pair_label, :lag_fraction]),
        :normalized_gain => finite_mean => :normalized_gain)
    group_all = combine(groupby(subject_all, [:state_space, :n_features, :n_delays, :group, :output_limb, :input_limb, :pair_label, :lag_fraction]),
        :normalized_gain => finite_mean => :mean,
        :normalized_gain => finite_sem => :sem,
        :subject => (x -> length(unique(x))) => :n_subjects)
    group_all.speed_bin .= "all"
    group_all.level .= "group"

    group_speed.lag_samples .= round.(Int, group_speed.lag_fraction .* 100)
    group_all.lag_samples .= round.(Int, group_all.lag_fraction .* 100)
    cols = [:level, :state_space, :n_features, :n_delays, :group, :speed_bin, :output_limb,
        :input_limb, :pair_label, :lag_samples, :lag_fraction, :mean, :sem, :n_subjects]
    out = vcat(group_speed[:, cols], group_all[:, cols])
    CSV.write(joinpath(RESULTS_DIR, "rodent_phase_lag_profiles.csv"), out)
    return model_df, subject_speed, out
end

function summarize_metrics(metrics)
    df = DataFrame(metrics)
    CSV.write(joinpath(RESULTS_DIR, "rodent_phase_coupling_model_metrics.csv"), df)
    metric_cols = names(df, r"(_mean_gain|_peak_gain|_peak_lag_fraction|_phase_selectivity|self_to_cross_ratio)$")

    subject_speed = combine(groupby(df, [:group, :subject, :speed_bin]),
        [Symbol(c) => finite_mean => Symbol(c) for c in metric_cols]...,
        :model_index => length => :n_models)
    CSV.write(joinpath(RESULTS_DIR, "rodent_phase_coupling_subject_metrics.csv"), subject_speed)

    rows = NamedTuple[]
    for sub in groupby(subject_speed, [:group, :speed_bin])
        source = df[(df.group .== sub.group[1]) .& (df.speed_bin .== sub.speed_bin[1]), :]
        for col in metric_cols
            vals = RodentDMDAnalysis.finite_vec(sub[!, Symbol(col)])
            push!(rows, (group=sub.group[1], speed_bin=sub.speed_bin[1], metric=col,
                n_models=nrow(source), n_subjects=length(unique(sub.subject)),
                mean=isempty(vals) ? NaN : mean(vals),
                sem=length(vals) > 1 ? std(vals) / sqrt(length(vals)) : 0.0))
        end
    end

    subject_all = combine(groupby(subject_speed, [:group, :subject]),
        [Symbol(c) => finite_mean => Symbol(c) for c in metric_cols]...)
    for sub in groupby(subject_all, :group)
        source = df[df.group .== sub.group[1], :]
        for col in metric_cols
            vals = RodentDMDAnalysis.finite_vec(sub[!, Symbol(col)])
            push!(rows, (group=sub.group[1], speed_bin="all", metric=col,
                n_models=nrow(source), n_subjects=length(unique(sub.subject)),
                mean=isempty(vals) ? NaN : mean(vals),
                sem=length(vals) > 1 ? std(vals) / sqrt(length(vals)) : 0.0))
        end
    end
    out = DataFrame(rows)
    CSV.write(joinpath(RESULTS_DIR, "rodent_phase_coupling_summary.csv"), out)
    return df, out
end

function make_phase_pair_heatmap(profiles; outfile="Figure 4.svg")
    RodentDMDAnalysis.pub_theme!()
    fig = Figure(size=(1360, 940))
    groups = GENOTYPE_GROUPS
    sub_all = profiles[(profiles.level .== "group") .& (profiles.speed_bin .== "all"), :]
    color_vals = RodentDMDAnalysis.finite_vec(sub_all.mean)
    crange = (0.0, isempty(color_vals) ? 1.0 : quantile(color_vals, 0.995))
    hm = nothing
    ordered_pairs = ["$o <- $i" for o in PAWS for i in PAWS]
    ymap = Dict(p => i for (i, p) in enumerate(ordered_pairs))
    for (i, g) in enumerate(groups)
        row = i <= 3 ? 1 : 2
        col = i <= 3 ? i : i - 3
        ax = Axis(fig[row, col], title=label_for(g), xlabel=row == 2 ? "output phase - input phase" : "",
            ylabel=col == 1 ? "directed limb pair" : "")
        sub = sub_all[sub_all.group .== g, :]
        sort!(sub, [:pair_label, :lag_fraction])
        xs = Float64.(sub.lag_fraction)
        ys = [ymap[String(p)] for p in sub.pair_label]
        hm = heatmap!(ax, xs, ys, Float64.(sub.mean); colormap=:magma, colorrange=crange)
        ax.yticks = (1:length(ordered_pairs), ordered_pairs)
        ax.yticklabelsize = 7
        ax.xticks = -0.5:0.25:0.5
        vlines!(ax, [0.0]; color=(:white, 0.45), linewidth=1.0, linestyle=:dash)
    end
    Colorbar(fig[:, 4], hm; label="mean normalized resolvent IO")
    save(joinpath(FIG_DIR, outfile), fig)
    println("Saved: ", joinpath(FIG_DIR, outfile))
    return fig
end

function make_category_profile_figure(profiles; outfile="Figure 4 phase lag profiles.svg")
    RodentDMDAnalysis.pub_theme!()
    cats = [:fore_homologs, :hind_homologs, :ipsilateral, :diagonal]
    line_groups = [
        ("Cntnap", ["Cntnap_neg", "Cntnap_het", "Cntnap_hom"]),
        ("TSC", ["TSC_neg", "TSC_het", "TSC_hom"]),
    ]
    fig = Figure(size=(1180, 720))
    sub_all = profiles[(profiles.level .== "group") .& (profiles.speed_bin .== "all"), :]
    for (row_i, (line_name, groups)) in enumerate(line_groups)
        for (col_i, cat) in enumerate(cats)
            ax = Axis(fig[row_i, col_i], title=row_i == 1 ? category_label(cat) : "",
                xlabel=row_i == 2 ? "output phase - input phase (stride fraction)" : "",
                ylabel=col_i == 1 ? "$line_name\nnormalized gain" : "")
            cat_pairs = category_pairs(cat)
            for g in groups
                gsub = sub_all[sub_all.group .== g, :]
                mask = [(String(o), String(i)) in cat_pairs for (o, i) in zip(gsub.output_limb, gsub.input_limb)]
                bylag = combine(groupby(gsub[mask, :], :lag_fraction), :mean => finite_mean => :gain)
                sort!(bylag, :lag_fraction)
                lines!(ax, bylag.lag_fraction, bylag.gain; color=color_for(g), linewidth=2.2, label=label_for(g))
            end
            vlines!(ax, [0.0]; color=(:black, 0.25), linewidth=1.0, linestyle=:dash)
            col_i == 4 && axislegend(ax; framevisible=false, position=:rt, labelsize=8)
        end
    end
    save(joinpath(FIG_DIR, outfile), fig)
    println("Saved: ", joinpath(FIG_DIR, outfile))
    return fig
end

function make_phase_metric_figure(summary; outfile="Figure 4 phase coupling metrics.svg")
    RodentDMDAnalysis.pub_theme!()
    metrics = [
        ("fore_homologs_peak_lag_fraction", "fore homolog peak phase"),
        ("hind_homologs_peak_lag_fraction", "hind homolog peak phase"),
        ("ipsilateral_peak_lag_fraction", "ipsilateral peak phase"),
        ("diagonal_peak_lag_fraction", "diagonal peak phase"),
        ("fore_homologs_phase_selectivity", "fore homolog phase selectivity"),
        ("ipsilateral_phase_selectivity", "ipsilateral phase selectivity"),
        ("diagonal_phase_selectivity", "diagonal phase selectivity"),
        ("self_to_cross_ratio", "self / cross gain"),
    ]
    fig = Figure(size=(1220, 860))
    for (i, (metric, title)) in enumerate(metrics)
        ax = Axis(fig[cld(i, 4), mod1(i, 4)], title=title, ylabel=occursin("lag", metric) ? "stride fraction" : "value")
        sub = summary[(summary.metric .== metric) .& (summary.speed_bin .== "all"), :]
        xs = 1:nrow(sub)
        barplot!(ax, xs, sub.mean; color=[color_for(g) for g in sub.group], width=0.72)
        for (x, y, e) in zip(xs, sub.mean, sub.sem)
            lines!(ax, [x, x], [y - e, y + e]; color=:black, linewidth=1.2)
        end
        ax.xticks = (collect(xs), [replace(label_for(g), " " => "\n") for g in sub.group])
        ax.xticklabelsize = 7
        occursin("lag", metric) && hlines!(ax, [0.0]; color=(:black, 0.2), linewidth=1.0, linestyle=:dash)
    end
    save(joinpath(FIG_DIR, outfile), fig)
    println("Saved: ", joinpath(FIG_DIR, outfile))
    return fig
end

function compute_phase_coupling(; frequency_multiple=1.0, force=false)
    mkpath(RESULTS_DIR)
    mkpath(FIG_DIR)
    df = RodentDMDAnalysis.load_metrics()
    cfg = deserialize(joinpath(RESULTS_DIR, "rodent_config.jls"))
    state_space = state_space_from_cache(df, cfg)
    profile_cache = joinpath(RESULTS_DIR, "rodent_phase_lag_profiles.csv")
    metrics_cache = joinpath(RESULTS_DIR, "rodent_phase_coupling_model_metrics.csv")
    summary_cache = joinpath(RESULTS_DIR, "rodent_phase_coupling_summary.csv")
    if !force && isfile(profile_cache) && isfile(metrics_cache) && isfile(summary_cache)
        metrics = CSV.read(metrics_cache, DataFrame)
        cache_ok = hasproperty(metrics, :frequency_multiple) &&
            all(isapprox.(Float64.(metrics.frequency_multiple), frequency_multiple; atol=1e-12))
        profiles = CSV.read(profile_cache, DataFrame)
        cache_ok &= hasproperty(profiles, :state_space) && all(String.(profiles.state_space) .== state_space)
        if cache_ok
            println("Using cached phase-coupling results for state_space=$state_space")
            return profiles, metrics, CSV.read(summary_cache, DataFrame)
        end
        println("Ignoring stale phase-coupling cache; current state_space=$state_space, frequency_multiple=$frequency_multiple")
    end
    cfg.phase_normalization == "phase" || @warn "Current rodent cache is not phase-normalized" phase_normalization=cfg.phase_normalization
    phase_steps = hasproperty(cfg, :phase_steps) ? Int(cfg.phase_steps) : 100
    results = deserialize(joinpath(RESULTS_DIR, "rodent_dmd_results.jls"))
    length(results) == nrow(df) || error("DMD result count $(length(results)) does not match metrics rows $(nrow(df))")

    model_profiles = NamedTuple[]
    metric_rows = NamedTuple[]
    println("Computing directed limb phase-lag resolvent profiles for $(nrow(df)) models at $(frequency_multiple)x stride frequency; state_space=$state_space")
    for i in ProgressBar(1:nrow(df))
        row = df[i, :]
        result = results[i]
        limb_groups = limb_feature_groups(state_space, result.n_features)
        freq = frequency_multiple * stride_frequency(row)
        H = phase_resolvent_matrix(result, freq; fs=RodentDMDAnalysis.FS)
        prof = directed_limb_lag_profiles(H, result, limb_groups, phase_steps)
        for p in prof
            push!(model_profiles, merge((group=row.group, subject=row.subject, speed_bin=row.speed_bin,
                model_index=row.model_index, mean_speed=row.mean_speed, frequency_multiple=frequency_multiple,
                forcing_frequency=freq, state_space=state_space, n_features=result.n_features,
                n_delays=result.hankel.n_delays), p))
        end
        push!(metric_rows, profile_metrics(prof, row, frequency_multiple))
    end
    _, _, profiles = subject_balanced_profiles(model_profiles)
    metric_df, summary = summarize_metrics(metric_rows)
    return profiles, metric_df, summary
end

function main(args=ARGS)
    force = "--force" in args
    multiple = begin
        i = findfirst(==("--frequency-multiple"), args)
        isnothing(i) ? 1.0 : parse(Float64, args[i + 1])
    end
    profiles, metric_df, summary = compute_phase_coupling(; frequency_multiple=multiple, force=force)
    profiles = CSV.read(joinpath(RESULTS_DIR, "rodent_phase_lag_profiles.csv"), DataFrame)
    make_phase_pair_heatmap(profiles)
    make_category_profile_figure(profiles)
    make_phase_metric_figure(summary)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
