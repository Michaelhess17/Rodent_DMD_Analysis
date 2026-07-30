#!/usr/bin/env julia

include(joinpath(@__DIR__, "..", "src", "RodentDMDAnalysis.jl"))

using .RodentDMDAnalysis
using CairoMakie
using CSV
using DataFrames
using DMDAnalysis
using Random
using Serialization
using Statistics

const ROOT = normpath(joinpath(@__DIR__, ".."))
const RESULTS_DIR = joinpath(ROOT, "results")
const FIG_DIR = joinpath(ROOT, "figures")

const LIMB_COLORS = [
    RGBf(0.00, 0.45, 0.70),
    RGBf(0.84, 0.37, 0.00),
    RGBf(0.00, 0.62, 0.45),
    RGBf(0.80, 0.47, 0.65),
]

function arg_value(args, name, default)
    i = findfirst(==(name), args)
    isnothing(i) && return default
    i == length(args) && error("Missing value after $name")
    return args[i + 1]
end

function load_current_config()
    path = joinpath(RESULTS_DIR, "rodent_config.jls")
    isfile(path) || error("Missing rodent_config.jls. Run scripts/00_fit_rodent_dmd_metrics.jl first.")
    return deserialize(path)
end

function load_training_sets(cfg; target_rows=nothing)
    manifest = RodentDMDAnalysis.rodent_manifest(; dataset=cfg.dataset)
    if !cfg.full
        manifest = RodentDMDAnalysis.sample_manifest(manifest;
            max_per_group=cfg.max_per_group,
            max_per_subject=cfg.max_per_subject,
            seed=cfg.seed)
    end
    if target_rows !== nothing
        targets = Set((String(r.group), String(r.subject)) for r in eachrow(target_rows))
        keep = [(String(r.group), String(r.subject)) in targets for r in eachrow(manifest)]
        manifest = manifest[keep, :]
    end
    bouts, meta = RodentDMDAnalysis.load_bouts(manifest;
        normalization=cfg.normalization,
        state_space=cfg.state_space,
        smooth_window=hasproperty(cfg, :smooth_window) ? cfg.smooth_window : 0,
        phase_normalization=hasproperty(cfg, :phase_normalization) ? cfg.phase_normalization : "none",
        phase_steps=hasproperty(cfg, :phase_steps) ? cfg.phase_steps : 100,
        progressbar=false,
        verbose=true,
        log_interval=25000)
    meta.state_space = fill(cfg.state_space, nrow(meta))
    RodentDMDAnalysis.add_speed_conditions!(meta;
        nbins=cfg.speed_conditions,
        mode=cfg.speed_mode)
    return bouts, meta
end

function candidate_rows(df)
    preferred = [
        ("Cntnap_neg", "medium"),
        ("TSC_hom", "medium"),
        ("Cntnap_hom", "medium"),
        ("TSC_neg", "medium"),
    ]
    rows = Int[]
    for (group, speed) in preferred
        sub = df[(df.group .== group) .& (df.speed_bin .== speed), :]
        nrow(sub) == 0 && continue
        idx = argmax(Float64.(sub.n_est_strides))
        push!(rows, Int(sub.model_index[idx]))
        length(rows) == 2 && return rows
    end
    order = sortperm(Float64.(df.n_est_strides); rev=true)
    return Int.(df.model_index[order[1:min(2, length(order))]])
end

function trial_for_model(bouts, meta, row)
    idx = findall((meta.group .== row.group) .&
                  (meta.subject .== row.subject) .&
                  (meta.speed_condition .== row.speed_bin))
    isempty(idx) && error("No bout found for $(row.group) $(row.subject) $(row.speed_bin)")
    lens = [size(bouts[i], 1) for i in idx]
    return bouts[idx[argmax(lens)]]
end

function choose_start(trial, n_delays, horizon; rng=MersenneTwister(20260729))
    available = size(trial, 1) - n_delays - horizon - 2
    available >= 1 || return 1
    lo = max(1, floor(Int, 0.10 * size(trial, 1)))
    hi = min(available, floor(Int, 0.55 * size(trial, 1)))
    hi <= lo && return max(1, min(lo, available))
    return rand(rng, lo:hi)
end

function prediction_window(result, trial; cycles=3.0)
    n_delays = result.hankel.n_delays
    f0 = DMDAnalysis.estimate_stride_frequency(trial; fs=RodentDMDAnalysis.FS, fmin=0.5, fmax=8.0)
    cycle_frames = isfinite(f0) && f0 > 0 ? RodentDMDAnalysis.FS / f0 : 45.0
    horizon = clamp(round(Int, cycles * cycle_frames), 90, 220)
    horizon = min(horizon, max(30, size(trial, 1) - n_delays - 3))
    start = choose_start(trial, n_delays, horizon)
    segment = trial[start:end, :]
    Hhat = DMDAnalysis.rollout(result, segment; n_steps=horizon + 1)
    block = ((n_delays - 1) * result.n_features + 1):(n_delays * result.n_features)
    pred = permutedims(Hhat[block, 1:horizon + 1])
    seed = trial[start:start + n_delays - 1, :]
    real = trial[start:start + n_delays + horizon - 1, :]
    return (
        real=real,
        seed=seed,
        pred=pred,
        n_delays=n_delays,
        horizon=horizon,
        f0=f0,
    )
end

function trace_indices_and_labels(state_space::AbstractString)
    labels = RodentDMDAnalysis.state_labels(state_space)
    if state_space == "centroid_geometry"
        return [2, 4, 6, 8], labels[[2, 4, 6, 8]]
    elseif length(labels) >= 4
        return collect(1:4), labels[1:4]
    else
        return collect(eachindex(labels)), labels
    end
end

function feature_scaled_seed(seed)
    Z = copy(seed)
    for j in axes(Z, 2)
        μ = mean(Z[:, j])
        σ = std(Z[:, j])
        if isfinite(σ) && σ > 0
            Z[:, j] .= (Z[:, j] .- μ) ./ σ
        else
            Z[:, j] .= 0
        end
    end
    return Z
end

function plot_example!(fig, rowpos, row, result, window; label)
    state_space = String(row.state_space)
    labels = RodentDMDAnalysis.state_labels(state_space)
    trace_idx, trace_labels = trace_indices_and_labels(state_space)
    t_real = (0:size(window.real, 1)-1) ./ RodentDMDAnalysis.FS
    t_pred = (window.n_delays-1:window.n_delays+window.horizon-1) ./ RodentDMDAnalysis.FS
    seed_end = (window.n_delays - 1) / RodentDMDAnalysis.FS
    vals = vcat(vec(window.real[:, trace_idx]), vec(window.pred[:, trace_idx]))
    vals = vals[isfinite.(vals)]
    ylo, yhi = extrema(vals)
    pad = max(0.15, 0.12 * (yhi - ylo))
    ylo -= pad
    yhi += pad
    ax = Axis(fig[rowpos, 1:3],
        title="$(RodentDMDAnalysis.label_for(row.group)) | $(row.subject) | $(row.speed_bin)",
        xlabel="Time from prediction window start (s)",
        ylabel=state_space == "centroid_geometry" ? "Centered limb angle (rad)" : "State value")
    RodentDMDAnalysis.panel_label!(fig, (rowpos, 1), label)
    poly!(ax, Point2f[(0.0, ylo), (seed_end, ylo), (seed_end, yhi), (0.0, yhi)];
        color=(RGBf(0.90, 0.72, 0.24), 0.22))
    vlines!(ax, [seed_end]; color=:black, linestyle=:dash, linewidth=1.5)
    for (k, j) in enumerate(trace_idx)
        c = LIMB_COLORS[k]
        lines!(ax, t_real, window.real[:, j]; color=(c, 0.75),
            linestyle=:dot, linewidth=2.0, label=k == 1 ? "measured" : nothing)
        lines!(ax, t_pred, window.pred[:, j]; color=c,
            linewidth=2.6, label=k == 1 ? "autonomous DMD" : nothing)
    end
    text!(ax, 0.02, 0.96;
        text="$(window.n_delays)-frame delay state",
        space=:relative,
        align=(:left, :top),
        fontsize=11,
        color=:black)
    axislegend(ax,
        [LineElement(color=:black, linestyle=:dot, linewidth=2.0),
         LineElement(color=:black, linestyle=:solid, linewidth=2.6),
         PolyElement(color=(RGBf(0.90, 0.72, 0.24), 0.35))],
        ["measured angles", "autonomous prediction", "delay state"];
        position=:rt,
        framevisible=false)
    ylims!(ax, ylo, yhi)

    ax2 = Axis(fig[rowpos, 4],
        title="Delay state",
        xlabel="Delay frame",
        ylabel="")
    hm = heatmap!(ax2, 1:window.n_delays, 1:length(labels), feature_scaled_seed(window.seed);
        colormap=:balance,
        colorrange=(-2.5, 2.5))
    ax2.yticks = (1:length(labels), labels)
    Colorbar(fig[rowpos, 5], hm; label="feature z-score")
    return ax
end

function make_figure0(; cycles=3.0)
    RodentDMDAnalysis.pub_theme!()
    cfg = load_current_config()
    df = CSV.read(joinpath(RESULTS_DIR, "rodent_subject_speed_metrics.csv"), DataFrame)
    results = deserialize(joinpath(RESULTS_DIR, "rodent_dmd_results.jls"))
    examples = candidate_rows(df)
    example_rows = df[in.(df.model_index, Ref(examples)), :]
    bouts, meta = load_training_sets(cfg; target_rows=example_rows)
    fig = Figure(size=(1300, 720))
    for (k, model_index) in enumerate(examples)
        row = df[df.model_index .== model_index, :][1, :]
        trial = trial_for_model(bouts, meta, row)
        window = prediction_window(results[model_index], trial; cycles=cycles)
        plot_example!(fig, k, row, results[model_index], window;
            label=string(Char(Int('A') + k - 1)))
    end
    Legend(fig[3, 1:5],
        [LineElement(color=LIMB_COLORS[i], linewidth=3) for i in eachindex(trace_indices_and_labels(String(df.state_space[1]))[1])],
        trace_indices_and_labels(String(df.state_space[1]))[2];
        orientation=:horizontal,
        framevisible=false,
        tellheight=true)
    resize_to_layout!(fig)
    svg = joinpath(FIG_DIR, "Figure 0.svg")
    save(svg, fig)
    println("Saved: ", svg)
    return fig
end

function main(args=ARGS)
    cycles = parse(Float64, arg_value(args, "--cycles", "3.0"))
    make_figure0(; cycles=cycles)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
