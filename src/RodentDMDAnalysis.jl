module RodentDMDAnalysis

using CairoMakie
using CSV
using DataFrames
using DMDAnalysis
using LinearAlgebra
using ProgressBars
using Random
using Serialization
using Statistics
using Dates

export ROOT, ORIGINAL_ROOT, DATA_DIR, RESULTS_DIR, FIG_DIR, FS,
    COORD_COLS, GROUP_ORDER, GROUP_LABELS,
    rodent_manifest, sample_manifest, load_bouts, add_speed_conditions!,
    combine_training_sets, fit_metrics_cache, load_metrics, make_main_figures

const ROOT = normpath(joinpath(@__DIR__, ".."))
const ORIGINAL_ROOT = normpath(get(ENV, "GAIT_ORIGINAL_REPO", joinpath(ROOT, "..", "Gait_NonNormality_Paper")))
const DATA_DIR = get(ENV, "RODENT_DATA_DIR", joinpath(ORIGINAL_ROOT, "data"))
const RESULTS_DIR = joinpath(ROOT, "results")
const FIG_DIR = joinpath(ROOT, "figures")

const FS = parse(Float64, get(ENV, "RODENT_FS", "100.0"))
const COORD_COLS = [:LFx, :LFy, :RFx, :RFy, :RHx, :RHy, :LHx, :LHy]
const COORD_LABELS = ["LFx", "LFy", "RFx", "RFy", "RHx", "RHy", "LHx", "LHy"]
const PAW_LABELS = ["LF", "RF", "RH", "LH"]
const STATE_LABELS = Dict(
    "raw_xy" => COORD_LABELS,
    "centroid_xy" => ["LFx-c", "LFy-c", "RFx-c", "RFy-c", "RHx-c", "RHy-c", "LHx-c", "LHy-c"],
    "body_aligned_xy" => ["LFx-body", "LFy-body", "RFx-body", "RFy-body", "RHx-body", "RHy-body", "LHx-body", "LHy-body"],
    "body_ap" => ["LF AP", "RF AP", "RH AP", "LH AP"],
    "centroid_geometry" => ["LF radius", "LF angle", "RF radius", "RF angle", "RH radius", "RH angle", "LH radius", "LH angle"],
    "centroid_angles" => ["LF angle", "RF angle", "RH angle", "LH angle"],
    "centroid_angle_sincos" => ["cos LF", "sin LF", "cos RF", "sin RF", "cos RH", "sin RH", "cos LH", "sin LH"],
    "centroid_geometry_sincos" => ["cos LF", "sin LF", "cos RF", "sin RF", "cos RH", "sin RH", "cos LH", "sin LH"],
)
const GROUP_ORDER = ["Cntnap_neg", "Cntnap_het", "Cntnap_hom",
                     "TSC_neg", "TSC_het", "TSC_hom", "new_preds"]
const GROUP_LABELS = Dict(
    "Cntnap_neg" => "Cntnap control",
    "Cntnap_het" => "Cntnap het",
    "Cntnap_hom" => "Cntnap hom/KO",
    "TSC_neg" => "TSC control",
    "TSC_het" => "TSC het",
    "TSC_hom" => "TSC hom/mut",
    "new_preds" => "New preds",
)
const GROUP_COLORS = Dict(
    "Cntnap_neg" => RGBf(0.00, 0.45, 0.70),
    "Cntnap_het" => RGBf(0.35, 0.70, 0.90),
    "Cntnap_hom" => RGBf(0.00, 0.62, 0.45),
    "TSC_neg" => RGBf(0.80, 0.47, 0.65),
    "TSC_het" => RGBf(0.90, 0.62, 0.00),
    "TSC_hom" => RGBf(0.84, 0.37, 0.00),
    "new_preds" => RGBf(0.25, 0.25, 0.25),
)

ensure_dirs!() = (mkpath(RESULTS_DIR); mkpath(FIG_DIR); nothing)
label_for(g) = get(GROUP_LABELS, String(g), String(g))
color_for(g) = get(GROUP_COLORS, String(g), RGBf(0.3, 0.3, 0.3))

function pub_theme!()
    CairoMakie.activate!(type="svg")
    set_theme!(Theme(
        fontsize=10,
        Axis=(xlabelsize=11, ylabelsize=11, titlesize=12,
              xgridvisible=false, ygridvisible=false),
        Lines=(linewidth=2.2,),
        Scatter=(markersize=7,),
    ))
end

function savefig(name, fig)
    ensure_dirs!()
    path = joinpath(FIG_DIR, name)
    save(path, fig)
    println("Saved: ", path)
    return path
end

function panel_label!(fig, pos, label)
    Label(fig[pos..., TopLeft()], label; fontsize=24, font=:bold,
        padding=(0, 8, 8, 0), halign=:right, tellwidth=false, tellheight=false)
end

function finite_vec(v)
    x = Float64.(collect(skipmissing(v)))
    return x[isfinite.(x)]
end

maybe_mean(v) = (vals = finite_vec(v); isempty(vals) ? NaN : mean(vals))
maybe_median(v) = (vals = finite_vec(v); isempty(vals) ? NaN : median(vals))

function sem(v)
    vals = finite_vec(v)
    length(vals) <= 1 && return 0.0
    return std(vals) / sqrt(length(vals))
end

timestamp() = Dates.format(now(), dateformat"yyyy-mm-dd HH:MM:SS")
logmsg(verbose, msg) = verbose && (println("[$(timestamp())] ", msg); flush(stdout))

function csv_paths(; data_dir=DATA_DIR, dataset="all")
    roots = String[]
    dataset in ("all", "asd") && push!(roots, joinpath(data_dir, "ASD_bouts_allAngles"))
    dataset in ("all", "new_preds") && push!(roots, joinpath(data_dir, "new_preds_bouts_allAngles"))
    paths = String[]
    for root in roots
        isdir(root) || continue
        for (dir, _, files) in walkdir(root)
            occursin("__MACOSX", dir) && continue
            for file in files
                endswith(file, ".csv") || continue
                startswith(file, "._") && continue
                push!(paths, joinpath(dir, file))
            end
        end
    end
    return sort(paths)
end

function path_metadata(path)
    base = basename(path)
    m = match(r"^(.*)\.predictions_([0-9]+)\.csv$", base)
    subject = isnothing(m) ? splitext(base)[1] : m.captures[1]
    bout = isnothing(m) ? missing : parse(Int, m.captures[2])
    parts = splitpath(path)
    dataset = occursin("ASD_bouts_allAngles", path) ? "ASD" : "new_preds"
    group = dataset == "ASD" ? parts[end - 1] : "new_preds"
    line = startswith(group, "Cntnap") ? "Cntnap" :
           startswith(group, "TSC") ? "TSC" : dataset
    genotype = occursin("_", group) ? split(group, "_"; limit=2)[2] : group
    return (path=path, dataset=dataset, group=group, line=line,
            genotype=genotype, subject=subject, bout=bout)
end

function rodent_manifest(; data_dir=DATA_DIR, dataset="all")
    return DataFrame([path_metadata(p) for p in csv_paths(; data_dir=data_dir, dataset=dataset)])
end

function sample_manifest(manifest; max_per_group=250, max_per_subject=30, seed=20260721)
    rng = MersenneTwister(seed)
    selected = Int[]
    for gdf in groupby(manifest, :group)
        group_idx = Int[]
        for sdf in groupby(gdf, :subject)
            idx = collect(parentindices(sdf)[1])
            shuffle!(rng, idx)
            append!(group_idx, idx[1:min(max_per_subject, length(idx))])
        end
        shuffle!(rng, group_idx)
        append!(selected, group_idx[1:min(max_per_group, length(group_idx))])
    end
    return manifest[sort(selected), :]
end

function normalizer(name::AbstractString)
    name == "none" && return NoNormalization()
    name == "center" && return Center()
    name == "zscore" && return ZScore()
    error("Unsupported normalization=$name")
end

function savgol_kernel(window::Int, order::Int)
    window == 0 && return Float64[1.0]
    isodd(window) || error("Savitzky-Golay window must be odd")
    window > order || error("Savitzky-Golay window must be larger than polynomial order")
    half = div(window, 2)
    x = collect(-half:half)
    V = hcat([Float64.(x) .^ p for p in 0:order]...)
    return vec((pinv(V)'[:, 1]))
end

reflect_index(i, n) = i < 1 ? 2 - i : (i > n ? 2n - i : i)

function smooth_columns_savgol(X::AbstractMatrix, window::Int, order::Int=3)
    window == 0 && return Matrix{Float64}(X)
    k = savgol_kernel(window, order)
    half = div(window, 2)
    Y = similar(Matrix{Float64}(X))
    @inbounds for j in axes(X, 2), i in axes(X, 1)
        s = 0.0
        for (kk, w) in enumerate(k)
            ii = reflect_index(i + kk - half - 1, size(X, 1))
            s += w * X[ii, j]
        end
        Y[i, j] = s
    end
    return Y
end

function canonical_state_space(state_space::AbstractString)
    state_space == "centroid_geometry_sincos" && return "centroid_angle_sincos"
    return state_space
end

state_labels(state_space::AbstractString) = get(STATE_LABELS, state_space,
    get(STATE_LABELS, canonical_state_space(state_space), COORD_LABELS))

function cached_state_space()
    path = joinpath(RESULTS_DIR, "rodent_config.jls")
    isfile(path) || return nothing
    try
        cfg = deserialize(path)
        :state_space in propertynames(cfg) && return String(cfg.state_space)
    catch err
        @warn "Could not read rodent config for state-space labels" path err
    end
    return nothing
end

function metric_state_space(df)
    if hasproperty(df, :state_space) && nrow(df) > 0
        return String(df.state_space[1])
    end
    cached = cached_state_space()
    isnothing(cached) || return cached
    return "raw_xy"
end

function split_paws(X)
    n = size(X, 1)
    xs = zeros(Float64, n, 4)
    ys = zeros(Float64, n, 4)
    for j in 1:4
        xs[:, j] .= X[:, 2j - 1]
        ys[:, j] .= X[:, 2j]
    end
    return xs, ys
end

function pack_paws(xs, ys)
    X = zeros(Float64, size(xs, 1), 8)
    for j in 1:4
        X[:, 2j - 1] .= xs[:, j]
        X[:, 2j] .= ys[:, j]
    end
    return X
end

function pack_angle_sincos(A)
    X = zeros(Float64, size(A, 1), 8)
    for j in 1:4
        X[:, 2j - 1] .= cos.(A[:, j])
        X[:, 2j] .= sin.(A[:, j])
    end
    return X
end

function centroid_xy(X)
    xs, ys = split_paws(X)
    cx = mean(xs; dims=2)
    cy = mean(ys; dims=2)
    return xs .- cx, ys .- cy, vec(cx), vec(cy)
end

function body_axis(xs, ys)
    fore_x = mean(xs[:, 1:2]; dims=2)
    fore_y = mean(ys[:, 1:2]; dims=2)
    hind_x = mean(xs[:, 3:4]; dims=2)
    hind_y = mean(ys[:, 3:4]; dims=2)
    dx = vec(fore_x .- hind_x)
    dy = vec(fore_y .- hind_y)
    len = sqrt.(dx .^ 2 .+ dy .^ 2)
    len .= ifelse.(len .<= 0, 1.0, len)
    return atan.(dy, dx), len
end

function unwrap_columns!(X)
    for j in axes(X, 2), i in 2:size(X, 1)
        d = X[i, j] - X[i - 1, j]
        if d > pi
            X[i:end, j] .-= 2pi
        elseif d < -pi
            X[i:end, j] .+= 2pi
        end
    end
    return X
end

function rodent_state_space(X; state_space="centroid_xy")
    state_space = canonical_state_space(state_space)
    state_space == "raw_xy" && return Matrix{Float64}(X)
    xs0, ys0, _, _ = centroid_xy(X)
    state_space == "centroid_xy" && return pack_paws(xs0, ys0)
    if state_space == "body_aligned_xy" || state_space == "body_ap"
        theta, body_length = body_axis(xs0, ys0)
        c = cos.(-theta)
        s = sin.(-theta)
        xs = xs0 .* c .- ys0 .* s
        ys = xs0 .* s .+ ys0 .* c
        state_space == "body_ap" && return xs ./ body_length
        return pack_paws(xs ./ body_length, ys ./ body_length)
    elseif state_space == "centroid_geometry"
        R = sqrt.(xs0 .^ 2 .+ ys0 .^ 2)
        A = atan.(ys0, xs0)
        unwrap_columns!(A)
        return pack_paws(R, A)
    elseif state_space == "centroid_angles"
        A = atan.(ys0, xs0)
        unwrap_columns!(A)
        return A
    elseif state_space == "centroid_angle_sincos"
        A = atan.(ys0, xs0)
        return pack_angle_sincos(A)
    end
    error("Unsupported state_space=$state_space. Use raw_xy, centroid_xy, body_aligned_xy, body_ap, centroid_geometry, centroid_angles, centroid_angle_sincos, or centroid_geometry_sincos.")
end

function phase_normalizer(name::AbstractString, steps::Int)
    name in ("none", "false", "off") && return NoPhaseNormalization()
    name in ("phase", "normalize", "true", "on") && return PhaseNormalize(steps=steps)
    error("Unsupported phase_normalization=$name")
end

function load_bout(path; normalization="center", state_space="centroid_xy", smooth_window=0,
                   phase_normalization="none", phase_steps=100)
    df = CSV.read(path, DataFrame)
    missing_cols = setdiff(COORD_COLS, Symbol.(names(df)))
    isempty(missing_cols) || error("Missing coordinate columns in $path: $missing_cols")
    X = reduce(hcat, [Float64[(ismissing(v) ? NaN : Float64(v)) for v in df[!, c]]
                      for c in COORD_COLS])
    ok = vec(all(isfinite, X; dims=2))
    X = X[ok, :]
    X = rodent_state_space(X; state_space=state_space)
    X = smooth_columns_savgol(X, smooth_window, 3)
    pre = Preprocessor(normalization=normalizer(normalization),
        phase=phase_normalizer(phase_normalization, phase_steps),
        state_mode=:position, fs=FS)
    X = only(DMDAnalysis.transform(pre, X).trials)
    speed_all = hasproperty(df, :speed) ? Float64[(ismissing(v) ? NaN : Float64(v)) for v in df.speed] :
                fill(NaN, nrow(df))
    frame_all = hasproperty(df, :frame) ? Float64[(ismissing(v) ? NaN : Float64(v)) for v in df.frame] :
                collect(1.0:nrow(df))
    return (data=X, speed=speed_all[ok], frame=frame_all[ok])
end

function load_bouts(manifest; normalization="center", state_space="centroid_xy", smooth_window=0,
                    phase_normalization="none", phase_steps=100, min_frames=35,
                    progressbar=true, verbose=true, log_interval=5000)
    trials = Matrix{Float64}[]
    rows = NamedTuple[]
    iter = progressbar ? ProgressBar(eachrow(manifest)) : eachrow(manifest)
    for (i, row) in enumerate(iter)
        try
            bout = load_bout(row.path; normalization=normalization, state_space=state_space,
                smooth_window=smooth_window, phase_normalization=phase_normalization,
                phase_steps=phase_steps)
            size(bout.data, 1) >= min_frames || continue
            sp = finite_vec(bout.speed)
            push!(trials, bout.data)
            push!(rows, merge(NamedTuple(row), (
                n_frames=size(bout.data, 1),
                mean_speed=isempty(sp) ? NaN : mean(sp),
                median_speed=isempty(sp) ? NaN : median(sp),
                p90_speed=isempty(sp) ? NaN : quantile(sp, 0.90),
            )))
        catch err
            @warn "Skipping rodent bout" path=row.path err
        end
        if !progressbar && verbose && (i == 1 || i == nrow(manifest) || i % log_interval == 0)
            logmsg(true, "loaded-bout scan $(i)/$(nrow(manifest)); kept=$(length(trials))")
        end
    end
    if isempty(rows)
        return trials, DataFrame(
            path=String[], dataset=String[], group=String[], line=String[],
            genotype=String[], subject=String[], bout=Union{Missing,Int}[],
            n_frames=Int[], mean_speed=Float64[], median_speed=Float64[],
            p90_speed=Float64[],
        )
    end
    return trials, DataFrame(rows)
end

function speed_labels(speeds; nbins=3)
    vals = finite_vec(speeds)
    isempty(vals) && return fill("unknown", length(speeds))
    labels = nbins == 3 ? ["slow", "medium", "fast"] : ["speed_$i" for i in 1:nbins]
    edges = quantile(vals, range(0, 1; length=nbins + 1))
    out = fill("unknown", length(speeds))
    for i in eachindex(speeds)
        v = Float64(speeds[i])
        isfinite(v) || continue
        out[i] = labels[clamp(searchsortedlast(edges, v), 1, nbins)]
    end
    return out
end

function add_speed_conditions!(meta; nbins=3, mode="subject")
    meta.speed_condition = fill("unknown", nrow(meta))
    if mode == "global"
        meta.speed_condition .= speed_labels(meta.mean_speed; nbins=nbins)
    elseif mode == "group"
        for sub in groupby(meta, :group)
            idx = parentindices(sub)[1]
            meta.speed_condition[idx] = speed_labels(sub.mean_speed; nbins=nbins)
        end
    elseif mode == "subject"
        for sub in groupby(meta, [:group, :subject])
            idx = parentindices(sub)[1]
            meta.speed_condition[idx] = speed_labels(sub.mean_speed; nbins=nbins)
        end
    else
        error("Unsupported speed-condition mode=$mode")
    end
    return meta
end

function estimated_stride_count(bouts)
    total = 0.0
    for b in bouts
        f0 = estimate_stride_frequency(b; fs=FS, fmin=0.5, fmax=8.0)
        isfinite(f0) && (total += size(b, 1) / FS * f0)
    end
    return total
end

function combine_training_sets(trials, meta; min_bouts=5, min_frames=0, min_strides=8.0)
    sets = Vector{Vector{Matrix{Float64}}}()
    rows = NamedTuple[]
    meta2 = hcat(meta, DataFrame(_trial_pos=1:nrow(meta)))
    for sub in groupby(meta2, [:group, :subject, :speed_condition])
        nrow(sub) >= min_bouts || continue
        idx = Int.(sub._trial_pos)
        n_frames = sum(Int.(sub.n_frames))
        n_frames >= min_frames || continue
        selected = trials[idx]
        n_strides = estimated_stride_count(selected)
        n_strides >= min_strides || continue
        speeds = finite_vec(sub.mean_speed)
        push!(sets, selected)
        push!(rows, (
            dataset=sub.dataset[1],
            group=sub.group[1],
            line=sub.line[1],
            genotype=sub.genotype[1],
            state_space=hasproperty(sub, :state_space) ? sub.state_space[1] : "unknown",
            subject=sub.subject[1],
            speed_condition=sub.speed_condition[1],
            n_bouts=length(idx),
            n_frames=n_frames,
            n_est_strides=n_strides,
            mean_speed=isempty(speeds) ? NaN : mean(speeds),
            median_speed=isempty(speeds) ? NaN : median(speeds),
            speed_min=isempty(speeds) ? NaN : minimum(speeds),
            speed_max=isempty(speeds) ? NaN : maximum(speeds),
        ))
    end
    return sets, DataFrame(rows)
end

function rank_rule(rank::Int, noise_floor)
    rank > 0 && return FixedRank(rank)
    return EnergyRank(noise_floor)
end

function fit_one_model(bouts, meta_row; n_delays=26, delay_interval=1, noise_floor=1e-6,
                       rank=0, stabilize=true, kmax=40, n_angles=720)
    h = HankelGenerator(n_delays=n_delays, delay_interval=delay_interval)
    model = DMD(rank=rank_rule(rank, noise_floor), stabilize=stabilize)
    result = DMDAnalysis.fit(model, DMDAnalysis.transform(h, bouts))
    f0s = [estimate_stride_frequency(b; fs=FS, fmin=0.5, fmax=8.0) for b in bouts]
    f0 = maybe_median(f0s)
    metricset = MetricSet(
        Henrici(),
        TransientGrowth(kmax=kmax),
        ResolventPeak(fs=FS, n_angles=n_angles, fmin=0.5),
        HarmonicMetrics(fs=FS, n_harmonics=5),
        IOGain(fs=FS),
    )
    row = analyze(result, metricset; metadata=merge(NamedTuple(meta_row), (stride_frequency=f0,)))
    pr = begin
        p = Float64.(result.singular_values).^2
        s = sum(p)
        s <= 0 ? NaN : (p ./= s; 1 / sum(p .^ 2))
    end
    qvals = finite_vec([row[Symbol("qfactor_h$h")] for h in 1:5])
    qvals = qvals[qvals .> 0]
    return merge(row, (
        participation_ratio=pr,
        qfactor_median=isempty(qvals) ? NaN : median(qvals),
    )), result
end

function cache_config(; dataset, full, max_per_group, max_per_subject, seed,
                      normalization, state_space, speed_mode, speed_conditions,
                      min_bouts, min_frames_per_model, min_strides_per_model,
                      n_delays, delay_interval, noise_floor, stabilize, n_angles,
                      rank, smooth_window, phase_normalization, phase_steps,
                      checkpoint_interval, progressbar)
    return (dataset=dataset, full=full, max_per_group=max_per_group,
        max_per_subject=max_per_subject, seed=seed, normalization=normalization,
        state_space=state_space, smooth_window=smooth_window,
        phase_normalization=phase_normalization, phase_steps=phase_steps,
        speed_mode=speed_mode, speed_conditions=speed_conditions,
        min_bouts=min_bouts, min_frames_per_model=min_frames_per_model,
        min_strides_per_model=min_strides_per_model,
        n_delays=n_delays, delay_interval=delay_interval,
        noise_floor=noise_floor, rank=rank, stabilize=stabilize, n_angles=n_angles,
        checkpoint_interval=checkpoint_interval, progressbar=progressbar,
        metric_version=7)
end

function fit_metrics_cache(; dataset="all", full=false, max_per_group=250,
                           max_per_subject=30, seed=20260721, normalization="center",
                           state_space="centroid_xy", smooth_window=0,
                           phase_normalization="none", phase_steps=100,
                           speed_mode="subject", speed_conditions=3, min_bouts=5,
                           min_frames_per_model=0, min_strides_per_model=8.0,
                           n_delays=26, delay_interval=1, noise_floor=1e-6,
                           rank=0, stabilize=true, n_angles=720, checkpoint_interval=25,
                           progressbar=true, verbose=true, log_interval=1,
                           force=false)
    ensure_dirs!()
    metric_path = joinpath(RESULTS_DIR, "rodent_subject_speed_metrics.csv")
    result_path = joinpath(RESULTS_DIR, "rodent_dmd_results.jls")
    config_path = joinpath(RESULTS_DIR, "rodent_config.jls")
    checkpoint_metric_path = joinpath(RESULTS_DIR, "rodent_subject_speed_metrics.checkpoint.csv")
    checkpoint_result_path = joinpath(RESULTS_DIR, "rodent_dmd_results.checkpoint.jls")
    checkpoint_config_path = joinpath(RESULTS_DIR, "rodent_config.checkpoint.jls")
    cfg = cache_config(; dataset, full, max_per_group, max_per_subject, seed,
        normalization, state_space, speed_mode, speed_conditions, min_bouts,
        min_frames_per_model, min_strides_per_model, n_delays, delay_interval,
        noise_floor, stabilize, n_angles, rank, smooth_window, phase_normalization,
        phase_steps, checkpoint_interval, progressbar)
    if !force && isfile(metric_path) && isfile(result_path) && isfile(config_path)
        try
            deserialize(config_path) == cfg && return CSV.read(metric_path, DataFrame)
        catch
        end
    end

    logmsg(verbose, "Rodent DMD fit starting")
    logmsg(verbose, "config=$(cfg)")
    logmsg(verbose, "Loading rodent manifest from $(DATA_DIR)")
    manifest = rodent_manifest(; dataset=dataset)
    logmsg(verbose, "found $(nrow(manifest)) CSV bouts")
    if !full
        manifest = sample_manifest(manifest; max_per_group=max_per_group,
            max_per_subject=max_per_subject, seed=seed)
        logmsg(verbose, "sampled $(nrow(manifest)) bouts")
    end
    logmsg(verbose, "Loading bouts with state_space=$(state_space), normalization=$(normalization), smooth_window=$(smooth_window), phase_normalization=$(phase_normalization), phase_steps=$(phase_steps)")
    bouts, meta = load_bouts(manifest; normalization=normalization, state_space=state_space,
        smooth_window=smooth_window, phase_normalization=phase_normalization,
        phase_steps=phase_steps, progressbar=progressbar, verbose=verbose)
    nrow(meta) > 0 || error("No rodent bouts were loaded. Check DATA_DIR=$(DATA_DIR), coordinate columns, and min_frames.")
    logmsg(verbose, "loaded $(length(bouts)) usable bouts after min frame and finite-coordinate filtering")
    meta.state_space = fill(state_space, nrow(meta))
    logmsg(verbose, "Assigning speed conditions with mode=$(speed_mode), bins=$(speed_conditions)")
    add_speed_conditions!(meta; nbins=speed_conditions, mode=speed_mode)
    logmsg(verbose, "Combining subject x speed-condition training sets")
    sets, model_meta = combine_training_sets(bouts, meta; min_bouts=min_bouts,
        min_frames=min_frames_per_model, min_strides=min_strides_per_model)
    length(sets) > 0 || error("No subject x speed-condition training sets met min_bouts=$min_bouts, min_frames_per_model=$min_frames_per_model, and min_strides_per_model=$min_strides_per_model. Increase sampling/full data or relax thresholds for diagnostics.")
    logmsg(verbose, "$(length(sets)) subject x speed-condition models from $(sum(model_meta.n_bouts)) bouts")
    logmsg(verbose, "median bouts/model=$(round(median(model_meta.n_bouts); digits=2)), median frames/model=$(round(median(model_meta.n_frames); digits=1)), median estimated strides/model=$(round(median(model_meta.n_est_strides); digits=1))")

    rows = NamedTuple[]
    results = DMDResult[]
    checkpoint!() = begin
        isempty(rows) && return nothing
        dfc = DataFrame(rows)
        dfc.speed_bin = String.(dfc.speed_condition)
        CSV.write(checkpoint_metric_path, dfc)
        serialize(checkpoint_result_path, results)
        serialize(checkpoint_config_path, cfg)
        logmsg(verbose, "checkpointed $(nrow(dfc)) fitted models to $(checkpoint_metric_path)")
        return nothing
    end
    iter = progressbar ? ProgressBar(eachindex(sets)) : eachindex(sets)
    total_fit_time = 0.0
    for i in iter
        meta_i = model_meta[i, :]
        dims = "$(size(first(sets[i]), 2)) vars x $(n_delays) delays"
        transitions = sum(max(size(b, 1) - (n_delays - 1) * delay_interval - 1, 0) for b in sets[i])
        should_log = verbose && (i == 1 || i == length(sets) || log_interval <= 1 || i % log_interval == 0)
        should_log && logmsg(true, "model $(i)/$(length(sets)) start: group=$(meta_i.group), subject=$(meta_i.subject), speed=$(meta_i.speed_condition), bouts=$(meta_i.n_bouts), frames=$(meta_i.n_frames), est_strides=$(round(meta_i.n_est_strides; digits=1)), transitions=$(transitions), state=$(dims)")
        t0 = time()
        try
            row, result = fit_one_model(sets[i], meta_i;
                n_delays=n_delays, delay_interval=delay_interval,
                noise_floor=noise_floor, rank=rank, stabilize=stabilize, n_angles=n_angles)
            elapsed = time() - t0
            total_fit_time += elapsed
            push!(rows, merge(row, (model_index=i,)))
            push!(results, result)
            should_log && logmsg(true, "model $(i)/$(length(sets)) done in $(round(elapsed; digits=2))s: fitted=$(length(rows)), rank=$(row.rank), henrici=$(round(row.henrici; digits=4)), log_tgmax=$(round(row.log_tgmax; digits=4))")
            checkpoint_interval > 0 && length(rows) % checkpoint_interval == 0 && checkpoint!()
        catch err
            elapsed = time() - t0
            @warn "Rodent model failed" i group=model_meta.group[i] subject=model_meta.subject[i] elapsed err
        end
    end
    checkpoint!()
    isempty(rows) && error("All rodent DMD models failed after filtering. Check warnings above and consider relaxing rank/state-space/metric settings.")
    logmsg(verbose, "Fitting loop complete: fitted $(length(rows))/$(length(sets)) models in $(round(total_fit_time / 60; digits=2)) minutes")
    df = DataFrame(rows)
    df.speed_bin = String.(df.speed_condition)
    logmsg(verbose, "Writing final metrics/results")
    CSV.write(metric_path, df)
    serialize(result_path, results)
    serialize(config_path, cfg)
    write_summaries(df)
    logmsg(verbose, "Rodent DMD fit complete: $(metric_path)")
    return df
end

load_metrics() = CSV.read(joinpath(RESULTS_DIR, "rodent_subject_speed_metrics.csv"), DataFrame)

function corr_or_nan(x, y)
    xx = Float64.(x); yy = Float64.(y)
    ok = isfinite.(xx) .& isfinite.(yy)
    count(ok) >= 3 || return NaN
    std(xx[ok]) == 0 || std(yy[ok]) == 0 ? NaN : cor(xx[ok], yy[ok])
end

function slope_or_nan(x, y)
    xx = Float64.(x); yy = Float64.(y)
    ok = isfinite.(xx) .& isfinite.(yy)
    count(ok) >= 3 || return NaN
    std(xx[ok]) == 0 && return NaN
    return cov(xx[ok], yy[ok]) / var(xx[ok])
end

function summarize_by(df, keys, metrics)
    rows = NamedTuple[]
    for sub in groupby(df, keys)
        base = NamedTuple(sub[1, keys])
        pairs = Pair{Symbol,Any}[:n_models => nrow(sub), :n_subjects => length(unique(sub.subject)),
                                 :n_bouts => sum(Int.(sub.n_bouts))]
        for m in metrics
            vals = finite_vec(sub[!, m])
            if isempty(vals)
                append!(pairs, [Symbol(m, "_mean") => NaN, Symbol(m, "_sem") => NaN,
                                Symbol(m, "_median") => NaN])
            else
                append!(pairs, [Symbol(m, "_mean") => mean(vals), Symbol(m, "_sem") => sem(vals),
                                Symbol(m, "_median") => median(vals)])
            end
        end
        push!(rows, merge(base, NamedTuple(pairs)))
    end
    return DataFrame(rows)
end

function write_summaries(df)
    metrics = [:mean_speed, :henrici, :log_tgmax, :participation_ratio,
        :rank, :stride_frequency, :resolvent_peak_gain, :qfactor_median,
        :harmonic_entropy, :h1_concentration, :spectral_decay_rate,
        :io_gain_max]
    CSV.write(joinpath(RESULTS_DIR, "rodent_group_summary.csv"), summarize_by(df, [:group], metrics))
    CSV.write(joinpath(RESULTS_DIR, "rodent_speed_bin_summary.csv"), summarize_by(df, [:group, :speed_bin], metrics))
    trend_rows = NamedTuple[]
    for sub in groupby(df, :group), m in [:henrici, :log_tgmax, :participation_ratio, :qfactor_median, :h1_concentration, :io_gain_max]
        push!(trend_rows, (group=sub.group[1], metric=String(m), n_models=nrow(sub),
            n_subjects=length(unique(sub.subject)),
            r_speed=corr_or_nan(sub.mean_speed, sub[!, m]),
            slope_per_speed=slope_or_nan(sub.mean_speed, sub[!, m]),
            r_henrici=m == :henrici ? NaN : corr_or_nan(sub.henrici, sub[!, m])))
    end
    CSV.write(joinpath(RESULTS_DIR, "rodent_speed_and_henrici_trends.csv"), DataFrame(trend_rows))
    io_rows = NamedTuple[]
    labels = state_labels(metric_state_space(df))
    for sub in groupby(df, [:group, :speed_bin])
        for j in 1:length(labels)
            push!(io_rows, (group=sub.group[1], speed_bin=sub.speed_bin[1],
                coordinate=labels[j],
                forcing_mean=maybe_mean(sub[!, Symbol("io_", j, "_from_", j)]),
                n_models=nrow(sub)))
        end
    end
    CSV.write(joinpath(RESULTS_DIR, "rodent_io_structure_summary.csv"), DataFrame(io_rows))
end

groups_present(df) = [g for g in GROUP_ORDER if g in unique(df.group)]

function binned_stats(df, metric; nbins=4)
    rows = NamedTuple[]
    for g in groups_present(df)
        sub = df[df.group .== g, :]
        ok = isfinite.(sub.mean_speed) .& isfinite.(sub[!, metric])
        count(ok) >= 3 || continue
        s = sub[ok, :]
        edges = quantile(Float64.(s.mean_speed), range(0, 1; length=nbins + 1))
        for b in 1:nbins
            idx = b == nbins ? findall((s.mean_speed .>= edges[b]) .& (s.mean_speed .<= edges[b + 1])) :
                                findall((s.mean_speed .>= edges[b]) .& (s.mean_speed .< edges[b + 1]))
            isempty(idx) && continue
            vals = Float64.(s[idx, metric])
            push!(rows, (group=g, speed=mean(Float64.(s.mean_speed[idx])),
                value=mean(vals), sem=length(vals) > 1 ? std(vals) / sqrt(length(vals)) : 0.0))
        end
    end
    return DataFrame(rows)
end

function speed_panel!(ax, df, metric; ylabel, title)
    ax.xlabel = "Mean bout speed"
    ax.ylabel = ylabel
    ax.title = title
    bs = binned_stats(df, metric)
    for g in groups_present(df)
        sub = bs[bs.group .== g, :]
        nrow(sub) == 0 && continue
        c = color_for(g)
        lines!(ax, sub.speed, sub.value; color=c, linewidth=2.8, label=label_for(g))
        scatter!(ax, sub.speed, sub.value; color=c, markersize=7)
        for r in eachrow(sub)
            lines!(ax, [r.speed, r.speed], [r.value - r.sem, r.value + r.sem]; color=c, linewidth=1.5)
        end
    end
end

function scatter_fit_panel!(ax, df, xcol, ycol; xlabel, ylabel, title)
    ax.xlabel = xlabel
    ax.ylabel = ylabel
    ax.title = title
    for g in groups_present(df)
        sub = df[df.group .== g, :]
        ok = isfinite.(sub[!, xcol]) .& isfinite.(sub[!, ycol])
        count(ok) >= 3 || continue
        c = color_for(g)
        x = Float64.(sub[ok, xcol]); y = Float64.(sub[ok, ycol])
        scatter!(ax, x, y; color=(c, 0.25), markersize=4)
        b = slope_or_nan(x, y)
        isfinite(b) || continue
        a = mean(y) - b * mean(x)
        xs = range(quantile(x, 0.05), quantile(x, 0.95); length=80)
        lines!(ax, xs, a .+ b .* xs; color=c, linewidth=2.2)
    end
end

function make_trend_figure(df)
    fig = Figure(size=(1080, 760))
    specs = [
        (:henrici, "Normalized Henrici departure", "Non-normality vs speed"),
        (:log_tgmax, "log10 max TG", "Transient growth vs speed"),
        (:participation_ratio, "Participation ratio", "Effective dimensionality"),
        (:qfactor_median, "Median harmonic Q", "Harmonic selectivity"),
        (:h1_concentration, "H1 / median H2-H5", "Fundamental harmonic strength"),
    ]
    for (i, (metric, ylabel, title)) in enumerate(specs)
        row = i <= 3 ? 1 : 2
        col = i <= 3 ? i : i - 3
        ax = Axis(fig[row, col])
        panel_label!(fig, (row, col), string(Char(Int('A') + i - 1)))
        speed_panel!(ax, df, metric; ylabel=ylabel, title=title)
    end
    axF = Axis(fig[2, 3])
    panel_label!(fig, (2, 3), "F")
    scatter_fit_panel!(axF, df, :log_tgmax, :henrici;
        xlabel="log10 max TG", ylabel="Normalized Henrici departure",
        title="Transient growth vs Henrici")
    Legend(fig[3, 1:3],
        [MarkerElement(color=color_for(g), marker=:circle, markersize=10) for g in groups_present(df)],
        [label_for(g) for g in groups_present(df)];
        orientation=:horizontal, framevisible=false, nbanks=2)
    savefig("Figure 1.svg", fig)
end

function io_matrix_from_row(row)
    n = hasproperty(row, :n_features) ? Int(row.n_features) : length(COORD_LABELS)
    M = fill(NaN, n, n)
    for i in 1:n, j in 1:n
        col = Symbol("io_", i, "_from_", j)
        hasproperty(row, col) && (M[i, j] = Float64(row[col]))
    end
    return M
end

function group_io_matrix(df, group)
    mats = Matrix{Float64}[]
    for row in eachrow(df[df.group .== group, :])
        M = io_matrix_from_row(row)
        any(.!isfinite.(M)) && continue
        m = maximum(M)
        m > 0 || continue
        push!(mats, M ./ m)
    end
    n = hasproperty(df, :n_features) && nrow(df) > 0 ? Int(first(df.n_features)) : length(COORD_LABELS)
    isempty(mats) && return fill(NaN, n, n)
    out = zeros(size(first(mats)))
    for M in mats
        out .+= M
    end
    return out ./ length(mats)
end

function make_io_figure(df)
    groups = groups_present(df)
    labels = state_labels(metric_state_space(df))
    ncols = min(4, length(groups))
    nrows = ceil(Int, length(groups) / ncols)
    fig = Figure(size=(1120, 320 * nrows))
    for (k, g) in enumerate(groups)
        row = div(k - 1, ncols) + 1
        col = mod(k - 1, ncols) + 1
        ax = Axis(fig[row, col], title=label_for(g),
            xticks=(1:length(labels), labels),
            yticks=(1:length(labels), labels),
            xticklabelrotation=pi / 4)
        panel_label!(fig, (row, col), string(Char(Int('A') + k - 1)))
        hm = heatmap!(ax, group_io_matrix(df, g); colormap=:viridis, colorrange=(0, 1))
        k == ncols && Colorbar(fig[row, ncols + 1], hm; label="Mean normalized gain")
    end
    savefig("Figure 2.svg", fig)
end

function make_main_figures()
    pub_theme!()
    df = load_metrics()
    make_trend_figure(df)
    make_io_figure(df)
end

end
