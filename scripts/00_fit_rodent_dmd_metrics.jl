#!/usr/bin/env julia

include(joinpath(@__DIR__, "..", "src", "RodentDMDAnalysis.jl"))
using .RodentDMDAnalysis

function arg_value(args, name, default)
    i = findfirst(==(name), args)
    isnothing(i) && return default
    i == length(args) && error("Missing value after $name")
    return args[i + 1]
end

has_flag(args, name) = any(==(name), args)

function fit_main(args=ARGS)
    RodentDMDAnalysis.fit_metrics_cache(
        dataset=arg_value(args, "--dataset", "all"),
        full=has_flag(args, "--full"),
        max_per_group=parse(Int, arg_value(args, "--max-per-group", "250")),
        max_per_subject=parse(Int, arg_value(args, "--max-per-subject", "30")),
        seed=parse(Int, arg_value(args, "--seed", "20260721")),
        normalization=arg_value(args, "--normalize", "center"),
        state_space=arg_value(args, "--state-space", "centroid_xy"),
        smooth_window=parse(Int, arg_value(args, "--smooth-window", "0")),
        phase_normalization=has_flag(args, "--phase-normalize") ? "phase" : arg_value(args, "--phase-normalization", "none"),
        phase_steps=parse(Int, arg_value(args, "--phase-steps", "100")),
        speed_mode=arg_value(args, "--speed-condition-mode", "subject"),
        speed_conditions=parse(Int, arg_value(args, "--speed-conditions", "3")),
        min_bouts=parse(Int, arg_value(args, "--min-bouts-per-model", "5")),
        min_frames_per_model=parse(Int, arg_value(args, "--min-frames-per-model", "0")),
        min_strides_per_model=parse(Float64, arg_value(args, "--min-strides-per-model", "8")),
        n_delays=parse(Int, arg_value(args, "--n-delays", "26")),
        delay_interval=parse(Int, arg_value(args, "--delay-interval", "1")),
        noise_floor=parse(Float64, arg_value(args, "--noise-floor", "1e-6")),
        rank=parse(Int, arg_value(args, "--rank", "0")),
        n_angles=parse(Int, arg_value(args, "--n-angles", "720")),
        checkpoint_interval=parse(Int, arg_value(args, "--checkpoint-interval", "25")),
        progressbar=!has_flag(args, "--quiet-progressbar"),
        verbose=!has_flag(args, "--quiet"),
        log_interval=parse(Int, arg_value(args, "--log-interval", "1")),
        force=has_flag(args, "--force"),
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    fit_main()
end
