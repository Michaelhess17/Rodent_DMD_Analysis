#!/usr/bin/env julia

include(joinpath(@__DIR__, "00_fit_rodent_dmd_metrics.jl"))
include(joinpath(@__DIR__, "01_make_rodent_figures.jl"))

fit_args = isempty(ARGS) ? [
    "--max-per-group", "80",
    "--max-per-subject", "18",
    "--n-delays", "6",
    "--n-angles", "48",
    "--min-bouts-per-model", "2",
    "--min-strides-per-model", "0",
] : ARGS

fit_main(fit_args)
figures_main()
