#!/usr/bin/env julia

include(joinpath(@__DIR__, "..", "src", "RodentDMDAnalysis.jl"))
using .RodentDMDAnalysis

function figures_main(args=ARGS)
    RodentDMDAnalysis.make_main_figures()
end

if abspath(PROGRAM_FILE) == @__FILE__
    figures_main()
end
