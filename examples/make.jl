import Pkg

Pkg.activate(@__DIR__)
Pkg.instantiate()

include(joinpath(@__DIR__, "explanation_diagram.jl"))
include(joinpath(@__DIR__, "plot_windows.jl"))