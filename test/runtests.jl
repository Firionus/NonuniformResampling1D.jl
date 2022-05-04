using NonuniformResampling1D
using Test
using Statistics
using Interpolations

@testset "NonuniformResampling1D Tests" begin
    include("range_utilities.test.jl")
    include("window_functions.test.jl")
    include("nuresample.test.jl")
end