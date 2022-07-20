using NonuniformResampling1D
using Test
using Statistics
using Interpolations
using Aqua

@testset "NonuniformResampling1D Tests" begin
    include("range_utilities.test.jl")
    include("window_functions.test.jl")
    include("nuresample.test.jl")

    @testset "Aqua.jl" begin
        Aqua.test_all(NonuniformResampling1D, ambiguities=(recursive=false))
    end
end