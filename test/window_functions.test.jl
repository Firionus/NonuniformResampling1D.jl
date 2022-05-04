include("fixtures.jl")

@testset "Window Functions" begin
    @testset "Rectangular Window" begin
        b = rect_window() # width = .5
        @test_throws Exception b(prevfloat(0.))
        @test b(0) == 1.
        @test b(prevfloat(0.5)) == 1.
        @test b(0.5) == 1.
        @test b(nextfloat(0.5)) == 0.
        @test b(100) == 0.
    end

    @testset "Triangular Window" begin
        b = tri_window() # width = 1.
        @test_throws Exception b(prevfloat(0.))
        @test b(0) == 1.
        @test b(.2) == .8
        @test b(.5) == .5
        @test b(1) == 0
    end

    @testset "Invalid width" begin
        @test_throws AssertionError rect_window(0)
        @test_throws AssertionError rect_window(0.)
    end

    @testset "API Example for Predefined Window Functions" begin
        nuresample(xin, yin, [3.3, 4.1], rect_window(.7), 
            required_points_per_slice=1, upsampling_function=rect_window(1.))
    end
end