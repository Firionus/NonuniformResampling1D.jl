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

    @testset "Lanczos Window" begin
        w = lanczos_window() # lobes=3, width=1
        @test w(0) == 1
        @test w(.5) > .5
        @test w(1) ≈ 0 atol=1e-15
        @test w(1.5) < -.1
        @test w(2) ≈ 0 atol=1e-15
        @test w(2.5) > .01
        @test w(3) ≈ 0 atol=1e-15

        w = lanczos_window(3, width = 2) # lobes=3, width=2
        @test w(0) == 1
        @test w(1) > .5
        @test w(2) ≈ 0 atol=1e-15
        @test w(3) < -.1
        @test w(4) ≈ 0 atol=1e-15
        @test w(5) > .01
        @test w(6) ≈ 0 atol=1e-15
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