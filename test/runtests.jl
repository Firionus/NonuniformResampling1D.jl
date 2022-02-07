using Regrid1D
using Test
using Statistics

@testset "Regrid1D Tests" begin

    @testset "find Utitilities" begin
        @testset "find_first_above_or_equal" begin
            x = 1.0:1.0:10.0
            @test Regrid1D.find_first_above_or_equal(.5, x) == 1
            @test Regrid1D.find_first_above_or_equal(0, x) == 1
            @test Regrid1D.find_first_above_or_equal(-.5, x) == 1
            @test Regrid1D.find_first_above_or_equal(1, x) == 1
            @test Regrid1D.find_first_above_or_equal(1.1, x) == 2
            @test Regrid1D.find_first_above_or_equal(5.4, x) == 6
            @test Regrid1D.find_first_above_or_equal(9.9, x) == 10
            @test Regrid1D.find_first_above_or_equal(10, x) == 10
            @test_throws Exception Regrid1D.find_first_above_or_equal(10.1, x)
        end

        @testset "find_last_below_or_equal" begin
            x = 1.0:1.0:10.0
            @test_throws Exception Regrid1D.find_last_below_or_equal(.9, x)
            @test Regrid1D.find_last_below_or_equal(1, x) == 1
            @test Regrid1D.find_last_below_or_equal(1.1, x) == 1
            @test Regrid1D.find_last_below_or_equal(5.4, x) == 5
            @test Regrid1D.find_last_below_or_equal(9.9, x) == 9
            @test Regrid1D.find_last_below_or_equal(10, x) == 10
            @test Regrid1D.find_last_below_or_equal(10.1, x) == 10
            @test Regrid1D.find_last_below_or_equal(30, x) == 10
        end
    end

    @testset "Basis Functions" begin
        @testset "Rectangular Basis" begin
            b = RectangularBasis() # width = .5
            @test_throws Exception basis_value(b, prevfloat(0.))
            @test basis_value(b, 0) == 1.
            @test basis_value(b, prevfloat(0.5)) == 1.
            @test basis_value(b, 0.5) == 1.
            @test basis_value(b, nextfloat(0.5)) == 0.
            @test basis_value(b, 100) == 0.
        end

        @testset "Triangular Basis" begin
            b = TriangularBasis() # width = 1.
            @test_throws Exception basis_value(b, prevfloat(0.))
            @test basis_value(b, 0) == 1.
            @test basis_value(b, .2) == .8
            @test basis_value(b, .5) == .5
            @test basis_value(b, 1) == 0
        end
    end

    @testset "Basic Example 1" begin
        xin = 1.0:1.0:7.0
        yin = [1.,2.,3.,4.,5.,6.,7.]
        xout = [3.2, 5.2]
        @test regrid(xin, yin, xout, required_input_points=1) == [3.5, 5.5]
    end

    @testset "Basic Example 2" begin
        xin = 1.0:1.0:9.0
        # xin=[1, 2, 3, 4, 5, 6, 7, 8, 9]
        yin = [0, 8, 2, 7, 4, 9, 7, 8, 4]
        xout = [3.2, 6] # Δ = 2.8, Δ/2 = 1.4
        # first point: from 1.8 to 4.6 -> 2 3 | 4
        # second points: from 4.6 to 7.4 -> 5 |6 7 (middle one should be counted to the higher interval)
        @test regrid(xin, yin, xout, required_input_points=1) == [
            mean(yin[2:4]), 
            mean(yin[5:7])
            ]
    end

    @testset "Basic Example 3" begin
        yin = [1.0, 212.47605778632067, 118.55093948572791, 273.88610509333614, 
        307.67343485176576, 257.1875847477985, 357.4538672443729, 484.87523154603514, 
        400.1796499550531, 398.29054790171745, 317.0589375578311, 501.18723362727246, 
        608.6533266004108, 428.7113682873226, 869.7311001362467]
        xin = range(0.0, step=0.18310546955466883, length=length(yin))
        xout = [1.0, 1.0599939505656388]
        output = regrid(xin, yin, xout)
        @test all(.!isnan.(output))
    end

end