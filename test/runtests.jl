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

    @testset "Basic Example 1" begin
        xin = 1.0:1.0:7.0
        yin = [1.,2.,3.,4.,5.,6.,7.]
        xout = [3.2, 5.2]
        @test regrid(xin, yin, xout) == [3.5, 5.5]
    end

    @testset "Basic Example 2" begin
        xin = 1.0:1.0:9.0
        # xin=[1, 2, 3, 4, 5, 6, 7, 8, 9]
        yin = [0, 8, 2, 7, 4, 9, 7, 8, 4]
        xout = [3.2, 6] # Δ = 2.8, Δ/2 = 1.4
        # first point: from 1.8 to 4.6 -> 2 3 | 4
        # second points: from 4.6 to 7.4 -> 5 |6 7 (middle one should be counted to the higher interval)
        @test regrid(xin, yin, xout) == [
            mean([mean(yin[2:3]), yin[4]]), 
            mean([mean(yin[5:5]), mean(yin[6:7])])
            ]
    end

end