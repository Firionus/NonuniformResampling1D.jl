@testset "Range Utitilities" begin
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