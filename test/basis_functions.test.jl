@testset "Basis Functions" begin
    @testset "Rectangular Basis" begin
        b = RectangularBasis() # width = .5
        @test_throws Exception b(prevfloat(0.))
        @test b(0) == 1.
        @test b(prevfloat(0.5)) == 1.
        @test b(0.5) == 1.
        @test b(nextfloat(0.5)) == 0.
        @test b(100) == 0.
    end

    @testset "Triangular Basis" begin
        b = TriangularBasis() # width = 1.
        @test_throws Exception b(prevfloat(0.))
        @test b(0) == 1.
        @test b(.2) == .8
        @test b(.5) == .5
        @test b(1) == 0
    end

    @testset "Invalid width" begin
        @test_throws AssertionError RectangularBasis(0)
        @test_throws AssertionError RectangularBasis(0.)
    end
end