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