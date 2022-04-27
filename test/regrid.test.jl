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

@testset "Balanced Average for Asymmetric Windows" begin
    yin = [9, 7, 8, 7, 1, 2, 4, 1, 9, 8, 5, 3]
    #      1--2-*3-*4--5--6--7-*8--9--10-11-12      (input axis and output points as "*")
    #            ll*rrrrrrrrrrr                     (point under test)
    xin = 1.:12.
    xout = [2.9, 3.8, 7.5] # middle point highly asymmetric
    output = regrid(xin, yin, xout, RectangularBasis(1.), required_input_points=1)
    @test output[2] == mean(yin[3:7]) # value without upsampling on one side is just a mean

    output_asymmetric_upsampling = regrid(xin, yin, xout, RectangularBasis(1.), 
        required_input_points=2, upsampling_basis=TriangularBasis(1.))

    # expected: left slice only contains 1 point, but 2 are required
    # so both slices should be upsampled with a step of 0.9/2 = 0.45
    # upsampling positions therefore are: [3.125, 3.575] (left) and 
    # [4.025, 4.475, 4.925, 5.375, 5.825, 6.275, 6.725, 7.175] (right)
    # the mean of all those points, linearly interpolated, should be the result
    itp = interpolate(yin, BSpline(Linear()))
    expected_positions = [3.125, 3.575, 4.025, 4.475, 4.925, 5.375, 5.825, 6.275, 6.725, 7.175]

    @test mean(itp(expected_positions)) ≈ output_asymmetric_upsampling[2]

    # relative change from upsampling should be smaller than 10 %
    @test abs(output_asymmetric_upsampling[2] - output[2])/output[2] < .1

    # compare with large oversampling
    big_oversampling = regrid(xin, yin, xout, RectangularBasis(1.), 
    required_input_points=9, upsampling_basis=TriangularBasis(1.))[2]
    @test abs(output_asymmetric_upsampling[2] - big_oversampling)/big_oversampling < .01
end

@testset "Large Oversampling Value Should Not Error" begin
    # TODO same example as above - consolidate?
    yin = [9, 7, 8, 7, 1, 2, 4, 1, 9, 8, 5, 3]
    xin = 1.:12.
    xout = [2.9, 3.8, 7.5] # middle point highly asymmetric

    result8 = regrid(xin, yin, xout, RectangularBasis(1.), 
    required_input_points=8, upsampling_basis=TriangularBasis(1.))[2]

    result32 = regrid(xin, yin, xout, RectangularBasis(1.), 
    required_input_points=32, upsampling_basis=TriangularBasis(1.))[2]

    # results should be very comparable
    @test abs(result8 - result32)/result32 < .001
end

@testset "Output Value Does Not Jump when Output Point Moves over Input Point" begin
    yin = [9, 9, 2, 1, 6, 4, 7, 3, 6, 4, 7]
    xin = 1.:11.
    # output point left of index 6
    output_left = regrid(xin, yin, [3.9, 5.9, 7.9], RectangularBasis(1.), required_input_points=1)
    # output point right of index 6
    output_right = regrid(xin, yin, [3.9, 6.1, 7.9], RectangularBasis(1.), required_input_points=1)
    @test output_left[2] == output_right[2]
end

@testset "Error When Not Enough Points at Beginning or End" begin
    # TODO same example as above - consolidate?
    yin = [9, 7, 8, 7, 1, 2, 4, 1, 9, 8, 5, 3]
    xin = 1.:12.

    # Just enough points at beginning
    regrid(xin, yin, [2., 3.9999999999999], RectangularBasis(1.), required_input_points=1)

    # not enough points (the point at index -1 would be needed but does not exist)
    @test_throws Exception regrid(xin, yin, [2., 4.], RectangularBasis(1.), required_input_points=1)

    # Just enough points at end
    regrid(xin, yin, [9.00000000000001, 11.], RectangularBasis(1.), required_input_points=1)

    # not enough points (the point at index 13 would be needed but does not exist)
    @test_throws Exception regrid(xin, yin, [9., 11.], RectangularBasis(1.), required_input_points=1)
end