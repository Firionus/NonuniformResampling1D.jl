# Example Default Data
xin = 1.:12.
yin = [9, 7, 8, 7, 1, 2, 4, 1, 9, 8, 5, 3]
#      1--2--3--4--5--6--7--8--9--10-11-12

@testset "Basic Example 1" begin
    xin = 1.0:1.0:7.0
    yin = [1.,2.,3.,4.,5.,6.,7.]
    xout = [3.2, 5.2]
    @test regrid(xin, yin, xout, required_points_per_slice=1) == [3.5, 5.5]
end

@testset "Basic Example 2" begin
    xin = 1.0:1.0:9.0
    # xin=[1, 2, 3, 4, 5, 6, 7, 8, 9]
    yin = [0, 8, 2, 7, 4, 9, 7, 8, 4]
    xout = [3.2, 6] # Δ = 2.8, Δ/2 = 1.4
    # first point: from 1.8 to 4.6 -> 2 3 | 4
    # second points: from 4.6 to 7.4 -> 5 |6 7 
    # (middle one should be counted to the higher interval)
    @test regrid(xin, yin, xout, required_points_per_slice=1) == [
        mean(yin[2:4]), 
        mean(yin[5:7])
        ]
end

@testset "Example that Should not Cause NaN Result" begin
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
    xout = [2.9, 3.8, 7.5] # middle point highly asymmetric
    # yin=[9, 7, 8, 7, 1, 2, 4, 1, 9, 8, 5, 3]
    #      1--2-*3-*4--5--6--7-*8--9--10-11-12  (output points as "*")
    #            ll*rrrrrrrrrrr                 (point under test)
    
    # without upsampling
    output = regrid(xin, yin, xout, RectangularBasis(1.), 
        required_points_per_slice=1)
    @test output[2] == mean(yin[3:7])

    # now with upsampling
    output_asymmetric_upsampling = regrid(xin, yin, xout, RectangularBasis(1.), 
        required_points_per_slice=2, upsampling_basis=TriangularBasis(1.))
    # expected: left slice only contains 1 point, but 2 are required
    # so both slices should be upsampled with a step of 0.9/2 = 0.45
    # upsampling positions therefore are: [3.125, 3.575] (left) and 
    # [4.025, 4.475, 4.925, 5.375, 5.825, 6.275, 6.725, 7.175] (right)
    # the mean of all those points, linearly interpolated, should be the result
    itp = interpolate(yin, BSpline(Linear()))
    expected_positions = [
        3.125, 3.575, 4.025, 4.475, 4.925, 5.375, 5.825, 6.275, 6.725, 7.175]
    @test mean(itp(expected_positions)) ≈ output_asymmetric_upsampling[2]

    # relative change from enabling upsampling should be smaller than 8.2 %
    # (value that is just achieved with implementation at time of writing test)
    rel_change = abs(output_asymmetric_upsampling[2] - output[2])/output[2]
    @test rel_change < .082
    
    # compare with large oversampling
    big_oversampling = regrid(xin, yin, xout, RectangularBasis(1.), 
    required_points_per_slice=9, upsampling_basis=TriangularBasis(1.))[2]
    # relative difference between 2x and 9x oversampling should be smaller 0.96 %
    # (value that is just achieved with implementation at time of writing test)
    rel_change = abs(output_asymmetric_upsampling[2] - big_oversampling)/
        big_oversampling
    @test rel_change < .0096
end

@testset "Large Oversampling Value Should Not Error" begin
    xout = [2.9, 3.8, 7.5] # middle point highly asymmetric

    result8 = regrid(xin, yin, xout, RectangularBasis(1.), 
    required_points_per_slice=8, upsampling_basis=TriangularBasis(1.))[2]

    result32 = regrid(xin, yin, xout, RectangularBasis(1.), 
    required_points_per_slice=32, upsampling_basis=TriangularBasis(1.))[2]

    # results should be very comparable
    @test abs(result8 - result32)/result32 < .001
end

@testset "Output Value Does Not Jump when Output Point Moves over Input Point" begin
    # output point left of index 6
    output_left = regrid(xin, yin, [3.9, 5.9, 7.9], RectangularBasis(1.), 
        required_points_per_slice=1)
    # output point right of index 6
    output_right = regrid(xin, yin, [3.9, 6.1, 7.9], RectangularBasis(1.), 
        required_points_per_slice=1)
    @test output_left[2] == output_right[2]
end

@testset "Error When not Enough Points at Beginning or End" begin
    # Just enough points at beginning
    regrid(xin, yin, [2., 3.9999999999999], RectangularBasis(1.), 
        required_points_per_slice=1)

    # not enough points (the point at index -1 would be needed but does not exist)
    @test_throws Exception regrid(xin, yin, [2., 4.], RectangularBasis(1.), 
        required_points_per_slice=1)

    # Just enough points at end
    regrid(xin, yin, [9.00000000000001, 11.], RectangularBasis(1.), 
        required_points_per_slice=1)

    # not enough points (the point at index 13 would be needed but does not exist)
    @test_throws Exception regrid(xin, yin, [9., 11.], RectangularBasis(1.), 
        required_points_per_slice=1)
end

@testset "xout not Monotonically Increasing Throws Error" begin
    # not monotonic
    @test_throws(
        ArgumentError("xout must be increasing everywhere. Violated between index 2 and 3: 6.5 >= 6.0"),
        regrid(xin, yin, [4.5, 6.5, 6., 8.5]))

    # monotonic, but not strict
    @test_throws(
        ArgumentError("xout must be increasing everywhere. Violated between index 2 and 3: 6.5 >= 6.5"),
        regrid(xin, yin, [4.5, 6.5, 6.5, 8.5]))

    # not monotonic at beginning
    @test_throws(
        ArgumentError("xout must be increasing everywhere. Violated between index 1 and 2: 4.5 >= 3.5"),
        regrid(xin, yin, [4.5, 3.5, 6., 8.5]))
    
    # not strict monotonic at beginning
    @test_throws(
        ArgumentError("xout must be increasing everywhere. Violated between index 1 and 2: 4.5 >= 4.5"),
        regrid(xin, yin, [4.5, 4.5, 6., 8.5]))
    
    # not monotonic at end
    @test_throws(
        ArgumentError("xout must be increasing everywhere. Violated between index 3 and 4: 6.0 >= 5.5"),
        regrid(xin, yin, [4.5, 5., 6., 5.5]))

    # not strict monotonic at end
    @test_throws(
        ArgumentError("xout must be increasing everywhere. Violated between index 3 and 4: 6.0 >= 6.0"),
        regrid(xin, yin, [4.5, 5., 6., 6.]))
end

@testset "Required Input Point per Slice" begin
    # The default for required_points_per_slice should be
    # required_points_per_unit = 4
    # required_points_per_slice = round(required_points_per_unit * smoothing_kernel.width)

    # So, with the default kernel of Rect(.5), at least 2 points have to be in the slice. 

    # 1 2 3 4 5 6 7 8 9 10 11 12
    #    *   | *   |
    # (.9:4.1) (4.1:7.3)
    # 2    2|1   2 -> no upsampling on first point, but on second point
    result = regrid(xin, yin, [2.5, 5.7], upsampling_basis=TriangularBasis(1.))
    @test result[1] == mean(yin[1:4])
    itp = interpolate(yin, BSpline(Linear()))
    expected_interpolation_positions = [
        4.1+1.6/2*.5, 4.1+1.6/2*1.5, 5.7+.5*1.6/2, 5.7+1.5*1.6/2
    ]
    @test result[2] == mean(itp(expected_interpolation_positions))

    # When using a larger kernel, more points need to be in the slice
    # e.g. Rect(.75) -> 3 points have to be in the slice
    # 1 2 3 4 5 6 7 8 9 10 11 12
    #      *     *    
    # (1.1:5.9) (4.3:9.1)
    # 2    2|2   3 -> upsampling on both points
    result = regrid(xin, yin, [3.5, 6.7], RectangularBasis(.75),
        upsampling_basis=TriangularBasis(1.))
    itp = interpolate(yin, BSpline(Linear()))
    s = 2.4/3 # step
    expected_interpolation_positions1 = [
        1.1+s*.5, 1.1+s*1.5, 1.1+s*2.5, 3.5+s*.5, 3.5+s*1.5, 3.5+s*2.5
    ]
    @test result[1] == mean(itp(expected_interpolation_positions1))
    expected_interpolation_positions2 = [
        4.3+s*.5, 4.3+s*1.5, 4.3+s*2.5, 6.7+s*.5, 6.7+s*1.5, 6.7+s*2.5
    ]
    @test result[2] ≈ mean(itp(expected_interpolation_positions2))

    # with "disabled" interpolation
    no_interpolation = regrid(xin, yin, [2.5, 5.7], required_points_per_slice = 1)
    @test no_interpolation[1] == mean(yin[1:4])
    @test no_interpolation[2] == mean(yin[5:7])
end

@testset "Strange required_points_per_slice Values" begin
    # disallow values smaller than 0, because slices without values in them do
    # not contribute anything, so the value of a point might become undefined
    # without upsampling, which then is never triggered. If high performance
    # like this is desired, required_points_per_slice=1 should be used with
    # nearest neighbor upsampling, i.e. upsampling_basis=Rect(.5)
    @test_throws AssertionError regrid(xin, yin, [4.5, 7.7], required_points_per_slice=0)
    @test_throws AssertionError regrid(xin, yin, [4.5, 7.7], required_points_per_slice=-1)
    # disallow float values, because they result in discontinuous values with
    # the upsampling algorithm, i.e. the outputs values with 4.5 are very
    # different from 4 or 5
    @test_throws TypeError regrid(xin, yin, [4.5, 7.7], required_points_per_slice=4.5)
end

# TODO test that basis=Rect(0.), upsampling_basis=Rect(.5) and required_points_per_slice=1
# results in nearest neighbor behavior