using SpecialFunctions

export hann_window, rect_window, kaiser_window, tri_window, lanczos_window

"""
    rect_window(width = 0.5)

Rectangular window which is 1 from 0 to `width` and 0 otherwise. 

Can be used to calculate a moving average, for example.
"""
rect_window(width = 0.5) = WindowFunction(x -> 1.0, width)

"""
    hann_window(width = 1)

Hann window (raised cosine) that is 1 at 0 and reaches 0 at `width`.
"""
hann_window(width = 1.0) = WindowFunction(
    x -> cos(pi * x / 2 / width)^2, 
    width
)

"""
    kaiser_window(width = 1.2; alpha = 3.5)

Kaiser window cut off at `width` with the shape parameter `alpha`.
"""
kaiser_window(width = 1.2; alpha = 3.5) = WindowFunction(
    x -> besseli(0, pi * alpha * sqrt(1 - (x / width)^2)) / besseli(0, pi * alpha),
    width
)

"""
    tri_window(width = 1)

Triangular window which is 1 at 0 and reaches 0 at `width`. 

Can be used to perform linear interpolation when used as upsampling_function
with `width=1`. 
"""
tri_window(width=1.) = WindowFunction(
    x -> 1 - x/width,
    width
)

"""
    lanczos_window(lobes = 3; width = 1)

Lanczos window which is 1 at 0 and 0 at 1*`width`, 2*`width`, 3*`width`, ...

The higher the number of lobes, the higher the accuracy of sinc approximation. 
"""
function lanczos_window(lobes::Int = 3; width=1.)
    @argcheck lobes > 0 "Lanczos window is only defined for at least one lobe"

    WindowFunction(
    x -> begin
        x == 0 && return 1.
        x = x/width
        lobes*sin(pi*x)*sin(pi*x/lobes)/(pi^2*x^2)
    end,
    lobes*width
    )
end