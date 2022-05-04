using SpecialFunctions

export hann_window, rect_window, kaiser_window, tri_window, lanczos_window

rect_window(width = 0.5) = WindowFunction(x -> 1.0, width)

hann_window(width = 1.0) = WindowFunction(
    x -> cos(pi * x / 2 / width)^2, 
    width
)

kaiser_window(width = 1.2; alpha = 3.5) = WindowFunction(
    x -> besseli(0, pi * alpha * sqrt(1 - (x / width)^2)) / besseli(0, pi * alpha),
    width
)

tri_window(width=1.) = WindowFunction(
    x -> 1 - x/width,
    width
)

function lanczos_window(lobes::Int = 3)
    @assert lobes > 0 "Lanczos Basis is only defined for at least one lobe"

    WindowFunction(
    x -> begin
        x == 1 && return 1.
        lobes*sin(pi*x)*sin(pi*x/lobes)/(pi^2*x^2)
    end,
    lobes
    )
end