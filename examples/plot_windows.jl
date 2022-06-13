using NonuniformResampling1D
using Plots

x = range(0, 1.1, length=300)

left_part = begin
    plot(x, rect_window(.5).(x), label="rect_window(0.5)")
    plot!(x, tri_window(1.).(x), label="tri_window(1.0)")
    plot!(x, hann_window(1.).(x), label="hann_window(1.0)")
    plot!(x, lanczos_window(1).(x), label="lanczos_window(1)")
    title!("Simple Windows")
end

x = range(0, 1.3, length=300)

right_part = begin
    plot(x, kaiser_window(1.2, alpha=1.5).(x), label="kaiser_window(1.2, alpha = 1.0)")
    plot!(x, kaiser_window(1.2, alpha=3.5).(x), label="kaiser_window(1.2, alpha = 3.5)")
    plot!(x, kaiser_window(1.2, alpha=5.5).(x), label="kaiser_window(1.2, alpha = 5.0)")
    title!("Kaiser Window")
    #xlabel!("Unit widths")
end

x = range(0, 3.1, length=300)

lanczos_part = begin
    plot(x, lanczos_window(3).(x), label="lanczos_window(3)")
    plot!(x, lanczos_window(2).(x), label="lanczos_window(2)")
    plot!(x, lanczos_window(2, width=1.5).(x), label="lanczos_window(2, width = 1.5)")
    title!("Lanczos Window")
end

plot(left_part, right_part, lanczos_part, layout=(1,3), link=:y,
size=(1200, 450), bottom_margin = 5Plots.mm)
xlabel!("Unit widths")

savefig("windows.svg")