using NonuniformResampling1D
using Plots
using Statistics

xin = 1:24
yin = [
    0.10221248850845677,
    -0.3017964617372699,
    0.6313812022052889,
    -0.469567525619184,
    -0.742550821690563,
    0.538340882717059,
    0.7645542960785265,
    -0.42847475695892334,
    -0.2030496164626363,
    -0.6961115394075117,
    0.9553667194266671,
    0.6773297980833843,
    0.6595211409431194,
    0.18762849994721864,
    0.9676546988790669,
    -0.6599170737336646,
    0.17851166014365583,
    0.7417238170162541,
    -0.7312025977875327,
    -0.4015347276663681,
    -0.685997771780074,
    -0.7558107010474338,
    0.12236047934097583,
    -0.8627888068515903]
xout = [4.6, 5.7, 8.2, 13.3, 17.3, 18.2, 19.1, 20, 20.8]
red_ind = 3
blue_ind = 6
red_color = :coral
blue_color = :dodgerblue

w = 1.35
win = hann_window(w)
N = 300
win_y = win.(range(0, w, length=N))
required_points_per_slice = 4

red_left_border = xout[red_ind]*(1-w)+w*xout[red_ind-1]
red_right_border = xout[red_ind]*(1-w)+w*xout[red_ind+1]
blue_left_border = xout[blue_ind]*(1-w)+w*xout[blue_ind-1]
blue_right_border = xout[blue_ind]*(1-w)+w*xout[blue_ind+1]

yout = nuresample(xin, yin, xout, win, required_points_per_slice=required_points_per_slice)

function vertical_lines()
    vline!([xout[red_ind]], c=red_color, linewidth=1)
    vline!([red_left_border, red_right_border], c=red_color, line=:dash)
    vline!([xout[blue_ind]], c=blue_color, linewidth=1)
    vline!([blue_left_border, blue_right_border], c=blue_color, line=:dash)
    # unit lines make the image too crowded
    #vline!([xout[red_ind-1], xout[red_ind+1], xout[blue_ind-1], xout[blue_ind+1]], c=:black, line=:dot)
    hline!([0], c=:gray, z_order=:back)
end

# blue upsampling + points
upsample_step = min(xout[blue_ind]-xout[blue_ind-1], xout[blue_ind+1]-xout[blue_ind])/
    required_points_per_slice*w
upsample_x = range(blue_left_border+upsample_step/2, step=upsample_step, 
stop=blue_right_border-upsample_step/2)
in_step = Float64(StepRangeLen(xin).step)
upsample_y = [NonuniformResampling1D.interpolate_point(
    StepRangeLen(xin), yin, up_x, in_step, in_step, lanczos_window()) for up_x in upsample_x]

# assume xout around blue is symmetrical
blue_weights = @. win(abs(upsample_x-xout[blue_ind])/(xout[blue_ind]-xout[blue_ind-1]))
blue_weighted = blue_weights .* upsample_y

plt_in = begin
    plot(xin, yin, m=:circle, c=:black, legend=false, xaxis=false)
    # red window
    red_x_left = range(stop=xout[red_ind], length=N, step=(xout[red_ind]-xout[red_ind-1])/N*w)
    plot!(red_x_left, reverse(win_y), c=red_color)
    red_x_right = range(start=xout[red_ind], length=N, step=(xout[red_ind+1]-xout[red_ind])/N*w)
    plot!(red_x_right, win_y, c=red_color)
    #blue window
    blue_x_left = range(stop=xout[blue_ind], length=N, step=(xout[blue_ind]-xout[blue_ind-1])/N*w)
    plot!(blue_x_left, reverse(win_y), c=blue_color)
    blue_x_right = range(start=xout[blue_ind], length=N, step=(xout[blue_ind+1]-xout[blue_ind])/N*w)
    plot!(blue_x_right, win_y, c=blue_color)
    plot!(upsample_x, upsample_y, c=:deepskyblue, m=:x, markersize=2)

    annotate!([
        (4, -1.1, ("adaptive, asymmetric\nwindow width", 9, red_color)),
        (21.5, .8, ("ad-hoc\nupsampling", 9, blue_color))
    ])

    vertical_lines()
    title!("Input", titlefontsize=13)
end

plt_windowed = begin
    red_start_ind = NonuniformResampling1D.find_first_above_or_equal(red_left_border, StepRangeLen(xin))
    red_middle_ind = NonuniformResampling1D.find_last_below(xout[red_ind], StepRangeLen(xin))
    red_stop_ind = NonuniformResampling1D.find_last_below_or_equal(red_right_border, StepRangeLen(xin))

    red_left_x = xin[red_start_ind:red_middle_ind]
    red_right_x = xin[red_middle_ind+1:red_stop_ind]

    red_left_weights = @. win(abs(red_left_x-xout[red_ind])/(xout[red_ind]-xout[red_ind-1]))
    red_right_weights = @. win(abs(red_right_x-xout[red_ind])/(xout[red_ind+1]-xout[red_ind]))

    red_y = yin[red_start_ind:red_stop_ind]
    red_weighted = [red_left_weights; red_right_weights] .* red_y
    plot(red_start_ind:red_stop_ind, red_weighted, 
        m=:circle, c=red_color, legend=false, markerstrokewidth=0)

    plot!(upsample_x, blue_weighted, m=:circle, c=blue_color, markerstrokewidth=0)

    vertical_lines()
    title!("Example Slices Windowed", titlefontsize=13)
end

plt_out = begin
    plot(xout, yout, m=:circle, c=:black, legend=false)
    scatter!([xout[red_ind]], [yout[red_ind]], c=red_color, markerstrokewidth=0)
    scatter!([xout[blue_ind]], [yout[blue_ind]], c=blue_color, markerstrokewidth=0)
    
    quiver!([xout[red_ind], xout[red_ind]], [-.5, -1.0], quiver=(-[(xout[red_ind]-xout[red_ind-1]), (xout[red_ind]-red_left_border)], [0,0]), c=:black)
    annotate!([
        (mean(xout[red_ind-1:red_ind]), -.7, ("Unit", 9)),
        (mean([xout[red_ind], red_left_border]), -1.2, ("Slice", 9))
    ])

    vertical_lines()
    title!("Output", titlefontsize=13)
end

plot(plt_in, plt_windowed, plt_out, layout=(3,1), link=:x, xticks=false, xaxis=false, yaxis=false)
ylims!(-1.1, 1.1)
