module nPendulumAnimation

using Dierckx: Spline1D
using Plots
using Printf

export GenPendulumAnimation

# from matplotlib: tab20 https://github.com/matplotlib/matplotlib/blob/a55909ae86049fc2e437c58dbfe7550b7990b687/lib/matplotlib/_cm.py#L1279
const tab20_data = [
    0.12156862745098039 0.4666666666666667  0.7058823529411765  ;  # 1f77b4
    0.6823529411764706  0.7803921568627451  0.9098039215686274  ;  # aec7e8
    1.0                 0.4980392156862745  0.054901960784313725;  # ff7f0e
    1.0                 0.7333333333333333  0.47058823529411764 ;  # ffbb78
    0.17254901960784313 0.6274509803921569  0.17254901960784313 ;  # 2ca02c
    0.596078431372549   0.8745098039215686  0.5411764705882353  ;  # 98df8a
    0.8392156862745098  0.15294117647058825 0.1568627450980392  ;  # d62728
    1.0                 0.596078431372549   0.5882352941176471  ;  # ff9896
    0.5803921568627451  0.403921568627451   0.7411764705882353  ;  # 9467bd
    0.7725490196078432  0.6901960784313725  0.8352941176470589  ;  # c5b0d5
    0.5490196078431373  0.33725490196078434 0.29411764705882354 ;  # 8c564b
    0.7686274509803922  0.611764705882353   0.5803921568627451  ;  # c49c94
    0.8901960784313725  0.4666666666666667  0.7607843137254902  ;  # e377c2
    0.9686274509803922  0.7137254901960784  0.8235294117647058  ;  # f7b6d2
    0.4980392156862745  0.4980392156862745  0.4980392156862745  ;  # 7f7f7f
    0.7803921568627451  0.7803921568627451  0.7803921568627451  ;  # c7c7c7
    0.7372549019607844  0.7411764705882353  0.13333333333333333 ;  # bcbd22
    0.8588235294117647  0.8588235294117647  0.5529411764705883  ;  # dbdb8d
    0.09019607843137255 0.7450980392156863  0.8117647058823529  ;  # 17becf
    0.6196078431372549  0.8549019607843137  0.8980392156862745;    # 9edae5
]

const colorSet = [ RGBA( tab20_data[i,1], tab20_data[i,2], tab20_data[i,3], mod(i,2)==1 ? 1.0 : 0.7 ) for i in 1:size(tab20_data,1) ]

function sampleZ( t, z, N, num_pts )

    # Evenly sample the data using interpolation
    T = LinRange(0, t[end], num_pts )
    Z = zeros(length(T), 2N)
    for i = 1:2N
        spl = Spline1D(t, z[:,i] )
        Z[:,i] = spl(T)
    end

    return T,Z
end


function GenPendulumAnimation( t, z, num_pts, params)
    N = params.N
    T,Z = sampleZ(t,z,N,num_pts)

    x = Z[:,1:2:2N]
    y = Z[:,2:2:2N]

    # find the bounding box
    lim_x = round(1.1*maximum( abs.(x[:]) ), digits=2)
    max_y = round(1.1*maximum([maximum(y[:]),0.05]), digits=2)
    min_y = round(1.1*minimum(y[:]), digits=2)

    # Set how long the traces are
    h = T[2] - T[1]
    n1 = round(Int, 2*1/h)

    anim = @animate for i = 1:length(T)

        # plot the pendulum bars
        plot( [0, x[i,1] ], [ 0, y[i,1] ],
            linewidth = 2,
            linecolor = colorSet[1],
            label=@sprintf("%d: m=%.2f, L=%.2f", 1, params.m[1], params.L[1] )
            )
        # plot the trace
        trace_idx = i < n1+1 ? (1:i) : (i-n1:i)
        scatter!( x[trace_idx,1], y[trace_idx,1],
            markercolor=colorSet[2],
            markersize = 3,
            markershape = :circle,
            markerstrokewidth = 0,
            label="" )

        # plot the mass
        scatter!( [x[i,1]], [y[i,1]],
            markercolor=colorSet[1],
            markersize = 4*params.m[1],
            markerstrokewidth = 1,
            label=""  )

        # repeat for all the other pendulums
        for j = 2:N
            plot!( x[i,j-1:j], y[i,j-1:j] ,
                linewidth = 2,
                linecolor = colorSet[2*j-1],
                label=@sprintf("%d: m=%.2f, L=%.2f", j, params.m[j], params.L[j] )
                )
            scatter!( x[trace_idx,j], y[trace_idx,j],
                markercolor=colorSet[2*j],
                markersize = 3,
                markershape = :circle,
                markerstrokewidth = 0,
                label="" )
            scatter!( [x[i,j]], [y[i,j]],
                markercolor=colorSet[2*j-1],
                markersize = 4*params.m[j],
                markerstrokewidth = 1,
                label="" )
        end

        # cleanup the plot
        plot!(aspect_ratio = 1, xlabel="x", ylabel="y",
            xlim=(-lim_x,lim_x), ylim=(min_y,max_y),
            legend=false )
    end

    return anim

end


end
