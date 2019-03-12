# IRK.jl
module IRK

using LinearAlgebra
using Printf

include("GaussLegendreRKRule.jl")

export IRKadapt, IRKfixed

function IRKadapt( f::Function, y0, t0::Float64, tf::Float64, h0::Float64;
    num_stages=3,
    step_tol=1e-4, step_max=4096, h_min = 1e-10, iter_tol=1e-4, iter_max=100)

    # build the Butcher Table
    butcherTable = GaussLegendreRKRule.buildRule(num_stages)

    # step 1
    y = Array{Float64,2}(undef, 1, length(y0))
    y[:] = y0

    # time
    time = Float64[]
    append!(time,t0)

    # time stepping
    i = 1
    h = h0

    while time[end] <= tf
        i += 1

        step_flag = 0
        step_iter = 0
        while step_flag == 0
            step_iter += 1

            z,sub_flag,sub_err = irk_step( y[i-1,:], time[i-1], h, f, iter_tol, iter_max, butcherTable)

            # did step work?
            if sub_err <= step_tol
                # yes it worked.  So accept the step, and move on.
                y = vcat( y, z' )
                append!(time, time[i-1] + h)

                # move to next time step
                step_flag = 1

            elseif step_iter >= iter_max
                step_flag = -1
                println("Failed to get error under control\n")
                return

            elseif h <= h_min
                step_flag = -2
                println("Step size below minimum\n")
                return

            elseif sub_flag != 1
                step_flag = -3
                @printf("substep had errors, sub_flag = %d\n", sub_flag)
                return

            end

            # Update the step size
            # based on https://en.wikipedia.org/wiki/Adaptive_stepsize
            # TODO replace with fancier scheme
            h = 0.9*h*min( max(step_tol/sub_err, 0.3), 2 )

        end

    end

    return (y,time)
end

function IRKfixed( f::Function, y0, t0::Float64, tf::Float64, h::Float64;
    num_stages=3, step_tol=1e-4, iter_tol=1e-4, iter_max=100)

    # build the Butcher Table
    butcherTable = GaussLegendreRKRule.buildRule(num_stages)

    # time
    time = t0:h:tf

    # compute number of steps to take
    num_steps = length(t0:h:tf)

    # step 1
    y = zeros(num_steps,length(y0))
    y[1,:] = y0

    # time stepping
    for i = 2:num_steps

        z,sub_flag,sub_err = irk_step( y[i-1,:], time[i-1], h, f, iter_tol, iter_max, butcherTable)

        # check to see if the step worked
        if sub_err ≥ step_tol
            @printf("Error: step size is too large, substep error est = %e\n",sub_err)
            return
        elseif sub_flag != 1
            @printf("substep failed, sub_flag = %d\n", sub_flag)
            return
        end

        # pack in the results, and move on
        y[i,:] = z
    end

    return (y,time)
end


function irk_step( xn, tn::Float64, h::Float64, f::Function,
    tol::Float64, iter_max::Int, butcherTable )

    n  = length(xn)
    s  = butcherTable.num_stages
    K  = zeros( Float64, n, s )
    Ks = zeros( Float64, n, s )

    # TODO replace fixed-point iterations with Broyden updates

    # Make an explicit guess
    for i = 1:s
        K[:,i] = f( tn+h*butcherTable.c[i], xn )
    end

    # fixed point iteration until convergence
    flag = 0
    iter = 0
    while flag == 0
        iter += 1

        for i = 1:s
            sum_ak = zeros(Float64, n)
            for j = 1:s
                sum_ak += butcherTable.A[i,j]*K[:,j]
            end
            Ks[:,i] = f( tn + h*butcherTable.c[i], xn + h*sum_ak )
        end

        rel_err = norm( (K-Ks)[:] )
        K = copy(Ks) # don't make a lazy copy

        if rel_err <= tol
            flag = 1
        elseif iter >= iter_max
            flag = -1
        end

    end

    # Order 2*s: compute x(n+1)
    sum_bk = zeros(Float64, n)
    for i = 1:s
        sum_bk += butcherTable.b[i]*K[:,i]
    end
    y = xn + h*sum_bk

    # Order 2*s-1: estimate
    sum_bk .= 0.0
    for i = 1:s
        sum_bk += butcherTable.btilde[i]*K[:,i]
    end
    z = xn + h*sum_bk

    step_err = norm( y - z )

    # local error estimate
    return ( y, flag, step_err )

end



end
