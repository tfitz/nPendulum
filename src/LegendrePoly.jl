module LegendrePoly

using LinearAlgebra

# setprecision(BigFloat, 512)

export P, getZeros

function P(x::BigFloat, n::Int)
    P0 = one(BigFloat)
    P1 = x

    if n == 0
        return (P0,0)
    elseif n == 1
        return (P1,1)
    else
        P = 0
        for I = 1:n-1
            i = BigInt(I)
            P = ( (2*i+1)*x*P1 - i*P0 )/(i+1)
            P0 = P1
            P1 = P
        end

        #calc dP/dx
        dP = (x*P1 - P0)*n/(x^2 - 1)

        return P, dP
    end
end
P(x::Float64, n::Int) = Float64.( P(BigFloat(x), n) )

function newton_legendraP(x0::BigFloat, n::Int64; epsilon=1e-20, iter_max = 3000)

    flag = 0
    iter = 0
    x1 = x0

    while flag == 0

        iter += 1

        p,dp = P(x1, n)

        xn = x1 - p/dp
        x1 = xn

        if abs(p) <= epsilon
            return xn

        elseif iter >= iter_max
            return NaN
        end

    end

end

function initialGuessOfZeros(n::Int)
    # This follows from the discussion in Bogaert (2014) doi:10.1137/140954969
    x = zeros(n)
    for k = 1:n
        a = pi*(k-1/4)
        j0 = a + 1/(8*a) - 124/3/(8*a)^3 + 120928/15/(8*a)^5 - 401743168/105/(8*a)^7 + 1071187749376/315/(8*a)^9
        v = 1/(n+1/2)
        α = v*j0
        θ = α + v^2*(α*cot(α) - 1)/(8*α)
        x[k] = cos(θ)
    end
    return x
end

getZeros(n::Int; tol = 1e-30, iter_max = 3000) = [ newton_legendraP( BigFloat(x0), n, epsilon=tol, iter_max=iter_max) for x0 in initialGuessOfZeros(n) ]

end
