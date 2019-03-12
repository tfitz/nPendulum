# UdwadiaKalaba.jl
module UdwadiaKalaba

using ForwardDiff
using LinearAlgebra

export ukEqnSS

function uk_qdd( q, qd, t::Float64, M::Array{Float64,2}, Qex::Array{Float64,1}, ϕ::Function, n_c::Integer )

    # build mass square root
    Ms = real(sqrt(M))

    # build holonomic constraint matrices
    if n_c == 1
        A = Array( ( ForwardDiff.gradient( ϕ, q) )' )
        b = - qd'*ForwardDiff.hessian(ϕ,q)*qd
    else
        A = ForwardDiff.jacobian( ϕ, q )
        b = zeros(n_c,1)
        for i = 1:n_c
            f = x -> ϕ(x)[i]
            H = ForwardDiff.hessian(f, q)
            b[i] = - qd'*H*qd
        end
    end

    # build contraint force
    Qc = Ms*((A/Ms)\( b - A*(M\Qex) ))

    qdd = M\(Qex + Qc)
    return qdd
end

function ukEqnSS( t::Float64, z, Mass::Function, Qex::Function, ϕ::Function, numHolonomicContraints::Integer, n::Integer )
    q   = z[1:n]
    qd  = z[n+1:2*n]
    M = Mass(q)
    qdd = uk_qdd( q, qd, t, M, Qex(q,t), ϕ, numHolonomicContraints )

    return vcat( qd, qdd )
end

end
