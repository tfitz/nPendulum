
module nPendulum

include("UdwadiaKalaba.jl")
include("IRK.jl")
using LinearAlgebra

export generateRandomParameters, generateRandomInitialConfig, solvePendulum

struct parameters
    N::Int
    m::Vector
    L::Vector
    g::Float64
end

function generateRandomParameters(N::Int; mass_range=(0.5, 3), length_range=( 0.2, 1), g=9.81 )
    mass_max = maximum(mass_range)
    mass_min = minimum(mass_range)
    m = (mass_max - mass_min)*rand(N) .+ mass_min

    L_max = maximum(length_range)
    L_min = minimum(length_range)
    L = (L_max - L_min)*rand(N) .+ L_min

    return parameters(N, m, L, g)
end

function generateRandomInitialConfig(params::parameters)
    N = params.N
    θ = 2π*rand(N)
    pos(θ,l) = l*(sin(θ)*[1,0] - cos(θ)*[0,1])

    z0 = pos( θ[1], params.L[1] )
    z = copy(z0)
    for i = 2:N
        y = pos( θ[i], params.L[i] ) + z0
        append!(z, y)
        z0 = copy( y )
    end

    append!(z, zeros(2*N) )

    return z
end

function constMass( p::parameters )
    # place the masses on the diagonal
    # diag([m1,m1,m2,m2,m3,m3, ... ])

    m=p.m
    M = zeros(2*p.N)
    for i=1:p.N
        j = 2*i-1
        M[j] = m[i]
        M[j+1] = m[i]
    end

    return diagm( 0 => M )
end

function extForceOfGravity( p::parameters )
    Qex = Array{Float64,1}()
    for i = 1:p.N
        append!(Qex, [0., -p.m[i]*p.g ])
    end
    return Qex
end

function constraintFunction( q, p::parameters )
    L = p.L
    N = p.N

    # extract x's and y's
    x = q[1:2:2N]
    y = q[2:2:2N]

    # build the constraints for 2:N.
    # This avoids issues that ForwardDiff.jacobian was having when
    # \phi=zeros(N) was used to initialize storage.
    ϕ = [ ( x[i] - x[i-1] )^2 + ( y[i] - y[i-1] )^2 - L[i]^2 for i in 2:N ]
    # add on #1
    append!(ϕ, (x[1] - 0)^2 + (y[1] - 0)^2 - L[1]^2 )

    return ϕ
end

function solvePendulum(p::parameters, z0::Array{Float64,1}, tf::Float64, h::Float64)

    Mass = constMass(p)
    Qex = extForceOfGravity( p )
    F(t,z) = UdwadiaKalaba.ukEqnSS( t, z, q -> Mass, (t,q)-> Qex,
                q->constraintFunction(q,p),
                p.N, 2*p.N )

    Z,t = IRK.IRKadapt( F, z0, 0.0, tf, h)

    return Z,t

end


end
