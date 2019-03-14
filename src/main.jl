using Pkg
Pkg.activate(".")

# include("nPendulum.jl")
# using .nPendulum

## Parameters
# number of pendulums:
# N = 2
# ending time:
# tf = 0.1

## Solve the system
# p   = generateRandomParameters(N)
# z0  = generateRandomInitialConfig(p)
# z,t = solvePendulum(p, z0, tf, 0.01 )


include("IRK.jl")
using .IRK
using LinearAlgebra
z0 = [1.]
h = 0.1
tn = 0 + h
f(t,x) = [ -x[1] ]
# num_stages = 3
# butcherTable = IRK.GaussLegendreRKRule.buildRule(num_stages)
# iter_max = 500
# tol = 1e-8
# n  = length(z0)
# s  = butcherTable.num_stages
# K  = zeros( Float64, n, s )
# # Make an explicit guess
# xn = z0
# for i = 1:s
#     K[:,i] = f( tn+h*butcherTable.c[i], xn )
# end
#
# # fixed point iteration until convergence
# # K, flag, rel_err = fixed_point_iteration(f, xn, tn, K, h, butcherTable, iter_max, tol)
# invJ0 = I
# K, invJ0, flag, rel_err,iter = IRK.broyden_iteration(f, invJ0, xn, tn, K, h, butcherTable, iter_max, tol)


z1,t1 = IRK.IRKadapt( f, z0, 0., 0.02, 0.01)
println(z1)
println(t1)

z2,t2 = IRK.IRKfixed( f, z0, 0., 0.02, 0.01)
println(z2)
println(t2)
