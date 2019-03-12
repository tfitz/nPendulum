using Pkg
Pkg.activate(".")

include("nPendulum.jl")
using .nPendulum

include("nPendulumAnimation.jl")
using .nPendulumAnimation

using Dierckx: Spline1D
using Plots

## Parameters
# number of pendulums:
N = 4
# ending time:
tf = 10.0
# number of frames in the animation
num_frames = 800
# frames per second in the animation
fps = 24


## Solve the system
p   = generateRandomParameters(N)
z0  = generateRandomInitialConfig(p)
z,t = solvePendulum(p, z0, tf, 0.01 )


## overall plot
T,Z = nPendulumAnimation.sampleZ(t,z, N, 400)
plt = plot( Z[:,1], Z[:,2] )
for i = 3:2:2N
    plot!(plt, Z[:,i], Z[:,i+1] )
end
plot!(plt, aspect_ratio=1)


## Generate the animation
anim = GenPendulumAnimation(t,z,num_frames,p)
gif(anim, "./anim.gif", fps = fps)
