using Pkg
Pkg.activate(".")

include("nPendulum.jl")
using .nPendulum
include("nPendulumAnimation.jl")
using .nPendulumAnimation
using Plots

## Parameters
# number of pendulums:
N = 8
# ending time:
tf = 1.
# frames per second in the animation
fps = 24
# number of frames in the animation
num_frames = tf*fps*2 |> round |> Int


## Solve the system
p   = generateRandomParameters(N)
z0  = generateRandomInitialConfig(p)

@elapsed z,t = solvePendulum(p,z0,tf,0.01, method=:broyden)

@elapsed z,t = solvePendulum(p,z0,tf,0.01, method=:fixedpoint)


## Generate the animation
# anim = GenPendulumAnimation(t,z,num_frames,p)
# gif(anim, "./anim.gif", fps = fps)
