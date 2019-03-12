# nPendulum
by T. Fitzgerald

## Purpose
For fun.  This was a project I made to get to know [Julia](https://julialang.org/) better, and interact with some new technologies.  This code sets up the equations of motion of $N$ simple pendulums using the [Udwadia-Kalaba equation](https://en.wikipedia.org/wiki/Udwadia%E2%80%93Kalaba_equation).  This set of $4N$ equations are integrated using a [Gauss-Legendre](https://en.wikipedia.org/wiki/Gauss%E2%80%93Legendre_method) implicit Runge-Kutta method.

It is not intended for production use.

## Organization of the code
Navigate to the `src` folder, and run `julia main.jl` to automatically make the animation.

- `LegendrePoly.jl` builds the [Legendre polynomials](https://en.wikipedia.org/wiki/Legendre_polynomials), and finds theirs roots using [arbitrary precision](https://docs.julialang.org/en/v1/manual/integers-and-floating-point-numbers/index.html#Arbitrary-Precision-Arithmetic-1).
- `GaussLengdreRKRule.jl` builds the coefficients of an $s$-stage Gauss-Legendre Runge-Kutta method, and populates these in a Butcher table.  The method to construct these coefficients is neatly laid out in [Rang (2014)](http://www.digibib.tu-bs.de/?docid=00055783).
- `IRK.jl` has both a fixed-step and variable time-step integrator for $s$-stage implicit Runge-Kutta methods.  It builds the Butcher table dynamically, just specify the number of stages $s$ and the order of the fixed-step method is $2s$.  The adaptive step method does suffer from numerical issues at high order (circa $s>15$) due to very large coefficients that arise in the error estimator.
- `UdwadiaKalaba.jl` builds the right-hand side of the [Udwadia-Kalaba equation](https://en.wikipedia.org/wiki/Udwadia%E2%80%93Kalaba_equation), and forms a function in state-space for use with the Runge-Kutta integrator.  The construction of the constraint matrix is performed by automatic differentiation using [ForwardDiff.jl](https://github.com/JuliaDiff/ForwardDiff.jl).
- `nPendulum.jl` sets up various functions needed by the Udwadia-Kalaba equation (such as the mass matrix, external forces, and the holonomic constraints) for $N$ simple pendulums.  It also has routines to get a random set of parameters, and a random set of initial conditions that are consistent with the constraints.
- `nPendulumAnimation.jl` takes the results of the simulation, and generates an animation of the pendulums.
- `main.jl` is the top-level script file that sets up and runs an example problem.

## TODO
- make documentation
- put **more better** annotations throughout the code
- update the solver inside the implicit Runge-Kutta step to use Broyden updates
- add [ProgressMeter.jl](https://github.com/timholy/ProgressMeter.jl) to both the IRK solver and the animation generator.
- repackage the code into proper modules for submission to the Pkg registry
