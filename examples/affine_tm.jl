using IntervalArithmetic
using AffineArithmetic
using ForwardDiff

 #=
 # Run taylor model constructor constructTM, and use affine forms
 # as initial conditions
 #
 # - order=4: affine z₀ bounds analytical solution up to t=0.3
 # - order=5: affine z₀ bounds analytical solution up to t=0.4
 #
 # Remark:
 # ForwardDiff stops working when using > 7 nested derivatives
 # taylor approximation from constructTM() does not work with order > 7
 # Gives the correct output with order = 7
=#
include("../Flowpipes.jl")

# problem z' = f(z) where f(z) = Az
f(z::Vector) = [2 4; 4 2] * z

# initial conditions
z0 = [Affine(5.0 ± 0.1); Affine(-1.0 ± 0.1)]
τ = 0.05

# analytical solution
z(t::Real) = 2*exp(6*t)*[1; 1] - 3*exp(-2*t)*[-1; 1]

# constructed taylor approximation
cT = constructTM(f; order=4)

# instance of taylor approximation
T(t::Real) = cT(t, 0.0, z0, z0)

# results
disp("getting taylor approximations")
r   = 0.0:τ:(τ*10)
av = T.(r)
disp("cleaning data")
iv  = [Interval(x[1]) for x in av]
uv  = [sup(x) for x in iv]
lv  = [inf(x) for x in iv]

disp("importing Plots...")
using Plots
pyplot()
disp("done")

# plot taylor approximation
plot(r, uv, label="T₁(t)", linecolor=:red, lw=1, xaxis="t", yaxis="z(t)")
plot!(r, lv, label="T₁(t)", linecolor=:red, lw=1)

# plot analytical solution
plot!(r, t -> (z(t))[1], label="z₁(t)", linecolor=:black, lw=1)

