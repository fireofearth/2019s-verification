using IntervalArithmetic
using AffineArithmetic
using ForwardDiff

 #=
 # Run taylor model constructor constructTM
 #
 # Remark:
 # ForwardDiff stops working when using > 7 nested derivatives
 # taylor approximation from constructTM() does not work with order > 7
 # Gives the correct output with order = 7
=#
include("../ODEIntegration.jl")

# problem z' = f(z) where f(z) = Az
f(z::Vector) = [2 4; 4 2] * z
# initial conditions
z0 = [5.0; -1.0]
# analytical solution
z(t::Real) = 2*exp(6*t)*[1; 1] - 3*exp(-2*t)*[-1; 1]

# constructed taylor approximation
cT = constructTM(f; order=7, T=Real)
# instance of taylor approximation
T(t::Real) = cT(t, 0.0, z0, z0)

disp("importing Plots...")
using Plots
pyplot()
disp("done")

a = 0 .. 1.0
r = range(inf(a), stop=sup(a), length=1000)
# plot taylor approximation
plot(r, t -> (T(t))[1], label="T₁(t)", lw=2, xaxis="t", yaxis="z(t)")
#plot!(t -> (T(t))[2], label="T₂(t)", lw=2)

# plot analytical solution
plot!(t -> (z(t))[1], label="z₁(t)", lw=5, markeralpha=0.3, linestyle=:dash)
#plot!(t -> (z(t))[2], label="z₂(t)", lw=5, markeralpha=0.3, linestyle=:dash)

