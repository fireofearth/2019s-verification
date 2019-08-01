using IntervalArithmetic
using AffineArithmetic
using ForwardDiff

 #=
 # Run taylor model constructor constructTM
 #
 # OMG it works beautifully.
=#
include("../ODEIntegration.jl")

# problem z' = f(z) where f(z) = Az
f(z::Vector) = [2 4; 4 2] * z
# initial conditions
z0 = [5.0; -1.0]
# analytical solution
z(t::Real) = 2*exp(6*t)*[1; 1] - 3*exp(-2*t)*[-1; 1]

# constructed taylor approximation
cT = constructTM(f; order=5, T=Real)
# instance of taylor approximation
T(t::Real) = cT(t, 0.0, z0, z0)

using Plots
pyplot()

a = 0 .. 1.0
r = range(inf(a), stop=sup(a), length=1000)
# plot taylor approximation
plot(r, t -> (T(t))[1], label="T₁(t)", linewidth=2, xaxis="t", yaxis="z(t)")
plot!(t -> (T(t))[2], label="T₂(t)", linewidth=2)

# plot analytical solution
#plot!(t -> (z(t))[1], label="z₁(t)", linewidth=5)
#plot!(t -> (z(t))[2], label="z₂(t)", linewidth=5)

