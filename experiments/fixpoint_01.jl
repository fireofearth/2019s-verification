using IntervalArithmetic
using AffineArithmetic
using ForwardDiff

 #=
 # This specific problem is failing and I have no idea why
 # it's almost as if F[z] = [z₀] + [0, τ]*[f]([z]) does not have a fixed point.
 # for this f
 #
 # UPDATE: it all depends on our choice of τ
 #
=#
include("../ODEIntegration.jl")

# problem z' = f(z) where f(z) = Az
f(z::Vector) = [11 -25; 4 -9] * z

# analytical solution
z(t::Real) = exp(t)*[5; 2] + exp(t)*[1; 0] + 2*t*exp(t)*[5; 2]

# time step
τ = 0.02

# initial conditions
z0 = [6.0 ± 0.1, 2.0 ± 0.1]
#z0 = [z(0)[1] ± 0.2; z(0)[2] ± 0.15]

disp("getting priori enclosures")
# priori enclosure
Fza = fixedPoint(f, z0, τ)
Fz  = Interval.(Fza)

disp("importing Plots...")
using Plots
pyplot()
disp("done")

a = 0 .. τ
r = range(inf(a), stop=sup(a), length=1000)
# plot analytical solution
plot(r, t -> (z(t))[1], label="z₁", linecolor=:black, lw=1, xaxis="t", yaxis="z(t)")
plot!(r, t -> (z(t))[2], label="z₂", linecolor=:grey, lw=1, xaxis="t", yaxis="z(t)")

plot!(t -> sup(Fz[1]), label="[z₁]", linecolor=:red, lw=1)
plot!(t -> inf(Fz[1]), label="[z₁]", linecolor=:red, lw=1)

plot!(t -> sup(Fz[2]), label="[z₂]", linecolor=:blue, lw=1)
plot!(t -> inf(Fz[2]), label="[z₂]", linecolor=:blue, lw=1)
