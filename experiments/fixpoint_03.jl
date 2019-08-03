using IntervalArithmetic
using AffineArithmetic
using ForwardDiff

 #=
=#
include("../ODEIntegration.jl")

 #=
 # Conversion of a critically damped oscillator y⃮⃮̲'' + 16y' + 64y = 0
 # to sytem of autonomous linear equations z' = A * z
=#
f(z::Vector) = [0 1; -64 -16] * z

# analytical solution representing [y; y']
z(t::Real) = exp(-8*t)*[1 + 28*t; 20 - 224*t]

# time step
τ = 0.05

# initial conditions
z0 = [1.0 ± 0.1, 20.0 ± 0.1]

# priori enclosure
disp("getting priori enclosures")
Fza = fixedPoint(f, z0, τ)
Fz  = Interval.(Fza)

disp("importing Plots...")
using Plots
pyplot()
disp("done")

r = range(0.0, stop=τ, length=100)

# plot analytical solution
p = plot(r, t -> (z(t))[1], label="z₁", linecolor=:black, lw=1, xaxis="t", yaxis="z(t)")
plot!(p, r, t -> (z(t))[2], label="z₂", linecolor=:grey, lw=1, xaxis="t", yaxis="z(t)")

plot!(p, r, t -> sup(Fz[1]), label="[z₁]", linecolor=:red, lw=1)
plot!(p, r, t -> inf(Fz[1]), label="[z₁]", linecolor=:red, lw=1)

plot!(p, r, t -> sup(Fz[2]), label="[z₂]", linecolor=:blue, lw=1)
plot!(p, r, t -> inf(Fz[2]), label="[z₂]", linecolor=:blue, lw=1)

p # show plot













