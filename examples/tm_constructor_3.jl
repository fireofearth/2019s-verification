using IntervalArithmetic
using AffineArithmetic
using ForwardDiff

 #=
 # Run taylor model constructor constructTM
 # 
 # It looks good.
=#
include("../ODEIntegration.jl")

 #=
 # Conversion of a critically damped oscillator y⃮⃮̲'' + 16y' + 64y = 0
 # to sytem of autonomous linear equations z' = A * z
=#
f(z::Vector) = [0 1; -64 -16] * z
# initial conditions
z0 = [1; 20]
# analytical solution representing [y; y']
z(t::Real) = exp(-8*t)*[1 + 28*t; 20 - 224*t]

# constructed taylor approximation
cT = constructTM(f; order=7, T=Real)
# instance of taylor approximation
T(t::Real) = cT(t, 0.0, z0, z0)

disp("importing Plots...")
using Plots
pyplot()
disp("done")

a = 0 .. 0.3
r = range(inf(a), stop=sup(a), length=1000)
# plot taylor approximation
plot(r, t -> (T(t))[1], label="T₁(t)", lw=2, xaxis="t", yaxis="z(t)")
plot!(t -> (T(t))[2], label="T₂(t)", lw=2)

# plot analytical solution
plot!(t -> (z(t))[1], label="z₁(t)", lw=5, markeralpha=0.3, linestyle=:dash)
plot!(t -> (z(t))[2], label="z₂(t)", lw=5, markeralpha=0.3, linestyle=:dash)
