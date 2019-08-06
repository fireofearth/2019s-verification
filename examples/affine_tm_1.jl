using IntervalArithmetic
using AffineArithmetic
using ForwardDiff

 #=
 # Run taylor model constructor constructTM, and use affine forms
 # as initial conditions
 #
 # - order=4: affine z₀ bounds analytical solution up to t=2
 # - order=5: affine z₀ bounds analytical solution up to t=0.4
 #
 # Remark:
 # ForwardDiff stops working when using > 7 nested derivatives
 # taylor approximation from constructTM() does not work with order > 7
 # Gives the correct output with order = 7
=#
include("../ODEIntegration.jl")

# problem z' = f(z) where f(z) = Az
f(z::Vector) = [11 -25; 4 -9] * z
#
# initial conditions
z0 = [Affine(6.0 ± 0.1); Affine(2.0 ± 0.1)]
τ = 0.05
d = 60

# analytical solution
z(t::Real) = exp(t)*[5; 2] + exp(t)*[1; 0] + 2*t*exp(t)*[5; 2]

# constructed taylor approximation
cT = constructTM(f; order=4)
# instance of taylor approximation
T(t::Real) = cT(t, 0.0, z0, z0)

# results
disp("getting taylor approximations")
r   = 0.0:τ:(τ*d)
av = T.(r)
disp("cleaning data")
iv  = [[Interval(x[1]) for x in av] [Interval(x[2]) for x in av]]
supiv  = [sup.(iv[:,1]) sup.(iv[:,2])]
infiv  = [inf.(iv[:,1]) inf.(iv[:,2])]

disp("importing Plots...")
using Plots
pyplot()
disp("done")

# plot taylor approximation
# plot analytical solution
p = plot(r, t -> (z(t))[1], label="z₁(t)", lw=1, linecolor=:black, xaxis="t", yaxis="z(t)")
plot!(p, r, t -> (z(t))[2], label="z₂(t)", lw=1, linecolor=:black)

# plot analytical solution
for i in 1:2
    ii = i == 1 ? "₁" : "₂"
    plot!(p, r, supiv[:,i], label="T$(ii)(t)", lw=1, linecolor=:red)
    plot!(p, r, infiv[:,i], label="T$(ii)(t)", lw=1, linecolor=:red)
end

p
