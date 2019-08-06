using IntervalArithmetic
using AffineArithmetic
using ForwardDiff

 #=
 # Run taylor model constructor constructTM
 #
=#
include("../ODEIntegration.jl")

 #=
 # Brusselator problem
 # z' = f(z)
=#
f(z::Vector) = [1 - (1.5 + 1)*z[1] + z[2]*z[1]^2; 1.5*z[1] - z[2]*z[1]^2]
# initial conditions
z0 = [1.0, 0.0]

# constructed taylor approximation
cT = constructTM(f; order=7, T=Real)
# instance of taylor approximation
T(t::Real) = cT(t, 0.0, z0, z0)

# results
disp("getting taylor approximations")
r = range(0.0, stop=0.5, length=50)
out = T.(r)
disp("cleaning data")
o1  = [o[1] for o in out]
o2  = [o[2] for o in out]

disp("importing Plots...")
using Plots
pyplot()
disp("done")

# plot taylor approximation
plot(r, o1, label="T₁(t)", lw=2, xaxis="t", yaxis="z(t)")
plot!(r, o2, label="T₂(t)", lw=2)



