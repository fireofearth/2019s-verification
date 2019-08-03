using IntervalArithmetic
using AffineArithmetic
using ForwardDiff

 #=
=#
include("../ODEIntegration.jl")

# problem z' = f(z) where f(z) = Az
f(z::Vector) = [11 -25; 4 -9] * z
#
# analytical solution
z(t::Real) = exp(t)*[5; 2] + exp(t)*[1; 0] + 2*t*exp(t)*[5; 2]

# time step
τ = 0.02
d = 10

# initial conditions
zInits = [[z(τ*Float64(i))[1] ± 0.2; z(τ*Float64(i))[2] ± 0.15] for i in 0:(d - 1)]

disp("getting priori enclosures")

# priori enclosure
Fza = map(z -> fixedPoint(f, z, τ), zInits)
Fz  = map(v -> Interval.(v), Fza)

disp("done")
disp("importing Plots...")
using Plots
pyplot()
disp("done")

r = range(0.0, stop=(τ*d), length=100)
# plot analytical solution
p = plot(r, t -> (z(t))[1], label="z₁", linecolor=:black, lw=1, xaxis="t", yaxis="z(t)")
plot!(p, r, t -> (z(t))[2], label="z₂", linecolor=:gold, lw=1, xaxis="t", yaxis="z(t)")

for i in 1:d
    plot!(p, r, t -> ((τ*Float64(i - 1) < t < τ*Float64(i)) ? sup(Fz[i][1]) : NaN), label="[z₁]", linecolor=:red, lw=1)
    plot!(p, r, t -> ((τ*Float64(i - 1) < t < τ*Float64(i)) ? inf(Fz[i][1]) : NaN), label="[z₁]", linecolor=:red, lw=1)

    plot!(p, r, t -> ((τ*Float64(i - 1) < t < τ*Float64(i)) ? sup(Fz[i][2]) : NaN), label="[z₁]", linecolor=:red, lw=1)
    plot!(p, r, t -> ((τ*Float64(i - 1) < t < τ*Float64(i)) ? inf(Fz[i][2]) : NaN), label="[z₁]", linecolor=:red, lw=1)
end

p # display the plot
