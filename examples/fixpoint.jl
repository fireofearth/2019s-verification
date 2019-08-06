using IntervalArithmetic
using AffineArithmetic
using ForwardDiff

 #=
 # Creating enclosures of ODE at regular time intervals [tᵢ, tᵢ + τ]
 # and verify that the analytic solution is bounded within these enclosures
 #
 # fixpoint() can fail when τ is too large
=#
include("../ODEIntegration.jl")

# problem z' = f(z) where f(z) = Az
f(z::Vector) = [2 4; 4 2] * z

# time step
τ = 0.1
d = 10

# initial conditions
zInits = [[5.00 ± 0.2; -1.0 ± 0.15],
          [6.10 ± 0.2; 1.18 ± 0.15],
          [8.65 ± 0.2; 4.62 ± 0.15],
          [13.7 ± 0.2; 10.4 ± 0.15],
          [23.3 ± 0.2; 20.6 ± 0.15],
          [41.2 ± 0.2; 39.0 ± 0.15],
          [74.1 ± 0.2; 72.2 ± 0.15],
          [134.1 ± 0.2; 132.6 ± 0.15],
          [243.6 ± 0.2; 242.4 ± 0.15],
          [443.3 ± 0.2; 442.3 ± 0.15]]

# analytical solution
z(t::Real) = 2*exp(6*t)*[1; 1] - 3*exp(-2*t)*[-1; 1]
# priori enclosure
Fza = map(z -> fixedPoint(f, z, τ), zInits)
Fz  = map(v -> Interval.(v), Fza)

#for i in 1:d
#    zi = z(τ*Float64(i - 1))
#    disp("$(zi[1]) $(zi[2])")
#end

disp("importing Plots...")
using Plots
pyplot()
disp("done")

r = range(0.0, stop=τ*Float64(d), length=100)
# plot analytical solution
p = plot(r, t -> (z(t))[1], label="z₁", linecolor=:black, lw=1, xaxis="t", yaxis="z(t)")

for i in 1:d
    plot!(p, t -> ((τ*Float64(i - 1) < t < τ*Float64(i)) ? sup(Fz[i][1]) : NaN), label="[z₁]", linecolor=:red, lw=1)
    plot!(p, t -> ((τ*Float64(i - 1) < t < τ*Float64(i)) ? inf(Fz[i][1]) : NaN), label="[z₁]", linecolor=:red, lw=1)
end
p # display the plot
