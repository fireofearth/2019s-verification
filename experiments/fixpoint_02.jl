using IntervalArithmetic
using AffineArithmetic
using ForwardDiff

 #=
 # very slow plotting for some reason?
=#
include("../ODEIntegration.jl")

 #=
 # Brusselator problem
 # z' = f(z)
=#
f(z::Vector) = [1 - (1.5 + 1)*z[1] + z[2]*z[1]^2; 1.5*z[1] - z[2]*z[1]^2]

# time step
τ = 0.1
d = 5

# initial conditions
z0 = [1.0; 0.0]
zInits = [[1.00 ± 0.05; 0.00 ± 0.05],
          [0.88 ± 0.08; 0.12 ± 0.08],
          [0.77 ± 0.08; 0.25 ± 0.08],
          [0.70 ± 0.08; 0.33 ± 0.08],
          [0.63 ± 0.08; 0.45 ± 0.08]]


# constructed taylor approximation
cT = constructTM(f; order=7, T=Real)
# instance of taylor approximation
T(t::Real) = cT(t, 0.0, z0, z0)

disp("getting priori enclosures")
# priori enclosure
Fza = map(x -> fixedPoint(f, x, τ), zInits)
Fz  = map(x -> Interval.(x), Fza)

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

# plot analytical solution
p = plot(r, o1, label="T₁", linecolor=:black, lw=1, xaxis="t", yaxis="z(t)")
plot!(p, r, o2, label="T₂", linecolor=:grey, lw=1)

for i in 1:d
    plot!(p, r, t -> (τ*Float64(i - 1) < t < τ*Float64(i) ? sup(Fz[i][1]) : NaN), label="[z₁]", linecolor=:red, lw=1)
    plot!(p, r, t -> (τ*Float64(i - 1) < t < τ*Float64(i) ? inf(Fz[i][1]) : NaN), label="[z₁]", linecolor=:red, lw=1)
    
    plot!(p, r, t -> (τ*Float64(i - 1) < t < τ*Float64(i) ? sup(Fz[i][2]) : NaN), label="[z₂]", linecolor=:blue, lw=1)
    plot!(p, r, t -> (τ*Float64(i - 1) < t < τ*Float64(i) ? inf(Fz[i][2]) : NaN), label="[z₂]", linecolor=:blue, lw=1)
end

p


