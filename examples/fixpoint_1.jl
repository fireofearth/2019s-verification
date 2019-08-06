using IntervalArithmetic
using AffineArithmetic
using ForwardDiff

 #=
 # Creating enclosures of ODE at regular time intervals [tᵢ, tᵢ + τ]
 # and verify that the analytic solution is bounded within these enclosures
 #
 # Picard–Lindelöf theorem can show that F[z] = [z₀] + [0, τ]*[f]([z]) has a fixed point
 # if τ is small enough(?)
 #
 # fixpoint() can fail when τ is too large (i.e. τ = 0.05)
=#
include("../ODEIntegration.jl")

# problem z' = f(z) where f(z) = Az
f(z::Vector) = [11 -25; 4 -9] * z

# analytical solution
z(t::Real) = exp(t)*[5; 2] + exp(t)*[1; 0] + 2*t*exp(t)*[5; 2]

# time variables
τ = 0.02
d = 40
r = 0.0:τ:(τ*d)
τₚ = 0.002
dₚ = Int(round(τ*d/τₚ, digits=0))
rₚ = 0.0:τₚ:(τₚ*dₚ)

# initial conditions
zInits = [ ]
for i in 0:d
    z₁ = z(τ*Float64(i))[1]
    z₂ = z(τ*Float64(i))[2]
    push!(zInits, [z₁ ± 0.1, z₂ ± 0.1])
end

# priori enclosure
disp("getting priori enclosures")
Fzas = map(z -> fixedPoint(f, z, τ), zInits)
Fzs  = hcat(Fzas...)
Fzs  = Interval.(Fzs)

disp("cleaning data")
for fn in (:inf, :sup), i in 1:2
    fname = Symbol("fixpt$fn$i")
    @eval function ($fname)(t::Real)
        ival = Fzs[$i, Int(floor(t/τ)) + 1]
        ($fn)(ival)
    end
end

disp("importing Plots...")
using Plots
pyplot()
disp("done")

# plot analytical solution
p = plot(rₚ, t -> (z(t))[1], label="z₁", linecolor=:black, lw=1, xaxis="t", yaxis="z(t)")
plot!(p, rₚ, t -> (z(t))[2], label="z₂", linecolor=:grey, lw=1, xaxis="t", yaxis="z(t)")

# plot enclosures
for i in 1:2
    plot!(p, rₚ, [fixptinf1.(rₚ), fixptsup1.(rₚ)], label="[z₁]", linecolor=:red, lw=1)
    plot!(p, rₚ, [fixptinf2.(rₚ), fixptsup2.(rₚ)], label="[z₂]", linecolor=:blue, lw=1)
end

p # show plot
