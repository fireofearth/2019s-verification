using IntervalArithmetic
using AffineArithmetic
using ForwardDiff

 #=
 # Creating enclosures of ODE at regular time intervals [tᵢ, tᵢ + τ]
 # and verify that the analytic solution is bounded within these enclosures
 #
 # fixpoint() can fail when τ is too large
=#
include("../Flowpipes.jl")

# problem z' = f(z) where f(z) = Az
f(z::Vector) = [2 4; 4 2] * z

# analytical solution
z(t::Real) = 2*exp(6*t)*[1; 1] - 3*exp(-2*t)*[-1; 1]
Jz(t::Real) = [0.5*(exp(6*t) + exp(-2*t)) 0.5*(exp(6*t) - exp(-2*t));
               0.5*(exp(6*t) - exp(-2*t)) 0.5*(exp(6*t) + exp(-2*t))]

# time variables
τ = 0.05
d = 40
r = 0.0:τ:(τ*d)
τₚ = 0.002
dₚ = Int(round(τ*d/τₚ, digits=0))
rₚ = 0.0:τₚ:(τₚ*dₚ)

# initial conditions
z0 = [5.0; -1.0]
J0 = [Affine(1.0 ± 0.1)  Affine(0.0 ± 0.1);
      Affine(0.0 ± 0.1)  Affine(1.0 ± 0.1)]

# initial conditions
inits = [ ]
for i in 0:d
    t = τ*Float64(i)
    z₁  = z(t)[1]
    z₂  = z(t)[2]
    J₁₁ = Jz(t)[1,1]
    J₁₂ = Jz(t)[1,2]
    J₂₁ = Jz(t)[2,1]
    J₂₂ = Jz(t)[2,2]
    push!(inits, (
        [J₁₁ ± 0.1  J₁₂ ± 0.1; J₂₁ ± 0.1 J₂₂ ± 0.1],
        [Affine(z₁ ± 0.1), Affine(z₂ ± 0.1)]
       ))
end

# priori enclosure
disp("getting priori enclosures")
Fzas = map(z -> fixedJacPoint(f, z[1], z[2], τ), inits)
Fzs  = cat(Fzas...; dims=3)
Fzs  = Interval.(Fzs)

disp("cleaning data")
for fn in (:inf, :sup), i in 1:2, j in 1:2
    fname = Symbol("fixpt$fn$i$j")
    @eval function ($fname)(t::Real)
        ival = Fzs[$i, $j, Int(floor(t/τ)) + 1]
        ($fn)(ival)
    end
end

disp("importing Plots...")
using Plots
pyplot()
disp("done")

# plot analytical solution
p = plot(rₚ, t -> (Jz(t))[1, 1], label="Jz₁₁", linecolor=:black, lw=1, xaxis="t", yaxis="z(t)")

# plot enclosures
plot!(p, rₚ, [fixptinf11.(rₚ), fixptsup11.(rₚ)], label="[Jz₁₁]", linecolor=:red, lw=1)







