using IntervalArithmetic
using AffineArithmetic
using ForwardDiff

 #=
 # Creating enclosures of ODE at regular time intervals [tᵢ, tᵢ + τ]
 # and verify that the analytic solution is bounded within these enclosures
 #
 # - τ ≤ 0.04 good fit
 # - τ = 0.05 bad fit
 #
 # fixJacPoint() can fail when τ is too large
=#
include("../ODEIntegration.jl")

 #=
 # Conversion of a critically damped oscillator y⃮⃮̲'' + 16y' + 64y = 0
 # to sytem of autonomous linear equations z' = A * z
=#
f(z::Vector) = [0 1; -64 -16] * z

# analytical solution representing [y; y']
z(t::Real)  = exp(-8*t)*[1 + 28*t; 20 - 224*t]
Jz(t::Real) = [exp(-8*t) + 8*t*exp(-8*t)  t*exp(-8*t);
               -64*t*exp(-8*t)            exp(-8*t) - 8*t*exp(-8*t)]

# time variables
τ = 0.04
d = 40
r = 0.0:τ:(τ*d)
τₚ = 0.002
dₚ = Int(round(τ*d/τₚ, digits=0))
rₚ = 0.0:τₚ:(τₚ*dₚ)

# initial conditions
z0 = [1.0; 20.0]
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

p = [plot(rₚ, x -> NaN, xaxis="t", yaxis="z(t)") plot(rₚ, x -> NaN, xaxis="t", yaxis="z(t)");
     plot(rₚ, x -> NaN, xaxis="t", yaxis="z(t)") plot(rₚ, x -> NaN, xaxis="t", yaxis="z(t)")]
for i in 1:2, j in 1:2
    fixptinf = Symbol("fixptinf$i$j")
    fixptsup = Symbol("fixptsup$i$j")
    fixptinf = @eval $fixptinf
    fixptsup = @eval $fixptsup
    ii = i == 1 ? "₁" : "₂"
    jj = j == 1 ? "₁" : "₂"
    
    # plot analytical solution
    plot!(p[i, j], rₚ, t -> (Jz(t))[i, j], 
          label="J(z)$ii$jj", linecolor=:red, lw=1)

    # plot priori enclosures
    plot!(p[i, j], rₚ, [fixptinf.(rₚ), fixptsup.(rₚ)], 
          label="[J(z)$ii$jj]", linecolor=:black, lw=1)
end

plot(p[1,1], p[1,2], p[2,1], p[2,2])





