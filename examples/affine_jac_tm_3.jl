using LinearAlgebra
using IntervalArithmetic
using AffineArithmetic
using ForwardDiff

 #=
 # Run the jacobian taylor model constructor constructJacTM
 #
 # Result:
 #  - initial conditions are different. I should double check
 #
 # Remark:
 # ForwardDiff stops working when using > 7 nested derivatives
 # taylor approximation of jacobian from constructJacTM() does not work with order > 6
 # Gives the correct output with order = 6
=#
include("../Flowpipes.jl")

 #=
 # Conversion of a critically damped oscillator y⃮⃮̲'' + 16y' + 64y = 0
 # to sytem of autonomous linear equations z' = A * z
=#
f(z::Vector) = [0 1; -64 -16] * z

# analytical solution representing [y; y']
z(t::Real)  = exp(-8*t)*[1 + 28*t; 20 - 224*t]
Jz(t::Real) = [exp(-8*t) + 8*t*exp(-8*t)  t*exp(-8*t);
               -64*t*exp(-8*t)            exp(-8*t) - 8*t*exp(-8*t)]

# initial conditions
z0 = [Affine(z(0.0)[1] ± 0.3); Affine(z(0.0)[2] ± 0.3)]
J0 = [Affine(1.0 ± 0.1)  Affine(0.0 ± 0.1);
      Affine(0.0 ± 0.1)  Affine(1.0 ± 0.1)]

# time step
τ = 0.02
d = 10
r = 0.0:τ:(τ*d)

# constructed taylor approximation
cJacT = constructJacTM(f; order=4)
# instance of taylor approximation with initial conditions
JacT(t::Real) = cJacT(t, 0.0, z0, J0, z0, J0)

disp("computing taylor approximation")
JacTo = JacT.(r)
JacTo = cat(JacTo...; dims=3)
JacTo = Interval.(JacTo)
Jzo   = Jz.(r)
Jzo = cat(Jzo...; dims=3)

disp("importing Plots...")
using Plots
pyplot()
disp("done")

p = plot(r, t -> NaN, xaxis="t", yaxis="∂zᵢ/∂z₀ⱼ (t)")
for i in 1:2, j in 1:2
    ii = i == 1 ? "₁" : "₂"
    jj = j == 1 ? "₁" : "₂"
    clr = (i, j) == (1, 1) ? (:red) :
        (i, j) == (1, 2) ? (:blue) :
        (i, j) == (2, 1) ? (:green) : (:orange)
    
    # plot taylor approximation
    plot!(p, r, [inf.(JacTo[i,j,:]) sup.(JacTo[i,j,:])], 
          label="J(T)$ii$jj", linecolor=clr, lw=1)

    # plot analytical solution
    plot!(p, r, Jzo[i,j,:], label="J(z)$ii$jj", linecolor=:black, lw=1)
end

p

