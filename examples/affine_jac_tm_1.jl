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
include("../ODEIntegration.jl")

# problem z' = f(z) where f(z) = Az
f(z::Vector) = [11 -25; 4 -9] * z

# analytical solution
# z(t::Real) = exp(t)*[5; 2] + exp(t)*[1; 0] + 2*t*exp(t)*[5; 2]
Jz(t::Real) = [exp(t)*[1; 0] + t*exp(t)*[10; 4] exp(t)*[5/2; 2] - t*exp(t)*[25; 10]]

# initial conditions
z0 = [Affine(6.0 ± 0.3); Affine(2.0 ± 0.3)]
J0 = [Affine(1.0 ± 0.1)  Affine(0.0 ± 0.1);
      Affine(0.0 ± 0.1)  Affine(1.0 ± 0.1)]

# time step
τ = 0.05
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

disp(JacTo)

disp("importing Plots...")
using Plots
pyplot()
disp("done")

p = plot(r, t -> NaN, xaxis="t", yaxis="∂zᵢ/∂z₀ⱼ (t)")
for i in 1:2, j in 1:2
    ii = i == 1 ? "₁" : "₂"
    jj = j == 1 ? "₁" : "₂"
    
    # plot taylor approximation
    plot!(p, r, [inf.(JacTo[i,j,:]) sup.(JacTo[i,j,:])], 
          label="J(T)$ii$jj", linecolor=:red, lw=1)

    # plot analytical solution
    plot!(p, r, Jzo[i,j,:], label="J(z)$ii$jj", linecolor=:black, lw=1)
end

p

