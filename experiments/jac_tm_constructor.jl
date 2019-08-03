using LinearAlgebra
using IntervalArithmetic
using AffineArithmetic
using ForwardDiff

 #=
 # Run the jacobian taylor model constructor constructJacTM
 #
 # Remark:
 # ForwardDiff stops working when using > 7 nested derivatives
 # taylor approximation of jacobian from constructJacTM() does not work with order > 6
 # Gives the correct output with order = 6
=#
include("../ODEIntegration.jl")

# problem z' = f(z) where f(z) = Az
f(z::Vector) = [2 4; 4 2] * z
# initial conditions
z0 = [5.0; -1.0]
J0 = Matrix{Float64}(I, 2, 2)
# analytical solution
# z(t::Real) = 2*exp(6*t)*[1; 1] - 3*exp(-2*t)*[-1; 1]
Jz(t::Real) = [0.5*(exp(6*t) + exp(-2*t)) 0.5*(exp(6*t) - exp(-2*t));
               0.5*(exp(6*t) - exp(-2*t)) 0.5*(exp(6*t) + exp(-2*t))]

# constructed taylor approximation
cJacT = constructJacTM(f; order=5, T=Real)
# instance of taylor approximation with initial conditions
JacT(t::Real) = cJacT(t, 0.0, z0, J0, z0, J0)

disp("importing Plots...")
using Plots
pyplot()
disp("done")

a = 0 .. 1.0
r = range(inf(a), stop=sup(a), length=1000)
# plot taylor approximation
plot(r, t -> (JacT(t))[1,1], label="J(T)₁₁", lw=2, xaxis="t", yaxis="∂zᵢ/∂z₀ⱼ (t)")
#plot!(t -> (JacT(t))[1,2], label="T(T)₁₂", lw=2)
#plot!(t -> (JacT(t))[2,1], label="T(T)₂₁", lw=2)
#plot!(t -> (JacT(t))[2,2], label="T(T)₂₂", lw=2)

# plot analytical solution
plot!(t -> (Jz(t))[1,1], label="J(z)₁₁", lw=5, markeralpha=0.3, linestyle=:dash)
#plot!(t -> (Jz(t))[1,2], label="J(z)₁₂", lw=5, markeralpha=0.3, linestyle=:dash)
#plot!(t -> (Jz(t))[2,1], label="J(z)₂₁", lw=5, markeralpha=0.3, linestyle=:dash)
#plot!(t -> (Jz(t))[2,2], label="J(z)₂₂", lw=5, markeralpha=0.3, linestyle=:dash)
