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

# initial conditions
z0 = [6.0; 2.0]
τ = 0.05
d = 40
dp = 10

J0 = Matrix{Float64}(I, 2, 2)
# analytical solution
# z(t::Real) = exp(t)*[5; 2] + exp(t)*[1; 0] + 2*t*exp(t)*[5; 2]
Jz(t::Real) = [exp(t)*[1; 0] + t*exp(t)*[10; 4] exp(t)*[5/2; 2] - t*exp(t)*[25; 10]]


# constructed taylor approximation
cJacT = constructJacTM(f; order=5)
# instance of taylor approximation with initial conditions
JacT(t::Real) = cJacT(t, 0.0, z0, J0, z0, J0)

# results
disp("getting taylor approximations")
r = 0.0:τ:(τ*d)
JacTMat = JacT.(r)
disp("cleaning data")
JaczMat = Jz.(r)

disp("importing Plots...")
using Plots
pyplot()
disp("done")

p1 = plot(r, x -> NaN, xaxis="t", yaxis="∂zᵢ/∂z₀ⱼ (t)")
for i in 1:2, j in 1:2
    ii = i == 1 ? "₁" : "₂"
    jj = j == 1 ? "₁" : "₂"

    # plot taylor approximation
    plot!(p1, r, (M -> getindex(M, i, j)).(JacTMat), label="J(T)$(ii)$(jj)",
          lw=1, linecolor=:red)

    # plot analytical solution
    plot!(p1, r, (M -> getindex(M, i, j)).(JaczMat), label="J(z)$(ii)$(jj)", 
          lw=1, linecolor=:black)
end

r = 0.0:τ:(τ*(dp - 1))
p2 = plot(r, x -> NaN, xaxis="t", yaxis="∂zᵢ/∂z₀ⱼ (t)")
for i in 1:2, j in 1:2
    ii = i == 1 ? "₁" : "₂"
    jj = j == 1 ? "₁" : "₂"

    # plot taylor approximation
    plot!(p2, r, (M -> getindex(M, i, j)).(JacTMat[1:dp]), label="J(T)$(ii)$(jj)",
          lw=1, linecolor=:red)

    # plot analytical solution
    plot!(p2, r, (M -> getindex(M, i, j)).(JaczMat[1:dp]), label="J(z)$(ii)$(jj)", 
          lw=1, linecolor=:black)
end

plot(p1, p2)
