using IntervalArithmetic

 #=
 # Run taylor model constructor constructTM
 #
 # Remark:
 # ForwardDiff stops working when using > 7 nested derivatives
 # taylor approximation from constructTM() does not work with order > 7
 # Gives the correct output with order = 7
=#
include("../Flowpipes.jl")

 #=
 # Conversion of a critically damped oscillator y⃮⃮̲'' + 16y' + 64y = 0
 # to sytem of autonomous linear equations z' = A * z
=#
f(z::Vector) = [0 1; -64 -16] * z

# initial conditions
# tₙ = 0.3. τ = 0.05 -> 6 steps
tspan = (0.0, 1.0)
τ     = 0.05
z₀ = [1.0 ± 1.0, -10.0 .. 12.0]

disp("getting inner and outer approximations")
st, sz, sia = approximate(f, tspan, τ, z₀)
disp("cleaning data")
sz1  = [inf.(sz[1,:]) sup.(sz[1,:])]
sz2  = [inf.(sz[2,:]) sup.(sz[2,:])]
sia1 = sia[1,:]
sia1 = [(x -> isnan(x) ? x : inf(x)).(sia1) (x -> isnan(x) ? x : sup(x)).(sia1)]
sia2 = sia[2,:]
sia2 = [(x -> isnan(x) ? x : inf(x)).(sia2) (x -> isnan(x) ? x : sup(x)).(sia2)]

disp("importing Plots...")
using Plots
pyplot()
disp("done")

p = plot(st, x -> NaN, xaxis="t", yaxis="z(t)")
plot!(p, st, sz1, label="[z₁](t)", lw=1, linecolor=:red)
plot!(p, st, sz2, label="[z₂](t)", lw=1, linecolor=:red)
plot!(p, st, sia1, label="]z₁[(t)", lw=1, linecolor=:blue)
plot!(p, st, sia2, label="]z₂[(t)", lw=1, linecolor=:blue)

p







