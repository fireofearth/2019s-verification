using IntervalArithmetic
using TaylorSeries
using ForwardDiff

include("../Flowpipes.jl")

 #=
 # Conversion of a critically damped oscillator y⃮⃮̲'' + 16y' + 64y = 0
 # to sytem of autonomous linear equations z' = A * z
=#
f(z::Vector) = [0 1; -64 -16] * z
# initial conditions
τ = 0.05
z0 = [1.; 20.]
# analytical solution representing [y; y']
z(t::Real) = exp(-8*t)*[1 + 28*t; 20 - 224*t]

ord = 5
vtns = set_variables("z", numvars=2, order=ord)

 #=
 # Task: figure out how to make TM constructor run faster.
=#

vftns = f(vtns)

let 
    vf = [ f ]
    for i in 2:ord
        fi = (x::Vector -> ForwardDiff.jacobian(vf[i - 1], x) * f(x))
        vf = vcat(vf, fi)
    end

    global acc = z0
    for i in 1:ord
        term = 1.0
        for l in 1:i
            term *= τ / l
        end
        acc += term * vf[i](z0)
    end
end

disp(z(τ))
disp(acc)
disp(vftns)
#disp(vftns(z0))
