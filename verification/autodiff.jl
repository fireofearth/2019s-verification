include("./helper.jl")

localModulePath = "/home/fireofearth/res/mitchell-ian/2019s-verification/verification/module"
if(!(localModulePath in LOAD_PATH))
    push!(LOAD_PATH, localModulePath)
end

using Random, IntervalArithmetic
using AffineArithmetic
using ForwardDiff

 
f(x) = 1. - x*x + x
df(x) = -x + 1.
a = AAF(2.1, [0.1], [1])

disp(f(a))
disp(ForwardDiff.derivative(f,a))

