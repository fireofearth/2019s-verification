include("../helper.jl")

using Random, IntervalArithmetic
using AffineArithmetic
using ForwardDiff


f(x::Real)  = x^2
df(x::Real) = 2*x

a1 = Affine(32.1, [0.1, -0.2, 1.5, -2.0], [1, 3, 4, 5])

disp(repr(df(a1)))
disp(repr(ForwardDiff.derivative(f, a1)))
