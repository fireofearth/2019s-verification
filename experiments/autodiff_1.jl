include("../helper.jl")

using AffineArithmetic
using DiffLogger
using ForwardDiff

f(x::Real)  = x^(-2)
df(x::Real) = -2*x^(-3)
# actual
df(x::Real) = (TestReal(1) * -(inv(x) * inv(x))) * (2 * inv(x)^1)
df(x::Real) = (TestReal(1) * -(inv(x) * inv(x))) * (2 * inv(x))

a1 = TestReal(2.0)
disp("begin")
disp(repr(df(a1)))
disp(repr(ForwardDiff.derivative(f, a1)))

#=
b1 = inv(a1)
b2 = b1 * b1
b3 = -(b2)
b4 = TestReal(1) * b3
b5 = b1^2 # unused
b6 = b1^1
b7 = 2 * b6
b8 = b4 * b7

df(a1) = (TestReal(1) * -(inv(a1) * inv(a1))) * (2 * inv(a1)^1)
df(a1) = (TestReal(1) * -abs2(inv(a1)) * (2 * inv(a1)^1)

[ Info: in function inv(p::Affine)
[ Info: in function *(a::Affine, p::Affine)
[ Info: in function *(a::Affine, p::Affine)
[ Info: in function ^(p::Affine, n::Int)
[ Info: in function ^(p::Affine, n::Int)
[ Info: in function *(a::Affine, p::Affine)
=#


















