include("./helper.jl")

using Random, IntervalArithmetic
using AffineArithmetic
using ForwardDiff

# try scalar-scalar function
f(x::Real) = 1.0 - x^2 + x
df(x::Real) = -2*x + 1.0
a = Affine(2.1, [0.1], [1])

#disp(f(a))
#disp(df(a))
#disp(ForwardDiff.derivative(f,a)) # OK

# try vector-scalar function
f(x::Vector) = x[1]*x[3] + 2.0*x[2]*x[1] - x[3]*x[2]
gf(x::Vector) = [x[3] + 2.0*x[2], 2.0*x[1] - x[3], x[1] - x[2]]

# verify that ForwardDiff.gradient gives the actual gradient
disp(gf([2 ,3 ,1]))
disp(ForwardDiff.gradient(f,[2 ,3 ,1])) # OK

ax = [Affine(2.1, [0.1], [1]), 
      Affine(3.0, [0.2], [2]),
      Affine(1.5, [0.1], [1])]

disp(f(ax))
disp(gf(ax))
disp(ForwardDiff.gradient(f, ax)) # OK




