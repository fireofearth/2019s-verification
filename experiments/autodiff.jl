include("../helper.jl")

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
#disp(gf([2 ,3 ,1]))
#disp(ForwardDiff.gradient(f,[2 ,3 ,1])) # OK

a1 = Affine(32.1, [0.1, -0.2, 1.5, -2.0], [1, 3, 4, 5])
a2 = Affine(27.3, [10.0, 0.5, 1.0], [1, 4, 6])
a3 = Affine(38.0,  [-3.33, 9.0, -1.5, 5.25], [2, 3, 5, 6])
ax = [a1, a2, a3]

#disp(f(ax))
#disp(gf(ax))
#disp(ForwardDiff.gradient(f, ax)) # OK

# try hessian
f(x::Vector) = x[1]*x[2]^2 - x[2]*x[1]^2
Hf(x::Vector) = [-2*x[2] (2*x[2] - 2*x[1]); (2*x[2] - 2*x[1]) 2*x[1]]

ax = [a1, a2]

#disp(f(ax))
#disp(Hf(ax))
#disp(compact(ForwardDiff.hessian(f, ax)))

# try inverting zero
f(x::Real) = 1/x

# try jacobian
f(x::Vector) = [x[1]*x[2]/x[3], x[1]*x[2]*x[3], x[1]^2 + x[2]^2 + x[3]^2]
Jf(x::Vector) = [x[2]/x[3]  x[1]/x[3]  -(x[1]*x[2] /x[3] /x[3]); x[2]*x[3]  x[1]*x[3]  x[1]*x[2]; 2*x[1]  2*x[2]  2*x[3]]

ax = [a1, a2, a3]

#disp("$(a3[0] - rad(a3)) $(a3[0] + rad(a3))")
#disp((a3[0] + rad(a3))*(a3[0] - rad(a3)))
#disp(f(ax))
#disp(Jf(ax))
#disp(compact(ForwardDiff.jacobian(f, ax)))

f(x::Vector) = x[1] /x[2]
gf(x::Vector) = [1 /x[2], -(x[1] /x[2] /x[2])]

disp(gf([1., 3.]))
disp(ForwardDiff.gradient(f, [1., 3.]))














