using Zygote
using TaylorSeries
#using DiffLogger
#using ForwardDiff
using AffineArithmetic

include("../helper.jl")

#=
Can Zygote can replace ForwardDiff for automatic differentiation?

Zygote seems very unstable.
=#

disp("using Taylor1 and Zygote with sin")

# original function sin
f(x) = sin(x)
# taylor expansion of sin
ftm = f(Taylor1(5))

# reference
#https://math.oregonstate.edu/home/programs/undergrad/CalculusQuestStudyGuides/SandS/PowerSeries/mclaurin_list.html
disp(ftm) # 1.0 t - 0.16666666666666666 t³ + 0.008333333333333333 t⁵ + 𝒪(t⁶)
print("\n")

df(x)  = cos(x)
zgf(x) = Zygote.gradient(f, x)

# Zygote.gradient gives the derivative of sin(x) for any real.
v = [π * x/6 for x in 0:6]
for ii in 1:7
    print("x = $(v[ii])\ndf(x) = $(df(v[ii]))\nZygote.gradient(f,x) = $(zgf(v[ii])[1])\n\n")
end

# Zygote.gradient can allow us to find the taylor series of a derivative.
dftm  = df(Taylor1(5))
zgftm =  Zygote.gradient(f, Taylor1(5))[1]

disp(dftm)
disp(zgftm)

##===============================
disp("using Taylor1 and Zygote and g(x) = x^3 + 1")

tm = Taylor1(5)
g(x) = x^3 + 1

gtm = g(tm)
disp(gtm)

dgtm  = df(Taylor1(5))
disp(dgtm)
zggtm =  Zygote.gradient(g, Taylor1(5))[1]
disp(zggtm)

##===============================

#=
disp("using Affine and Zygote")

g(x::Affine) = x^3 + 1

a = Affine(2.0, [0.1, -0.2], [1, 3])
ga = g(a)
disp(ga)

dg(x) = 3*x^2
zgg(x) = Zygote.gradient(g, x)

dga  = dg(a)
disp(dga)
zgga = zgg(a)[1]
disp(zgga)
=#








