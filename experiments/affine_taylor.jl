
using TaylorSeries
using ForwardDiff
using AffineArithmetic

include("../helper.jl")

# original function sin
f(x) = sin(x)
# taylor expansion of sin
ftm = f(Taylor1(5))

# reference
#https://math.oregonstate.edu/home/programs/undergrad/CalculusQuestStudyGuides/SandS/PowerSeries/mclaurin_list.html
disp(ftm) # 1.0 t - 0.16666666666666666 t³ + 0.008333333333333333 t⁵ + 𝒪(t⁶)

df(x) = ForwardDiff.derivative(f, x)
dftm  = df(Taylor1(5))

disp(dftm)


