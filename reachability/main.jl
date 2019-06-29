include("ODECommon.jl")
include("ODEIntegration.jl")

# example using Brusselator
f = ODEFunc(BRUSSELATOR)
x = SVector(0.5, 0.5)
print(f(x))
