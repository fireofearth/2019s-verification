#using Test, Random
#using IntervalArithmetic
#using TaylorSeries
#using ModalIntervalArithmetic
#using AffineArithmetic
#using ForwardDiff
#using Logging

#include("../ODEIntegration.jl")

include("ModalIntervalTest.jl")

include("ForwardDiffTest.jl")

include("AffineTest.jl")

include("AffineDiffTest.jl")
