using Test, Random, Dates
using IntervalArithmetic
using ModalIntervalArithmetic
using AffineArithmetic
using Logging

@testset "modal interval function extensions" begin
    @testset "affines" begin
        a = Affine(12.0, [0.1, -0.1, 0.02], [1, 2, 3])
        X = ModalInterval(a)
        @test inf(X) == min(a)
        @test sup(X) == max(a)
    end
end
