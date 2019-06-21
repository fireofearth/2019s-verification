using Test, IntervalArithmetic

include("interval.jl")

@testset "modal interval getters" begin
    X = ModalInterval(1,2) # proper
    Y = ModalInterval(2,1) # improper
    dX = dual(X)
    @test inf(X) == 1
    @test sup(X) == 2
    @test inf(dX) == 2
    @test sup(dX) == 1
    @test inf(Y) == 2
    @test sup(Y) == 1
    @test dX == Y
end

@testset "modal interval arithmetic" begin
    A = ModalInterval(1,2)
    B = ModalInterval(3,4)
    @test A + B == ModalInterval(1+3,2+4)
    @test -A == ModalInterval(-2,-1)
    @test A - B == ModalInterval(1-4,2-3)
end



