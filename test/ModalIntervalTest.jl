using Test, Random
using IntervalArithmetic
using ModalIntervalArithmetic
using Logging

 #=
 # Modal Interval testing
=#

@testset "modal interval getters" begin
    X = ModalInterval(1,2) # proper
    Y = ModalInterval(2,1) # improper
    Z = ModalInterval(2,2) # real_number
    dX = dual(X)
    dY = dual(Y)
    dZ = dual(Z)
    @test inf(X)  == 1 && sup(X) == 2
    @test inf(dX) == 2 && sup(dX) == 1
    @test inf(Y)  == 2 && sup(Y) == 1
    @test inf(dY) == 1 && sup(dY) == 2
    @test inf(Z)  == sup(Z)  == 2
    @test inf(dZ) == sup(dZ) == 2
    @test dX == Y && dY == X
end

@testset "modal interval comparison" begin
    X  = ModalInterval(2,3)
    Y  = ModalInterval(1,4)
    dX = ModalInterval(3,2)
    dY = ModalInterval(4,1)
    # ([2,3]',∃) = [2,3] ⊆ [1,4] = ([1,4]',∃)
    # => [2,3]' ⊆ [1,4]'
    @test  X ⊆ Y  &&  Y ⊇ X
    # ([2,3]',∀) = [3,2] ⊇ [4,1] = ([1,4]',∀)
    # => [2,3]' ⊆ [1,4]'
    @test dX ⊇ dY && dY ⊆ dX 
end

#@testset "modal interval additive group" begin
#    A = ModalInterval(1,2)
#    B = ModalInterval(3,4)
#    @test A + B == ModalInterval(1+3,2+4)
#    @test -A == ModalInterval(-2,-1)
#    @test A - B == ModalInterval(1-4,2-3)
#    @test A * B == ModalInterval(1 * 3, 2 * 4)
#end

