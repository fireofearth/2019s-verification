using Test, Random
using IntervalArithmetic
using ModalIntervalArithmetic
using Logging

 #=
 # Modal Interval testing
=#

MIN = typemin(Int)
MAX = typemax(Int)
EPS = eps(Float64)

@testset "modal interval general" begin
    xₘ = rand((MIN + 100):(MAX - 100))
    x₁ = rand(MIN:xₘ)
    x₂ = rand((x₁ + 1):MAX)
    @test x₁ < x₂

    @testset "proper interval" begin
        X = ModalInterval(x₁, x₂)
        @test inf(X)  == x₁ && sup(X) == x₂
        @test mod(X)  == PROPER
        @test Interval(X) == Interval(x₁, x₂)
        @test dual(X) == ModalInterval(x₂, x₁)
        @test prop(X) == X
    end

    @testset "improper integral" begin
        X = ModalInterval(x₂, x₁)
        @test inf(X)  == x₂ && sup(X) == x₁
        @test mod(X)  == IMPROPER
        @test dual(X) == prop(X)
        @test_throws DomainError Interval(X)
    end

    @testset "real number integral" begin
        X = ModalInterval(x₁, x₁)
        @test inf(X)  == x₁ && sup(X) == x₁
        @test mod(X)  == REAL_NUMBER
        @test isreal(X)
        @test Interval(X) == Interval(x₁)
    end

    @testset "conversion" begin
        X = Interval(x₁, x₂)
        @test one(ModalInterval{Int}) == ModalInterval{Int}(1, 1)
        @test zero(ModalInterval{Int}) == ModalInterval{Int}(0, 0)
        @test convert(ModalInterval{Float64}, X) == ModalInterval{Float64}(x₁, x₂)
        @test convert(ModalInterval{Float64}, x₁) == ModalInterval{Float64}(x₁)
    end
end

@testset "modal interval arithmetic" begin
    q₁ = rand(-Inf .. Inf)
    q₃ = rand(q₁ .. Inf)
    q₃ = rand(q₁ .. q₂)
    disp(q₁)
    disp(q₂)
    disp(q₃)
    x₁ = rand(-Inf .. q₁)
    x₂ = rand(q₁ .. q₂)
    x₃ = rand(q₂ .. q₃)
    x₄ = rand(q₃ .. Inf)
    @test x₁ ≤ x₂ ≤ x₃ ≤ x₄

    @testset "comparison" begin
        X = ModalInterval(x₁, x₄)
        Y = ModalInterval(x₂, x₃)
        disp(X)
        disp(Y)
        @test X ⊆ Y
        @test dual(Y) ⊆ dual(X)
    end
end

#@testset "modal interval comparison" begin
#    X  = ModalInterval(2,3)
#    Y  = ModalInterval(1,4)
#    dX = ModalInterval(3,2)
#    dY = ModalInterval(4,1)
#    # ([2,3]',∃) = [2,3] ⊆ [1,4] = ([1,4]',∃)
#    # => [2,3]' ⊆ [1,4]'
#    @test  X ⊆ Y  &&  Y ⊇ X
#    # ([2,3]',∀) = [3,2] ⊇ [4,1] = ([1,4]',∀)
#    # => [2,3]' ⊆ [1,4]'
#    @test dX ⊇ dY && dY ⊆ dX 
#end

#@testset "modal interval additive group" begin
#    A = ModalInterval(1,2)
#    B = ModalInterval(3,4)
#    @test A + B == ModalInterval(1+3,2+4)
#    @test -A == ModalInterval(-2,-1)
#    @test A - B == ModalInterval(1-4,2-3)
#    @test A * B == ModalInterval(1 * 3, 2 * 4)
#end

