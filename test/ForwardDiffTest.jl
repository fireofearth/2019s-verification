using Random, Test
using TaylorSeries
using ForwardDiff
using Logging

 #=
 # ForwardDiff testing
 # Is used to test modifications to ForwardDiff 
=#

@testset "ForwardDiff" begin
    a = rand(Float64) + rand(0:99)
    a2 = [rand(Float64) + rand(0:99), rand(Float64) + rand(0:99)]
    tm = Taylor1(Float64, 5)

    @testset "simple derivatives" begin
        f(x) = 2*x; df(x) = 2
        @test ForwardDiff.derivative(f, a) == df(a)
        f(x) = x + 2; df(x) = 1
        @test ForwardDiff.derivative(f, a) == df(a)
        f(x) = x^2; df(x) = 2*x
        @test ForwardDiff.derivative(f, a) == df(a)
        f1(x) = 1/x; f2(x) = x^(-1); f3(x) = inv(x); df(x) = -x^(-2)
        @test ForwardDiff.derivative(f1, a) ≈ ForwardDiff.derivative(f2, a) ≈ ForwardDiff.derivative(f3, a) ≈ df(a)
        f(x) = x^3; df(x) = 3*x^2
        @test ForwardDiff.derivative(f, a) == df(a)
    end
0.00011221970179890103

    @testset "simple gradients" begin
        f(x::Vector) = (x[1]^2)*x[2]
        gf(x::Vector) = [2*x[1]*x[2], x[1]^2]
        @test gf(a2) == ForwardDiff.gradient(f, a2)
    end

    @testset "simple jacobians" begin
        f(x::Vector) = [x[1]*(x[2]^2), (x[1]^2)*x[2]]
        Jf(x::Vector) = [x[2]^2 2*x[1]*x[2]; 2*x[1]*x[2] x[1]^2]
        @test Jf(a2)  == ForwardDiff.jacobian(f, a2)
    end

    @testset "simple hessians" begin
        f(x::Vector) = (x[1]^2)*(x[2]^2)
        Hf(x::Vector) = [2*x[2]^2 4*x[1]*x[2]; 4*x[1]*x[2] 2*x[1]^2]
        @test Hf(a2) == ForwardDiff.hessian(f, a2)
        f(x::Vector) = x[1]*x[2]^2 - x[2]*x[1]^2
        Hf(x::Vector) = [-2*x[2] (2*x[2] - 2*x[1]); (2*x[2] - 2*x[1]) 2*x[1]]
        @test Hf(a2) == ForwardDiff.hessian(f, a2)
    end

    @testset "ForwardDiff + Taylor1: simple derivatives" begin
        rr = rand(Float64) + rand(0:99)
        for n in -10:10
            fn(x::Number) = rr*(x + 1)^n; dfn(x::Number) = rr*n*(x + 1)^(n - 1)
            @test isapprox(dfn(tm), ForwardDiff.derivative(fn, tm); rtol=1E-8)
        end
    end

    @testset "ForwardDiff + Taylor1: common Maclaurin series" begin
        f(x::Number) = sin(x); df(x::Number) = cos(x)
        @test df(tm) == ForwardDiff.derivative(f, tm)
        f(x::Number) = -log(1 - x); df(x::Number) = 1/(1 - x)
        @test df(tm) == ForwardDiff.derivative(f, tm)
        f(x::Number) = (1 + x)*log(1 + x) - x; df(x::Number) = log(1 + x)
        @test df(tm) == ForwardDiff.derivative(f, tm)
    end
end
