using Test, Random
using IntervalArithmetic
using AffineArithmetic
using ForwardDiff
using Logging

 #=
 # Affine Arithmetic ForwardDiff
 # Tests compatibility with ForwardDiff
 # 
 # Remark: ForwardDiff almost works as is; must make affines compact by removing zeros.
 # Remark: ForwardDiff gives unpredictable derivatives when f contains division.
=#

@testset "affine arithmetic ForwardDiff" begin
    centers = [32.1, 27.3, 58.0, 32.1]
    devs    = [[0.1, -0.2, 1.5, -2.0],
               [10.0, 0.5, 1.0], 
               [-3.33, 9.0, -1.5, 5.25],
               [0.1]]
    inds    = [[1, 3, 4, 5],
               [1, 4, 6],
               [2, 3, 5, 6],
               [1]]
    a1     = Affine(centers[1], devs[1], inds[1])
    a2     = Affine(centers[2], devs[2], inds[2])
    a3     = Affine(centers[3], devs[3], inds[3])
    va     = [a1, a2, a3]
    
    @testset "derivative simple" begin
        f(x::Number)  = 1.0 - x^2 + x
        df(x::Number) = -2*x + 1.0
        for a in va
            @test sameForm(df(a), ForwardDiff.derivative(f, a))
        end
    end

    @testset "derivative of 1/x" begin
        f(x::Number)  = 1/x
        df(x::Number) = -(1 /x /x)
        for a in va
            @test sameForm(df(a), ForwardDiff.derivative(f, a))
        end
    end

    @testset "derivative of const / x" begin
        rr = rand(Float64) + rand(0:99)
        f(x::Number) = rr/x
        df(x::Number) = -(rr /x /x)
        for a in va
            @test sameForm(df(a), ForwardDiff.derivative(f, a))
        end
    end

     #=
     # Remark: `df(x::Number) = 1` also works when using `==`
    =#
    @testset "derivative of x^n" begin
        rr = rand(Float64) + rand(0:99)
        f(x::Number)  = rr*x
        df(x::Number) = Affine(rr)
        for a in va
            @test sameForm(df(a), ForwardDiff.derivative(f, a))
        end
        for n in 1:10
            fn(x::Number)  = rr * x^n
            dfn(x::Number) = rr* n * (x^(n-1))
            for a in va
                @test sameForm(dfn(a), ForwardDiff.derivative(fn, a))
            end
        end
    end

     #=
     # FAIL when f(x::Number) = 1 but success when f(x::Number)  = Affine(1)
     # Remark: `df(x::Number) = 0` also works when using `==`
    =#
    @testset "derivative of constants" begin
        rr = rand(Float64) + rand(0:99)
        f(x::Number)  = Affine(rr)
        df(x::Number) = Affine(0)
        for a in va
            @test sameForm(df(a), ForwardDiff.derivative(f, a))
        end
    end

     #=
     # For higher order inverses, there is no clear closed form solution for
     # derivatives of affines. Instead we will check that the affine forms do indeed act as
     # bounds of real numbers.
    =#
    @testset "derivative of x^-n" begin
        for n in 1:10
            fn(x::Number) = x^(-n)
            dfn(x::Number) = -n*x^(-n-1)
            x = 1.2
            a = Affine(1.2, [0.0001, -0.01, 0.001], [1, 2, 3])
            @test dfn(x) ⊆ Interval(ForwardDiff.derivative(fn, a))
        end
    end

    @testset "gradient" begin
        f(x::Vector)  = x[1]*x[3] + 2.0*x[2]*x[1] - x[3]*x[2]
        gf(x::Vector) = [x[3] + 2.0*x[2], 2.0*x[1] - x[3], x[1] - x[2]]
        ax = [a1, a2, a3]
        @test sameForm(gf(ax), compact(ForwardDiff.gradient(f, ax)))
    end

    @testset "hessian" begin
        f(x::Vector) = (x[1]^2)*(x[2]^2)
        Hf(x::Vector) = [2*x[2]^2 4*x[1]*x[2]; 4*x[1]*x[2] 2*x[1]^2]
        ax = [a1, a2]
        @test sameForm(Hf(ax), compact(ForwardDiff.hessian(f, ax)))
    end
end

