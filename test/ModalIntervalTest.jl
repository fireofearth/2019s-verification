using Test, Random, Dates
using IntervalArithmetic
using ModalIntervalArithmetic
using Logging

 #=
 # Modal Interval testing
=#

Random.seed!(Dates.value(Dates.now()))
MIN = -10^4
MAX =  10^4

 #=
 # Helper functions
 #
 # Notation:
 # R := {[x₁, x₂] : x₁ ≤ x₂}
 # P := {[x₁, x₂] : x₁, x₂ ≥ 0}
 # D := {[x₁, x₂] : x₁ ≤ 0 ≤ x₂ }
 # NP := -P, DZ = dual(Z), DR = dual(R)
=# 
random_interval()  = ModalInterval(rand(MIN .. MAX), rand(MIN .. MAX))
random_p_interval()  = ModalInterval(rand(0.0 .. MAX), rand(0.0 .. MAX))
random_np_interval() = -random_p_interval()
random_z_interval()  = ModalInterval(rand(MIN .. 0.0), rand(0.0 .. MAX))
random_dz_interval() = dual(random_z_interval())
random_r_interval()  = rand(MIN .. (MAX - 100)) |> x₁ -> (x₁, rand(x₁ .. MAX)) |> v -> ModalInterval(v...)
random_dr_interval() = dual(random_r_interval())

@testset "modal interval general" begin
    xₘ = rand((MIN + 100):(MAX - 100))
    x₁ = rand(MIN:xₘ)
    x₂ = rand((xₘ + 1):MAX)
    @test x₁ < x₂

    @testset "proper interval" begin
        X = ModalInterval(x₁, x₂)
        @test inf(X)  == x₁ && sup(X) == x₂
        @test mid(X) == (x₁ + x₂) / 2
        @test mod(X)  == PROPER
        @test !isreal(X)
        @test isproper(X)
        @test !isimproper(X)
        @test Interval(X) == Interval(x₁, x₂)
        @test dual(X) == ModalInterval(x₂, x₁)
        @test prop(X) == X
    end

    @testset "improper integral" begin
        X = ModalInterval(x₂, x₁)
        @test inf(X)  == x₂ && sup(X) == x₁
        @test mid(X) == (x₁ + x₂) / 2
        @test mod(X)  == IMPROPER
        @test !isreal(X)
        @test !isproper(X)
        @test isimproper(X)
        @test dual(X) == prop(X)
        @test_throws DomainError Interval(X)
    end

    @testset "real number integral" begin
        X = ModalInterval(x₁, x₁)
        @test inf(X)  == x₁ && sup(X) == x₁
        @test mid(X) == Float64(x₁)
        @test mod(X)  == REAL_NUMBER
        @test isreal(X)
        @test isproper(X)
        @test isimproper(X)
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

@testset "modal interval set operations" begin
    q₂ = rand(MIN .. MAX)
    q₃ = rand(q₂ .. MAX)
    q₁ = rand(MIN .. q₂)
    x₁ = rand(MIN .. q₁)
    x₂ = rand(q₁ .. q₂)
    x₃ = rand(q₂ .. q₃)
    x₄ = rand(q₃ .. MAX)
    @test x₁ ≤ x₂ ≤ x₃ ≤ x₄

    @testset "comparison hardcode" begin
        X = ModalInterval(1, 4)
        Y = ModalInterval(2, 3)
        Y ⊆ X
        @test Y ⊂ X
        # ([2,3]',∃) = [2,3] ⊆ [1,4] = ([1,4]',∃)
        @test Y ⊆ X
        @test !(Y ⊇ X)
        @test Y ⊉ X
        # ([2,3]',∀) = [3,2] ⊇ [4,1] = ([1,4]',∀)
        @test dual(Y) ⊇ dual(X)
        @test !(dual(Y) ⊆ dual(X))
        @test dual(Y) ⊈ dual(X)
        @test !(X ⊆ dual(Y))
        @test !(Y ⊆ dual(Y))
        W = ModalInterval(1, 3)
        Z = ModalInterval(2, 4)
        # ([1,3]',∀) = [3,1] ⊆ [2,4] = ([2,4]',∃)
        @test dual(W) ⊆ Z
        @test !(dual(W) ⊇ Z)
        # ([2,4]',∀) = [4,2] ⊆ [1,3] = ([1,3]',∃)
        @test dual(Z) ⊆ W
        @test !(dual(Z) ⊇ W)
        # ([2,2]',∃) = [2,2] ⊆ [2,2] = ([2,2]',∀)
        @test ModalInterval(2, 2) ⊆ ModalInterval(2, 2)
    end

    @testset "comparison gen." begin
        X = ModalInterval(x₁, x₄)
        Y = ModalInterval(x₂, x₃)
        @test X == ModalInterval(x₁, x₄)
        @test Y ⊆ X
        @test dual(X) ⊆ dual(Y)
        # ...
    end
end

@testset "test the test itself" begin
    @testset "test interval cases" begin
        X = random_p_interval()
        @test inf(X) ≥ 0 && sup(X) ≥ 0
        X = random_np_interval()
        @test inf(X) ≤ 0 && sup(X) ≤ 0
        X = random_z_interval()
        @test inf(X) ≤ 0 && sup(X) ≥ 0
        X = random_dz_interval()
        @test inf(X) ≥ 0 && sup(X) ≤ 0
        X = random_r_interval()
        @test inf(X) ≤ sup(X)
        X = random_dr_interval()
        @test inf(X) ≥ sup(X)
    end
end

@testset "modal interval arithmetic operations" begin
    x₁ = rand(MIN .. MAX)
    x₂ = rand(MIN .. MAX)
    x₃ = rand(MIN .. MAX)
    x₄ = rand(MIN .. MAX)

    @testset "additive group gen." begin
        X = ModalInterval(x₁, x₂)
        Y = ModalInterval(x₃, x₄)
        @test X + Y == ModalInterval(x₁ + x₃, x₂ + x₄)
        @test -X == ModalInterval(-x₂, -x₁)
        @test X - Y == ModalInterval(x₁ - x₄, x₂ - x₃)
        @test X - dual(X) == zero(ModalInterval{Float64})
    end

    @testset "additive group with real numbers gen." begin
        X = ModalInterval(x₁, x₂)
        @test X + x₃ == ModalInterval(x₁ + x₃, x₂ + x₃)
        @test x₃ + X == ModalInterval(x₁ + x₃, x₂ + x₃)
        @test X - x₃ == ModalInterval(x₁ - x₃, x₂ - x₃)
        @test x₃ - X == ModalInterval(x₃ - x₂, x₃ - x₁)
    end

     #=
     # Using the Kaucher × table we want to test
     # NP := -P, DZ = dual(Z)
     # P * P, NP * NP, Z * Z, DZ * DZ
     # P * NP
    =#
    @testset "multiplication gen." begin
        # P * {P, NP, Z, DZ}
        X = random_p_interval()
        Y = random_p_interval()
        @test X * Y == ModalInterval(inf(X)*inf(Y), sup(X)*sup(Y))
        Y = random_z_interval()
        @test X * Y == Y * X == ModalInterval(sup(X)*inf(Y), sup(X)*sup(Y))
        Y = random_np_interval()
        @test X * Y == Y * X == ModalInterval(sup(X)*inf(Y), inf(X)*sup(Y))
        Y = random_dz_interval()
        @test X * Y == Y * X == ModalInterval(inf(X)*inf(Y), inf(X)*sup(Y))
        # Z * {Z, NP, DZ}
        X = random_z_interval()
        Y = random_z_interval()
        @test X * Y == ModalInterval(min(inf(X)*sup(Y),sup(X)*inf(Y)), 
                                              max(inf(X)*inf(Y), sup(X)*sup(Y)))
        Y = random_np_interval()
        @test X * Y == Y * X == ModalInterval(sup(X)*inf(Y), inf(X)*inf(Y))
        Y = random_dz_interval()
        @test X * Y == Y * X == zero(typeof(X))
        # NP * {NP, DZ}
        X = random_np_interval()
        Y = random_np_interval()
        @test X * Y == ModalInterval(sup(X)*sup(Y), inf(X)*inf(Y))
        Y = random_dz_interval()
        @test X * Y == Y * X == ModalInterval(sup(X)*sup(Y), sup(X)*inf(Y))
        # DZ * DZ
        X = random_dz_interval()
        Y = random_dz_interval()
        @test X * Y == ModalInterval(max(inf(X)*inf(Y), sup(X)*sup(Y)),
                                     min(inf(X)*sup(Y), sup(X)*inf(Y)))
    end

    @testset "interval modules" begin
        # row column
        X11 = random_interval()
        X12 = random_interval()
        X21 = random_interval()
        X22 = random_interval()
        Y11 = random_interval()
        Y12 = random_interval()
        Y21 = random_interval()
        Y22 = random_interval()
        Xvec = [X11; X21]
        Yvec = [Y11; Y21]
        Xmat = [X11 X12; X21 X22]
        Ymat = [Y11 Y12; Y21 Y22]

        @testset "vector eval" begin
            @test isproper([random_r_interval(); 
                            random_r_interval();
                            random_r_interval()])
            @test isimproper([random_dr_interval(); 
                              random_dr_interval();
                              random_dr_interval()])
        end

        @testset "additive module" begin
            @test Xvec + Yvec == [X11 + Y11; X21 + Y21]
            @test Xvec - Yvec == [X11 - Y11; X21 - Y21]
            @test Xvec .+ Y11 == [X11 + Y11; X21 + Y11]
            @test Xvec .- Y11 == [X11 - Y11; X21 - Y11]
            @test Xvec .+ x₁ == [X11 + x₁; X21 + x₁]
            @test Xvec .- x₁ == [X11 - x₁; X21 - x₁]
        end

         #=
         # Xmat * Yvec produces types ModalInterval{Int} ?
         #
        =#
        @testset "vector-matrix and matrix-matrix mult" begin
            @test Xmat * Yvec == [X11*Y11 + X12*Y21; X21*Y11 + X22*Y21]
            @test Xmat * Ymat == [X11*Y11 + X12*Y21 X11*Y12 + X12*Y22;
                                  X21*Y11 + X22*Y21 X21*Y12 + X22*Y22]
        end
    end
end

