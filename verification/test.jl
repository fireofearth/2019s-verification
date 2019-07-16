include("./helper.jl")

localModulePath = "/home/fireofearth/res/mitchell-ian/2019s-verification/verification/module"
if(!(localModulePath in LOAD_PATH))
    push!(LOAD_PATH, localModulePath)
end

using Test, Random, IntervalArithmetic
using ModalIntervalArithmetic
using AffineArithmetic

@testset "modal interval getters" begin
    X = ModalInterval(1,2) # proper
    Y = ModalInterval(2,1) # improper
    Z = ModalInterval(2,2) # real_number
    dX = dual(X)
    dY = dual(Y)
    dZ = dual(Z)
    @test inf(X) == 1 && sup(X) == 2
    @test inf(dX) == 2 && sup(dX) == 1
    @test inf(Y) == 2 && sup(Y) == 1
    @test inf(dY) == 1 && sup(dY) == 2
    @test inf(Z) == sup(Z) == 2
    @test inf(dZ) == sup(dZ) == 2
    @test dX == Y && dY == X
end

@testset "modal interval comparison" begin
    X = ModalInterval(2,3)
    Y = ModalInterval(1,4)
    dX = ModalInterval(3,2)
    dY = ModalInterval(4,1)
    # ([2,3]',∃) = [2,3] ⊆ [1,4] = ([1,4]',∃)
    # => [2,3]' ⊆ [1,4]'
    @test  X ⊆ Y  &&  Y ⊇ X
    # ([2,3]',∀) = [3,2] ⊇ [4,1] = ([1,4]',∀)
    # => [2,3]' ⊆ [1,4]'
    @test dX ⊇ dY && dY ⊆ dX 
end

@testset "modal interval additive group" begin
    A = ModalInterval(1,2)
    B = ModalInterval(3,4)
    @test A + B == ModalInterval(1+3,2+4)
    @test -A == ModalInterval(-2,-1)
    @test A - B == ModalInterval(1-4,2-3)
    @test A * B == ModalInterval(1 * 3, 2 * 4)
end

@testset "modal interval multiplicative group" begin
    #A = ModalInterval()
end

 #=
 # Affine arithmetic testing
 # 
 # TODO: constructors, getters, utility functionality
 # TODO: hardcoded tests to check boundary cases
 # TODO: tests to compare solutions with aaflib
 # TODO: improve random number generators
 # TODO: 
=#

@testset "affine arithmetic" begin
    @testset "basic" begin
        center = 3.14
        dev    = [0.75, 0.01]
        ind    = [1, 3]
        a = AAF(center, dev, ind)
        aOne = one(a)
        aZero = zero(a)
        @test a*aOne == a == aOne*a
        @test a + aZero == a == aZero + a
    end

    @testset "unitary" begin
        center1 = 100.001
        dev1    = [2.0, 1.0, -5.0, 4.0]
        ind1    = [1  , 3  ,  4,   6]
        a1 = AAF(center1, dev1, ind1)
        @test a1[0] == center1
        for ii in 1:4
            @test a1[ii] == dev1[ii]
            @test a1.indexes[ii] == ind1[ii]
        end
        @test rad(a1) == 2.0 + 1.0 + 5.0 + 4.0
        @test length(a1) == 4
        @test Interval(a1) == Interval(a1[0] - rad(a1), a1[0] + rad(a1))
    end

    @testset "equality" begin
        center1 = 1.321
        dev1    = [2.1, 0.0, 1.5, -5.3, 4.2]
        ind1    = [1  , 2,   5  ,  7,   8]
        a1 = AAF(center1, dev1, ind1)
        center2 = 1.321
        dev2    = [2.1, 0.0, 1.5, -5.3, 4.2]
        ind2    = [1  , 2,   5  ,  7,   8]
        a2 = AAF(center2, dev2, ind2)

        @test a1 == a2
        ind1[2] = 3
        a1 = AAF(center1, dev1, ind1)
        @test a1 != a2
        ind1[2] = 2
        a1 = AAF(center1, dev1, ind1)
        @test a1 == a2
        dev1[3] = 0.5
        a1 = AAF(center1, dev1, ind1)
        @test a1 != a2
        dev1[3] = 1.5
        a1 = AAF(center1, dev1, ind1)
        @test a1 == a2
    end

    @testset "equality boundary conds" begin
        center1 = 1.321
        dev1    = [2.1, 0.0, 1.5, -5.3, 4.0]
        ind1    = [1  , 2,   5  ,  7,   8]
        a1 = AAF(center1, dev1, ind1)
        center2 = 1.321
        dev2    = [2.1, 0.0, 1.5, -5.3, 4.2]
        ind2    = [1  , 2,   5  ,  7,   8]
        a2 = AAF(center2, dev2, ind2)

        @test a1 != a2
        dev1    = [2.1, 0.0, 1.5, -5.3, 4.2]
        ind1    = [1  , 2,   5  ,  7,   9]
        a1 = AAF(center1, dev1, ind1)
        dev2    = [2.1, 0.0, 1.5, -5.3, 4.2]
        ind2    = [1  , 2,   5  ,  7,   8]
        a2 = AAF(center2, dev2, ind2)
        @test a1 != a2
        dev1    = [1.1, 0.0, 1.5, -5.3, 4.2]
        ind1    = [1  , 2,   5  ,  7,   8]
        a1 = AAF(center1, dev1, ind1)
        dev2    = [2.1, 0.0, 1.5, -5.3, 4.2]
        ind2    = [1  , 2,   5  ,  7,   8]
        a2 = AAF(center2, dev2, ind2)
        @test a1 != a2
        dev1    = [2.1, 0.0, 1.5, -5.3, 4.2]
        ind1    = [2  , 3,   5  ,  7,   8]
        a1 = AAF(center1, dev1, ind1)
        dev2    = [2.1, 0.0, 1.5, -5.3, 4.2]
        ind2    = [1  , 3,   5  ,  7,   8]
        a2 = AAF(center2, dev2, ind2)
        @test a1 != a2
        center1 = 0.321
        dev1    = [2.1, 0.0, 1.5, -5.3, 4.2]
        ind1    = [1  , 2,   5  ,  7,   8]
        a1 = AAF(center1, dev1, ind1)
        center2 = 1.321
        dev2    = [2.1, 0.0, 1.5, -5.3, 4.2]
        ind2    = [1  , 2,   5  ,  7,   8]
        a2 = AAF(center2, dev2, ind2)
        @test a1 != a2
    end
end

 #=
 # Computes the new noise coefficients, given two collections of noise
 #
 # Example: given function f
 #  ind1 = [ 2,   4,   6, ...]
 #  dev1 = [10.2, 5.4, 3.1,..] corresponding to 10.2ϵ₂ + 5.4ϵ₄ + 3.1ϵ₆ + ...
 #
 #  ind2 = [ 3,   4,    6, ...]
 #  dev2 = [72.0, 2.1, 27.5...] corresponding to 72.0ϵ₃ + 2.1ϵ₄ + 27.5ϵ₆ + ...
 #
 #  returns [(2, f(10.2, 0)), (3, f(0, 72.0)), (4, f(5.4,2.1)), (6, f(3.1, 27.5)),...]
 #  corresponding f(10.2, 0))ϵ₂ + f(0, 72.0)ϵ₃ + f(5.4,2.1)ϵ₄ + f(3.1, 27.5)ϵ₆ + ...
 #
 # Returns:
 #  An array, with each index ii a tuple with first and second entries corresponding 
 #  to index ii and deviation ii respectively
=#
function getGroupSolCoeffs(
    f::Function, ind1::Vector, ind2::Vector ,dev1::Vector, dev2::Vector
)
    mm   = max(last(ind1), last(ind2))
    rg   = Array(1:mm)
    addSol = Vector{Union{Nothing,AbstractFloat}}(nothing, mm)
    for ii in rg
        idx1 = findfirst(x -> x == ii, ind1)
        idx2 = findfirst(x -> x == ii, ind2)
        if(idx1 != nothing && idx2 != nothing)
            addSol[ii] = f(dev1[idx1], dev2[idx2])
        elseif(idx1 != nothing)
            addSol[ii] = f(dev1[idx1], 0)
        elseif(idx2 != nothing)
            addSol[ii] = f(0, dev2[idx2])
        end
    end
    return sort(
        filter(x -> x[2] != nothing, tuple.(rg,addSol)), 
        by = y -> y[1]
    )
end

@testset "affine arithmetic automatic" begin
    center1 = rand(Float64) + rand(0:99)
    n1 = rand(6:10)
    dev1 = rand(Float64,n1) .+ rand(-9:9, n1)
    ind1 = sort(shuffle(Array(1:20))[1:n1])
    a1 = AAF(center1,dev1,ind1)

    center2 = rand(Float64) + rand(0:99)
    n2 = rand(6:10)
    dev2 = rand(Float64,n2) .+ rand(-9:9, n2)
    ind2 = sort(shuffle(Array(1:20))[1:n2])
    a2 = AAF(center2,dev2,ind2)

    rn1 = rand(Float64) + rand(0:99)

    @testset "affine-affine addition" begin
        aAddn = a1 + a2
        @test aAddn[0] ≈ center1 + center2
        addnSol = getGroupSolCoeffs(+,ind1,ind2,dev1,dev2)
        @test length(aAddn) == length(aAddn.indexes) == length(addnSol)
        addnComp = [(aAddn.indexes[ii], dev) for (ii, dev) in enumerate(aAddn.deviations)]
        for ((compIdx, compDev), (solIdx, solDev)) in tuple.(addnComp,addnSol)
            @test compIdx == solIdx
            @test compDev ≈ solDev
        end
    end

    @testset "affine-affine substraction" begin
        aSubt = a1 - a2
        @test aSubt[0] ≈ center1 - center2
        subtSol = getGroupSolCoeffs(-,ind1,ind2,dev1,dev2)
        @test length(aSubt) == length(aSubt.indexes) == length(subtSol)
        subtComp = [(aSubt.indexes[ii], dev) for (ii, dev) in enumerate(aSubt.deviations)]
        for ((compIdx, compDev), (solIdx, solDev)) in tuple.(subtComp,subtSol)
            @test compIdx == solIdx
            @test compDev ≈ solDev
        end
    end

    @testset "affine + constant" begin
        aAddcr = a1 + rn1
        @test aAddcr[0] == center1 + rn1
        @test length(aAddcr) == length(aAddcr.indexes) == n1
        for ii in 1:n1
            @test aAddcr[ii] ≈ dev1[ii]
            @test aAddcr.indexes[ii] == ind1[ii]
        end
    end

    @testset "constant + affine" begin
        aAddcl = rn1 + a1
        @test aAddcl[0] == rn1 + center1
        @test length(aAddcl) == length(aAddcl.indexes) == n1
        for ii in 1:n1
            @test aAddcl[ii] ≈ dev1[ii]
            @test aAddcl.indexes[ii] == ind1[ii]
        end
    end

    @testset "-affine" begin
        aNeg = -a2
        @test aNeg[0] == -center2
        @test length(aNeg) == length(aNeg.indexes) == n2
        for ii in 1:n2
            @test aNeg[ii] ≈ -dev2[ii]
            @test aNeg.indexes[ii]≈ ind2[ii]
        end
    end

    @testset "affine - constant" begin
        aSubtcr = a1 - rn1
        @test aSubtcr[0] == center1 - rn1
        @test length(aSubtcr) == length(aSubtcr.indexes) == n1
        for ii in 1:n1
            @test aSubtcr[ii] ≈ dev1[ii]
            @test aSubtcr.indexes[ii] == ind1[ii]
        end
    end

    @testset "constant - affine" begin
        aSubcl = rn1 - a1
        @test aSubcl[0] == rn1 - center1
        @test length(aSubcl) == length(aSubcl.indexes) == n1
        for ii in 1:n1
            @test aSubcl[ii] ≈ dev1[ii]
            @test aSubcl.indexes[ii] == ind1[ii]
        end
    end

    @testset "affine * constant" begin
        aMultcr = a1 * rn1
        @test aMultcr[0] == center1 * rn1
        @test length(aMultcr) == length(aMultcr.indexes) == n1
        for ii in 1:n1
            @test aMultcr[ii] ≈ dev1[ii] * rn1
            @test aMultcr.indexes[ii] == ind1[ii]
        end
    end

    @testset "constant * affine" begin
        aMultcl = rn1 * a1
        @test aMultcl[0] == rn1 * center1
        @test length(aMultcl) == length(aMultcl.indexes) == n1
        for ii in 1:n1
            @test aMultcl[ii] ≈ rn1 * dev1[ii]
            @test aMultcl.indexes[ii] == ind1[ii]
        end
    end

    # xy = x₀ŷ₀ + ½∑ᴺᵢxᵢyᵢ + ∑ᴺᵢ(xᵢy₀+yᵢx₀)ϵᵢ + [(∑ᴺᵢ|xᵢ|)(∑ᴺᵢ|yᵢ|) - ½∑ᴺᵢ|xᵢyᵢ|]μₖ
    # TODO: fix the tests
    @testset "affine * affine" begin
        aMult = a1 * a2
        #@test aMult[0] ≈ center1*center2 + 0.5*(dev1'*dev2)
        multSol = getGroupSolCoeffs((x,y) -> (x*center2 + y*center1) ,ind1,ind2,dev1,dev2)
        @test length(aMult) == length(aMult.indexes) == length(multSol) + 1
        for (ii, (idx, dev)) in enumerate(multSol)
            @test aMult[ii] ≈ dev
            @test aMult.indexes[ii]  == idx
        end
        kk = length(multSol) + 1
        #@test aMult[kk] ≈ (sum(abs.(dev1)) * sum(abs.(dev2))) - 0.5*sum(abs.(dev1 .* dev2))
        @test aMult.indexes[kk]   == getLastAAFIndex()
    end
end






