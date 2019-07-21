include("./helper.jl")

disp(msg) = print("$(msg)\n")
debug() = print("DEBUG\n")

using Test, Random, IntervalArithmetic
using ModalIntervalArithmetic
using AffineArithmetic
using ForwardDiff

 #=
 # Modal Interval testing
 #
 # TODO:
=#

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
 # Specifications:
 # - algebriac identities one, zero, isone, iszero
 # - unitary operations and constructors
 # - detailed constructors
 #
 # - when comparing two affine forms we require the indexes to match thus we call 
 # resetLastAffineIndex() before generating them. This is especially true when testing for
 # ForwardDiff
 #
 # TODO: constructors, getters, utility functionality
 # TODO: hardcoded tests to check boundary cases
 # TODO: tests to compare solutions with aaflib
 # TODO: improve random number generators
 # TODO: tests to check compatibility with ForwardDiff
 #
 # TODO: it is not possible to recreate which indexes are assigned to noise symbols when
 # using ForwardDiff. Now using `sameForm()` instead of `==`
=#

 #=
 # Comparable to `==` except we do not check indexes. Used for testing only.
=#
function sameForm(a::Affine, p::Affine; tol::Float64=affineTOL)
    if(length(a) != length(p))
        return false
    end

    if(abs(a[0]) < 1 && abs(p[0]) < 1)
        if(abs(a[0] - p[0]) > tol)
            return false
        end
    else
        if(abs((a[0] - p[0]) / (a[0] + p[0])) > tol)
            return false
        end
    end

    for i in 1:length(a)
        if(abs(a[i]) < 1 && abs(p[i]) < 1)
            if(abs(a[i] - p[i]) > tol)
                return false
            end
        else
            if(abs(a[i] - p[i]) / 
                   (abs(a[i]) + abs(p[i])) > tol)
                return false
            end
        end
    end
    
    return true
end

 #=
 # Affine Arithmetic Common
 # All functionality except for elementary functions and binary operations
=#
@testset "affine arithmetic common" begin
    @testset "algebriac identities" begin
        center = 3.14
        dev    = [0.75, 0.01]
        ind    = [1, 3]
        a = Affine(center, dev, ind)
        aOne = one(a)
        aZero = zero(a)
        @test a*aOne == a == aOne*a
        @test a + aZero == a == aZero + a
    end

    @testset "constructors; unitary" begin
        center1 = 100.001
        dev1    = [2.0, 1.0, -5.0, 4.0]
        ind1    = [1  , 3  ,  4,   6]
        a1 = Affine(center1, dev1, ind1)
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
        a1 = Affine(center1, dev1, ind1)
        center2 = 1.321
        dev2    = [2.1, 0.0, 1.5, -5.3, 4.2]
        ind2    = [1  , 2,   5  ,  7,   8]
        a2 = Affine(center2, dev2, ind2)

        @test a1 == a2
        ind1[2] = 3
        a1 = Affine(center1, dev1, ind1)
        @test a1 != a2
        ind1[2] = 2
        a1 = Affine(center1, dev1, ind1)
        @test a1 == a2
        dev1[3] = 0.5
        a1 = Affine(center1, dev1, ind1)
        @test a1 != a2
        dev1[3] = 1.5
        a1 = Affine(center1, dev1, ind1)
        @test a1 == a2
    end

    @testset "equality boundary conds" begin
        center1 = 1.321
        dev1    = [2.1, 0.0, 1.5, -5.3, 4.0]
        ind1    = [1  , 2,   5  ,  7,   8]
        a1 = Affine(center1, dev1, ind1)
        center2 = 1.321
        dev2    = [2.1, 0.0, 1.5, -5.3, 4.2]
        ind2    = [1  , 2,   5  ,  7,   8]
        a2 = Affine(center2, dev2, ind2)

        @test a1 != a2
        dev1    = [2.1, 0.0, 1.5, -5.3, 4.2]
        ind1    = [1  , 2,   5  ,  7,   9]
        a1 = Affine(center1, dev1, ind1)
        dev2    = [2.1, 0.0, 1.5, -5.3, 4.2]
        ind2    = [1  , 2,   5  ,  7,   8]
        a2 = Affine(center2, dev2, ind2)
        @test a1 != a2
        dev1    = [1.1, 0.0, 1.5, -5.3, 4.2]
        ind1    = [1  , 2,   5  ,  7,   8]
        a1 = Affine(center1, dev1, ind1)
        dev2    = [2.1, 0.0, 1.5, -5.3, 4.2]
        ind2    = [1  , 2,   5  ,  7,   8]
        a2 = Affine(center2, dev2, ind2)
        @test a1 != a2
        dev1    = [2.1, 0.0, 1.5, -5.3, 4.2]
        ind1    = [2  , 3,   5  ,  7,   8]
        a1 = Affine(center1, dev1, ind1)
        dev2    = [2.1, 0.0, 1.5, -5.3, 4.2]
        ind2    = [1  , 3,   5  ,  7,   8]
        a2 = Affine(center2, dev2, ind2)
        @test a1 != a2
        center1 = 0.321
        dev1    = [2.1, 0.0, 1.5, -5.3, 4.2]
        ind1    = [1  , 2,   5  ,  7,   8]
        a1 = Affine(center1, dev1, ind1)
        center2 = 1.321
        dev2    = [2.1, 0.0, 1.5, -5.3, 4.2]
        ind2    = [1  , 2,   5  ,  7,   8]
        a2 = Affine(center2, dev2, ind2)
        @test a1 != a2
    end

    @testset "compact" begin
        center = 12.0
        dev = [0.0, 2.0, 3.0, 0.0, 1.0, 0.0]
        ind = [1,   2,   3,   5,   8,   9]
        a = Affine(center, dev, ind)
        @test compact(a) == Affine(center, [2.0, 3.0, 1.0], [2, 3, 8])
        a = Affine(center, [0.0], [1])
        @test compact(a) == Affine(center)
        a = Affine(center)
        @test compact(a) == Affine(center)
    end
end

 #=
 # Verify Affine-Constant operations where constant can be floats or ints.
=#
@testset "affine arithmetic constant ops" begin
    center = 12.0
    dev    = [3.0, 6.0, 9.0]
    ind    = [1,    3,    4]
    a = Affine(center, dev, ind)

    @testset "addition / subtraction" begin
        @test a + 2.0 == Affine(14.0, dev, ind)
        @test a + 2   == Affine(14.0, dev, ind)
        @test 2.0 + a == Affine(14.0, dev, ind)
        @test 2   + a == Affine(14.0, dev, ind)
        @test a - 2.0 == Affine(10.0, dev, ind)
        @test a - 2   == Affine(10.0, dev, ind)
        @test 2.0 - a == Affine(-10.0, dev, ind)
        @test 2   - a == Affine(-10.0, dev, ind)
    end

    @testset "multiplication" begin
        @test a * 2.0 == Affine(24.0, [6.0, 12.0, 18.0], ind)
        @test a * 2   == Affine(24.0, [6.0, 12.0, 18.0], ind)
        @test 2.0 * a == Affine(24.0, [6.0, 12.0, 18.0], ind)
        @test 2   * a == Affine(24.0, [6.0, 12.0, 18.0], ind)
    end

    @testset "division" begin
        @test a / 3.0 == Affine(4.0, [1.0, 2.0, 3.0], ind)
        @test a / 3   == Affine(4.0, [1.0, 2.0, 3.0], ind)
    end
end

@testset "affine arithmetic hardcode ops" begin
    center  = 26.10
    dev     = [2.11, -3.03, 4.59, 1.0, -10.0]
    ind     = [1,     3,    5,    8,    10]

    # RINO uses CHEBYSHEV
    @testset "inverse" begin
        resetLastAffineIndex()
        nCenter = 0.06305953707868751
        nDev    = [-0.00839043, 0.0120488, -0.0182522, -0.00397651, 0.0397651, 0.0407272]
        nInd    = [ 1,          3,          5,          8,         10,        11]
        a = Affine(center)
        @test inv(a) == Affine(1 / center)
        a = Affine(center, dev, ind)
        @test isapprox(inv(a), Affine(nCenter, nDev, nInd); tol=10E-8)
        resetLastAffineIndex()
        a = Affine(center, dev, ind)
        @test isapprox(a^(-1), Affine(nCenter, nDev, nInd); tol=10E-8)
    end

    @testset "inverse equivalencies" begin
        resetLastAffineIndex()
        a     = Affine(center, dev, ind)
        inva1 = inv(a)
        resetLastAffineIndex()
        a     = Affine(center, dev, ind)
        inva2 = 1/a
        resetLastAffineIndex()
        a     = Affine(center, dev, ind)
        inva3 = a^(-1)
        @test inva1 == inva2 == inva3
    end

    @testset "power" begin
        resetLastAffineIndex()
        nCenter = 896.07645
        nDev    = [110.142, -158.166, 239.598, 52.2, -522.0, 214.86645]
        nInd    = [  1,        3,       5,      8,     10,    11]
        a = Affine(center, dev, ind)
        @test a^1 == a
        @test isapprox(a^2, Affine(nCenter, nDev, nInd); tol=10E-8)
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
    a1 = Affine(center1,dev1,ind1)

    center2 = rand(Float64) + rand(0:99)
    n2 = rand(6:10)
    dev2 = rand(Float64,n2) .+ rand(-9:9, n2)
    ind2 = sort(shuffle(Array(1:20))[1:n2])
    a2 = Affine(center2,dev2,ind2)

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
        @test aMult.indexes[kk] == getLastAffineIndex()
    end
end

 #=
 # Affine Arithmetic ForwardDiff
 # Tests compatibility with ForwardDiff
 # 
 # Remark: ForwardDiff almost works as is; must make affines compact by removing zeros
=#
@testset "affine arithmetic ForwardDiff" begin
    centers = [32.1, 27.3, 58.0]
    devs    = [[0.1, -0.2, 1.5, -2.0],
               [10.0, 0.5, 1.0], 
               [-3.33, 9.0, -1.5, 5.25]]
    inds    = [[1, 3, 4, 5],
              [1, 4, 6],
              [2, 3, 5, 6]]
    a1 = Affine(32.1, [0.1, -0.2, 1.5, -2.0], [1, 3, 4, 5])
    a2 = Affine(27.3, [10.0, 0.5, 1.0], [1, 4, 6])
    a3 = Affine(4.0,  [-3.33, 9.0, -1.5, 5.25], [2, 3, 5, 6])
    
    @testset "derivative simple" begin
        f(x::Real)  = 1.0 - x^2 + x
        df(x::Real) = -2*x + 1.0
        @test ForwardDiff.derivative(f,a1) == df(a1)
    end

    @testset "derivative of 1/x" begin
        f(x::Real)  = 1/x
        df(x::Real) = -(1 /x /x)
        resetLastAffineIndex()
        a1     = Affine(centers[1], devs[1], inds[1])
        actual = df(a1)
        resetLastAffineIndex()
        a1     = Affine(centers[1], devs[1], inds[1])
        res    = ForwardDiff.derivative(f, a1)
        @test res == actual
        resetLastAffineIndex()
        a2     = Affine(centers[2], devs[2], inds[2])
        actual = df(a2)
        resetLastAffineIndex()
        a2     = Affine(centers[2], devs[2], inds[2])
        res    = ForwardDiff.derivative(f, a2)
        @test res == actual
        resetLastAffineIndex()
        a3     = Affine(centers[3], devs[3], inds[3])
        actual = df(a3)
        resetLastAffineIndex()
        a3     = Affine(centers[3], devs[3], inds[3])
        res    = ForwardDiff.derivative(f, a3)
        @test res == actual
    end

    @testset "derivative of x^2" begin
        f(x::Real)  = x^2
        df(x::Real) = 2*x
        a1          = Affine(centers[1], devs[1], inds[1])
        @test df(a1) == ForwardDiff.derivative(f, a1)
    end

     #=
     # Remark: `df(x::Real) = 1` also works
    =#
    @testset "derivative of x^n" begin
        f(x::Real)  = x
        df(x::Real) = Affine(1)
        resetLastAffineIndex()
        a1          = Affine(centers[1], devs[1], inds[1])
        @test sameForm(df(a1), ForwardDiff.derivative(f, a1))
        for n in 1:5
            fn(x::Real)  = x^n
            dfn(x::Real) = n * (x^(n-1))
            @test sameForm(dfn(a1), ForwardDiff.derivative(fn, a1))
        end
    end

    @testset "gradient" begin
        f(x::Vector)  = x[1]*x[3] + 2.0*x[2]*x[1] - x[3]*x[2]
        gf(x::Vector) = [x[3] + 2.0*x[2], 2.0*x[1] - x[3], x[1] - x[2]]
        ax = [a1, a2, a3]
        res = ForwardDiff.gradient(f, ax)
        @test compact(res) == gf(ax)
    end

    @testset "hessian" begin
        f(x::Vector) = x[1]*x[2]^2 - x[2]*x[1]^2
        Hf(x::Vector) = [-2*x[2] (2*x[2] - 2*x[1]); (2*x[2] - 2*x[1]) 2*x[1]]
        ax = [a1, a2]
        res = ForwardDiff.hessian(f, ax)
        @test compact(res) == Hf(ax)
    end

#=
    @testset "jacobian" begin
        f(x::Vector) = [(x[1]*x[2]) /x[3], (x[1]*x[2])*x[3], x[1]^2 + x[2]^2 + x[3]^2]
        Jf(x::Vector) = [x[2]/x[3]  x[1]/x[3]  -((x[1]*x[2]) /x[3] /x[3]); 
                         x[2]*x[3]  x[1]*x[3]  x[1]*x[2]; 
                         2*x[1]     2*x[2]     2*x[3]]
        resetLastAffineIndex()
        a1     = Affine(centers[1], devs[1], inds[1])
        a2     = Affine(centers[2], devs[2], inds[2])
        a3     = Affine(centers[3], devs[3], inds[3])
        ax = [a1, a2, a3]
        actual = Jf(ax)
        resetLastAffineIndex()
        a1     = Affine(centers[1], devs[1], inds[1])
        a2     = Affine(centers[2], devs[2], inds[2])
        a3     = Affine(centers[3], devs[3], inds[3])
        ax     = [a1, a2, a3]
        res    = ForwardDiff.jacobian(f, ax)
        #@test res == actual
        @test compact(res) == actual
    end
    =#
end

#=
@testset "affine arithmetic ForwardDiff 2" begin
    centers = [180.1, 205.25, 58.0]
    devs    = [[0.1, -0.2, 1.5, -2.0],
               [10.0, 0.5, 1.0], 
               [-3.33, 9.0, -1.5, 5.25]]
    inds    = [[1, 3, 4, 5],
              [1, 4, 6],
              [2, 3, 5, 6]]

     #=
     # First coord OK
    =#
    @testset "gradient of f(x, y) = x / y" begin
        f(x::Vector) = x[1]/x[2]
        gf(x::Vector) = [inv(x[2]), -x[1]*inv(x[2])]
        resetLastAffineIndex()
        a1     = Affine(centers[1], devs[1], inds[1])
        a2     = Affine(centers[2], devs[2], inds[2])
        ax     = [a1, a2]
        actual = gf(ax)
        disp(actual)
        resetLastAffineIndex()
        a1     = Affine(centers[1], devs[1], inds[1])
        a2     = Affine(centers[2], devs[2], inds[2])
        ax     = [a1, a2]
        res    = ForwardDiff.gradient(f, ax)
        @test compact(res) == actual
    end
end
=#

#=
Affine(0.00487977, [-0.000238122, -1.19061e-5, -2.38122e-5, 7.67756e-6], [1, 4, 6, 9]), 
Affine(-0.00430362, [0.000419778, 4.77924e-6, -1.47359e-5, 4.77924e-5, 4.22167e-5, 1.35495e-5, -2.88369e-5, -1.105e-5], [1, 3, 4, 5, 6, 10, 11, 12])]

# 1st coord OK
Affine(0.00487977, [-0.000238122, -1.19061e-5, -2.38122e-5, 7.67756e-6], [1, 4, 6, 7])
Affine(-0.00428858, [0.000416162, 4.76244e-6, -1.47911e-5, 4.76244e-5, 4.18543e-5, -6.7474e-6, -5.17665e-6, -6.74733e-6, -1.82591e-5], [1, 3, 4, 5, 6, 8, 9, 10, 11])]
=#


#=
# Jacobian: First round
Affine(0.497556, [0.182576, 0.0303035, -0.0819013, 0.00912878, 0.0136502, -0.0295182, 0.0293768, 0.084641], [1, 2, 3, 4, 5, 6, 10, 14]) 
Affine(0.586068, [0.00182576, 0.0356316, -0.099953, 0.0273863, -0.0204649, -0.0561759, 0.0345419, 0.0282574], [1, 2, 3, 4, 5, 6, 10, 15]) 
Affine(-0.396579, [-0.147186, -0.0700039, 0.191682, -0.0259154, -0.00670904, 0.095772, -0.019698, 0.063477, -0.213541, -0.397266], [1, 2, 3, 4, 5, 6, 7, 11, 12, 13]); 

Affine(1586.03, [580.0, -90.909, 245.7, 29.0, -40.95, 201.325, 216.795], [1, 2, 3, 4, 5, 6, 19]) 
Affine(1861.8, [5.8, -106.893, 277.3, 87.0, -164.15, 168.525, 72.504], [1, 2, 3, 4, 5, 6, 20]) 
Affine(876.705, [323.73, -5.46, 57.0, -54.6, 32.1, 43.325], [1, 3, 4, 5, 6, 17]); 

Affine(64.2, [0.2, -0.4, 3.0, -4.0], [1, 3, 4, 5]) 
Affine(54.6, [20.0, 1.0, 2.0], [1, 4, 6]) 
Affine(116.0, [-6.66, 18.0, -3.0, 10.5], [2, 3, 5, 6])

==============

Affine(0.497556, [0.182576, 0.0303035, -0.0819013, 0.00912878, 0.0136502, -0.0295182, 0.0293768, 0.084641], [1, 2, 3, 4, 5, 6, 7, 8]) 
# center
Affine(0.585568, [0.00182576, 0.0356316, -0.099953, 0.0273863, -0.0204649, -0.0561759, 0.0345419, 0.0277574], [1, 2, 3, 4, 5, 6, 9, 10]) 
# everything
Affine(-0.292556, [-0.107912, -0.0355038, 0.0977762, -0.0190003, 0.00220758, 0.0452743, -0.0144419, -0.0172241, -0.0695718, -0.0171939, -0.130834], [1, 2, 3, 4, 5, 6, 11, 12, 13, 14, 15]); 

Affine(1586.03, [580.0, -90.909, 245.7, 29.0, -40.95, 201.325, 216.795], [1, 2, 3, 4, 5, 6, 16]) 
# center
Affine(1863.3, [5.8, -106.893, 277.3, 87.0, -164.15, 168.525, 71.004], [1, 2, 3, 4, 5, 6, 17]) 
Affine(876.705, [323.73, -5.46, 57.0, -54.6, 32.1, 43.325], [1, 3, 4, 5, 6, 18]); 

Affine(64.2, [0.2, -0.4, 3.0, -4.0], [1, 3, 4, 5]) 
Affine(54.6, [20.0, 1.0, 2.0], [1, 4, 6]) 
Affine(116.0, [-6.66, 18.0, -3.0, 10.5], [2, 3, 5, 6])]
=#

#=
# Jacobian: Second round
Affine(0.497556, [0.182576, 0.0303035, -0.0819013, 0.00912878, 0.0136502, -0.0295182, 0.0293768, 0.084641], [1, 2, 3, 4, 5, 6, 10, 14]) 
Affine(0.586068, [0.00182576, 0.0356316, -0.099953, 0.0273863, -0.0204649, -0.0561759, 0.0345419, 0.0282574], [1, 2, 3, 4, 5, 6, 10, 15]) 
Affine(-0.396579, [-0.147186, -0.0700039, 0.191682, -0.0259154, -0.00670904, 0.095772, -0.019698, 0.063477, -0.213541, -0.397266], [1, 2, 3, 4, 5, 6, 7, 11, 12, 13]); 

Affine(1586.03, [580.0, -90.909, 245.7, 29.0, -40.95, 201.325, 216.795], [1, 2, 3, 4, 5, 6, 19]) 
Affine(1861.8, [5.8, -106.893, 277.3, 87.0, -164.15, 168.525, 72.504], [1, 2, 3, 4, 5, 6, 20]) 
Affine(876.705, [323.73, -5.46, 57.0, -54.6, 32.1, 43.325], [1, 3, 4, 5, 6, 17]); 

Affine(64.2, [0.2, -0.4, 3.0, -4.0], [1, 3, 4, 5]) 
Affine(54.6, [20.0, 1.0, 2.0], [1, 4, 6]) 
Affine(116.0, [-6.66, 18.0, -3.0, 10.5], [2, 3, 5, 6])

==============================

Affine(0.497556, [0.182576, 0.0303035, -0.0819013, 0.00912878, 0.0136502, -0.0295182, 0.0293768, 0.084641], [1, 2, 3, 4, 5, 6, 7, 8]) 
Affine(0.585568, [0.00182576, 0.0356316, -0.099953, 0.0273863, -0.0204649, -0.0561759, 0.0345419, 0.0277574], [1, 2, 3, 4, 5, 6, 9, 10]) 
Affine(-0.292556, [-0.107912, -0.0355038, 0.0977762, -0.0190003, 0.00220758, 0.0452743, -0.0144419, -0.0172241, -0.0695718, -0.0171939, -0.130834], [1, 2, 3, 4, 5, 6, 11, 12, 13, 14, 15]); 

Affine(1586.03, [580.0, -90.909, 245.7, 29.0, -40.95, 201.325, 216.795], [1, 2, 3, 4, 5, 6, 16]) 
Affine(1863.3, [5.8, -106.893, 277.3, 87.0, -164.15, 168.525, 71.004], [1, 2, 3, 4, 5, 6, 17]) 
Affine(876.705, [323.73, -5.46, 57.0, -54.6, 32.1, 43.325], [1, 3, 4, 5, 6, 18]); 

Affine(64.2, [0.2, -0.4, 3.0, -4.0], [1, 3, 4, 5]) 
Affine(54.6, [20.0, 1.0, 2.0], [1, 4, 6])
Affine(116.0, [-6.66, 18.0, -3.0, 10.5], [2, 3, 5, 6])
=#


